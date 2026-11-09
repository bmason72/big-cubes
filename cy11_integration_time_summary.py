#!/usr/bin/env python3
"""Summarize observed Cycle-11 science visibility integration times by MOUS.

The preferred inputs are ``*_targets.ms__listobs.txt`` files.  When a MOUS
does not have those pipeline-transformed products, the full-MS listobs is used
and restricted to target-intent ``#FULL_RES`` spectral windows.  This avoids
mixing science visibility dump times with WVR, SQLD, channel-averaged, or
atmospheric products.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from config import ConfigClassifier
from cy11_geometry import parse_listobs_geometry
from spw_tabulate import extract_eb_tuning_rows, parse_listobs


TARGET_SUFFIX = "_targets.ms__listobs.txt"
TARGET_LINE_SUFFIX = "_targets_line.ms__listobs.txt"
FULL_SUFFIX = ".ms__listobs.txt"
RB_RE = re.compile(r"#ALMA_RB_(?P<band>\d{2})#")


@dataclass(frozen=True)
class PprMetadata:
    project_code: str = ""
    sb_uid: str = ""
    procedure: str = ""


def format_number(value: float) -> str:
    return f"{value:.9g}"


def format_bandwidth(total_bw_hz: float) -> str:
    if total_bw_hz >= 1.0e9:
        return f"{format_number(total_bw_hz / 1.0e9)}GHz"
    if total_bw_hz >= 1.0e6:
        return f"{format_number(total_bw_hz / 1.0e6)}MHz"
    if total_bw_hz >= 1.0e3:
        return f"{format_number(total_bw_hz / 1.0e3)}kHz"
    return f"{format_number(total_bw_hz)}Hz"


def format_counts(counts: Counter[float]) -> str:
    return ";".join(f"{format_number(value)}:{counts[value]}" for value in sorted(counts))


def receiver_band(spw_name: str) -> str:
    match = RB_RE.search(spw_name)
    return str(int(match.group("band"))) if match else ""


def is_full_listobs(path: Path) -> bool:
    return (
        path.name.endswith(FULL_SUFFIX)
        and not path.name.endswith(TARGET_SUFFIX)
        and not path.name.endswith(TARGET_LINE_SUFFIX)
    )


def choose_listobs_files(mous_dir: Path) -> tuple[list[Path], str]:
    listobs_dir = mous_dir / "listobs"
    target_files = sorted(listobs_dir.glob(f"*{TARGET_SUFFIX}"))
    if target_files:
        return target_files, "targets_ms"
    return sorted(path for path in listobs_dir.glob(f"*{FULL_SUFFIX}") if is_full_listobs(path)), "full_ms_fallback"


def parse_ppr_metadata(mous_dir: Path) -> PprMetadata:
    paths = sorted(mous_dir.glob("*.pprequest.xml"))
    if not paths:
        return PprMetadata()
    try:
        root = ET.parse(paths[0]).getroot()
    except (ET.ParseError, OSError):
        return PprMetadata()

    def text_at(path: str) -> str:
        elem = root.find(path)
        return (elem.text or "").strip() if elem is not None else ""

    sb_ref = root.find(".//SchedBlockRef")
    return PprMetadata(
        project_code=text_at(".//ProposalCode"),
        sb_uid=(sb_ref.attrib.get("entityId", "").strip() if sb_ref is not None else ""),
        procedure=text_at(".//ProcedureTitle"),
    )


def load_pipeline_l80(mous_dir: Path, mous_uid: str) -> dict[str, float]:
    paths = sorted(mous_dir.glob("pipeline_stats_*.json"))
    if not paths:
        return {}
    try:
        payload = json.loads(paths[0].read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}
    mous_payload = payload.get(mous_uid, {})
    result: dict[str, float] = {}
    for eb_name, eb_payload in mous_payload.get("EB", {}).items():
        try:
            result[eb_name.removesuffix(".ms")] = float(eb_payload["L80"]["value"])
        except (KeyError, TypeError, ValueError):
            continue
    return result


def classify_array(path: Path, procedure: str) -> tuple[str, int, float]:
    geometry = parse_listobs_geometry(path)
    if procedure.lower().startswith("hsd_"):
        return "TP", geometry.n_ant, geometry.l80_m
    if geometry.n_ant <= 0:
        return "unknown", geometry.n_ant, geometry.l80_m
    return geometry.array_type, geometry.n_ant, geometry.l80_m


def estimated_config(array_type: str, l80_m: float) -> str:
    if array_type == "TP":
        return "TP"
    if array_type not in {"7m", "12m"}:
        return "unknown"
    return ConfigClassifier().classify(array_type, l80_m)


def spw_signature(rows: list[dict[str, Any]]) -> str:
    tuning_groups: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        tuning_groups[int(row["tuning_group_index"])].append(row)
    parts: list[str] = []
    for index in sorted(tuning_groups):
        group = sorted(tuning_groups[index], key=lambda row: (float(row["center_freq_hz"]), int(row["spw_id"])))
        tokens = [
            f"B{receiver_band(str(row['spw_name'])) or '?'}:{int(row['nchan'])}ch@{format_bandwidth(float(row['total_bw_hz']))}"
            for row in group
        ]
        parts.append(f"TG{index:02d}[{','.join(tokens)}]")
    return ";".join(parts)


def finite_median(values: Iterable[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return statistics.median(finite) if finite else math.nan


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_index_statuses(path: Path) -> dict[str, set[str]]:
    statuses: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            uid = row.get("mous_uid", "").strip()
            if uid:
                statuses[uid].add(row.get("status", "").strip())
    return statuses


def build_catalog(
    input_dir: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[str]]:
    eb_spw_rows: list[dict[str, Any]] = []
    errors: list[str] = []

    for mous_dir in sorted(path for path in input_dir.glob("uid___*") if path.is_dir()):
        mous_uid = mous_dir.name
        files, input_kind = choose_listobs_files(mous_dir)
        if not files:
            continue
        ppr = parse_ppr_metadata(mous_dir)
        l80_by_eb = load_pipeline_l80(mous_dir, mous_uid)

        for listobs_path in files:
            parsed = parse_listobs(listobs_path, mous_uid)
            if any(issue.severity == "error" for issue in parsed.issues):
                errors.extend(f"{mous_uid},{listobs_path},{issue.message}" for issue in parsed.issues)
                continue
            _, tuning_spws = extract_eb_tuning_rows(parsed)
            if input_kind == "full_ms_fallback":
                tuning_spws = [row for row in tuning_spws if "#FULL_RES" in str(row["spw_name"])]
            if not tuning_spws:
                errors.append(f"{mous_uid},{listobs_path},no science visibility SPWs")
                continue

            array_type, n_ant, geometry_l80 = classify_array(listobs_path, ppr.procedure)
            l80_m = l80_by_eb.get(parsed.eb_uid, geometry_l80)
            l80_source = "pipeline_stats" if parsed.eb_uid in l80_by_eb else "listobs_geometry"
            config = estimated_config(array_type, l80_m)
            interval_values_by_spw: dict[int, set[float]] = defaultdict(set)
            for scan in parsed.scans:
                if not any("OBSERVE_TARGET" in intent for intent in scan.scan_intents):
                    continue
                intervals = scan.intervals_s
                if len(intervals) == 1 and len(scan.spw_ids) > 1:
                    intervals = intervals * len(scan.spw_ids)
                for spw_id, interval_s in zip(scan.spw_ids, intervals):
                    interval_values_by_spw[spw_id].add(interval_s)

            eb_signature = spw_signature(tuning_spws)
            for row in tuning_spws:
                spw_id = int(row["spw_id"])
                interval_values = sorted(interval_values_by_spw.get(spw_id, set()))
                eb_spw_rows.append(
                    {
                        "mous_uid": mous_uid,
                        "project_code": ppr.project_code,
                        "sb_uid": ppr.sb_uid,
                        "procedure": ppr.procedure,
                        "eb_uid": parsed.eb_uid,
                        "input_kind": input_kind,
                        "listobs_path": str(listobs_path),
                        "array": array_type,
                        "n_ant": n_ant,
                        "band": receiver_band(str(row["spw_name"])),
                        "l80_m": l80_m,
                        "l80_source": l80_source,
                        "estimated_config": config,
                        "tuning_group": row["tuning_group_id"],
                        "target_names": row["target_names"],
                        "spw_id": spw_id,
                        "nchan": int(row["nchan"]),
                        "total_bw_hz": float(row["total_bw_hz"]),
                        "median_interval_s": float(row["median_interval_s"]),
                        "interval_values_s": ";".join(format_number(value) for value in interval_values),
                        "mixed_interval_within_eb_spw": len(interval_values) > 1,
                        "eb_spectral_signature": eb_signature,
                    }
                )

    by_mous: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in eb_spw_rows:
        by_mous[str(row["mous_uid"])].append(row)

    mous_rows: list[dict[str, Any]] = []
    for mous_uid, rows in sorted(by_mous.items()):
        eb_uids = sorted({str(row["eb_uid"]) for row in rows})
        intervals = [float(row["median_interval_s"]) for row in rows]
        interval_counts = Counter(intervals)
        mode_count = max(interval_counts.values())
        mode_s = min(value for value, count in interval_counts.items() if count == mode_count)
        l80_values = [float(next(row["l80_m"] for row in rows if row["eb_uid"] == eb_uid)) for eb_uid in eb_uids]
        configs = sorted({str(row["estimated_config"]) for row in rows})
        bands = sorted({str(row["band"]) for row in rows if row["band"]}, key=int)
        arrays = sorted({str(row["array"]) for row in rows})
        target_names = sorted(
            {
                target
                for row in rows
                for target in str(row["target_names"]).split(";")
                if target
            }
        )
        signatures = Counter(
            str(next(row["eb_spectral_signature"] for row in rows if row["eb_uid"] == eb_uid))
            for eb_uid in eb_uids
        )
        signature_text = " || ".join(
            f"{count}EB:{signature}" for signature, count in sorted(signatures.items())
        )
        first = rows[0]
        mous_rows.append(
            {
                "mous_uid": mous_uid,
                "project_code": first["project_code"],
                "sb_uid": first["sb_uid"],
                "procedure": first["procedure"],
                "target_names": ";".join(target_names),
                "array": ";".join(arrays),
                "bands": ";".join(bands),
                "n_eb": len(eb_uids),
                "input_kinds": ";".join(sorted({str(row["input_kind"]) for row in rows})),
                "integration_time_mode_s": mode_s,
                "integration_time_median_s": finite_median(intervals),
                "integration_time_min_s": min(intervals),
                "integration_time_max_s": max(intervals),
                "integration_time_counts_eb_spw": format_counts(interval_counts),
                "mixed_integration_time": len(interval_counts) > 1,
                "l80_min_m": min(l80_values),
                "l80_median_m": finite_median(l80_values),
                "l80_max_m": max(l80_values),
                "l80_sources": ";".join(sorted({str(row["l80_source"]) for row in rows})),
                "estimated_configurations": ";".join(configs),
                "estimated_config_at_max_l80": estimated_config(str(first["array"]), max(l80_values)),
                "configuration_basis": "Cycle-13 nominal L80, log-nearest",
                "spectral_signature": signature_text,
                "spectral_signature_consistent_across_ebs": len(signatures) == 1,
                "n_science_spw_in_representative_eb": sum(1 for row in rows if row["eb_uid"] == eb_uids[0]),
            }
        )
    return mous_rows, eb_spw_rows, errors


def build_distribution_rows(mous_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    counts: Counter[tuple[str, str, str, float]] = Counter()
    totals: Counter[str] = Counter()
    for row in mous_rows:
        array = str(row["array"])
        totals[array] += 1
        for band in str(row["bands"]).split(";"):
            if band:
                counts[(array, band, str(row["estimated_config_at_max_l80"]), float(row["integration_time_mode_s"]))] += 1
    return [
        {
            "array": array,
            "band": band,
            "estimated_config_at_max_l80": config,
            "integration_time_s": interval,
            "n_mous": count,
            "fraction_of_array_mous": count / totals[array],
        }
        for (array, band, config, interval), count in sorted(counts.items())
    ]


def compact_spectral_signature(rows: list[dict[str, Any]]) -> str:
    tuning_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        tuning_groups[str(row["tuning_group"])].append(row)

    parts: list[str] = []
    for tuning_group in sorted(tuning_groups):
        spws = sorted(
            tuning_groups[tuning_group],
            key=lambda row: (int(row["spw_id"]), int(row["nchan"]), float(row["total_bw_hz"])),
        )
        signature = ",".join(
            f"{int(row['nchan'])}ch@{format_bandwidth(float(row['total_bw_hz']))}"
            for row in spws
        )
        parts.append(f"{tuning_group}[{signature}]")
    return ";".join(parts)


def build_compact_rows(eb_spw_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Build one concise row per MOUS and science receiver band."""
    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in eb_spw_rows:
        groups[(str(row["mous_uid"]), str(row["band"]))].append(row)

    compact_rows: list[dict[str, Any]] = []
    for (mous_uid, band), rows in sorted(groups.items(), key=lambda item: (item[0][0], int(item[0][1]))):
        intervals = sorted({float(row["median_interval_s"]) for row in rows})
        arrays = sorted({str(row["array"]) for row in rows})
        if len(intervals) != 1:
            raise ValueError(f"{mous_uid} Band {band} has multiple science integration times: {intervals}")
        if len(arrays) != 1:
            raise ValueError(f"{mous_uid} Band {band} has multiple array classifications: {arrays}")

        by_eb: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in rows:
            by_eb[str(row["eb_uid"])].append(row)
        signatures = sorted({compact_spectral_signature(eb_rows) for eb_rows in by_eb.values()})
        l80_values = [float(eb_rows[0]["l80_m"]) for eb_rows in by_eb.values()]
        l80_m = finite_median(l80_values)
        array_type = arrays[0]
        compact_rows.append(
            {
                "mous_uid": mous_uid,
                "integration_time_s": intervals[0],
                "array": array_type,
                "band": band,
                "spectral_signature": " || ".join(signatures),
                "l80_m": round(l80_m, 3),
                "estimated_config": estimated_config(array_type, l80_m),
            }
        )
    return compact_rows


def write_summary(path: Path, mous_rows: list[dict[str, Any]], coverage_rows: list[dict[str, Any]], errors: list[str]) -> None:
    array_counts = Counter(str(row["array"]) for row in mous_rows)
    tint_counts: dict[str, Counter[float]] = defaultdict(Counter)
    for row in mous_rows:
        tint_counts[str(row["array"])][float(row["integration_time_mode_s"])] += 1
    geometry_fallbacks = sum("listobs_geometry" in str(row["l80_sources"]) for row in mous_rows)
    multi_config = sum(";" in str(row["estimated_configurations"]) for row in mous_rows)
    inconsistent_setup = sum(not bool(row["spectral_signature_consistent_across_ebs"]) for row in mous_rows)

    lines = [
        "# Cycle-11 Science Integration Times",
        "",
        "Integration time means the listobs `Average Interval(s)` for target-intent science visibility SPWs. It is not scan duration or total on-source time. Values are the rounded values printed by CASA listobs.",
        "",
        "## Coverage",
        "",
        f"- Indexed MOUS: {len(coverage_rows)}",
        f"- MOUS summarized: {len(mous_rows)}",
        f"- EBs summarized: {sum(int(row['n_eb']) for row in mous_rows)}",
        f"- Arrays: {', '.join(f'{key}={value}' for key, value in sorted(array_counts.items()))}",
        f"- TP MOUS: {array_counts.get('TP', 0)} (none are present in this hifa-only extraction)",
        f"- MOUS using listobs-geometry L80 fallback for at least one EB: {geometry_fallbacks}",
        f"- MOUS spanning multiple estimated configurations: {multi_config}",
        f"- MOUS with differing exact setup signatures across EBs: {inconsistent_setup}",
        f"- Parse errors: {len(errors)}",
        "",
        "## MOUS Distribution",
        "",
        "| Array | Integration time (s) | MOUS | Fraction within array |",
        "|---|---:|---:|---:|",
    ]
    for array in sorted(tint_counts):
        total = array_counts[array]
        for interval, count in sorted(tint_counts[array].items()):
            lines.append(f"| {array} | {format_number(interval)} | {count} | {count / total:.3%} |")
    lines.extend(
        [
            "",
            "## Conventions",
            "",
            "- `*_targets.ms` is preferred; full-MS fallback rows retain only target-intent `#FULL_RES` SPWs.",
            "- Receiver band comes from the `ALMA_RB_XX` SPW name token because frequency-only band assignment is ambiguous in overlap regions.",
            "- Pipeline-statistics L80 is primary; antenna-coordinate geometry is the fallback.",
            "- Estimated 12m configurations use the nearest Cycle-13 nominal L80 in log space. They are estimates, not recovered Cycle-11 scheduling configuration names.",
            "- The local PPRs contain SB UIDs, but their human-readable `SBTitle` values are all `Undefined`; this local-only product therefore omits SB names.",
            "- `integration_time_counts_eb_spw` is encoded as `seconds:number_of_EB-SPW_rows`.",
            "",
            "## Products",
            "",
            "- `mous_integration_times.csv`: requested one-row-per-MOUS summary.",
            "- `eb_spw_integration_times.csv`: audit table retaining EB, tuning, and SPW provenance.",
            "- `integration_time_distribution.csv`: array/band/config distribution.",
            "- `coverage.csv`: disposition of every MOUS in the extraction index.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


MOUS_FIELDS = [
    "mous_uid", "project_code", "sb_uid", "procedure", "target_names", "array", "bands", "n_eb", "input_kinds",
    "integration_time_mode_s", "integration_time_median_s", "integration_time_min_s", "integration_time_max_s",
    "integration_time_counts_eb_spw", "mixed_integration_time", "l80_min_m", "l80_median_m", "l80_max_m",
    "l80_sources", "estimated_configurations", "estimated_config_at_max_l80", "configuration_basis",
    "spectral_signature", "spectral_signature_consistent_across_ebs", "n_science_spw_in_representative_eb",
]

EB_SPW_FIELDS = [
    "mous_uid", "project_code", "sb_uid", "procedure", "eb_uid", "input_kind", "listobs_path",
    "array", "n_ant", "band", "l80_m", "l80_source", "estimated_config", "tuning_group", "target_names", "spw_id", "nchan",
    "total_bw_hz", "median_interval_s", "interval_values_s", "mixed_interval_within_eb_spw", "eb_spectral_signature",
]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, nargs="?", default=Path("ryanCy11"))
    parser.add_argument("output_dir", type=Path, nargs="?", default=Path("ryanCy11out/integration_times"))
    args = parser.parse_args(argv)

    mous_rows, eb_spw_rows, errors = build_catalog(args.input_dir)
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "mous_integration_times.csv", mous_rows, MOUS_FIELDS)
    write_csv(
        output_dir / "cy11_mous_science_setups.csv",
        build_compact_rows(eb_spw_rows),
        ["mous_uid", "integration_time_s", "array", "band", "spectral_signature", "l80_m", "estimated_config"],
    )
    write_csv(output_dir / "eb_spw_integration_times.csv", eb_spw_rows, EB_SPW_FIELDS)
    write_csv(
        output_dir / "integration_time_distribution.csv",
        build_distribution_rows(mous_rows),
        ["array", "band", "estimated_config_at_max_l80", "integration_time_s", "n_mous", "fraction_of_array_mous"],
    )

    statuses = load_index_statuses(args.input_dir / "listobs_extraction_index.csv")
    included = {str(row["mous_uid"]) for row in mous_rows}
    coverage_rows = [
        {
            "mous_uid": uid,
            "included": uid in included,
            "index_statuses": ";".join(sorted(values)),
            "disposition": "summarized" if uid in included else "no usable listobs",
        }
        for uid, values in sorted(statuses.items())
    ]
    write_csv(output_dir / "coverage.csv", coverage_rows, ["mous_uid", "included", "index_statuses", "disposition"])
    (output_dir / "parse_errors.txt").write_text("\n".join(errors) + ("\n" if errors else ""), encoding="utf-8")
    write_summary(output_dir / "summary.md", mous_rows, coverage_rows, errors)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
