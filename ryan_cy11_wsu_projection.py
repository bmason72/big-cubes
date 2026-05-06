#!/usr/bin/env python3
"""Project WSU channel/data-rate growth for the ryanCy11 sample.

This analysis consumes the EB/SPW tabulation produced by ``spw_tabulate.py``
and evaluates several WSU projection methods against the current data:

1. ``memo_uniform_binned``:
   reproduce the memo procedure: use the finest current SPW resolution per
   MOUS, map it into the 5 stepped2 velocity bins, then assume every WSU SPW
   uses the floor of that bin.
2. ``memo_uniform_exact``:
   same as above but keep the current finest requested resolution instead of
   snapping it to the bin floor.
3. ``distributed_binned``:
   preserve the current within-MOUS SPW resolution distribution, but map each
   current SPW into the 5 stepped2 bins before projecting to WSU.
4. ``distributed_exact``:
   preserve the current within-MOUS SPW resolution distribution using the
   current requested resolution directly (still quantized to realizable TALON
   channels).

The script writes per-EB projected rates/factors plus aggregate summaries.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from statistics import median
from typing import Any

from config import DEFAULT_CONFIG
from wsu_projection import (
    alma_band_from_freq_hz,
    build_templates_from_rows,
    normalize_mous_uid,
    project_nchan_agg_for_templates,
    projected_spw_equivalents,
)


ANTENNA_RE = re.compile(r"^\s*Antennas:\s*(?P<nant>\d+):")
SCAN_ROW_RE = re.compile(
    r"^\s*(?:(?P<date>\d{2}-[A-Za-z]{3}-\d{4})/)?"
    r"(?P<start>\d{2}:\d{2}:\d{2}\.\d+)\s*-\s*"
    r"(?P<end>\d{2}:\d{2}:\d{2}\.\d+)\s+"
    r"(?P<scan>\d+)\s+"
    r"(?P<field_id>\d+)\s+"
    r"(?P<field_name>\S+)\s+"
    r"(?P<nrows>\d+)\s+"
    r"(?P<spw_ids>\[[^\]]*\])\s+"
    r"(?P<intervals>\[[^\]]*\])\s+"
    r"(?P<intents>\[[^\]]*\])\s*$"
)
CORR_TOKEN_RE = re.compile(r"[A-Z]{2}")

METHODS = (
    "memo_uniform_binned",
    "memo_uniform_exact",
    "distributed_binned",
    "distributed_exact",
)
MILESTONES = ("M1", "M4", "M5")


@dataclass
class EbRuntime:
    listobs_path: str
    actual_antenna_count: int | None
    array_type: str
    science_interval_median_s: float
    science_interval_by_spw_s: dict[int, float]
    science_time_s: float


def parse_int_list(text: str) -> list[int]:
    inner = text.strip()[1:-1].strip()
    if not inner:
        return []
    return [int(item.strip()) for item in inner.split(",") if item.strip()]


def parse_float_list(text: str) -> list[float]:
    inner = text.strip()[1:-1].strip()
    if not inner:
        return []
    return [float(item.strip()) for item in inner.split(",") if item.strip()]


def parse_token_list(text: str) -> list[str]:
    inner = text.strip()[1:-1].strip()
    if not inner:
        return []
    return [item.strip() for item in inner.split(",") if item.strip()]


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def finite(values: list[float]) -> list[float]:
    return [value for value in values if math.isfinite(value)]


def mean(values: list[float]) -> float:
    vals = finite(values)
    return sum(vals) / len(vals) if vals else math.nan


def resolve_local_path(path_str: str, search_roots: list[Path]) -> Path:
    raw = Path(path_str)
    if raw.is_absolute() and raw.exists():
        return raw
    for root in search_roots:
        candidate = (root / raw).resolve()
        if candidate.exists():
            return candidate
    return raw


def parse_time_of_day(value: str) -> float:
    stamp = datetime.strptime(value, "%H:%M:%S.%f")
    return stamp.hour * 3600.0 + stamp.minute * 60.0 + stamp.second + stamp.microsecond / 1.0e6


def parse_listobs_runtime(path: Path) -> EbRuntime:
    actual_antenna_count: int | None = None
    intervals_by_spw: dict[int, list[float]] = defaultdict(list)
    all_science_intervals: list[float] = []
    science_time_s = 0.0
    in_scans = False

    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.rstrip()
        stripped = line.strip()

        antenna_match = ANTENNA_RE.match(line)
        if antenna_match:
            actual_antenna_count = int(antenna_match.group("nant"))

        if stripped.startswith("Date        Timerange"):
            in_scans = True
            continue
        if in_scans and "(nRows =" in stripped:
            in_scans = False
            continue
        if not in_scans:
            continue

        match = SCAN_ROW_RE.match(line)
        if not match:
            continue
        intents = parse_token_list(match.group("intents"))
        if not any("OBSERVE_TARGET" in intent for intent in intents):
            continue

        start_s = parse_time_of_day(match.group("start"))
        end_s = parse_time_of_day(match.group("end"))
        if end_s < start_s:
            end_s += 24.0 * 3600.0
        science_time_s += end_s - start_s

        spw_ids = parse_int_list(match.group("spw_ids"))
        intervals = parse_float_list(match.group("intervals"))
        if len(intervals) == 1 and len(spw_ids) > 1:
            intervals = intervals * len(spw_ids)
        for spw_id, interval in zip(spw_ids, intervals):
            intervals_by_spw[spw_id].append(interval)
            all_science_intervals.append(interval)

    if actual_antenna_count is None:
        array_type = "unknown"
    elif actual_antenna_count <= 20:
        array_type = "7m"
    else:
        array_type = "12m"

    return EbRuntime(
        listobs_path=str(path),
        actual_antenna_count=actual_antenna_count,
        array_type=array_type,
        science_interval_median_s=float(median(all_science_intervals)) if all_science_intervals else math.nan,
        science_interval_by_spw_s={
            spw_id: float(median(values)) for spw_id, values in intervals_by_spw.items()
        },
        science_time_s=science_time_s,
    )


def count_corr_products(corr_products: str) -> int:
    tokens = CORR_TOKEN_RE.findall(corr_products)
    if tokens:
        return len(tokens)
    return len([part for part in corr_products.split() if part])


def projected_tint_s(array_type: str, current_interval_s: float, milestone: str) -> float:
    if array_type == "7m":
        return 9.984 if milestone in {"M1", "M4", "M5"} else 10.1
    if milestone == "M1":
        return 3.072 if current_interval_s <= 4.0 else 6.144
    return 3.072


def projected_nant(array_type: str, milestone: str) -> int:
    if array_type == "7m":
        return DEFAULT_CONFIG.arrays["7m"].nant_typical
    return DEFAULT_CONFIG.arrays["12m"].nant_typical


def datarate_gbps(nant: int, terms: list[tuple[float, float]], tint_s: float) -> float:
    if not math.isfinite(tint_s) or tint_s <= 0:
        return math.nan
    nbase = nant * (nant - 1.0) / 2.0
    cross_auto = 2.0 * 2.0 * 1.0 * nbase + 4.0 * nant
    return sum(cross_auto * channel_pol / tint / 1.0e9 for channel_pol, tint in terms)


def current_eb_rate_gbps(rows: list[dict[str, str]], runtime: EbRuntime) -> float:
    if runtime.actual_antenna_count is None:
        return math.nan
    terms: list[tuple[float, float]] = []
    for row in rows:
        spw_id = int(row["spw_id"])
        tint_s = runtime.science_interval_by_spw_s.get(spw_id, runtime.science_interval_median_s)
        nchan = int(row["nchan"])
        ncorr = count_corr_products(row["corr_products"])
        terms.append((nchan * ncorr, tint_s))
    return datarate_gbps(runtime.actual_antenna_count, terms, 1.0)


def projected_channel_pol_per_spw(
    method: str,
    milestone: str,
    templates,
    corr_counts: list[int],
) -> float:
    projected_nchan_agg = project_nchan_agg_for_templates(
        mode=method,
        milestone=milestone,
        templates=templates,
        projected_nspw=float(len(templates)),
    )
    mean_corr = sum(corr_counts) / len(corr_counts)
    return projected_nchan_agg * mean_corr / len(templates)


def summarize(values: list[float]) -> dict[str, float]:
    clean = finite(values)
    if not clean:
        return {"mean": math.nan, "median": math.nan}
    return {"mean": sum(clean) / len(clean), "median": float(median(clean))}


def format_factor(value: float) -> str:
    return "nan" if not math.isfinite(value) else f"{value:.3f}x"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tabulation_dir", type=Path, help="Directory containing spw_level_table.csv from spw_tabulate.py")
    parser.add_argument("output_dir", type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    spw_rows_all = read_csv_rows(args.tabulation_dir / "spw_level_table.csv")
    eb_unique: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in spw_rows_all:
        key = (row["mous_uid"], row["eb_uid"], row["spw_id"])
        eb_unique.setdefault(key, row)
    eb_spw_rows = sorted(eb_unique.values(), key=lambda row: (row["mous_uid"], row["eb_uid"], int(row["spw_id"])))

    mous_templates = build_templates_from_rows(eb_spw_rows)
    mous_corr_counts: dict[str, list[int]] = defaultdict(list)
    seen_corr_keys: set[tuple[str, int]] = set()
    for row in eb_spw_rows:
        mous_uid = normalize_mous_uid(row["mous_uid"])
        spw_id = int(row["spw_id"])
        key = (mous_uid, spw_id)
        if key in seen_corr_keys:
            continue
        seen_corr_keys.add(key)
        mous_corr_counts[mous_uid].append(count_corr_products(row["corr_products"]))
    eb_groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in eb_spw_rows:
        eb_groups[(row["mous_uid"], row["eb_uid"])].append(row)

    search_roots = [Path.cwd(), args.tabulation_dir.parent, Path(__file__).resolve().parent]
    runtime_cache: dict[str, EbRuntime] = {}
    eb_rows: list[dict[str, Any]] = []

    for (mous_uid, eb_uid), rows in sorted(eb_groups.items()):
        resolved_path = resolve_local_path(rows[0]["listobs_path"], search_roots)
        cache_key = str(resolved_path)
        if cache_key not in runtime_cache:
            runtime_cache[cache_key] = parse_listobs_runtime(resolved_path)
        runtime = runtime_cache[cache_key]
        normalized_mous_uid = normalize_mous_uid(mous_uid)
        templates = mous_templates[normalized_mous_uid]
        corr_counts = mous_corr_counts[normalized_mous_uid]
        band = templates[0].band
        current_rate = current_eb_rate_gbps(rows, runtime)
        current_volume = current_rate * runtime.science_time_s if math.isfinite(current_rate) else math.nan

        record: dict[str, Any] = {
            "mous_uid": mous_uid,
            "eb_uid": eb_uid,
            "array_type": runtime.array_type,
            "current_science_time_s": runtime.science_time_s,
            "current_interval_median_s": runtime.science_interval_median_s,
            "current_rate_gbps": current_rate,
            "current_volume_gb": current_volume,
            "listobs_path": runtime.listobs_path,
            "band": band,
        }

        for milestone in MILESTONES:
            eq_nspw = projected_spw_equivalents(milestone, band)
            tint_s = projected_tint_s(runtime.array_type, runtime.science_interval_median_s, milestone)
            nant = projected_nant(runtime.array_type, milestone)
            for method in METHODS:
                channel_pol_per_spw = projected_channel_pol_per_spw(
                    method, milestone, templates, corr_counts
                )
                total_channel_pol = eq_nspw * channel_pol_per_spw
                projected_rate = datarate_gbps(nant, [(total_channel_pol, tint_s)], 1.0)
                projected_volume = projected_rate * runtime.science_time_s if math.isfinite(projected_rate) else math.nan
                factor = projected_rate / current_rate if math.isfinite(projected_rate) and math.isfinite(current_rate) and current_rate > 0 else math.nan
                prefix = f"{milestone.lower()}_{method}"
                record[f"{prefix}_eq_nspw"] = eq_nspw
                record[f"{prefix}_rate_gbps"] = projected_rate
                record[f"{prefix}_volume_gb"] = projected_volume
                record[f"{prefix}_factor"] = factor
        eb_rows.append(record)

    eb_fieldnames = list(eb_rows[0].keys()) if eb_rows else []
    write_csv(args.output_dir / "eb_projection_factors.csv", eb_rows, eb_fieldnames)

    sidecar_rows = []
    for mous_uid, templates in sorted(mous_templates.items()):
        for template in templates:
            sidecar_rows.append(
                {
                    "mous_uid": mous_uid,
                    "spw_id": template.spw_id,
                    "center_freq_hz": template.center_freq_hz,
                    "velocity_resolution_kms": template.velocity_resolution_kms,
                    "band": template.band,
                }
            )
    write_csv(
        args.output_dir / "mous_spw_templates.csv",
        sidecar_rows,
        ["mous_uid", "spw_id", "center_freq_hz", "velocity_resolution_kms", "band"],
    )

    summary_rows: list[dict[str, Any]] = []
    by_array_rows: list[dict[str, Any]] = []
    for milestone in MILESTONES:
        for method in METHODS:
            prefix = f"{milestone.lower()}_{method}"
            factors = [float(row[f"{prefix}_factor"]) for row in eb_rows]
            projected_volumes = [float(row[f"{prefix}_volume_gb"]) for row in eb_rows]
            current_volumes = [float(row["current_volume_gb"]) for row in eb_rows]
            total_current = sum(value for value in current_volumes if math.isfinite(value))
            total_projected = sum(value for value in projected_volumes if math.isfinite(value))
            total_factor = total_projected / total_current if total_current > 0 else math.nan
            stats = summarize(factors)
            summary_rows.append(
                {
                    "milestone": milestone,
                    "method": method,
                    "n_eb": len(eb_rows),
                    "sample_total_factor": total_factor,
                    "mean_eb_factor": stats["mean"],
                    "median_eb_factor": stats["median"],
                    "current_total_volume_gb": total_current,
                    "projected_total_volume_gb": total_projected,
                }
            )
            for array_type in ("12m", "7m"):
                sub = [row for row in eb_rows if row["array_type"] == array_type]
                if not sub:
                    continue
                sub_factors = [float(row[f"{prefix}_factor"]) for row in sub]
                sub_proj = [float(row[f"{prefix}_volume_gb"]) for row in sub]
                sub_cur = [float(row["current_volume_gb"]) for row in sub]
                sub_current_total = sum(value for value in sub_cur if math.isfinite(value))
                sub_projected_total = sum(value for value in sub_proj if math.isfinite(value))
                sub_total_factor = sub_projected_total / sub_current_total if sub_current_total > 0 else math.nan
                sub_stats = summarize(sub_factors)
                by_array_rows.append(
                    {
                        "milestone": milestone,
                        "method": method,
                        "array_type": array_type,
                        "n_eb": len(sub),
                        "sample_total_factor": sub_total_factor,
                        "mean_eb_factor": sub_stats["mean"],
                        "median_eb_factor": sub_stats["median"],
                    }
                )

    write_csv(
        args.output_dir / "summary_overall.csv",
        summary_rows,
        [
            "milestone",
            "method",
            "n_eb",
            "sample_total_factor",
            "mean_eb_factor",
            "median_eb_factor",
            "current_total_volume_gb",
            "projected_total_volume_gb",
        ],
    )
    write_csv(
        args.output_dir / "summary_by_array.csv",
        by_array_rows,
        [
            "milestone",
            "method",
            "array_type",
            "n_eb",
            "sample_total_factor",
            "mean_eb_factor",
            "median_eb_factor",
        ],
    )

    rows_by_key = {(row["milestone"], row["method"]): row for row in summary_rows}
    summary_lines = [
        "# ryanCy11 WSU Projection Summary",
        "",
        "## Methods",
        "",
        "- `memo_uniform_binned`: memo procedure using finest current SPW per MOUS and the 5 stepped2 bins.",
        "- `memo_uniform_exact`: same, but keep the current finest requested resolution instead of snapping to the bin floor.",
        "- `distributed_binned`: preserve the current within-MOUS SPW resolution distribution, but map each SPW into the 5 stepped2 bins.",
        "- `distributed_exact`: preserve the current within-MOUS SPW resolution distribution using the current requested resolution directly.",
        "",
        "## Overall Results",
        "",
    ]
    for milestone in MILESTONES:
        summary_lines.append(f"### {milestone}")
        summary_lines.append("")
        for method in METHODS:
            row = rows_by_key[(milestone, method)]
            summary_lines.append(
                "- "
                + f"{method}: sample-total factor {format_factor(float(row['sample_total_factor']))}, "
                + f"mean EB factor {format_factor(float(row['mean_eb_factor']))}, "
                + f"median EB factor {format_factor(float(row['median_eb_factor']))}"
            )
        memo = rows_by_key[(milestone, "memo_uniform_binned")]
        uniform_exact = rows_by_key[(milestone, "memo_uniform_exact")]
        dist_binned = rows_by_key[(milestone, "distributed_binned")]
        dist_exact = rows_by_key[(milestone, "distributed_exact")]
        summary_lines.extend(
            [
                "",
                "- Change from keeping exact finest resolution instead of bin-flooring: "
                + format_factor(float(uniform_exact["sample_total_factor"]) / float(memo["sample_total_factor"]))
                + " relative to memo_uniform_binned.",
                "- Change from preserving the SPW resolution distribution at binned resolution: "
                + format_factor(float(dist_binned["sample_total_factor"]) / float(memo["sample_total_factor"]))
                + " relative to memo_uniform_binned.",
                "- Additional change from using exact per-SPW resolutions on top of preserving the distribution: "
                + format_factor(float(dist_exact["sample_total_factor"]) / float(dist_binned["sample_total_factor"]))
                + " relative to distributed_binned.",
                "",
            ]
        )

    (args.output_dir / "summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    (args.output_dir / "run_metadata.json").write_text(
        json.dumps(
            {
                "tabulation_dir": str(args.tabulation_dir),
                "n_eb": len(eb_rows),
                "methods": list(METHODS),
                "milestones": list(MILESTONES),
                "mous_spw_sidecar": str(args.output_dir / "mous_spw_templates.csv"),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
