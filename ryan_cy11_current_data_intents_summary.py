#!/usr/bin/env python3
"""Tabulate current scan-level, intent-aware data properties for the ryanCy11 sidecar sample."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from statistics import median
from typing import Any

C_KMS = 299792.458
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
SPW_ROW_RE = re.compile(
    r"^\s*(?P<spw_id>\d+)\s+"
    r"(?P<name>\S+)\s+"
    r"(?P<nchan>\d+)\s+"
    r"(?P<frame>\S+)\s+"
    r"(?P<ch0_mhz>[-+]?\d+(?:\.\d+)?)\s+"
    r"(?P<chan_width_khz>[-+]?\d+(?:\.\d+)?)\s+"
    r"(?P<tot_bw_khz>[-+]?\d+(?:\.\d+)?)\s+"
    r"(?P<center_mhz>[-+]?\d+(?:\.\d+)?)\s+"
    r"(?P<bbc_num>\d+)\s+"
    r"(?P<corrs>\S.*\S|\S)\s*$"
)
TARGET_SUFFIXES = (
    "_targets.ms__listobs.txt",
    "_targets_line.ms__listobs.txt",
)
BUCKETS = (
    "science",
    "phase",
    "bandpass",
    "check_source",
    "pointing",
    "polarization",
    "diffgain_reference",
    "diffgain_on_source",
    "atmosphere",
    "other",
)
PLOT_BUCKETS = [bucket for bucket in BUCKETS if bucket != "other"]


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


def count_corr_products(corr_products: str) -> int:
    tokens = CORR_TOKEN_RE.findall(corr_products)
    if tokens:
        return len(tokens)
    return len([part for part in corr_products.split() if part])


def datarate_gbps(nant: int, terms: list[tuple[float, float]], tint_s: float) -> float:
    if not math.isfinite(tint_s) or tint_s <= 0:
        return math.nan
    nbase = nant * (nant - 1.0) / 2.0
    cross_auto = 2.0 * 2.0 * 1.0 * nbase + 4.0 * nant
    return sum(cross_auto * channel_pol / tint / 1.0e9 for channel_pol, tint in terms)


@dataclass
class FullSpwInfo:
    spw_id: int
    spw_name: str
    nchan: int
    frame: str
    ch0_mhz: float
    chan_width_khz: float
    total_bw_khz: float
    center_mhz: float
    corrs: str

    @property
    def center_freq_hz(self) -> float:
        return self.center_mhz * 1.0e6

    @property
    def abs_chan_width_hz(self) -> float:
        return abs(self.chan_width_khz) * 1.0e3

    @property
    def velocity_resolution_kms(self) -> float:
        if self.center_freq_hz <= 0:
            return math.nan
        return C_KMS * self.abs_chan_width_hz / self.center_freq_hz

    @property
    def is_sqld(self) -> bool:
        return "#SQLD" in self.spw_name or self.spw_name.endswith("SQLD")

    @property
    def is_ch_avg(self) -> bool:
        return "#CH_AVG" in self.spw_name or self.spw_name.endswith("CH_AVG")

    def signature_token(self) -> str:
        return (
            f"{self.spw_name}:{self.nchan}:{count_corr_products(self.corrs)}:"
            f"{self.center_mhz:.6f}:{self.velocity_resolution_kms:.9f}"
        )


@dataclass
class FullScan:
    scan_id: int
    field_id: int
    field_name: str
    duration_s: float
    spw_ids: list[int]
    intervals_s: list[float]
    intents: list[str]
    raw_intent_signature: str
    bucket: str


@dataclass
class FullEbData:
    listobs_path: str
    actual_antenna_count: int | None
    array_type: str
    spws: dict[int, FullSpwInfo]
    scans: list[FullScan]


def classify_bucket(intents: list[str]) -> str:
    intent_set = set(intents)
    if "OBSERVE_TARGET#ON_SOURCE" in intent_set:
        return "science"
    if "CALIBRATE_BANDPASS#ON_SOURCE" in intent_set or "CALIBRATE_FLUX#ON_SOURCE" in intent_set:
        return "bandpass"
    if "CALIBRATE_PHASE#ON_SOURCE" in intent_set:
        return "phase"
    if "OBSERVE_CHECK_SOURCE#ON_SOURCE" in intent_set:
        return "check_source"
    if "CALIBRATE_POINTING#ON_SOURCE" in intent_set:
        return "pointing"
    if "CALIBRATE_POLARIZATION#ON_SOURCE" in intent_set:
        return "polarization"
    if "CALIBRATE_DIFFGAIN#ON_SOURCE" in intent_set:
        return "diffgain_on_source"
    if "CALIBRATE_DIFFGAIN#REFERENCE" in intent_set:
        return "diffgain_reference"
    if any(intent.startswith("CALIBRATE_ATMOSPHERE#") for intent in intents):
        return "atmosphere"
    return "other"


def matched_full_listobs_path(path_str: str, search_roots: list[Path]) -> Path:
    target_path = resolve_local_path(path_str, search_roots)
    for suffix in TARGET_SUFFIXES:
        if target_path.name.endswith(suffix):
            full_name = target_path.name[: -len(suffix)] + ".ms__listobs.txt"
            candidate = target_path.with_name(full_name)
            if candidate.exists():
                return candidate
    return target_path


def percentile(values: list[float], frac: float) -> float:
    clean = sorted(value for value in values if math.isfinite(value))
    if not clean:
        return math.nan
    if len(clean) == 1:
        return clean[0]
    pos = frac * (len(clean) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return clean[lo]
    weight = pos - lo
    return clean[lo] * (1.0 - weight) + clean[hi] * weight


def fraction_summary(values: list[float]) -> dict[str, float]:
    clean = [value for value in values if math.isfinite(value)]
    if not clean:
        return {
            "mean": math.nan,
            "median": math.nan,
            "q25": math.nan,
            "q75": math.nan,
            "iqr": math.nan,
        }
    q25 = percentile(clean, 0.25)
    q75 = percentile(clean, 0.75)
    return {
        "mean": sum(clean) / len(clean),
        "median": float(median(clean)),
        "q25": q25,
        "q75": q75,
        "iqr": q75 - q25,
    }


def parse_full_listobs(path: Path) -> FullEbData:
    actual_antenna_count: int | None = None
    spws: dict[int, FullSpwInfo] = {}
    scans: list[FullScan] = []
    in_scans = False
    in_spws = False

    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.rstrip()
        stripped = line.strip()

        antenna_match = ANTENNA_RE.match(line)
        if antenna_match:
            actual_antenna_count = int(antenna_match.group("nant"))

        if stripped.startswith("Date        Timerange"):
            in_scans = True
            in_spws = False
            continue
        if stripped.startswith("Spectral Windows:"):
            in_spws = True
            in_scans = False
            continue
        if stripped.startswith("Sources:"):
            in_spws = False
            in_scans = False
            continue

        if in_scans:
            if "(nRows =" in stripped:
                in_scans = False
                continue
            match = SCAN_ROW_RE.match(line)
            if not match:
                continue
            start_s = parse_time_of_day(match.group("start"))
            end_s = parse_time_of_day(match.group("end"))
            if end_s < start_s:
                end_s += 24.0 * 3600.0
            spw_ids = parse_int_list(match.group("spw_ids"))
            intervals = parse_float_list(match.group("intervals"))
            if len(intervals) == 1 and len(spw_ids) > 1:
                intervals = intervals * len(spw_ids)
            intents = parse_token_list(match.group("intents"))
            scans.append(
                FullScan(
                    scan_id=int(match.group("scan")),
                    field_id=int(match.group("field_id")),
                    field_name=match.group("field_name"),
                    duration_s=end_s - start_s,
                    spw_ids=spw_ids,
                    intervals_s=intervals,
                    intents=intents,
                    raw_intent_signature=",".join(intents),
                    bucket=classify_bucket(intents),
                )
            )
            continue

        if in_spws:
            if stripped.startswith("SpwID") or not stripped:
                continue
            match = SPW_ROW_RE.match(line)
            if not match:
                continue
            spw = FullSpwInfo(
                spw_id=int(match.group("spw_id")),
                spw_name=match.group("name"),
                nchan=int(match.group("nchan")),
                frame=match.group("frame"),
                ch0_mhz=float(match.group("ch0_mhz")),
                chan_width_khz=float(match.group("chan_width_khz")),
                total_bw_khz=float(match.group("tot_bw_khz")),
                center_mhz=float(match.group("center_mhz")),
                corrs=match.group("corrs").strip(),
            )
            spws[spw.spw_id] = spw

    if actual_antenna_count is None:
        array_type = "unknown"
    elif actual_antenna_count <= 20:
        array_type = "7m"
    else:
        array_type = "12m"

    return FullEbData(
        listobs_path=str(path),
        actual_antenna_count=actual_antenna_count,
        array_type=array_type,
        spws=spws,
        scans=scans,
    )


def scan_metrics(scan: FullScan, eb: FullEbData) -> dict[str, Any]:
    retained_tokens: list[str] = []
    retained_velres: list[float] = []
    retained_terms: list[tuple[float, float]] = []
    retained_spw_ids: list[int] = []
    excluded_sqld = 0
    excluded_ch_avg = 0
    missing_spw = 0
    retained_nchan_sum = 0

    intervals = scan.intervals_s
    if len(intervals) == 1 and len(scan.spw_ids) > 1:
        intervals = intervals * len(scan.spw_ids)

    for spw_id, interval_s in zip(scan.spw_ids, intervals):
        spw = eb.spws.get(spw_id)
        if spw is None:
            missing_spw += 1
            continue
        if spw.is_sqld:
            excluded_sqld += 1
            continue
        if spw.is_ch_avg:
            excluded_ch_avg += 1
            continue
        retained_spw_ids.append(spw_id)
        retained_tokens.append(spw.signature_token())
        retained_velres.append(spw.velocity_resolution_kms)
        retained_nchan_sum += spw.nchan
        retained_terms.append((spw.nchan * count_corr_products(spw.corrs), interval_s))

    if eb.actual_antenna_count is None:
        rate_gbps = math.nan
    elif retained_terms:
        rate_gbps = datarate_gbps(eb.actual_antenna_count, retained_terms, 1.0)
    else:
        rate_gbps = 0.0

    volume_gb = rate_gbps * scan.duration_s if math.isfinite(rate_gbps) else math.nan
    retained_signature = "|".join(sorted(retained_tokens)) if retained_tokens else "<no-retained-spw>"
    median_velres = float(median(retained_velres)) if retained_velres else math.nan
    return {
        "retained_spw_ids": ",".join(str(spw_id) for spw_id in retained_spw_ids),
        "retained_spw_count": len(retained_spw_ids),
        "excluded_sqld_count": excluded_sqld,
        "excluded_ch_avg_count": excluded_ch_avg,
        "missing_spw_count": missing_spw,
        "retained_nchan_sum": retained_nchan_sum,
        "median_retained_velocity_resolution_kms": median_velres,
        "retained_signature": retained_signature,
        "rate_gbps": rate_gbps,
        "volume_gb": volume_gb,
    }


def bucket_zero_record(key_fields: dict[str, Any], bucket: str) -> dict[str, Any]:
    row = dict(key_fields)
    row.update(
        {
            "bucket": bucket,
            "scan_count": 0,
            "time_s": 0.0,
            "time_fraction": 0.0,
            "volume_gb": 0.0,
            "volume_fraction": 0.0,
            "signature_count": 0,
        }
    )
    return row


def aggregate_bucket_rows(
    base_rows: list[dict[str, Any]],
    group_fields: list[str],
) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in base_rows:
        key = tuple(row[field] for field in group_fields)
        grouped[key].append(row)

    out_rows: list[dict[str, Any]] = []
    for key, rows in sorted(grouped.items()):
        key_fields = {field: value for field, value in zip(group_fields, key)}
        total_time_s = sum(float(row["duration_s"]) for row in rows if math.isfinite(float(row["duration_s"])))
        total_volume_gb = sum(float(row["volume_gb"]) for row in rows if math.isfinite(float(row["volume_gb"])))
        by_bucket: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in rows:
            by_bucket[str(row["bucket"])].append(row)
        for bucket in BUCKETS:
            bucket_rows = by_bucket.get(bucket, [])
            if not bucket_rows:
                out_rows.append(bucket_zero_record(key_fields, bucket))
                continue
            time_s = sum(float(row["duration_s"]) for row in bucket_rows if math.isfinite(float(row["duration_s"])))
            volume_gb = sum(float(row["volume_gb"]) for row in bucket_rows if math.isfinite(float(row["volume_gb"])))
            signatures = {str(row["retained_signature"]) for row in bucket_rows}
            record = dict(key_fields)
            record.update(
                {
                    "bucket": bucket,
                    "scan_count": len(bucket_rows),
                    "time_s": time_s,
                    "time_fraction": time_s / total_time_s if total_time_s > 0 else math.nan,
                    "volume_gb": volume_gb,
                    "volume_fraction": volume_gb / total_volume_gb if total_volume_gb > 0 else math.nan,
                    "signature_count": len(signatures),
                }
            )
            out_rows.append(record)
    return out_rows


def build_fraction_summary_rows(rows: list[dict[str, Any]], level_name: str) -> list[dict[str, Any]]:
    out_rows: list[dict[str, Any]] = []
    for bucket in BUCKETS:
        bucket_values = [float(row["volume_fraction"]) for row in rows if row["bucket"] == bucket]
        stats = fraction_summary(bucket_values)
        out_rows.append(
            {
                "level": level_name,
                "bucket": bucket,
                "n": len(bucket_values),
                "mean_volume_fraction": stats["mean"],
                "median_volume_fraction": stats["median"],
                "q25_volume_fraction": stats["q25"],
                "q75_volume_fraction": stats["q75"],
                "iqr_volume_fraction": stats["iqr"],
            }
        )
    return out_rows


def render_markdown_table(rows: list[dict[str, Any]], headers: list[tuple[str, str]]) -> list[str]:
    lines = ["| " + " | ".join(label for _, label in headers) + " |"]
    lines.append("| " + " | ".join("---" for _ in headers) + " |")
    for row in rows:
        formatted: list[str] = []
        for key, _label in headers:
            value = row[key]
            if isinstance(value, float):
                formatted.append("nan" if not math.isfinite(value) else f"{value:.4f}")
            else:
                formatted.append(str(value))
        lines.append("| " + " | ".join(formatted) + " |")
    return lines


def write_fraction_plots(
    eb_rows: list[dict[str, Any]],
    mous_rows: list[dict[str, Any]],
    output_dir: Path,
) -> dict[str, Any]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        return {"plots_written": False, "reason": f"matplotlib unavailable: {exc}"}

    plot_dir = output_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    for level_name, rows in (("eb", eb_rows), ("mous", mous_rows)):
        data = [
            [float(row["volume_fraction"]) for row in rows if row["bucket"] == bucket]
            for bucket in PLOT_BUCKETS
        ]
        fig, ax = plt.subplots(figsize=(11, 5))
        ax.boxplot(data, labels=PLOT_BUCKETS, showfliers=False)
        ax.set_ylim(0.0, 1.0)
        ax.set_ylabel("Current volume fraction")
        ax.set_title(f"{level_name.upper()} current volume fraction by bucket")
        ax.tick_params(axis="x", rotation=25)
        fig.tight_layout()
        fig.savefig(plot_dir / f"{level_name}_bucket_volume_fraction_boxplot.png", dpi=150)
        plt.close(fig)

    return {"plots_written": True, "plot_dir": str(plot_dir)}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tabulation_dir", type=Path, help="Directory containing spw_level_table.csv")
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--skip-plots", action="store_true", help="Do not attempt to generate matplotlib plots")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    spw_rows = read_csv_rows(args.tabulation_dir / "spw_level_table.csv")
    eb_sample_rows: dict[tuple[str, str], dict[str, str]] = {}
    for row in spw_rows:
        eb_sample_rows.setdefault((row["mous_uid"], row["eb_uid"]), row)

    search_roots = [Path.cwd(), args.tabulation_dir.parent, Path(__file__).resolve().parent]
    full_cache: dict[str, FullEbData] = {}
    sample_scan_rows: list[dict[str, Any]] = []
    sample_metadata_rows: list[dict[str, Any]] = []
    missing_full_paths: list[dict[str, Any]] = []

    for (mous_uid, eb_uid), row in sorted(eb_sample_rows.items()):
        full_path = matched_full_listobs_path(row["listobs_path"], search_roots)
        if not full_path.exists():
            missing_full_paths.append(
                {"mous_uid": mous_uid, "eb_uid": eb_uid, "expected_full_listobs_path": str(full_path)}
            )
            continue
        cache_key = str(full_path)
        if cache_key not in full_cache:
            full_cache[cache_key] = parse_full_listobs(full_path)
        full_data = full_cache[cache_key]
        sample_metadata_rows.append(
            {
                "mous_uid": mous_uid,
                "eb_uid": eb_uid,
                "array_type": full_data.array_type,
                "target_listobs_path": row["listobs_path"],
                "full_listobs_path": str(full_path),
                "actual_antenna_count": full_data.actual_antenna_count,
                "n_scans": len(full_data.scans),
                "n_spws": len(full_data.spws),
            }
        )
        for scan in full_data.scans:
            metrics = scan_metrics(scan, full_data)
            sample_scan_rows.append(
                {
                    "mous_uid": mous_uid,
                    "eb_uid": eb_uid,
                    "array_type": full_data.array_type,
                    "target_listobs_path": row["listobs_path"],
                    "full_listobs_path": str(full_path),
                    "scan_id": scan.scan_id,
                    "field_id": scan.field_id,
                    "field_name": scan.field_name,
                    "bucket": scan.bucket,
                    "raw_intent_signature": scan.raw_intent_signature,
                    "duration_s": scan.duration_s,
                    **metrics,
                }
            )

    write_csv(
        args.output_dir / "sample_ebs.csv",
        sample_metadata_rows,
        [
            "mous_uid",
            "eb_uid",
            "array_type",
            "target_listobs_path",
            "full_listobs_path",
            "actual_antenna_count",
            "n_scans",
            "n_spws",
        ],
    )
    write_csv(
        args.output_dir / "missing_full_listobs.csv",
        missing_full_paths,
        ["mous_uid", "eb_uid", "expected_full_listobs_path"],
    )
    write_csv(
        args.output_dir / "scan_level_table.csv",
        sample_scan_rows,
        [
            "mous_uid",
            "eb_uid",
            "array_type",
            "target_listobs_path",
            "full_listobs_path",
            "scan_id",
            "field_id",
            "field_name",
            "bucket",
            "raw_intent_signature",
            "duration_s",
            "retained_spw_ids",
            "retained_spw_count",
            "excluded_sqld_count",
            "excluded_ch_avg_count",
            "missing_spw_count",
            "retained_nchan_sum",
            "median_retained_velocity_resolution_kms",
            "retained_signature",
            "rate_gbps",
            "volume_gb",
        ],
    )

    raw_grouped: dict[tuple[str, str, str, str, str], dict[str, Any]] = {}
    for row in sample_scan_rows:
        key = (
            str(row["mous_uid"]),
            str(row["eb_uid"]),
            str(row["array_type"]),
            str(row["bucket"]),
            str(row["raw_intent_signature"]),
        )
        record = raw_grouped.setdefault(
            key,
            {
                "mous_uid": row["mous_uid"],
                "eb_uid": row["eb_uid"],
                "array_type": row["array_type"],
                "bucket": row["bucket"],
                "raw_intent_signature": row["raw_intent_signature"],
                "scan_count": 0,
                "time_s": 0.0,
                "volume_gb": 0.0,
            },
        )
        record["scan_count"] += 1
        record["time_s"] += float(row["duration_s"])
        if math.isfinite(float(row["volume_gb"])):
            record["volume_gb"] += float(row["volume_gb"])

    raw_intent_rows = sorted(raw_grouped.values(), key=lambda item: (item["mous_uid"], item["eb_uid"], item["bucket"], item["raw_intent_signature"]))
    write_csv(
        args.output_dir / "raw_intent_by_eb.csv",
        raw_intent_rows,
        ["mous_uid", "eb_uid", "array_type", "bucket", "raw_intent_signature", "scan_count", "time_s", "volume_gb"],
    )

    bucket_by_eb_rows = aggregate_bucket_rows(sample_scan_rows, ["mous_uid", "eb_uid", "array_type"])
    bucket_by_mous_rows = aggregate_bucket_rows(sample_scan_rows, ["mous_uid"])
    write_csv(
        args.output_dir / "bucket_by_eb.csv",
        bucket_by_eb_rows,
        ["mous_uid", "eb_uid", "array_type", "bucket", "scan_count", "time_s", "time_fraction", "volume_gb", "volume_fraction", "signature_count"],
    )
    write_csv(
        args.output_dir / "bucket_by_mous.csv",
        bucket_by_mous_rows,
        ["mous_uid", "bucket", "scan_count", "time_s", "time_fraction", "volume_gb", "volume_fraction", "signature_count"],
    )

    fraction_summary_rows = build_fraction_summary_rows(bucket_by_eb_rows, "eb")
    fraction_summary_rows.extend(build_fraction_summary_rows(bucket_by_mous_rows, "mous"))
    write_csv(
        args.output_dir / "bucket_fraction_summary.csv",
        fraction_summary_rows,
        ["level", "bucket", "n", "mean_volume_fraction", "median_volume_fraction", "q25_volume_fraction", "q75_volume_fraction", "iqr_volume_fraction"],
    )

    signature_grouped: dict[tuple[str, str], dict[str, Any]] = {}
    signature_ebs: dict[tuple[str, str], set[str]] = defaultdict(set)
    signature_mous: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in sample_scan_rows:
        key = (str(row["bucket"]), str(row["retained_signature"]))
        record = signature_grouped.setdefault(
            key,
            {
                "bucket": row["bucket"],
                "retained_signature": row["retained_signature"],
                "scan_count": 0,
                "time_s": 0.0,
                "volume_gb": 0.0,
            },
        )
        record["scan_count"] += 1
        record["time_s"] += float(row["duration_s"])
        if math.isfinite(float(row["volume_gb"])):
            record["volume_gb"] += float(row["volume_gb"])
        signature_ebs[key].add(f"{row['mous_uid']}::{row['eb_uid']}")
        signature_mous[key].add(str(row["mous_uid"]))

    signature_rows = []
    for key, record in sorted(signature_grouped.items(), key=lambda item: (item[1]["bucket"], -float(item[1]["volume_gb"]))):
        signature_rows.append(
            {
                **record,
                "n_eb": len(signature_ebs[key]),
                "n_mous": len(signature_mous[key]),
            }
        )
    write_csv(
        args.output_dir / "bucket_signature_summary.csv",
        signature_rows,
        ["bucket", "retained_signature", "scan_count", "n_eb", "n_mous", "time_s", "volume_gb"],
    )

    sanity_rows: list[dict[str, Any]] = []
    bucket_rows_by_eb: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in bucket_by_eb_rows:
        bucket_rows_by_eb[(str(row["mous_uid"]), str(row["eb_uid"]))].append(row)
    for (mous_uid, eb_uid), rows in sorted(bucket_rows_by_eb.items()):
        time_fraction_sum = sum(float(row["time_fraction"]) for row in rows if math.isfinite(float(row["time_fraction"])))
        volume_fraction_sum = sum(float(row["volume_fraction"]) for row in rows if math.isfinite(float(row["volume_fraction"])))
        total_time_s = sum(float(row["time_s"]) for row in rows if math.isfinite(float(row["time_s"])))
        total_volume_gb = sum(float(row["volume_gb"]) for row in rows if math.isfinite(float(row["volume_gb"])))
        sanity_rows.append(
            {
                "mous_uid": mous_uid,
                "eb_uid": eb_uid,
                "time_fraction_sum": time_fraction_sum,
                "volume_fraction_sum": volume_fraction_sum,
                "time_fraction_error": abs(time_fraction_sum - 1.0) if total_time_s > 0 else math.nan,
                "volume_fraction_error": abs(volume_fraction_sum - 1.0) if total_volume_gb > 0 else math.nan,
                "total_time_s": total_time_s,
                "total_volume_gb": total_volume_gb,
            }
        )
    write_csv(
        args.output_dir / "sanity_checks.csv",
        sanity_rows,
        ["mous_uid", "eb_uid", "time_fraction_sum", "volume_fraction_sum", "time_fraction_error", "volume_fraction_error", "total_time_s", "total_volume_gb"],
    )

    if args.skip_plots:
        plot_status = {"plots_written": False, "reason": "skipped by flag"}
    else:
        plot_status = write_fraction_plots(bucket_by_eb_rows, bucket_by_mous_rows, args.output_dir)

    summary_by_level: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in fraction_summary_rows:
        summary_by_level[str(row["level"])].append(row)
    diffgain_eb_rows = [
        row for row in bucket_by_eb_rows
        if row["bucket"] in {"diffgain_reference", "diffgain_on_source"} and float(row["volume_gb"]) > 0.0
    ]
    diffgain_mous = sorted({str(row["mous_uid"]) for row in diffgain_eb_rows})
    diffgain_eb = sorted({str(row["eb_uid"]) for row in diffgain_eb_rows})

    summary_lines = [
        "# Current Data Property Summary",
        "",
        f"- Sample EBs analyzed: {len(sample_metadata_rows)}",
        f"- Sample MOUS analyzed: {len(sorted({row['mous_uid'] for row in sample_metadata_rows}))}",
        f"- Missing matched full-MS listobs: {len(missing_full_paths)}",
        f"- Plot status: {'written' if plot_status.get('plots_written') else plot_status.get('reason', 'unknown')}",
        "",
        "## EB-Level Volume Fractions",
        "",
    ]
    summary_lines.extend(
        render_markdown_table(
            summary_by_level["eb"],
            [
                ("bucket", "Bucket"),
                ("mean_volume_fraction", "Mean"),
                ("median_volume_fraction", "Median"),
                ("q25_volume_fraction", "Q25"),
                ("q75_volume_fraction", "Q75"),
                ("iqr_volume_fraction", "IQR"),
            ],
        )
    )
    summary_lines.extend(
        [
            "",
            "## MOUS-Level Volume Fractions",
            "",
        ]
    )
    summary_lines.extend(
        render_markdown_table(
            summary_by_level["mous"],
            [
                ("bucket", "Bucket"),
                ("mean_volume_fraction", "Mean"),
                ("median_volume_fraction", "Median"),
                ("q25_volume_fraction", "Q25"),
                ("q75_volume_fraction", "Q75"),
                ("iqr_volume_fraction", "IQR"),
            ],
        )
    )
    max_time_error = max(
        (float(row["time_fraction_error"]) for row in sanity_rows if math.isfinite(float(row["time_fraction_error"]))),
        default=math.nan,
    )
    max_volume_error = max(
        (float(row["volume_fraction_error"]) for row in sanity_rows if math.isfinite(float(row["volume_fraction_error"]))),
        default=math.nan,
    )
    summary_lines.extend(
        [
            "",
            "## Sanity Checks",
            "",
            f"- Max EB time-fraction closure error: {'nan' if not math.isfinite(max_time_error) else f'{max_time_error:.6g}'}",
            f"- Max EB volume-fraction closure error: {'nan' if not math.isfinite(max_volume_error) else f'{max_volume_error:.6g}'}",
            "",
            "## B2B / DIFFGAIN",
            "",
            f"- EBs with non-zero DIFFGAIN retained volume: {len(diffgain_eb)}",
            f"- MOUS with non-zero DIFFGAIN retained volume: {len(diffgain_mous)}",
            "",
            "## Key Artifacts",
            "",
            f"- `scan_level_table.csv`",
            f"- `bucket_by_eb.csv`",
            f"- `bucket_by_mous.csv`",
            f"- `bucket_fraction_summary.csv`",
            f"- `bucket_signature_summary.csv`",
            f"- `raw_intent_by_eb.csv`",
            f"- `sanity_checks.csv`",
        ]
    )

    (args.output_dir / "summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    (args.output_dir / "run_metadata.json").write_text(
        json.dumps(
            {
                "tabulation_dir": str(args.tabulation_dir),
                "n_sample_eb": len(sample_metadata_rows),
                "n_sample_mous": len(sorted({row["mous_uid"] for row in sample_metadata_rows})),
                "n_missing_full_listobs": len(missing_full_paths),
                "buckets": list(BUCKETS),
                "plots": plot_status,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
