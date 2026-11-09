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
from ryan_cy11_current_data_intents_summary import matched_full_listobs_path, parse_full_listobs, scan_metrics
from wsu_projection import (
    MousSpwTemplate,
    TuningGroupTemplate,
    alma_band_from_freq_hz,
    build_templates_from_rows,
    build_tuning_groups_from_rows,
    normalize_mous_uid,
    project_nchan_agg_for_templates,
    project_nchan_agg_for_tuning_groups,
    projected_nchan_per_spw,
    projected_spw_equivalents,
    stepped2_floor_velocity,
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
INTENT_PROJECTION_BUCKETS = ("science", "bandpass", "phase", "check_source")
CALIBRATOR_CAP_BUCKETS = ("bandpass", "phase", "check_source")


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
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
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
    seen_scan_ids: set[int] = set()

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
        scan_id = int(match.group("scan"))
        if scan_id in seen_scan_ids:
            continue
        seen_scan_ids.add(scan_id)

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


def projected_channel_pol_agg_for_templates(
    method: str,
    milestone: str,
    templates: list[MousSpwTemplate],
    corr_counts: list[int],
    projected_nspw: float,
    velocity_cap_kms: float | None = None,
) -> float:
    if not templates:
        return math.nan
    if len(templates) != len(corr_counts):
        raise ValueError("templates and corr_counts length mismatch")

    if method.startswith("memo_uniform"):
        finest = min(template.velocity_resolution_kms for template in templates)
        target = finest if method.endswith("exact") else stepped2_floor_velocity(finest)
        if velocity_cap_kms is not None:
            target = max(target, velocity_cap_kms)
        mean_freq = sum(template.center_freq_hz for template in templates) / len(templates)
        # the mean of SPWs from different bands can fall in an inter-band gap
        band = min(templates, key=lambda t: abs(t.center_freq_hz - mean_freq)).band
        nchan = projected_nchan_per_spw(
            requested_vel_kms=target,
            original_vel_kms=finest,
            center_freq_hz=mean_freq,
            band=band,
            milestone=milestone,
        )
        mean_corr = sum(corr_counts) / len(corr_counts)
        return projected_nspw * nchan * mean_corr

    projected_channel_pol = []
    for template, corr_count in zip(templates, corr_counts):
        target = (
            template.velocity_resolution_kms
            if method.endswith("exact")
            else stepped2_floor_velocity(template.velocity_resolution_kms)
        )
        if velocity_cap_kms is not None:
            target = max(target, velocity_cap_kms)
        projected_channel_pol.append(
            projected_nchan_per_spw(
                requested_vel_kms=target,
                original_vel_kms=template.velocity_resolution_kms,
                center_freq_hz=template.center_freq_hz,
                band=template.band,
                milestone=milestone,
            )
            * corr_count
        )
    return projected_nspw * (sum(projected_channel_pol) / len(projected_channel_pol))


def build_eb_tuning_groups(rows: list[dict[str, str]]) -> list[TuningGroupTemplate]:
    groups_by_mous = build_tuning_groups_from_rows(rows)
    if not rows:
        return []
    mous_uid = normalize_mous_uid(rows[0]["mous_uid"])
    return groups_by_mous.get(mous_uid, [])


def current_eb_metrics_from_tuning_rows(
    tuning_summary_rows: list[dict[str, str]],
    tuning_spw_rows: list[dict[str, str]],
    runtime: EbRuntime,
) -> tuple[float, float, float]:
    """Return current mean rate, total volume, and total science time for one EB."""
    if runtime.actual_antenna_count is None:
        return math.nan, math.nan, math.nan

    summary_by_group = {
        str(row["tuning_group_id"]): row for row in tuning_summary_rows
    }
    rows_by_group: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in tuning_spw_rows:
        rows_by_group[str(row["tuning_group_id"])].append(row)

    total_time_s = 0.0
    total_volume_gb = 0.0
    for tuning_group_id, rows in rows_by_group.items():
        summary_row = summary_by_group[tuning_group_id]
        terms: list[tuple[float, float]] = []
        for row in rows:
            interval_s = float(row["median_interval_s"])
            if not math.isfinite(interval_s) or interval_s <= 0.0:
                interval_s = runtime.science_interval_median_s
            terms.append(
                (
                    int(row["nchan"]) * count_corr_products(row["corr_products"]),
                    interval_s,
                )
            )
        group_rate = datarate_gbps(runtime.actual_antenna_count, terms, 1.0)
        group_time_s = float(summary_row["total_scan_time_s"])
        total_time_s += group_time_s
        if math.isfinite(group_rate):
            total_volume_gb += group_rate * group_time_s
    if total_time_s <= 0.0:
        return math.nan, math.nan, math.nan
    return total_volume_gb / total_time_s, total_volume_gb, total_time_s


def current_eb_nchan_agg_from_tuning_rows(tuning_summary_rows: list[dict[str, str]]) -> float:
    if not tuning_summary_rows:
        return math.nan
    return float(tuning_summary_rows[0]["eb_current_nchan_agg"])


def retained_scan_templates(scan, full_eb) -> tuple[list[MousSpwTemplate], list[int]]:
    templates: list[MousSpwTemplate] = []
    corr_counts: list[int] = []
    for spw_id in scan.spw_ids:
        spw = full_eb.spws.get(spw_id)
        if spw is None or spw.is_sqld or spw.is_ch_avg:
            continue
        try:
            band = alma_band_from_freq_hz(spw.center_freq_hz)
        except ValueError:
            continue
        templates.append(
            MousSpwTemplate(
                mous_uid="scan",
                spw_id=spw.spw_id,
                center_freq_hz=spw.center_freq_hz,
                velocity_resolution_kms=spw.velocity_resolution_kms,
                band=band,
            )
        )
        corr_counts.append(count_corr_products(spw.corrs))
    return templates, corr_counts


def cap_ratio_for_scan(
    method: str,
    milestone: str,
    scan,
    full_eb,
    velocity_cap_kms: float,
) -> float:
    templates, corr_counts = retained_scan_templates(scan, full_eb)
    if not templates:
        return 1.0
    uncapped = projected_channel_pol_agg_for_templates(
        method=method,
        milestone=milestone,
        templates=templates,
        corr_counts=corr_counts,
        projected_nspw=1.0,
    )
    capped = projected_channel_pol_agg_for_templates(
        method=method,
        milestone=milestone,
        templates=templates,
        corr_counts=corr_counts,
        projected_nspw=1.0,
        velocity_cap_kms=velocity_cap_kms,
    )
    if not math.isfinite(uncapped) or uncapped <= 0.0 or not math.isfinite(capped):
        return 1.0
    return capped / uncapped


def summarize(values: list[float]) -> dict[str, float]:
    clean = finite(values)
    if not clean:
        return {"mean": math.nan, "median": math.nan}
    return {"mean": sum(clean) / len(clean), "median": float(median(clean))}


def format_factor(value: float) -> str:
    return "nan" if not math.isfinite(value) else f"{value:.3f}x"


def latex_factor(value: float) -> str:
    return "nan" if not math.isfinite(value) else f"{value:.3g}x"


def write_latex_factor_table(
    path: Path,
    rows: list[dict[str, Any]],
    caption: str,
    label: str,
) -> None:
    lines = [
        r"\begin{table}[h]",
        r"\centering",
        rf"\caption{{{caption}}}",
        rf"\label{{{label}}}",
        r"\begin{tabular}{llccc}",
        r"\toprule",
        r"Milestone & Method & Total & Mean EB & Median EB \\",
        r"\midrule",
    ]
    for milestone in MILESTONES:
        for i, method in enumerate(METHODS):
            row = next(item for item in rows if item["milestone"] == milestone and item["method"] == method)
            lines.append(
                f"{milestone} & {method.replace('_', r'\_')} & "
                f"{latex_factor(float(row['sample_total_factor']))} & "
                f"{latex_factor(float(row['mean_eb_factor']))} & "
                f"{latex_factor(float(row['median_eb_factor']))} \\\\"
            )
        if milestone != MILESTONES[-1]:
            lines.append(r"\addlinespace")
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\end{table}",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def build_intent_aware_rows(
    eb_rows: list[dict[str, Any]],
    eb_groups: dict[tuple[str, str], list[dict[str, str]]],
    search_roots: list[Path],
    calibrator_cap_kms: float,
    apply_calibrator_cap: bool,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    full_cache: dict[str, Any] = {}
    detail_rows: list[dict[str, Any]] = []
    eb_intent_rows: list[dict[str, Any]] = []
    missing_rows: list[dict[str, Any]] = []
    base_by_eb = {(row["mous_uid"], row["eb_uid"]): row for row in eb_rows}

    for (mous_uid, eb_uid), science_rows in sorted(eb_groups.items()):
        base_row = base_by_eb[(mous_uid, eb_uid)]
        full_path = matched_full_listobs_path(science_rows[0]["listobs_path"], search_roots)
        if not full_path.exists():
            missing_rows.append(
                {"mous_uid": mous_uid, "eb_uid": eb_uid, "expected_full_listobs_path": str(full_path)}
            )
            continue
        cache_key = str(full_path)
        if cache_key not in full_cache:
            full_cache[cache_key] = parse_full_listobs(full_path)
        full_eb = full_cache[cache_key]

        current_volume_by_bucket = defaultdict(float)
        current_time_by_bucket = defaultdict(float)
        scan_count_by_bucket = defaultdict(int)
        cap_ratio_weight_sums = defaultdict(float)
        cap_ratio_weighted = defaultdict(float)

        for scan in full_eb.scans:
            metrics = scan_metrics(scan, full_eb)
            current_volume = float(metrics["volume_gb"])
            current_time = float(scan.duration_s)
            if not math.isfinite(current_volume):
                current_volume = 0.0
            if not math.isfinite(current_time):
                current_time = 0.0
            bucket = scan.bucket
            current_volume_by_bucket[bucket] += current_volume
            current_time_by_bucket[bucket] += current_time
            scan_count_by_bucket[bucket] += 1
            detail_rows.append(
                {
                    "mous_uid": mous_uid,
                    "eb_uid": eb_uid,
                    "array_type": full_eb.array_type,
                    "full_listobs_path": str(full_path),
                    "scan_id": scan.scan_id,
                    "bucket": bucket,
                    "duration_s": current_time,
                    "current_volume_gb": current_volume,
                    "retained_spw_count": metrics["retained_spw_count"],
                    "median_retained_velocity_resolution_kms": metrics["median_retained_velocity_resolution_kms"],
                    "raw_intent_signature": scan.raw_intent_signature,
                    "retained_signature": metrics["retained_signature"],
                }
            )
            if (
                not apply_calibrator_cap
                or current_volume <= 0.0
                or bucket not in CALIBRATOR_CAP_BUCKETS
            ):
                continue
            for milestone in MILESTONES:
                for method in METHODS:
                    ratio = cap_ratio_for_scan(
                        method=method,
                        milestone=milestone,
                        scan=scan,
                        full_eb=full_eb,
                        velocity_cap_kms=calibrator_cap_kms,
                    )
                    cap_ratio_weight_sums[(milestone, method, bucket)] += current_volume
                    cap_ratio_weighted[(milestone, method, bucket)] += current_volume * ratio

        record: dict[str, Any] = {
            "mous_uid": mous_uid,
            "eb_uid": eb_uid,
            "array_type": full_eb.array_type,
            "full_listobs_path": str(full_path),
        }
        total_current_all = sum(current_volume_by_bucket.values())
        selected_current = sum(current_volume_by_bucket[bucket] for bucket in INTENT_PROJECTION_BUCKETS)
        record["current_total_volume_all_buckets_gb"] = total_current_all
        record["current_selected_volume_gb"] = selected_current
        record["current_excluded_volume_gb"] = total_current_all - selected_current
        for bucket in INTENT_PROJECTION_BUCKETS + ("atmosphere", "other", "pointing", "polarization", "diffgain_reference", "diffgain_on_source"):
            record[f"current_{bucket}_volume_gb"] = current_volume_by_bucket[bucket]
            record[f"current_{bucket}_time_s"] = current_time_by_bucket[bucket]
            record[f"current_{bucket}_scan_count"] = scan_count_by_bucket[bucket]

        for milestone in MILESTONES:
            for method in METHODS:
                prefix = f"{milestone.lower()}_{method}"
                science_factor = float(base_row[f"{prefix}_factor"])
                uncapped_total = 0.0
                capped_total = 0.0
                for bucket in INTENT_PROJECTION_BUCKETS:
                    current_bucket_volume = current_volume_by_bucket[bucket]
                    projected_uncapped = current_bucket_volume * science_factor if math.isfinite(science_factor) else math.nan
                    ratio = 1.0
                    if bucket in CALIBRATOR_CAP_BUCKETS:
                        weight_sum = cap_ratio_weight_sums[(milestone, method, bucket)]
                        if weight_sum > 0.0:
                            ratio = cap_ratio_weighted[(milestone, method, bucket)] / weight_sum
                    projected_capped = projected_uncapped
                    if apply_calibrator_cap and bucket in CALIBRATOR_CAP_BUCKETS and math.isfinite(projected_uncapped):
                        projected_capped = projected_uncapped * ratio
                    record[f"{prefix}_{bucket}_projected_uncapped_volume_gb"] = projected_uncapped
                    record[f"{prefix}_{bucket}_projected_volume_gb"] = projected_capped
                    record[f"{prefix}_{bucket}_cap_ratio"] = ratio
                    if math.isfinite(projected_uncapped):
                        uncapped_total += projected_uncapped
                    else:
                        uncapped_total = math.nan
                    if math.isfinite(projected_capped):
                        capped_total += projected_capped
                    else:
                        capped_total = math.nan
                record[f"{prefix}_selected_projected_uncapped_volume_gb"] = uncapped_total
                record[f"{prefix}_selected_projected_volume_gb"] = capped_total
                record[f"{prefix}_selected_factor_uncapped"] = (
                    uncapped_total / selected_current
                    if selected_current > 0.0 and math.isfinite(uncapped_total)
                    else math.nan
                )
                record[f"{prefix}_selected_factor"] = (
                    capped_total / selected_current
                    if selected_current > 0.0 and math.isfinite(capped_total)
                    else math.nan
                )
                record[f"{prefix}_cap_delta_volume_gb"] = capped_total - uncapped_total
        eb_intent_rows.append(record)

    return eb_intent_rows, detail_rows, missing_rows


def build_intent_aware_summary_rows(eb_intent_rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    overall_rows: list[dict[str, Any]] = []
    by_array_rows: list[dict[str, Any]] = []
    for milestone in MILESTONES:
        for method in METHODS:
            prefix = f"{milestone.lower()}_{method}"
            factors = [float(row[f"{prefix}_selected_factor"]) for row in eb_intent_rows]
            factors_uncapped = [float(row[f"{prefix}_selected_factor_uncapped"]) for row in eb_intent_rows]
            projected = [float(row[f"{prefix}_selected_projected_volume_gb"]) for row in eb_intent_rows]
            projected_uncapped = [float(row[f"{prefix}_selected_projected_uncapped_volume_gb"]) for row in eb_intent_rows]
            current = [float(row["current_selected_volume_gb"]) for row in eb_intent_rows]
            total_current = sum(value for value in current if math.isfinite(value))
            total_projected = sum(value for value in projected if math.isfinite(value))
            total_projected_uncapped = sum(value for value in projected_uncapped if math.isfinite(value))
            overall_rows.append(
                {
                    "milestone": milestone,
                    "method": method,
                    "n_eb": len(eb_intent_rows),
                    "current_selected_total_volume_gb": total_current,
                    "projected_selected_total_uncapped_volume_gb": total_projected_uncapped,
                    "projected_selected_total_volume_gb": total_projected,
                    "sample_total_factor_uncapped": total_projected_uncapped / total_current if total_current > 0.0 else math.nan,
                    "sample_total_factor": total_projected / total_current if total_current > 0.0 else math.nan,
                    "mean_eb_factor_uncapped": summarize(factors_uncapped)["mean"],
                    "median_eb_factor_uncapped": summarize(factors_uncapped)["median"],
                    "mean_eb_factor": summarize(factors)["mean"],
                    "median_eb_factor": summarize(factors)["median"],
                    "cap_delta_total_volume_gb": total_projected - total_projected_uncapped,
                }
            )
            for array_type in ("12m", "7m"):
                sub = [row for row in eb_intent_rows if row["array_type"] == array_type]
                if not sub:
                    continue
                sub_factors = [float(row[f"{prefix}_selected_factor"]) for row in sub]
                sub_factors_uncapped = [float(row[f"{prefix}_selected_factor_uncapped"]) for row in sub]
                sub_projected = [float(row[f"{prefix}_selected_projected_volume_gb"]) for row in sub]
                sub_projected_uncapped = [float(row[f"{prefix}_selected_projected_uncapped_volume_gb"]) for row in sub]
                sub_current = [float(row["current_selected_volume_gb"]) for row in sub]
                sub_total_current = sum(value for value in sub_current if math.isfinite(value))
                sub_total_projected = sum(value for value in sub_projected if math.isfinite(value))
                sub_total_projected_uncapped = sum(value for value in sub_projected_uncapped if math.isfinite(value))
                by_array_rows.append(
                    {
                        "milestone": milestone,
                        "method": method,
                        "array_type": array_type,
                        "n_eb": len(sub),
                        "sample_total_factor_uncapped": (
                            sub_total_projected_uncapped / sub_total_current if sub_total_current > 0.0 else math.nan
                        ),
                        "sample_total_factor": (
                            sub_total_projected / sub_total_current if sub_total_current > 0.0 else math.nan
                        ),
                        "mean_eb_factor_uncapped": summarize(sub_factors_uncapped)["mean"],
                        "median_eb_factor_uncapped": summarize(sub_factors_uncapped)["median"],
                        "mean_eb_factor": summarize(sub_factors)["mean"],
                        "median_eb_factor": summarize(sub_factors)["median"],
                        "cap_delta_total_volume_gb": sub_total_projected - sub_total_projected_uncapped,
                    }
                )
    return overall_rows, by_array_rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tabulation_dir", type=Path, help="Directory containing spw_level_table.csv from spw_tabulate.py")
    parser.add_argument("output_dir", type=Path)
    parser.add_argument(
        "--intent-aware",
        action="store_true",
        help="Add an intent-aware EB-level projection for science, bandpass, phase, and check-source scans",
    )
    parser.add_argument(
        "--apply-calibrator-cap",
        action="store_true",
        help="Apply a spectral-resolution cap to the intent-aware bandpass/phase/check-source buckets",
    )
    parser.add_argument(
        "--calibrator-cap-kms",
        type=float,
        default=1.0,
        help="Velocity-resolution cap in km/s for calibrator scans when --apply-calibrator-cap is set",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.calibrator_cap_kms <= 0.0:
        raise SystemExit("--calibrator-cap-kms must be > 0")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    intent_aware_requested = args.intent_aware or args.apply_calibrator_cap

    spw_rows_all = read_csv_rows(args.tabulation_dir / "spw_level_table.csv")
    eb_unique: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in spw_rows_all:
        key = (row["mous_uid"], row["eb_uid"], row["spw_id"])
        eb_unique.setdefault(key, row)
    eb_spw_rows = sorted(eb_unique.values(), key=lambda row: (row["mous_uid"], row["eb_uid"], int(row["spw_id"])))
    eb_groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in eb_spw_rows:
        eb_groups[(row["mous_uid"], row["eb_uid"])].append(row)

    eb_tuning_summary_all = read_csv_rows(args.tabulation_dir / "eb_tuning_summary.csv")
    eb_tuning_spw_all = read_csv_rows(args.tabulation_dir / "eb_tuning_spw_table.csv")
    mous_tuning_template_rows_all = read_csv_rows(args.tabulation_dir / "mous_tuning_spw_templates.csv")
    eb_tuning_summary_groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    eb_tuning_spw_groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in eb_tuning_summary_all:
        eb_tuning_summary_groups[(row["mous_uid"], row["eb_uid"])].append(row)
    for row in eb_tuning_spw_all:
        eb_tuning_spw_groups[(row["mous_uid"], row["eb_uid"])].append(row)

    search_roots = [Path.cwd(), args.tabulation_dir.parent, Path(__file__).resolve().parent]
    runtime_cache: dict[str, EbRuntime] = {}
    eb_rows: list[dict[str, Any]] = []

    for (mous_uid, eb_uid), rows in sorted(eb_groups.items()):
        tuning_summary_rows = eb_tuning_summary_groups.get((mous_uid, eb_uid), [])
        tuning_spw_rows = eb_tuning_spw_groups.get((mous_uid, eb_uid), [])
        if not tuning_summary_rows or not tuning_spw_rows:
            continue

        resolved_path = resolve_local_path(tuning_summary_rows[0]["listobs_path"], search_roots)
        cache_key = str(resolved_path)
        if cache_key not in runtime_cache:
            runtime_cache[cache_key] = parse_listobs_runtime(resolved_path)
        runtime = runtime_cache[cache_key]
        tuning_groups = build_eb_tuning_groups(tuning_spw_rows)
        if not tuning_groups:
            continue
        all_templates = [template for group in tuning_groups for template in group.spw_templates]
        corr_counts = [
            count_corr_products(row["corr_products"])
            for row in sorted(
                tuning_spw_rows,
                key=lambda row: (int(row["tuning_group_index"]), int(row["spw_id"])),
            )
        ]
        band = all_templates[0].band
        current_rate, current_volume, current_science_time_s = current_eb_metrics_from_tuning_rows(
            tuning_summary_rows,
            tuning_spw_rows,
            runtime,
        )
        current_nchan_agg = current_eb_nchan_agg_from_tuning_rows(tuning_summary_rows)
        current_tuning_groups = int(tuning_summary_rows[0]["n_tuning_groups_in_eb"])
        current_total_science_bw_hz = float(tuning_summary_rows[0]["eb_total_science_bw_hz"])
        is_spectral_scan = str(tuning_summary_rows[0]["is_spectral_scan"]).lower() == "true"

        record: dict[str, Any] = {
            "mous_uid": mous_uid,
            "eb_uid": eb_uid,
            "array_type": runtime.array_type,
            "current_science_time_s": current_science_time_s,
            "current_interval_median_s": runtime.science_interval_median_s,
            "current_rate_gbps": current_rate,
            "current_volume_gb": current_volume,
            "current_nchan_agg": current_nchan_agg,
            "current_tuning_groups": current_tuning_groups,
            "current_total_science_bw_hz": current_total_science_bw_hz,
            "is_spectral_scan": is_spectral_scan,
            "listobs_path": runtime.listobs_path,
            "band": band,
        }

        for milestone in MILESTONES:
            eq_nspw = projected_spw_equivalents(milestone, band)
            milestone_bw_hz = DEFAULT_CONFIG.milestones[milestone].bandwidth_for(band) * 1.0e9
            projected_tuning_groups = max(1, math.ceil(current_total_science_bw_hz / milestone_bw_hz))
            projected_nspw_total = projected_tuning_groups * eq_nspw
            tint_s = projected_tint_s(runtime.array_type, runtime.science_interval_median_s, milestone)
            nant = projected_nant(runtime.array_type, milestone)
            for method in METHODS:
                projected_nchan_agg = project_nchan_agg_for_tuning_groups(
                    method,
                    milestone,
                    tuning_groups,
                    projected_nspw_total=projected_nspw_total,
                )
                total_channel_pol = projected_channel_pol_agg_for_templates(
                    method,
                    milestone,
                    all_templates,
                    corr_counts,
                    projected_nspw_total,
                )
                mean_channel_pol_per_tuning = total_channel_pol / projected_tuning_groups
                projected_rate = datarate_gbps(nant, [(mean_channel_pol_per_tuning, tint_s)], 1.0)
                projected_volume = projected_rate * current_science_time_s if math.isfinite(projected_rate) else math.nan
                factor = projected_rate / current_rate if math.isfinite(projected_rate) and math.isfinite(current_rate) and current_rate > 0 else math.nan
                nchan_factor = (
                    projected_nchan_agg / current_nchan_agg
                    if math.isfinite(projected_nchan_agg) and math.isfinite(current_nchan_agg) and current_nchan_agg > 0
                    else math.nan
                )
                prefix = f"{milestone.lower()}_{method}"
                record[f"{prefix}_eq_nspw_per_tuning"] = eq_nspw
                record[f"{prefix}_projected_tuning_groups"] = projected_tuning_groups
                record[f"{prefix}_projected_total_nspw"] = projected_nspw_total
                record[f"{prefix}_nchan_agg"] = projected_nchan_agg
                record[f"{prefix}_nchan_factor"] = nchan_factor
                record[f"{prefix}_rate_gbps"] = projected_rate
                record[f"{prefix}_volume_gb"] = projected_volume
                record[f"{prefix}_factor"] = factor
        eb_rows.append(record)

    eb_fieldnames = list(eb_rows[0].keys()) if eb_rows else []
    write_csv(args.output_dir / "eb_projection_factors.csv", eb_rows, eb_fieldnames)

    sidecar_rows = []
    for row in mous_tuning_template_rows_all:
        out = dict(row)
        out["band"] = alma_band_from_freq_hz(float(row["center_freq_hz"]))
        sidecar_rows.append(out)
    write_csv(
        args.output_dir / "mous_spw_templates.csv",
        sidecar_rows,
        [
            "mous_uid",
            "representative_eb_uid",
            "tuning_group_id",
            "tuning_group_index",
            "tuning_total_bw_hz",
            "spw_id",
            "center_freq_hz",
            "velocity_resolution_kms",
            "band",
            "n_tuning_groups_in_mous",
            "is_spectral_scan_mous",
            "mous_total_science_bw_hz",
            "mous_current_nchan_agg",
            "eb_tuning_consistency_flag",
            "eb_tuning_consistency_notes",
        ],
    )
    write_csv(
        args.output_dir / "mous_tuning_spw_templates.csv",
        sidecar_rows,
        [
            "mous_uid",
            "representative_eb_uid",
            "tuning_group_id",
            "tuning_group_index",
            "tuning_total_bw_hz",
            "spw_id",
            "center_freq_hz",
            "velocity_resolution_kms",
            "band",
            "n_tuning_groups_in_mous",
            "is_spectral_scan_mous",
            "mous_total_science_bw_hz",
            "mous_current_nchan_agg",
            "eb_tuning_consistency_flag",
            "eb_tuning_consistency_notes",
        ],
    )

    summary_rows: list[dict[str, Any]] = []
    by_array_rows: list[dict[str, Any]] = []
    nchan_summary_rows: list[dict[str, Any]] = []
    nchan_by_array_rows: list[dict[str, Any]] = []
    for milestone in MILESTONES:
        for method in METHODS:
            prefix = f"{milestone.lower()}_{method}"
            factors = [float(row[f"{prefix}_factor"]) for row in eb_rows]
            projected_volumes = [float(row[f"{prefix}_volume_gb"]) for row in eb_rows]
            current_volumes = [float(row["current_volume_gb"]) for row in eb_rows]
            nchan_factors = [float(row[f"{prefix}_nchan_factor"]) for row in eb_rows]
            projected_nchan = [float(row[f"{prefix}_nchan_agg"]) for row in eb_rows]
            current_nchan = [float(row["current_nchan_agg"]) for row in eb_rows]
            total_current = sum(value for value in current_volumes if math.isfinite(value))
            total_projected = sum(value for value in projected_volumes if math.isfinite(value))
            total_factor = total_projected / total_current if total_current > 0 else math.nan
            total_current_nchan = sum(value for value in current_nchan if math.isfinite(value))
            total_projected_nchan = sum(value for value in projected_nchan if math.isfinite(value))
            total_nchan_factor = total_projected_nchan / total_current_nchan if total_current_nchan > 0 else math.nan
            stats = summarize(factors)
            nchan_stats = summarize(nchan_factors)
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
            nchan_summary_rows.append(
                {
                    "milestone": milestone,
                    "method": method,
                    "n_eb": len(eb_rows),
                    "sample_total_factor": total_nchan_factor,
                    "mean_eb_factor": nchan_stats["mean"],
                    "median_eb_factor": nchan_stats["median"],
                    "current_total_nchan_agg": total_current_nchan,
                    "projected_total_nchan_agg": total_projected_nchan,
                }
            )
            for array_type in ("12m", "7m"):
                sub = [row for row in eb_rows if row["array_type"] == array_type]
                if not sub:
                    continue
                sub_factors = [float(row[f"{prefix}_factor"]) for row in sub]
                sub_proj = [float(row[f"{prefix}_volume_gb"]) for row in sub]
                sub_cur = [float(row["current_volume_gb"]) for row in sub]
                sub_nchan_factors = [float(row[f"{prefix}_nchan_factor"]) for row in sub]
                sub_proj_nchan = [float(row[f"{prefix}_nchan_agg"]) for row in sub]
                sub_cur_nchan = [float(row["current_nchan_agg"]) for row in sub]
                sub_current_total = sum(value for value in sub_cur if math.isfinite(value))
                sub_projected_total = sum(value for value in sub_proj if math.isfinite(value))
                sub_total_factor = sub_projected_total / sub_current_total if sub_current_total > 0 else math.nan
                sub_current_nchan_total = sum(value for value in sub_cur_nchan if math.isfinite(value))
                sub_projected_nchan_total = sum(value for value in sub_proj_nchan if math.isfinite(value))
                sub_total_nchan_factor = (
                    sub_projected_nchan_total / sub_current_nchan_total
                    if sub_current_nchan_total > 0
                    else math.nan
                )
                sub_stats = summarize(sub_factors)
                sub_nchan_stats = summarize(sub_nchan_factors)
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
                nchan_by_array_rows.append(
                    {
                        "milestone": milestone,
                        "method": method,
                        "array_type": array_type,
                        "n_eb": len(sub),
                        "sample_total_factor": sub_total_nchan_factor,
                        "mean_eb_factor": sub_nchan_stats["mean"],
                        "median_eb_factor": sub_nchan_stats["median"],
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
    write_csv(
        args.output_dir / "summary_overall_nchan.csv",
        nchan_summary_rows,
        [
            "milestone",
            "method",
            "n_eb",
            "sample_total_factor",
            "mean_eb_factor",
            "median_eb_factor",
            "current_total_nchan_agg",
            "projected_total_nchan_agg",
        ],
    )
    write_csv(
        args.output_dir / "summary_by_array_nchan.csv",
        nchan_by_array_rows,
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
    write_latex_factor_table(
        args.output_dir / "summary_overall_nchan.tex",
        nchan_summary_rows,
        "Science aggregate channel count, relative to current, for different WSU projection scenarios.",
        "tab:sidecarproj_nchan",
    )

    intent_summary_lines: list[str] = []
    intent_overall_rows: list[dict[str, Any]] = []
    if intent_aware_requested:
        eb_intent_rows, intent_detail_rows, intent_missing_rows = build_intent_aware_rows(
            eb_rows=eb_rows,
            eb_groups=eb_groups,
            search_roots=search_roots,
            calibrator_cap_kms=args.calibrator_cap_kms,
            apply_calibrator_cap=args.apply_calibrator_cap,
        )
        write_csv(
            args.output_dir / "eb_intent_projection_factors.csv",
            eb_intent_rows,
            list(eb_intent_rows[0].keys()) if eb_intent_rows else ["mous_uid", "eb_uid"],
        )
        write_csv(
            args.output_dir / "intent_aware_scan_details.csv",
            intent_detail_rows,
            list(intent_detail_rows[0].keys()) if intent_detail_rows else [
                "mous_uid",
                "eb_uid",
                "array_type",
                "full_listobs_path",
                "scan_id",
                "bucket",
                "duration_s",
                "current_volume_gb",
                "retained_spw_count",
                "median_retained_velocity_resolution_kms",
                "raw_intent_signature",
                "retained_signature",
            ],
        )
        write_csv(
            args.output_dir / "intent_aware_missing_full_listobs.csv",
            intent_missing_rows,
            ["mous_uid", "eb_uid", "expected_full_listobs_path"],
        )
        intent_overall_rows, intent_by_array_rows = build_intent_aware_summary_rows(eb_intent_rows)
        write_csv(
            args.output_dir / "summary_overall_intent_aware.csv",
            intent_overall_rows,
            [
                "milestone",
                "method",
                "n_eb",
                "current_selected_total_volume_gb",
                "projected_selected_total_uncapped_volume_gb",
                "projected_selected_total_volume_gb",
                "sample_total_factor_uncapped",
                "sample_total_factor",
                "mean_eb_factor_uncapped",
                "median_eb_factor_uncapped",
                "mean_eb_factor",
                "median_eb_factor",
                "cap_delta_total_volume_gb",
            ],
        )
        write_csv(
            args.output_dir / "summary_by_array_intent_aware.csv",
            intent_by_array_rows,
            [
                "milestone",
                "method",
                "array_type",
                "n_eb",
                "sample_total_factor_uncapped",
                "sample_total_factor",
                "mean_eb_factor_uncapped",
                "median_eb_factor_uncapped",
                "mean_eb_factor",
                "median_eb_factor",
                "cap_delta_total_volume_gb",
            ],
        )
        intent_rows_by_key = {(row["milestone"], row["method"]): row for row in intent_overall_rows}
        intent_summary_lines.extend(
            [
                "## Intent-Aware Projection",
                "",
                "This optional analysis projects science, bandpass, phase, and check-source scans separately at the EB level.",
                "It anchors the uncapped calibrator buckets to the same science growth factor as the legacy projection, then applies any requested calibrator-only spectral-resolution cap as an EB-level adjustment.",
                "ATMOSPHERE and other buckets are tabulated in the current-data summary but excluded from this projection path for now.",
                "",
                f"- Calibrator cap enabled: `{args.apply_calibrator_cap}`",
                f"- Calibrator cap value: `{args.calibrator_cap_kms:.3f}` km/s",
                "",
            ]
        )
        for milestone in MILESTONES:
            intent_summary_lines.append(f"### {milestone}")
            intent_summary_lines.append("")
            for method in METHODS:
                row = intent_rows_by_key[(milestone, method)]
                intent_summary_lines.append(
                    "- "
                    + f"{method}: uncapped sample-total factor {format_factor(float(row['sample_total_factor_uncapped']))}, "
                    + f"capped sample-total factor {format_factor(float(row['sample_total_factor']))}, "
                    + f"cap delta {float(row['cap_delta_total_volume_gb']):.2f} GB"
                )
            intent_summary_lines.append("")

    rows_by_key = {(row["milestone"], row["method"]): row for row in summary_rows}
    nchan_rows_by_key = {(row["milestone"], row["method"]): row for row in nchan_summary_rows}
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

    summary_lines.extend(
        [
            "## Aggregate Nchan Results",
            "",
            "These factors are based on tuning-aware aggregate science channel counts over the whole EB, not raw visibility volume.",
            "",
        ]
    )
    for milestone in MILESTONES:
        summary_lines.append(f"### {milestone} Aggregate Nchan")
        summary_lines.append("")
        for method in METHODS:
            row = nchan_rows_by_key[(milestone, method)]
            summary_lines.append(
                "- "
                + f"{method}: sample-total factor {format_factor(float(row['sample_total_factor']))}, "
                + f"mean EB factor {format_factor(float(row['mean_eb_factor']))}, "
                + f"median EB factor {format_factor(float(row['median_eb_factor']))}"
            )
        summary_lines.append("")

    if intent_summary_lines:
        summary_lines.extend(intent_summary_lines)

    (args.output_dir / "summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    (args.output_dir / "run_metadata.json").write_text(
        json.dumps(
            {
                "tabulation_dir": str(args.tabulation_dir),
                "n_eb": len(eb_rows),
                "methods": list(METHODS),
                "milestones": list(MILESTONES),
                "mous_spw_sidecar": str(args.output_dir / "mous_tuning_spw_templates.csv"),
                "intent_aware": intent_aware_requested,
                "apply_calibrator_cap": args.apply_calibrator_cap,
                "calibrator_cap_kms": args.calibrator_cap_kms,
                "intent_projection_buckets": list(INTENT_PROJECTION_BUCKETS),
                "calibrator_cap_buckets": list(CALIBRATOR_CAP_BUCKETS),
                "nchan_summary": str(args.output_dir / "summary_overall_nchan.csv"),
                "nchan_summary_latex": str(args.output_dir / "summary_overall_nchan.tex"),
                "intent_aware_summary": (
                    str(args.output_dir / "summary_overall_intent_aware.csv") if intent_aware_requested else None
                ),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
