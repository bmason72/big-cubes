#!/usr/bin/env python3
"""Summarize ALMA spectral setup signatures from spw_tabulate outputs.

This script consumes the CSV outputs from ``spw_tabulate.py`` and builds a
setup-level summary tuned for visibility data-volume proxy analysis:

- derive bandwidth and spectral-resolution bins from unique EB/SPW rows,
- assign each MOUS a sorted setup signature based on MOUS-combined SPWs,
- assign each EB to the MOUS-combined setup while preserving EB-level SPW and
  data-rate information,
- summarize setup prevalence by MOUS, EB, and EB-SPW counts,
- compute line-bandwidth fraction statistics per setup,
- compute mean EB data rate per setup,
- generate checkable histogram plots with the derived cut positions.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
ARRAY_7M_MAX_ANTENNAS = 20
FIDUCIAL_NANT_12M = 43
FIDUCIAL_NANT_7M = 10
BYTES_PER_CORRELATION = 4.0
DEFAULT_BANDWIDTH_LOW_CUT_MHZ = 90.0
DEFAULT_BANDWIDTH_HIGH_CUT_MHZ = 1200.0


@dataclass
class EbRuntimeMetadata:
    listobs_path: str
    actual_antenna_count: int | None
    fiducial_nant: int | None
    array_type: str
    science_interval_median_s: float
    science_interval_by_spw_s: dict[int, float]


def percentile(values: list[float], frac: float) -> float:
    if not values:
        return math.nan
    if len(values) == 1:
        return values[0]
    idx = (len(values) - 1) * frac
    lo = values[int(math.floor(idx))]
    hi = values[int(math.ceil(idx))]
    return lo + (hi - lo) * (idx - math.floor(idx))


def midpoint(a: float, b: float) -> float:
    return (a + b) / 2.0


def choose_gap_cut(values: list[float], frac: float) -> dict[str, float]:
    """Choose a quantile cut, shifting off spikes so equal values are unsplit."""
    sorted_vals = sorted(values)
    raw_quantile = percentile(sorted_vals, frac)
    unique = sorted(set(sorted_vals))

    if not unique:
        return {
            "raw_quantile": math.nan,
            "cut": math.nan,
            "lower_value": math.nan,
            "upper_value": math.nan,
        }
    if len(unique) == 1:
        only = unique[0]
        return {
            "raw_quantile": only,
            "cut": only,
            "lower_value": only,
            "upper_value": only,
        }

    pos = 0
    while pos < len(unique) and unique[pos] < raw_quantile:
        pos += 1

    if pos == 0:
        lower, upper = unique[0], unique[1]
    elif pos == len(unique):
        lower, upper = unique[-2], unique[-1]
    elif unique[pos] == raw_quantile:
        lower = unique[pos]
        upper = unique[pos + 1] if pos + 1 < len(unique) else unique[pos]
        if upper == lower and pos > 0:
            lower = unique[pos - 1]
    else:
        lower, upper = unique[pos - 1], unique[pos]

    return {
        "raw_quantile": raw_quantile,
        "cut": midpoint(lower, upper),
        "lower_value": lower,
        "upper_value": upper,
    }


def classify(value: float, low_cut: float, high_cut: float, labels: tuple[str, str, str]) -> str:
    if value <= low_cut:
        return labels[0]
    if value <= high_cut:
        return labels[1]
    return labels[2]


def finite_values(values: list[float]) -> list[float]:
    return [value for value in values if math.isfinite(value)]


def mean(values: list[float]) -> float:
    finite = finite_values(values)
    return sum(finite) / len(finite) if finite else math.nan


def population_stddev(values: list[float]) -> float:
    finite = finite_values(values)
    if not finite:
        return math.nan
    mu = mean(finite)
    return math.sqrt(sum((value - mu) ** 2 for value in finite) / len(finite))


def finite_median(values: list[float]) -> float:
    finite = finite_values(values)
    return float(median(finite)) if finite else math.nan


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


def dedupe_eb_spw_rows(spw_level_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    unique: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in spw_level_rows:
        key = (row["mous_uid"], row["eb_uid"], row["spw_id"])
        unique.setdefault(key, row)
    return sorted(unique.values(), key=lambda row: (row["mous_uid"], row["eb_uid"], int(row["spw_id"])))


def derive_cut_metadata(
    eb_spw_rows: list[dict[str, str]],
    *,
    bandwidth_cut_mode: str,
    bandwidth_low_cut_mhz: float,
    bandwidth_high_cut_mhz: float,
) -> dict[str, Any]:
    bandwidths_ghz = [float(row["total_bw_hz"]) / 1.0e9 for row in eb_spw_rows]
    resolutions_kms = [float(row["velocity_resolution_kms"]) for row in eb_spw_rows]

    if bandwidth_cut_mode == "fixed":
        bw_low = {
            "raw_quantile": math.nan,
            "cut": bandwidth_low_cut_mhz / 1.0e3,
            "lower_value": math.nan,
            "upper_value": math.nan,
        }
        bw_high = {
            "raw_quantile": math.nan,
            "cut": bandwidth_high_cut_mhz / 1.0e3,
            "lower_value": math.nan,
            "upper_value": math.nan,
        }
    elif bandwidth_cut_mode == "percentile":
        bw_low = choose_gap_cut(bandwidths_ghz, 1.0 / 3.0)
        bw_high = choose_gap_cut(bandwidths_ghz, 2.0 / 3.0)
    else:
        raise ValueError(f"unsupported bandwidth_cut_mode={bandwidth_cut_mode}")

    res_low = choose_gap_cut(resolutions_kms, 1.0 / 3.0)
    res_high = choose_gap_cut(resolutions_kms, 2.0 / 3.0)

    return {
        "bandwidth": {
            "mode": bandwidth_cut_mode,
            "q33": bw_low["raw_quantile"],
            "q66": bw_high["raw_quantile"],
            "low_cut": bw_low["cut"],
            "high_cut": bw_high["cut"],
            "low_gap": [bw_low["lower_value"], bw_low["upper_value"]],
            "high_gap": [bw_high["lower_value"], bw_high["upper_value"]],
        },
        "resolution": {
            "q33": res_low["raw_quantile"],
            "q66": res_high["raw_quantile"],
            "low_cut": res_low["cut"],
            "high_cut": res_high["cut"],
            "low_gap": [res_low["lower_value"], res_low["upper_value"]],
            "high_gap": [res_high["lower_value"], res_high["upper_value"]],
        },
        "n_unique_eb_spw": len(eb_spw_rows),
    }


def signature_from_rows(
    rows: list[dict[str, str]],
    bandwidth_cuts: dict[str, float],
    resolution_cuts: dict[str, float],
) -> tuple[str, list[str]]:
    tokens: list[str] = []
    for row in sorted(rows, key=lambda item: (float(item["center_freq_hz"]), int(item["spw_id"]))):
        bw_ghz = float(row["total_bw_hz"]) / 1.0e9
        res_kms = float(row["velocity_resolution_kms"])
        bw_code = classify(bw_ghz, bandwidth_cuts["low_cut"], bandwidth_cuts["high_cut"], ("N", "M", "W"))
        res_code = classify(res_kms, resolution_cuts["low_cut"], resolution_cuts["high_cut"], ("F", "M", "C"))
        tokens.append(f"{bw_code}{res_code}")
    sorted_tokens = sorted(tokens)
    return "[" + ",".join(sorted_tokens) + "]", sorted_tokens


def build_mous_assignments(
    mous_spw_rows: list[dict[str, str]],
    bandwidth_cuts: dict[str, float],
    resolution_cuts: dict[str, float],
) -> list[dict[str, Any]]:
    mous_groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in mous_spw_rows:
        mous_groups[row["mous_uid"]].append(row)

    assignments: list[dict[str, Any]] = []
    for mous_uid, group in sorted(mous_groups.items()):
        signature, tokens = signature_from_rows(group, bandwidth_cuts, resolution_cuts)
        n_ebs = max(int(row["n_ebs_parsed"]) for row in group)
        total_spw_bw_hz = 0.0
        total_line_bw_hz = 0.0
        for row in group:
            total_bw_hz = float(row["total_bw_hz"])
            total_spw_bw_hz += total_bw_hz
            total_line_bw_hz += total_bw_hz * float(row["mean_summed_line_fraction_of_spw"])
        line_fraction_percent = 100.0 * total_line_bw_hz / total_spw_bw_hz if total_spw_bw_hz else math.nan
        assignments.append(
            {
                "mous_uid": mous_uid,
                "n_ebs": n_ebs,
                "n_spws": len(tokens),
                "setup_signature": signature,
                "setup_tokens": " ".join(tokens),
                "line_fraction_percent": line_fraction_percent,
            }
        )
    return assignments


def parse_eb_runtime_metadata(listobs_path: Path) -> EbRuntimeMetadata:
    actual_antenna_count: int | None = None
    science_intervals: list[float] = []
    intervals_by_spw: dict[int, list[float]] = defaultdict(list)
    in_scans = False

    for raw in listobs_path.read_text(encoding="utf-8", errors="ignore").splitlines():
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
        spw_ids = parse_int_list(match.group("spw_ids"))
        intervals = parse_float_list(match.group("intervals"))
        if len(intervals) == 1 and len(spw_ids) > 1:
            intervals = intervals * len(spw_ids)
        for spw_id, interval in zip(spw_ids, intervals):
            intervals_by_spw[spw_id].append(interval)
            science_intervals.append(interval)

    if actual_antenna_count is None:
        array_type = "unknown"
        fiducial_nant = None
    elif actual_antenna_count <= ARRAY_7M_MAX_ANTENNAS:
        array_type = "7m"
        fiducial_nant = FIDUCIAL_NANT_7M
    else:
        array_type = "12m"
        fiducial_nant = FIDUCIAL_NANT_12M

    return EbRuntimeMetadata(
        listobs_path=str(listobs_path),
        actual_antenna_count=actual_antenna_count,
        fiducial_nant=fiducial_nant,
        array_type=array_type,
        science_interval_median_s=finite_median(science_intervals),
        science_interval_by_spw_s={spw_id: finite_median(vals) for spw_id, vals in intervals_by_spw.items()},
    )


def count_correlation_products(corr_products: str) -> int:
    tokens = CORR_TOKEN_RE.findall(corr_products)
    if tokens:
        return len(tokens)
    parts = [part for part in corr_products.split() if part]
    return len(parts)


def compute_eb_data_rate_gbps(rows: list[dict[str, str]], metadata: EbRuntimeMetadata) -> float:
    if metadata.fiducial_nant is None:
        return math.nan
    nbase = metadata.fiducial_nant * (metadata.fiducial_nant - 1) / 2.0
    sum_bytes_per_second = 0.0
    for row in rows:
        spw_id = int(row["spw_id"])
        tint_s = metadata.science_interval_by_spw_s.get(spw_id, metadata.science_interval_median_s)
        if not math.isfinite(tint_s) or tint_s <= 0:
            return math.nan
        nchan = int(row["nchan"])
        ncorr = count_correlation_products(row.get("corr_products", ""))
        bytes_per_visibility = BYTES_PER_CORRELATION * ncorr
        sum_bytes_per_second += nchan * bytes_per_visibility / tint_s
    return nbase * sum_bytes_per_second / 1.0e9


def build_eb_assignments(
    eb_spw_rows: list[dict[str, str]],
    mous_assignments: list[dict[str, Any]],
    tabulation_dir: Path,
) -> list[dict[str, Any]]:
    mous_signature_map = {row["mous_uid"]: row["setup_signature"] for row in mous_assignments}
    eb_groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in eb_spw_rows:
        eb_groups[(row["mous_uid"], row["eb_uid"])].append(row)

    search_roots = [Path.cwd(), tabulation_dir.parent, Path(__file__).resolve().parent]
    runtime_cache: dict[str, EbRuntimeMetadata] = {}
    assignments: list[dict[str, Any]] = []

    for (mous_uid, eb_uid), group in sorted(eb_groups.items()):
        listobs_path_str = group[0]["listobs_path"]
        resolved_path = resolve_local_path(listobs_path_str, search_roots)
        cache_key = str(resolved_path)
        if cache_key not in runtime_cache:
            runtime_cache[cache_key] = parse_eb_runtime_metadata(resolved_path)
        metadata = runtime_cache[cache_key]
        assignments.append(
            {
                "mous_uid": mous_uid,
                "eb_uid": eb_uid,
                "setup_signature": mous_signature_map[mous_uid],
                "actual_n_spws": len(group),
                "actual_total_nchan": sum(int(row["nchan"]) for row in group),
                "science_interval_median_s": metadata.science_interval_median_s,
                "actual_antenna_count": metadata.actual_antenna_count if metadata.actual_antenna_count is not None else "",
                "fiducial_nant": metadata.fiducial_nant if metadata.fiducial_nant is not None else "",
                "array_type": metadata.array_type,
                "data_rate_gbps": compute_eb_data_rate_gbps(group, metadata),
                "listobs_path": metadata.listobs_path,
            }
        )
    return assignments


def summarize_setups(
    mous_assignments: list[dict[str, Any]],
    eb_assignments: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    total_mous = len(mous_assignments)
    total_ebs = len(eb_assignments)
    total_eb_spws = sum(int(row["actual_n_spws"]) for row in eb_assignments)

    mous_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    eb_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in mous_assignments:
        mous_groups[str(row["setup_signature"])].append(row)
    for row in eb_assignments:
        eb_groups[str(row["setup_signature"])].append(row)

    signatures = sorted(set(mous_groups) | set(eb_groups))
    summaries: list[dict[str, Any]] = []
    for signature in signatures:
        mous_members = mous_groups.get(signature, [])
        eb_members = eb_groups.get(signature, [])
        line_fracs = [float(member["line_fraction_percent"]) for member in mous_members]
        data_rates = [float(member["data_rate_gbps"]) for member in eb_members]
        mous_count = len(mous_members)
        eb_count = len(eb_members)
        eb_spw_count = sum(int(member["actual_n_spws"]) for member in eb_members)
        summaries.append(
            {
                "setup_signature": signature,
                "mous_count": mous_count,
                "mous_percent": 100.0 * mous_count / total_mous if total_mous else math.nan,
                "eb_count": eb_count,
                "eb_percent": 100.0 * eb_count / total_ebs if total_ebs else math.nan,
                "eb_spw_count": eb_spw_count,
                "eb_spw_percent": 100.0 * eb_spw_count / total_eb_spws if total_eb_spws else math.nan,
                "line_stat_mous_count": sum(1 for value in line_fracs if math.isfinite(value)),
                "mean_line_fraction_percent": mean(line_fracs),
                "median_line_fraction_percent": finite_median(line_fracs),
                "stddev_line_fraction_percent": population_stddev(line_fracs),
                "mean_eb_data_rate_gbps": mean(data_rates),
            }
        )

    summaries.sort(
        key=lambda row: (
            -int(row["eb_count"]),
            -int(row["eb_spw_count"]),
            -int(row["mous_count"]),
            str(row["setup_signature"]),
        )
    )
    return summaries


def select_top_rows_by_coverage(setup_rows: list[dict[str, Any]], coverage_fraction: float) -> list[dict[str, Any]]:
    total_ebs = sum(int(row["eb_count"]) for row in setup_rows)
    if total_ebs <= 0:
        return []
    target = coverage_fraction * total_ebs
    selected: list[dict[str, Any]] = []
    cumulative = 0
    for row in setup_rows:
        selected.append(row)
        cumulative += int(row["eb_count"])
        if cumulative >= target:
            break
    return selected


def discrete_bin_edges(values: list[float]) -> list[float]:
    unique = sorted(set(values))
    if len(unique) == 1:
        delta = abs(unique[0]) * 0.1 if unique[0] else 0.5
        return [unique[0] - delta, unique[0] + delta]
    mids = [midpoint(unique[idx], unique[idx + 1]) for idx in range(len(unique) - 1)]
    first_edge = unique[0] - (mids[0] - unique[0])
    last_edge = unique[-1] + (unique[-1] - mids[-1])
    return [first_edge] + mids + [last_edge]


def log_hist_edges(values: list[float], min_bins: int = 60, max_bins: int = 120) -> list[float]:
    finite = sorted(value for value in values if value > 0 and math.isfinite(value))
    if len(finite) < 2:
        only = finite[0] if finite else 1.0
        return [only * 0.8, only * 1.2]
    logs = [math.log10(value) for value in finite]
    q1 = percentile(logs, 0.25)
    q3 = percentile(logs, 0.75)
    iqr = q3 - q1
    span = max(logs) - min(logs)
    if iqr <= 0 or span <= 0:
        nbins = min_bins
    else:
        width = 2.0 * iqr * (len(logs) ** (-1.0 / 3.0))
        nbins = int(math.ceil(span / width)) if width > 0 else min_bins
        nbins = max(min_bins, min(max_bins, nbins))
    lo = min(logs)
    hi = max(logs)
    return [10 ** (lo + (hi - lo) * idx / nbins) for idx in range(nbins + 1)]


def plot_histogram(
    values: list[float],
    low_cut: float,
    high_cut: float,
    xlabel: str,
    title: str,
    path: Path,
    bin_edges: list[float],
    *,
    log_x: bool = False,
    xticks: list[float] | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    ax.hist(values, bins=bin_edges, color="#4C78A8", edgecolor="white", alpha=0.9)
    ax.axvline(low_cut, color="#E45756", linestyle="--", linewidth=2, label=f"lower cut = {low_cut:.4g}")
    ax.axvline(high_cut, color="#72B7B2", linestyle="--", linewidth=2, label=f"upper cut = {high_cut:.4g}")
    if log_x:
        ax.set_xscale("log")
    if xticks:
        ax.set_xticks(xticks)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("EB-SPW count")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_summary_markdown(
    path: Path,
    cut_metadata: dict[str, Any],
    setup_rows: list[dict[str, Any]],
    coverage_fraction: float,
    selected_rows: list[dict[str, Any]],
    bandwidth_counts: Counter[str],
    resolution_counts: Counter[str],
) -> None:
    total_mous = sum(int(row["mous_count"]) for row in setup_rows)
    total_ebs = sum(int(row["eb_count"]) for row in setup_rows)
    total_eb_spws = sum(int(row["eb_spw_count"]) for row in setup_rows)
    top_mous = sum(int(row["mous_count"]) for row in selected_rows)
    top_ebs = sum(int(row["eb_count"]) for row in selected_rows)
    top_eb_spws = sum(int(row["eb_spw_count"]) for row in selected_rows)

    bw = cut_metadata["bandwidth"]
    res = cut_metadata["resolution"]
    if bw["mode"] == "fixed":
        bandwidth_cut_line = (
            f"- Bandwidth cuts (GHz, fixed): low={bw['low_cut']:.6f}, "
            f"high={bw['high_cut']:.6f}"
        )
    else:
        bandwidth_cut_line = (
            f"- Bandwidth cuts (GHz, percentile): q33={bw['q33']:.6f}, q66={bw['q66']:.6f}, "
            f"applied at {bw['low_cut']:.6f} between {bw['low_gap'][0]:.6f} and {bw['low_gap'][1]:.6f}, "
            f"and {bw['high_cut']:.6f} between {bw['high_gap'][0]:.6f} and {bw['high_gap'][1]:.6f}"
        )

    lines = [
        "# Spectral Setup Signature Summary",
        "",
        "## Inputs",
        "",
        f"- Unique EB-SPW rows used for binning: {cut_metadata['n_unique_eb_spw']}",
        f"- Total MOUS assigned a setup: {total_mous}",
        f"- Total EBs represented by those MOUS: {total_ebs}",
        f"- Total EB-SPWs represented by those EBs: {total_eb_spws}",
        "",
        "## Cut Definitions",
        "",
        bandwidth_cut_line,
        (
            f"- Resolution cuts (km/s): q33={res['q33']:.6f}, q66={res['q66']:.6f}, "
            f"applied at {res['low_cut']:.6f} between {res['low_gap'][0]:.6f} and {res['low_gap'][1]:.6f}, "
            f"and {res['high_cut']:.6f} between {res['high_gap'][0]:.6f} and {res['high_gap'][1]:.6f}"
        ),
        "- Token convention: first letter is bandwidth (`N`, `M`, `W`), second letter is resolution (`F`, `M`, `C`).",
        "- EB data rate uses fiducial Nant = 43 for 12m and 10 for 7m, with 4 bytes per complex correlation product and the science-scan median integration time per SPW.",
        "",
        "## Category Counts",
        "",
        f"- Bandwidth bins over EB-SPW rows: {dict(sorted(bandwidth_counts.items()))}",
        f"- Resolution bins over EB-SPW rows: {dict(sorted(resolution_counts.items()))}",
        "",
        f"## Top {len(selected_rows)} Setup Signatures (covering {coverage_fraction:.0%} of EBs target)",
        "",
    ]
    for row in selected_rows:
        lines.append(
            "- "
            + f"{row['setup_signature']}: "
            + f"{row['mous_count']} MOUS ({float(row['mous_percent']):.2f}%), "
            + f"{row['eb_count']} EBs ({float(row['eb_percent']):.2f}%), "
            + f"{row['eb_spw_count']} EB-SPWs ({float(row['eb_spw_percent']):.2f}%), "
            + f"mean EB data rate = {float(row['mean_eb_data_rate_gbps']):.4f} GB/s, "
            + f"line fraction mean/median/stddev = "
            + f"{float(row['mean_line_fraction_percent']):.2f}% / "
            + f"{float(row['median_line_fraction_percent']):.2f}% / "
            + f"{float(row['stddev_line_fraction_percent']):.2f}%"
        )
    lines.extend(
        [
            "",
            f"## Coverage Of Top {len(selected_rows)} Setups",
            "",
            f"- Top {len(selected_rows)} setups cover {top_mous} / {total_mous} MOUS ({100.0 * top_mous / total_mous:.2f}%).",
            f"- Top {len(selected_rows)} setups cover {top_ebs} / {total_ebs} EBs ({100.0 * top_ebs / total_ebs:.2f}%).",
            f"- Top {len(selected_rows)} setups cover {top_eb_spws} / {total_eb_spws} EB-SPWs ({100.0 * top_eb_spws / total_eb_spws:.2f}%).",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tabulation_dir", type=Path, help="Directory containing spw_level_table.csv and mous_spw_summary.csv")
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--coverage-fraction", type=float, default=0.85)
    parser.add_argument("--bandwidth-cut-mode", choices=("fixed", "percentile"), default="fixed")
    parser.add_argument("--bandwidth-low-cut-mhz", type=float, default=DEFAULT_BANDWIDTH_LOW_CUT_MHZ)
    parser.add_argument("--bandwidth-high-cut-mhz", type=float, default=DEFAULT_BANDWIDTH_HIGH_CUT_MHZ)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    spw_level_rows = read_csv_rows(args.tabulation_dir / "spw_level_table.csv")
    mous_spw_rows = read_csv_rows(args.tabulation_dir / "mous_spw_summary.csv")

    eb_spw_rows = dedupe_eb_spw_rows(spw_level_rows)
    cut_metadata = derive_cut_metadata(
        eb_spw_rows,
        bandwidth_cut_mode=args.bandwidth_cut_mode,
        bandwidth_low_cut_mhz=args.bandwidth_low_cut_mhz,
        bandwidth_high_cut_mhz=args.bandwidth_high_cut_mhz,
    )
    bandwidth_cuts = cut_metadata["bandwidth"]
    resolution_cuts = cut_metadata["resolution"]

    mous_assignments = build_mous_assignments(mous_spw_rows, bandwidth_cuts, resolution_cuts)
    eb_assignments = build_eb_assignments(eb_spw_rows, mous_assignments, args.tabulation_dir)
    setup_rows = summarize_setups(mous_assignments, eb_assignments)
    selected_rows = select_top_rows_by_coverage(setup_rows, args.coverage_fraction)

    bandwidth_values_ghz = [float(row["total_bw_hz"]) / 1.0e9 for row in eb_spw_rows]
    resolution_values_kms = [float(row["velocity_resolution_kms"]) for row in eb_spw_rows]
    bandwidth_labels = Counter(
        classify(value, bandwidth_cuts["low_cut"], bandwidth_cuts["high_cut"], ("narrow", "medium", "wide"))
        for value in bandwidth_values_ghz
    )
    resolution_labels = Counter(
        classify(value, resolution_cuts["low_cut"], resolution_cuts["high_cut"], ("fine", "mid", "coarse"))
        for value in resolution_values_kms
    )

    write_csv(
        args.output_dir / "mous_setup_assignments.csv",
        mous_assignments,
        ["mous_uid", "n_ebs", "n_spws", "setup_signature", "setup_tokens", "line_fraction_percent"],
    )
    write_csv(
        args.output_dir / "eb_setup_assignments.csv",
        eb_assignments,
        [
            "mous_uid",
            "eb_uid",
            "setup_signature",
            "actual_n_spws",
            "actual_total_nchan",
            "science_interval_median_s",
            "actual_antenna_count",
            "fiducial_nant",
            "array_type",
            "data_rate_gbps",
            "listobs_path",
        ],
    )
    write_csv(
        args.output_dir / "setup_summary.csv",
        setup_rows,
        [
            "setup_signature",
            "mous_count",
            "mous_percent",
            "eb_count",
            "eb_percent",
            "eb_spw_count",
            "eb_spw_percent",
            "line_stat_mous_count",
            "mean_line_fraction_percent",
            "median_line_fraction_percent",
            "stddev_line_fraction_percent",
            "mean_eb_data_rate_gbps",
        ],
    )

    cut_metadata["coverage_fraction"] = args.coverage_fraction
    cut_metadata["selected_setup_count"] = len(selected_rows)
    (args.output_dir / "cut_metadata.json").write_text(json.dumps(cut_metadata, indent=2) + "\n", encoding="utf-8")

    plot_histogram(
        bandwidth_values_ghz,
        bandwidth_cuts["low_cut"],
        bandwidth_cuts["high_cut"],
        "Total bandwidth (GHz)",
        "EB-SPW Bandwidth Distribution",
        args.output_dir / "bandwidth_hist_eb_spw.png",
        discrete_bin_edges(bandwidth_values_ghz),
        xticks=sorted(set(bandwidth_values_ghz)),
    )
    plot_histogram(
        resolution_values_kms,
        resolution_cuts["low_cut"],
        resolution_cuts["high_cut"],
        "Velocity resolution (km/s)",
        "EB-SPW Resolution Distribution",
        args.output_dir / "resolution_hist_eb_spw.png",
        log_hist_edges(resolution_values_kms),
        log_x=True,
    )

    write_summary_markdown(
        args.output_dir / "setup_summary.md",
        cut_metadata,
        setup_rows,
        args.coverage_fraction,
        selected_rows,
        bandwidth_labels,
        resolution_labels,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
