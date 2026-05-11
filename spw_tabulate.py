#!/usr/bin/env python3
"""Tabulate ALMA spectral-window configurations from listobs + cont.dat files.

This script is intentionally stand-alone: it uses only the Python standard
library and writes CSV/Markdown outputs suitable for downstream inspection or
for folding into the big-cubes projection workflow.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median
from typing import Any, Iterable


C_KMS = 299792.458
MOUS_RE = re.compile(r"(uid___A001_X[0-9A-Za-z]+_X[0-9A-Za-z]+)")
EB_RE = re.compile(r"(uid___A002_X[0-9A-Za-z]+_X[0-9A-Za-z]+)")
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
FIELD_ROW_RE = re.compile(r"^\s*(?P<field_id>\d+)\s+\S+\s+(?P<field_name>\S+)\s+")
CONTDAT_FIELD_RE = re.compile(r"^Field:\s*(?P<name>.+?)\s*$")
CONTDAT_SPW_RE = re.compile(r"^SpectralWindow:\s*(?P<spw_id>\d+)\s+(?P<name>\S+)\s*$")
CONTDAT_RANGE_RE = re.compile(
    r"^\s*(?P<low>\d+(?:\.\d+)?)~(?P<high>\d+(?:\.\d+)?)GHz(?:\s+\S+)?\s*$"
)
TARGET_LISTOBS_SUFFIX = "_targets.ms__listobs.txt"


@dataclass
class ParseIssue:
    mous_uid: str
    file_type: str
    path: str
    severity: str
    message: str
    eb_uid: str = ""
    target_name: str = ""
    spw_id: str = ""


@dataclass
class SpwInfo:
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
    def abs_chan_width_hz(self) -> float:
        return abs(self.chan_width_khz) * 1.0e3

    @property
    def total_bw_hz(self) -> float:
        return abs(self.total_bw_khz) * 1.0e3

    @property
    def center_freq_hz(self) -> float:
        return self.center_mhz * 1.0e6

    @property
    def freq_min_hz(self) -> float:
        first_center_hz = self.ch0_mhz * 1.0e6
        last_center_hz = first_center_hz + (self.nchan - 1) * self.chan_width_khz * 1.0e3
        half_width_hz = self.abs_chan_width_hz / 2.0
        return min(first_center_hz, last_center_hz) - half_width_hz

    @property
    def freq_max_hz(self) -> float:
        first_center_hz = self.ch0_mhz * 1.0e6
        last_center_hz = first_center_hz + (self.nchan - 1) * self.chan_width_khz * 1.0e3
        half_width_hz = self.abs_chan_width_hz / 2.0
        return max(first_center_hz, last_center_hz) + half_width_hz

    @property
    def velocity_resolution_kms(self) -> float:
        if self.center_freq_hz <= 0:
            return math.nan
        return C_KMS * self.abs_chan_width_hz / self.center_freq_hz


@dataclass
class ScanInfo:
    scan_id: int
    field_id: int
    field_name: str
    spw_ids: list[int]
    scan_intents: list[str]


@dataclass
class ListobsParseResult:
    mous_uid: str
    eb_uid: str
    path: Path
    spws: dict[int, SpwInfo] = field(default_factory=dict)
    fields: dict[int, str] = field(default_factory=dict)
    scans: list[ScanInfo] = field(default_factory=list)
    issues: list[ParseIssue] = field(default_factory=list)


@dataclass
class ContdatEntry:
    mous_uid: str
    path: Path
    target_name: str
    spw_id: int
    spw_name: str
    flags: list[str]
    continuum_ranges_hz: list[tuple[float, float]]


def parse_uid(value: str, regex: re.Pattern[str]) -> str:
    match = regex.search(value)
    return match.group(1) if match else ""


def parse_int_list(text: str) -> list[int]:
    inner = text.strip()[1:-1].strip()
    if not inner:
        return []
    return [int(item.strip()) for item in inner.split(",") if item.strip()]


def parse_token_list(text: str) -> list[str]:
    inner = text.strip()[1:-1].strip()
    if not inner:
        return []
    return [item.strip() for item in inner.split(",") if item.strip()]


def is_target_listobs_path(path: Path) -> bool:
    return path.name.endswith(TARGET_LISTOBS_SUFFIX)


def merge_intervals(intervals: Iterable[tuple[float, float]]) -> list[tuple[float, float]]:
    cleaned = sorted((min(a, b), max(a, b)) for a, b in intervals if math.isfinite(a) and math.isfinite(b) and a != b)
    if not cleaned:
        return []
    merged: list[tuple[float, float]] = [cleaned[0]]
    for start, end in cleaned[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def complement_intervals(
    continuum: list[tuple[float, float]],
    low_hz: float,
    high_hz: float,
) -> list[tuple[float, float]]:
    merged = merge_intervals(
        (max(low_hz, a), min(high_hz, b))
        for a, b in continuum
        if b > low_hz and a < high_hz
    )
    if not merged:
        return [(low_hz, high_hz)]
    out: list[tuple[float, float]] = []
    cursor = low_hz
    for start, end in merged:
        if start > cursor:
            out.append((cursor, start))
        cursor = max(cursor, end)
    if cursor < high_hz:
        out.append((cursor, high_hz))
    return [(a, b) for a, b in out if b > a]


def buffer_intervals(
    intervals: list[tuple[float, float]],
    buffer_hz: float,
    low_hz: float,
    high_hz: float,
) -> list[tuple[float, float]]:
    buffered = [
        (max(low_hz, start - buffer_hz), min(high_hz, end + buffer_hz))
        for start, end in intervals
    ]
    return merge_intervals(buffered)


def interval_total_width(intervals: Iterable[tuple[float, float]]) -> float:
    return sum(max(0.0, end - start) for start, end in intervals)


def parse_listobs(path: Path, mous_uid: str) -> ListobsParseResult:
    eb_uid = parse_uid(path.name, EB_RE)
    result = ListobsParseResult(mous_uid=mous_uid, eb_uid=eb_uid, path=path)
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()

    in_scans = False
    in_fields = False
    in_spws = False

    for raw in lines:
        line = raw.rstrip()
        stripped = line.strip()
        if stripped.startswith("Date        Timerange"):
            in_scans = True
            in_fields = False
            in_spws = False
            continue
        if stripped.startswith("Fields:"):
            in_fields = True
            in_scans = False
            in_spws = False
            continue
        if stripped.startswith("Spectral Windows:"):
            in_spws = True
            in_fields = False
            in_scans = False
            continue
        if stripped.startswith("Sources:"):
            in_scans = False
            in_fields = False
            in_spws = False
            continue

        if in_scans:
            if "(nRows =" in stripped:
                in_scans = False
                continue
            match = SCAN_ROW_RE.match(line)
            if not match:
                if stripped:
                    result.issues.append(
                        ParseIssue(
                            mous_uid=mous_uid,
                            file_type="listobs",
                            path=str(path),
                            severity="warning",
                            eb_uid=eb_uid,
                            message=f"Unparsed scan row: {stripped}",
                        )
                    )
                continue
            result.scans.append(
                ScanInfo(
                    scan_id=int(match.group("scan")),
                    field_id=int(match.group("field_id")),
                    field_name=match.group("field_name"),
                    spw_ids=parse_int_list(match.group("spw_ids")),
                    scan_intents=parse_token_list(match.group("intents")),
                )
            )
            continue

        if in_fields:
            if stripped.startswith("ID   Code"):
                continue
            if not stripped:
                continue
            match = FIELD_ROW_RE.match(line)
            if match:
                result.fields[int(match.group("field_id"))] = match.group("field_name")
            continue

        if in_spws:
            if stripped.startswith("SpwID"):
                continue
            if not stripped:
                continue
            if stripped.startswith("Sources:"):
                in_spws = False
                continue
            match = SPW_ROW_RE.match(line)
            if not match:
                result.issues.append(
                    ParseIssue(
                        mous_uid=mous_uid,
                        file_type="listobs",
                        path=str(path),
                        severity="warning",
                        eb_uid=eb_uid,
                        message=f"Unparsed spectral-window row: {stripped}",
                    )
                )
                continue
            spw = SpwInfo(
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
            result.spws[spw.spw_id] = spw

    if not result.spws:
        result.issues.append(
            ParseIssue(
                mous_uid=mous_uid,
                file_type="listobs",
                path=str(path),
                severity="error",
                eb_uid=eb_uid,
                message="No spectral-window table parsed",
            )
        )
    if not result.scans:
        result.issues.append(
            ParseIssue(
                mous_uid=mous_uid,
                file_type="listobs",
                path=str(path),
                severity="error",
                eb_uid=eb_uid,
                message="No scan table parsed",
            )
        )
    return result


def extract_spw_level_rows(result: ListobsParseResult) -> list[dict[str, Any]]:
    spw_usage: dict[tuple[int, str], set[int]] = defaultdict(set)
    warnings_by_target: dict[tuple[int, str], list[str]] = defaultdict(list)

    for scan in result.scans:
        if not any("OBSERVE_TARGET" in intent for intent in scan.scan_intents):
            continue
        field_name = result.fields.get(scan.field_id, scan.field_name)
        key = (scan.field_id, field_name)
        for spw_id in scan.spw_ids:
            spw_usage[key].add(spw_id)
            if spw_id not in result.spws:
                warnings_by_target[key].append(f"scan references missing spw_id={spw_id}")

    rows: list[dict[str, Any]] = []
    for (field_id, field_name), spw_ids in sorted(spw_usage.items()):
        field_warnings = warnings_by_target.get((field_id, field_name), [])
        for spw_id in sorted(spw_ids):
            spw = result.spws.get(spw_id)
            if spw is None:
                continue
            rows.append(
                {
                    "mous_uid": result.mous_uid,
                    "eb_uid": result.eb_uid,
                    "listobs_path": str(result.path),
                    "target_name": field_name,
                    "field_id": field_id,
                    "spw_id": spw.spw_id,
                    "spw_name": spw.spw_name,
                    "nchan": spw.nchan,
                    "center_freq_hz": spw.center_freq_hz,
                    "freq_min_hz": spw.freq_min_hz,
                    "freq_max_hz": spw.freq_max_hz,
                    "total_bw_hz": spw.total_bw_hz,
                    "abs_chan_width_hz": spw.abs_chan_width_hz,
                    "velocity_resolution_kms": spw.velocity_resolution_kms,
                    "corr_products": spw.corrs,
                    "science_target_flag": True,
                    "parse_warning": "; ".join(sorted(set(field_warnings))),
                }
            )
    return rows


def parse_contdat(path: Path, mous_uid: str) -> tuple[list[ContdatEntry], list[ParseIssue]]:
    issues: list[ParseIssue] = []
    entries: list[ContdatEntry] = []
    current_target = ""
    current_spw_id: int | None = None
    current_spw_name = ""
    current_flags: list[str] = []
    current_ranges_hz: list[tuple[float, float]] = []

    def flush_current() -> None:
        nonlocal current_spw_id, current_spw_name, current_flags, current_ranges_hz
        if current_target and current_spw_id is not None:
            entries.append(
                ContdatEntry(
                    mous_uid=mous_uid,
                    path=path,
                    target_name=current_target,
                    spw_id=current_spw_id,
                    spw_name=current_spw_name,
                    flags=list(current_flags),
                    continuum_ranges_hz=merge_intervals(current_ranges_hz),
                )
            )
        current_spw_id = None
        current_spw_name = ""
        current_flags = []
        current_ranges_hz = []

    for line_no, raw in enumerate(path.read_text(encoding="utf-8", errors="ignore").splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        match = CONTDAT_FIELD_RE.match(line)
        if match:
            flush_current()
            current_target = match.group("name")
            continue
        match = CONTDAT_SPW_RE.match(line)
        if match:
            flush_current()
            current_spw_id = int(match.group("spw_id"))
            current_spw_name = match.group("name")
            continue
        if line.startswith("Flags:"):
            current_flags.append(line.split(":", 1)[1].strip())
            continue
        match = CONTDAT_RANGE_RE.match(line)
        if match:
            low_hz = float(match.group("low")) * 1.0e9
            high_hz = float(match.group("high")) * 1.0e9
            current_ranges_hz.append((low_hz, high_hz))
            continue
        issues.append(
            ParseIssue(
                mous_uid=mous_uid,
                file_type="contdat",
                path=str(path),
                severity="warning",
                message=f"Unparsed cont.dat line {line_no}: {line}",
            )
        )
    flush_current()
    return entries, issues


def safe_mean(values: Iterable[float]) -> float:
    finite = [value for value in values if value is not None and math.isfinite(value)]
    if not finite:
        return math.nan
    return sum(finite) / len(finite)


def safe_median(values: Iterable[float]) -> float:
    finite = [value for value in values if value is not None and math.isfinite(value)]
    if not finite:
        return math.nan
    return float(median(finite))


def build_spw_lookup(spw_rows: list[dict[str, Any]]) -> tuple[dict[tuple[str, int], dict[str, Any]], dict[int, dict[str, Any]]]:
    by_target: dict[tuple[str, int], dict[str, Any]] = {}
    by_spw: dict[int, dict[str, Any]] = {}
    for row in spw_rows:
        key = (str(row["target_name"]), int(row["spw_id"]))
        by_target.setdefault(key, row)
        by_spw.setdefault(int(row["spw_id"]), row)
    return by_target, by_spw


def summarize_contdat_entry(
    entry: ContdatEntry,
    spw_row: dict[str, Any],
    issues: list[ParseIssue],
) -> dict[str, Any]:
    total_bw_hz = float(spw_row["total_bw_hz"])
    spw_min_hz = float(spw_row["freq_min_hz"])
    spw_max_hz = float(spw_row["freq_max_hz"])
    chan_width_hz = float(spw_row["abs_chan_width_hz"])
    continuum_ranges = merge_intervals(entry.continuum_ranges_hz)
    line_ranges = complement_intervals(continuum_ranges, spw_min_hz, spw_max_hz)
    buffer_hz = 10.0 * chan_width_hz

    if not line_ranges:
        line_range_min_hz = math.nan
        line_range_max_hz = math.nan
        line_range_width_hz = 0.0
        buffered_line_range_min_hz = math.nan
        buffered_line_range_max_hz = math.nan
        buffered_line_range_width_hz = 0.0
        summed_line_width_hz = 0.0
        buffered_summed_line_width_hz = 0.0
        buffered_line_ranges: list[tuple[float, float]] = []
        warning = "no apparent line range after continuum complement"
    else:
        line_range_min_hz = min(start for start, _ in line_ranges)
        line_range_max_hz = max(end for _, end in line_ranges)
        line_range_width_hz = line_range_max_hz - line_range_min_hz
        buffered_line_range_min_hz = max(spw_min_hz, line_range_min_hz - buffer_hz)
        buffered_line_range_max_hz = min(spw_max_hz, line_range_max_hz + buffer_hz)
        buffered_line_range_width_hz = buffered_line_range_max_hz - buffered_line_range_min_hz
        summed_line_width_hz = interval_total_width(line_ranges)
        buffered_line_ranges = buffer_intervals(line_ranges, buffer_hz, spw_min_hz, spw_max_hz)
        buffered_summed_line_width_hz = interval_total_width(buffered_line_ranges)
        warning = ""

    known_flags = {"ALLCONT", "LOWSPREAD", "LOWBANDWIDTH"}
    unknown_flags = [flag for flag in entry.flags if flag not in known_flags]
    if unknown_flags:
        warning = "; ".join(filter(None, [warning, f"unrecognized flags={','.join(sorted(unknown_flags))}"]))

    if warning:
        issues.append(
            ParseIssue(
                mous_uid=entry.mous_uid,
                file_type="contdat",
                path=str(entry.path),
                severity="warning",
                target_name=entry.target_name,
                spw_id=str(entry.spw_id),
                message=warning,
            )
        )

    return {
        "mous_uid": entry.mous_uid,
        "eb_uid": "",
        "contdat_path": str(entry.path),
        "target_name": entry.target_name,
        "spw_id": entry.spw_id,
        "continuum_ranges_hz": json.dumps(continuum_ranges),
        "disjoint_line_ranges_hz": json.dumps(line_ranges),
        "buffered_disjoint_line_ranges_hz": json.dumps(buffered_line_ranges),
        "line_range_min_hz": line_range_min_hz,
        "line_range_max_hz": line_range_max_hz,
        "buffered_line_range_min_hz": buffered_line_range_min_hz,
        "buffered_line_range_max_hz": buffered_line_range_max_hz,
        "line_range_width_hz": line_range_width_hz,
        "buffered_line_range_width_hz": buffered_line_range_width_hz,
        "summed_line_width_hz": summed_line_width_hz,
        "buffered_summed_line_width_hz": buffered_summed_line_width_hz,
        "line_range_fraction_of_spw": line_range_width_hz / total_bw_hz if total_bw_hz else math.nan,
        "buffered_line_range_fraction_of_spw": buffered_line_range_width_hz / total_bw_hz if total_bw_hz else math.nan,
        "summed_line_fraction_of_spw": summed_line_width_hz / total_bw_hz if total_bw_hz else math.nan,
        "buffered_summed_line_fraction_of_spw": buffered_summed_line_width_hz / total_bw_hz if total_bw_hz else math.nan,
        "contdat_parse_warning": warning,
    }


def aggregate_mous_spw_summary(
    mous_uid: str,
    spw_rows: list[dict[str, Any]],
    cont_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    cont_by_spw: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in cont_rows:
        cont_by_spw[int(row["spw_id"])].append(row)

    grouped: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in spw_rows:
        grouped[int(row["spw_id"])].append(row)

    out: list[dict[str, Any]] = []
    for spw_id, rows in sorted(grouped.items()):
        unique_eb_rows: dict[str, dict[str, Any]] = {}
        for row in rows:
            unique_eb_rows.setdefault(str(row["eb_uid"]), row)
        eb_rows = list(unique_eb_rows.values())
        rep = sorted(eb_rows, key=lambda item: str(item["eb_uid"]))[0]

        nchan_values = {int(item["nchan"]) for item in eb_rows}
        bw_values = [float(item["total_bw_hz"]) for item in eb_rows]
        chan_values = [float(item["abs_chan_width_hz"]) for item in eb_rows]
        ctr_values = [float(item["center_freq_hz"]) for item in eb_rows]

        eb_consistency_flag = "consistent"
        notes: list[str] = []
        if len(nchan_values) > 1:
            eb_consistency_flag = "inconsistent"
            notes.append(f"nchan values={sorted(nchan_values)}")
        if max(bw_values) - min(bw_values) > max(1.0, 1.0e-6 * max(bw_values)):
            eb_consistency_flag = "inconsistent"
            notes.append("total bandwidth differs across EBs")
        if max(chan_values) - min(chan_values) > max(1.0, 1.0e-6 * max(chan_values)):
            eb_consistency_flag = "inconsistent"
            notes.append("channel width differs across EBs")
        if max(ctr_values) - min(ctr_values) > max(1.0, 1.0e-6 * max(ctr_values)):
            if eb_consistency_flag == "consistent":
                eb_consistency_flag = "minor_differences"
            notes.append("center frequency differs slightly across EBs")

        cont_matches = cont_by_spw.get(spw_id, [])
        out.append(
            {
                "mous_uid": mous_uid,
                "representative_eb_uid": rep["eb_uid"],
                "spw_id": spw_id,
                "nchan": round(safe_median(float(item["nchan"]) for item in eb_rows)),
                "center_freq_hz": safe_median(float(item["center_freq_hz"]) for item in eb_rows),
                "freq_min_hz": safe_median(float(item["freq_min_hz"]) for item in eb_rows),
                "freq_max_hz": safe_median(float(item["freq_max_hz"]) for item in eb_rows),
                "total_bw_hz": safe_median(float(item["total_bw_hz"]) for item in eb_rows),
                "abs_chan_width_hz": safe_median(float(item["abs_chan_width_hz"]) for item in eb_rows),
                "velocity_resolution_kms": safe_median(float(item["velocity_resolution_kms"]) for item in eb_rows),
                "mean_line_range_fraction_of_spw": safe_mean(float(item["line_range_fraction_of_spw"]) for item in cont_matches),
                "mean_buffered_line_range_fraction_of_spw": safe_mean(float(item["buffered_line_range_fraction_of_spw"]) for item in cont_matches),
                "mean_summed_line_fraction_of_spw": safe_mean(float(item["summed_line_fraction_of_spw"]) for item in cont_matches),
                "mean_buffered_summed_line_fraction_of_spw": safe_mean(float(item["buffered_summed_line_fraction_of_spw"]) for item in cont_matches),
                "n_targets_with_contdat": len({str(item["target_name"]) for item in cont_matches}),
                "n_ebs_parsed": len(unique_eb_rows),
                "eb_consistency_flag": eb_consistency_flag,
                "eb_consistency_notes": "; ".join(notes),
            }
        )
    return out


def load_metadata(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    out: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            mous_uid = (row.get("member_ous_uid") or row.get("mous_uid") or "").strip()
            if not mous_uid:
                continue
            out[mous_uid] = {key: (value or "").strip() for key, value in row.items()}
    return out


def quantile_edges(values: list[float]) -> tuple[float, float]:
    finite = sorted(value for value in values if math.isfinite(value))
    if len(finite) < 3:
        return (math.nan, math.nan)
    def pick(frac: float) -> float:
        idx = int(round((len(finite) - 1) * frac))
        return finite[max(0, min(idx, len(finite) - 1))]
    return pick(1.0 / 3.0), pick(2.0 / 3.0)


def classify_value(value: float, low_edge: float, high_edge: float, labels: tuple[str, str, str]) -> str:
    if not math.isfinite(value):
        return "unknown"
    if math.isfinite(low_edge) and value <= low_edge:
        return labels[0]
    if math.isfinite(high_edge) and value <= high_edge:
        return labels[1]
    return labels[2]


def annotate_setups(
    mous_spw_rows: list[dict[str, Any]],
    metadata: dict[str, dict[str, str]],
    binning_mode: str,
    velocity_edges_kms: tuple[float, float],
    bandwidth_edges_hz: tuple[float, float],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if binning_mode == "tertile":
        velocity_edges_kms = quantile_edges([float(row["velocity_resolution_kms"]) for row in mous_spw_rows])
        bandwidth_edges_hz = quantile_edges([float(row["total_bw_hz"]) for row in mous_spw_rows])

    by_mous: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in mous_spw_rows:
        by_mous[str(row["mous_uid"])].append(row)

    setup_rows: list[dict[str, Any]] = []
    signature_rows: list[dict[str, Any]] = []
    for mous_uid, rows in sorted(by_mous.items()):
        metadata_row = metadata.get(mous_uid, {})
        classes: list[str] = []
        for row in rows:
            bw_class = classify_value(float(row["total_bw_hz"]), bandwidth_edges_hz[0], bandwidth_edges_hz[1], ("narrow", "medium", "wide"))
            vel_class = classify_value(float(row["velocity_resolution_kms"]), velocity_edges_kms[0], velocity_edges_kms[1], ("fine", "medium", "coarse"))
            classes.append(f"{vel_class}-{bw_class}")
        class_counts = Counter(classes)
        canonical_signature = " + ".join(
            f"{count}x_{label}" for label, count in sorted(class_counts.items())
        )
        consistency_flags = {str(row["eb_consistency_flag"]) for row in rows}
        eb_consistency_flag = "consistent" if consistency_flags == {"consistent"} else ",".join(sorted(consistency_flags))
        total_science_bw_hz = sum(float(row["total_bw_hz"]) for row in rows if math.isfinite(float(row["total_bw_hz"])))
        freq_min = min(float(row["freq_min_hz"]) for row in rows if math.isfinite(float(row["freq_min_hz"])))
        freq_max = max(float(row["freq_max_hz"]) for row in rows if math.isfinite(float(row["freq_max_hz"])))
        setup_row = {
            "mous_uid": mous_uid,
            "representative_eb_uid": rows[0]["representative_eb_uid"],
            "n_science_spw": len(rows),
            "total_science_bw_hz": total_science_bw_hz,
            "min_channel_width_hz": min(float(row["abs_chan_width_hz"]) for row in rows),
            "median_channel_width_hz": safe_median(float(row["abs_chan_width_hz"]) for row in rows),
            "max_channel_width_hz": max(float(row["abs_chan_width_hz"]) for row in rows),
            "min_velocity_resolution_kms": min(float(row["velocity_resolution_kms"]) for row in rows),
            "median_velocity_resolution_kms": safe_median(float(row["velocity_resolution_kms"]) for row in rows),
            "max_velocity_resolution_kms": max(float(row["velocity_resolution_kms"]) for row in rows),
            "min_spw_bandwidth_hz": min(float(row["total_bw_hz"]) for row in rows),
            "median_spw_bandwidth_hz": safe_median(float(row["total_bw_hz"]) for row in rows),
            "max_spw_bandwidth_hz": max(float(row["total_bw_hz"]) for row in rows),
            "frequency_span_hz": freq_max - freq_min,
            "mean_line_range_fraction_of_spw": safe_mean(float(row["mean_line_range_fraction_of_spw"]) for row in rows),
            "mean_buffered_line_range_fraction_of_spw": safe_mean(float(row["mean_buffered_line_range_fraction_of_spw"]) for row in rows),
            "mean_summed_line_fraction_of_spw": safe_mean(float(row["mean_summed_line_fraction_of_spw"]) for row in rows),
            "mean_buffered_summed_line_fraction_of_spw": safe_mean(float(row["mean_buffered_summed_line_fraction_of_spw"]) for row in rows),
            "number_of_distinct_spw_bandwidth_classes": len({classify_value(float(row["total_bw_hz"]), bandwidth_edges_hz[0], bandwidth_edges_hz[1], ("narrow", "medium", "wide")) for row in rows}),
            "number_of_distinct_velocity_resolution_classes": len({classify_value(float(row["velocity_resolution_kms"]), velocity_edges_kms[0], velocity_edges_kms[1], ("fine", "medium", "coarse")) for row in rows}),
            "canonical_setup_signature": canonical_signature,
            "setup_family": canonical_signature,
            "eb_consistency_flag": eb_consistency_flag,
            "science_category": metadata_row.get("science_category") or metadata_row.get("scientific_category") or "",
            "band": metadata_row.get("band", ""),
            "array_type": metadata_row.get("array_type") or metadata_row.get("config_proxy") or "",
            "project_code": metadata_row.get("project_code", ""),
        }
        setup_rows.append(setup_row)
        signature_rows.append(
            {
                "mous_uid": mous_uid,
                "canonical_setup_signature": canonical_signature,
                "science_category": setup_row["science_category"],
                "band": setup_row["band"],
                "array_type": setup_row["array_type"],
                "project_code": setup_row["project_code"],
                "n_science_spw": len(rows),
            }
        )
    return setup_rows, signature_rows


def summarize_signatures(
    signature_rows: list[dict[str, Any]],
    group_field: str | None = None,
) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in signature_rows:
        group_value = "overall" if group_field is None else str(row.get(group_field) or "unknown")
        grouped[(group_value, str(row["canonical_setup_signature"]))].append(row)

    out: list[dict[str, Any]] = []
    for (group_value, signature), rows in sorted(grouped.items()):
        out.append(
            {
                "group_value": group_value,
                "canonical_setup_signature": signature,
                "n_mous": len(rows),
                "mean_n_science_spw": safe_mean(float(row["n_science_spw"]) for row in rows),
            }
        )
    return out


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_summary(
    path: Path,
    spw_rows: list[dict[str, Any]],
    cont_rows: list[dict[str, Any]],
    mous_spw_rows: list[dict[str, Any]],
    setup_rows: list[dict[str, Any]],
    issues: list[ParseIssue],
) -> None:
    mous_count = len({row["mous_uid"] for row in setup_rows})
    eb_count = len({row["eb_uid"] for row in spw_rows})
    signature_counts = Counter(str(row["canonical_setup_signature"]) for row in setup_rows)
    nspw_counts = Counter(int(row["n_science_spw"]) for row in setup_rows)
    issue_counts = Counter((issue.severity, issue.file_type) for issue in issues)

    lines = [
        "# Spectral Setup Summary",
        "",
        f"- MOUS processed: {mous_count}",
        f"- EBs parsed: {eb_count}",
        f"- SPW-level rows: {len(spw_rows)}",
        f"- cont.dat rows: {len(cont_rows)}",
        f"- MOUS/SPW summary rows: {len(mous_spw_rows)}",
        f"- setup-level rows: {len(setup_rows)}",
        "",
        "## NSPW Distribution",
        "",
    ]
    for nspw, count in sorted(nspw_counts.items()):
        lines.append(f"- {nspw} science SPWs: {count} MOUS")
    lines.extend(["", "## Top Setup Signatures", ""])
    for signature, count in signature_counts.most_common(20):
        lines.append(f"- {signature}: {count}")
    lines.extend(["", "## Parse Issues", ""])
    if not issue_counts:
        lines.append("- none")
    else:
        for (severity, file_type), count in sorted(issue_counts.items()):
            lines.append(f"- {severity}/{file_type}: {count}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def process_mous_dir(mous_dir: Path, issues: list[ParseIssue]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    mous_uid = parse_uid(mous_dir.name, MOUS_RE) or mous_dir.name
    listobs_files = sorted(
        path for path in (mous_dir / "listobs").glob("*.txt")
        if is_target_listobs_path(path)
    )
    spw_rows: list[dict[str, Any]] = []
    for listobs_file in listobs_files:
        result = parse_listobs(listobs_file, mous_uid)
        issues.extend(result.issues)
        spw_rows.extend(extract_spw_level_rows(result))

    cont_rows: list[dict[str, Any]] = []
    cont_path = mous_dir / "cont.dat"
    if cont_path.exists():
        entries, cont_issues = parse_contdat(cont_path, mous_uid)
        issues.extend(cont_issues)
        by_target, by_spw = build_spw_lookup(spw_rows)
        for entry in entries:
            spw_row = by_target.get((entry.target_name, entry.spw_id))
            if spw_row is None:
                spw_row = by_spw.get(entry.spw_id)
            if spw_row is None:
                issues.append(
                    ParseIssue(
                        mous_uid=mous_uid,
                        file_type="contdat",
                        path=str(cont_path),
                        severity="warning",
                        target_name=entry.target_name,
                        spw_id=str(entry.spw_id),
                        message="Could not match cont.dat SPW to parsed listobs science SPW",
                    )
                )
                cont_rows.append(
                    {
                        "mous_uid": mous_uid,
                        "eb_uid": "",
                        "contdat_path": str(cont_path),
                        "target_name": entry.target_name,
                        "spw_id": entry.spw_id,
                        "continuum_ranges_hz": json.dumps(merge_intervals(entry.continuum_ranges_hz)),
                        "disjoint_line_ranges_hz": json.dumps([]),
                        "buffered_disjoint_line_ranges_hz": json.dumps([]),
                        "line_range_min_hz": math.nan,
                        "line_range_max_hz": math.nan,
                        "buffered_line_range_min_hz": math.nan,
                        "buffered_line_range_max_hz": math.nan,
                        "line_range_width_hz": math.nan,
                        "buffered_line_range_width_hz": math.nan,
                        "summed_line_width_hz": math.nan,
                        "buffered_summed_line_width_hz": math.nan,
                        "line_range_fraction_of_spw": math.nan,
                        "buffered_line_range_fraction_of_spw": math.nan,
                        "summed_line_fraction_of_spw": math.nan,
                        "buffered_summed_line_fraction_of_spw": math.nan,
                        "contdat_parse_warning": "missing matching science spw metadata",
                    }
                )
                continue
            cont_rows.append(summarize_contdat_entry(entry, spw_row, issues))
    else:
        issues.append(
            ParseIssue(
                mous_uid=mous_uid,
                file_type="contdat",
                path=str(cont_path),
                severity="warning",
                message="cont.dat not present",
            )
        )

    return spw_rows, cont_rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, help="Root directory containing one MOUS directory per entry")
    parser.add_argument("output_dir", type=Path, help="Directory to write output CSV/summary files")
    parser.add_argument("--metadata-csv", type=Path, default=None, help="Optional metadata CSV keyed by member_ous_uid")
    parser.add_argument("--binning-mode", choices=["fixed", "tertile"], default="fixed")
    parser.add_argument("--velocity-bin-edges-kms", default="0.5,2.0", help="Comma-separated low,high edges for fine/medium/coarse")
    parser.add_argument("--bandwidth-bin-edges-ghz", default="0.25,1.0", help="Comma-separated low,high edges for narrow/medium/wide")
    parser.add_argument("--max-mous", type=int, default=None, help="Optional limit for smoke runs")
    args = parser.parse_args(argv)

    metadata = load_metadata(args.metadata_csv)
    issues: list[ParseIssue] = []
    spw_rows: list[dict[str, Any]] = []
    cont_rows: list[dict[str, Any]] = []

    mous_dirs = sorted(path for path in args.input_dir.iterdir() if path.is_dir() and path.name.startswith("uid"))
    if args.max_mous is not None:
        mous_dirs = mous_dirs[: args.max_mous]

    for mous_dir in mous_dirs:
        mous_spw_rows, mous_cont_rows = process_mous_dir(mous_dir, issues)
        spw_rows.extend(mous_spw_rows)
        cont_rows.extend(mous_cont_rows)

    mous_spw_summary_rows: list[dict[str, Any]] = []
    spw_by_mous: dict[str, list[dict[str, Any]]] = defaultdict(list)
    cont_by_mous: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in spw_rows:
        spw_by_mous[str(row["mous_uid"])].append(row)
    for row in cont_rows:
        cont_by_mous[str(row["mous_uid"])].append(row)
    for mous_uid in sorted(spw_by_mous):
        mous_spw_summary_rows.extend(
            aggregate_mous_spw_summary(mous_uid, spw_by_mous[mous_uid], cont_by_mous.get(mous_uid, []))
        )

    vel_edges_parts = [part.strip() for part in args.velocity_bin_edges_kms.split(",") if part.strip()]
    bw_edges_parts = [part.strip() for part in args.bandwidth_bin_edges_ghz.split(",") if part.strip()]
    velocity_edges_kms = (
        float(vel_edges_parts[0]) if len(vel_edges_parts) > 0 else math.nan,
        float(vel_edges_parts[1]) if len(vel_edges_parts) > 1 else math.nan,
    )
    bandwidth_edges_hz = (
        float(bw_edges_parts[0]) * 1.0e9 if len(bw_edges_parts) > 0 else math.nan,
        float(bw_edges_parts[1]) * 1.0e9 if len(bw_edges_parts) > 1 else math.nan,
    )
    setup_rows, signature_rows = annotate_setups(
        mous_spw_summary_rows,
        metadata,
        args.binning_mode,
        velocity_edges_kms,
        bandwidth_edges_hz,
    )

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    write_csv(
        output_dir / "spw_level_table.csv",
        spw_rows,
        [
            "mous_uid",
            "eb_uid",
            "listobs_path",
            "target_name",
            "field_id",
            "spw_id",
            "spw_name",
            "nchan",
            "center_freq_hz",
            "freq_min_hz",
            "freq_max_hz",
            "total_bw_hz",
            "abs_chan_width_hz",
            "velocity_resolution_kms",
            "corr_products",
            "science_target_flag",
            "parse_warning",
        ],
    )
    write_csv(
        output_dir / "contdat_line_ranges.csv",
        cont_rows,
        [
            "mous_uid",
            "eb_uid",
            "contdat_path",
            "target_name",
            "spw_id",
            "continuum_ranges_hz",
            "disjoint_line_ranges_hz",
            "buffered_disjoint_line_ranges_hz",
            "line_range_min_hz",
            "line_range_max_hz",
            "buffered_line_range_min_hz",
            "buffered_line_range_max_hz",
            "line_range_width_hz",
            "buffered_line_range_width_hz",
            "summed_line_width_hz",
            "buffered_summed_line_width_hz",
            "line_range_fraction_of_spw",
            "buffered_line_range_fraction_of_spw",
            "summed_line_fraction_of_spw",
            "buffered_summed_line_fraction_of_spw",
            "contdat_parse_warning",
        ],
    )
    write_csv(
        output_dir / "mous_spw_summary.csv",
        mous_spw_summary_rows,
        [
            "mous_uid",
            "representative_eb_uid",
            "spw_id",
            "nchan",
            "center_freq_hz",
            "freq_min_hz",
            "freq_max_hz",
            "total_bw_hz",
            "abs_chan_width_hz",
            "velocity_resolution_kms",
            "mean_line_range_fraction_of_spw",
            "mean_buffered_line_range_fraction_of_spw",
            "mean_summed_line_fraction_of_spw",
            "mean_buffered_summed_line_fraction_of_spw",
            "n_targets_with_contdat",
            "n_ebs_parsed",
            "eb_consistency_flag",
            "eb_consistency_notes",
        ],
    )
    write_csv(
        output_dir / "setup_level_table.csv",
        setup_rows,
        [
            "mous_uid",
            "representative_eb_uid",
            "n_science_spw",
            "total_science_bw_hz",
            "min_channel_width_hz",
            "median_channel_width_hz",
            "max_channel_width_hz",
            "min_velocity_resolution_kms",
            "median_velocity_resolution_kms",
            "max_velocity_resolution_kms",
            "min_spw_bandwidth_hz",
            "median_spw_bandwidth_hz",
            "max_spw_bandwidth_hz",
            "frequency_span_hz",
            "mean_line_range_fraction_of_spw",
            "mean_buffered_line_range_fraction_of_spw",
            "mean_summed_line_fraction_of_spw",
            "mean_buffered_summed_line_fraction_of_spw",
            "number_of_distinct_spw_bandwidth_classes",
            "number_of_distinct_velocity_resolution_classes",
            "canonical_setup_signature",
            "setup_family",
            "eb_consistency_flag",
            "science_category",
            "band",
            "array_type",
            "project_code",
        ],
    )
    write_csv(
        output_dir / "signature_summary_overall.csv",
        summarize_signatures(signature_rows),
        ["group_value", "canonical_setup_signature", "n_mous", "mean_n_science_spw"],
    )
    if metadata:
        write_csv(
            output_dir / "signature_summary_by_science_category.csv",
            summarize_signatures(signature_rows, "science_category"),
            ["group_value", "canonical_setup_signature", "n_mous", "mean_n_science_spw"],
        )
        write_csv(
            output_dir / "signature_summary_by_band.csv",
            summarize_signatures(signature_rows, "band"),
            ["group_value", "canonical_setup_signature", "n_mous", "mean_n_science_spw"],
        )
        write_csv(
            output_dir / "signature_summary_by_array_or_config.csv",
            summarize_signatures(signature_rows, "array_type"),
            ["group_value", "canonical_setup_signature", "n_mous", "mean_n_science_spw"],
        )
    write_csv(
        output_dir / "parse_errors.csv",
        [issue.__dict__ for issue in issues],
        ["mous_uid", "file_type", "path", "severity", "message", "eb_uid", "target_name", "spw_id"],
    )
    write_summary(
        output_dir / "summary_overall.md",
        spw_rows,
        cont_rows,
        mous_spw_summary_rows,
        setup_rows,
        issues,
    )

    print(f"Wrote outputs to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
