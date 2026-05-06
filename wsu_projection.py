"""Shared WSU nchan projection helpers.

This module centralizes the TALON/stepped2 channelization logic so the
main recompute path and sample-driven analyses can share one implementation.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from config import DEFAULT_CONFIG


C_KMS = 299792.458

ALMA_BANDS_GHZ = {
    1: (35.0, 50.0),
    2: (67.0, 90.0),
    3: (84.0, 116.0),
    4: (125.0, 163.0),
    5: (163.0, 211.0),
    6: (211.0, 275.0),
    7: (275.0, 373.0),
    8: (385.0, 500.0),
    9: (602.0, 720.0),
    10: (787.0, 950.0),
}


@dataclass(frozen=True)
class MousSpwTemplate:
    mous_uid: str
    spw_id: int
    center_freq_hz: float
    velocity_resolution_kms: float
    band: int


def normalize_mous_uid(value: str) -> str:
    """Normalize uid:// and uid___ MOUS identifiers to one form."""
    text = str(value).strip()
    if text.startswith("uid://"):
        return text.replace("uid://", "uid___").replace("/", "_")
    return text


def alma_band_from_freq_hz(freq_hz: float) -> int:
    freq_ghz = freq_hz / 1.0e9
    for band, (lo, hi) in ALMA_BANDS_GHZ.items():
        if lo <= freq_ghz <= hi:
            return band
    raise ValueError(f"frequency {freq_ghz:.3f} GHz is outside ALMA bands 1-10")


def stepped2_floor_velocity(current_vel_kms: float) -> float:
    if current_vel_kms > 10.0:
        return 10.0
    if current_vel_kms > 2.0:
        return 2.0
    if current_vel_kms > 0.5:
        return 0.5
    if current_vel_kms > 0.1:
        return 0.1
    return current_vel_kms


def calc_talon_specwidth_khz(
    requested_vel_kms: float,
    original_vel_kms: float,
    center_freq_hz: float,
    band: int,
    milestone: str,
) -> tuple[float, float]:
    """Return realizable TALON channel width and chanavg for a request."""
    requested_specwidth_khz = (requested_vel_kms / C_KMS) * center_freq_hz / 1.0e3
    chanavg = math.floor(requested_specwidth_khz / DEFAULT_CONFIG.talon_channel_kHz)
    if chanavg < 1:
        chanavg = 1

    min_lookup = (
        DEFAULT_CONFIG.chanavg_min_initial if milestone == "M1"
        else DEFAULT_CONFIG.chanavg_min
    )
    if 0.095 < original_vel_kms <= 0.5 and chanavg < min_lookup[band]:
        chanavg = int(min_lookup[band])

    specwidth_khz = DEFAULT_CONFIG.talon_channel_kHz * chanavg
    return specwidth_khz, float(chanavg)


def projected_nchan_per_spw(
    requested_vel_kms: float,
    original_vel_kms: float,
    center_freq_hz: float,
    band: int,
    milestone: str,
) -> float:
    specwidth_khz, _ = calc_talon_specwidth_khz(
        requested_vel_kms=requested_vel_kms,
        original_vel_kms=original_vel_kms,
        center_freq_hz=center_freq_hz,
        band=band,
        milestone=milestone,
    )
    return math.floor((DEFAULT_CONFIG.spw_bandwidth_GHz * 1.0e6) / specwidth_khz)


def projected_spw_equivalents(milestone: str, band: int) -> float:
    """Return the milestone-specific WSU SPW-equivalent count for a band."""
    bandwidth_ghz = DEFAULT_CONFIG.milestones[milestone].bandwidth_for(band)
    if milestone == "M1":
        return math.ceil(bandwidth_ghz / DEFAULT_CONFIG.spw_bandwidth_GHz)
    return bandwidth_ghz / DEFAULT_CONFIG.spw_bandwidth_GHz


def project_nchan_agg_for_templates(
    mode: str,
    milestone: str,
    templates: list[MousSpwTemplate],
    projected_nspw: float,
) -> float:
    """Project aggregate WSU channel count for one MOUS.

    Supported modes:
    - ``memo_uniform_binned``
    - ``memo_uniform_exact``
    - ``distributed_binned``
    - ``distributed_exact``
    """
    if not templates:
        return float("nan")
    if mode.startswith("memo_uniform"):
        finest = min(template.velocity_resolution_kms for template in templates)
        target = finest if mode.endswith("exact") else stepped2_floor_velocity(finest)
        mean_freq = sum(template.center_freq_hz for template in templates) / len(templates)
        band = alma_band_from_freq_hz(mean_freq)
        nchan = projected_nchan_per_spw(
            requested_vel_kms=target,
            original_vel_kms=finest,
            center_freq_hz=mean_freq,
            band=band,
            milestone=milestone,
        )
        return projected_nspw * nchan

    projected_per_spw = []
    for template in templates:
        target = (
            template.velocity_resolution_kms
            if mode.endswith("exact")
            else stepped2_floor_velocity(template.velocity_resolution_kms)
        )
        projected_per_spw.append(
            projected_nchan_per_spw(
                requested_vel_kms=target,
                original_vel_kms=template.velocity_resolution_kms,
                center_freq_hz=template.center_freq_hz,
                band=template.band,
                milestone=milestone,
            )
        )
    return projected_nspw * (sum(projected_per_spw) / len(projected_per_spw))


def build_templates_from_rows(
    rows: Iterable[dict[str, str]],
    mous_field: str = "mous_uid",
    spw_field: str = "spw_id",
    freq_field: str = "center_freq_hz",
    vel_field: str = "velocity_resolution_kms",
) -> dict[str, list[MousSpwTemplate]]:
    """Aggregate potentially repeated EB/SPW rows into one template per MOUS/SPW."""
    grouped: dict[tuple[str, int], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        mous_uid = normalize_mous_uid(row[mous_field])
        spw_id = int(row[spw_field])
        grouped[(mous_uid, spw_id)].append(row)

    out: dict[str, list[MousSpwTemplate]] = defaultdict(list)
    for (mous_uid, spw_id), group in grouped.items():
        center_freq = sum(float(row[freq_field]) for row in group) / len(group)
        velres = sum(float(row[vel_field]) for row in group) / len(group)
        out[mous_uid].append(
            MousSpwTemplate(
                mous_uid=mous_uid,
                spw_id=spw_id,
                center_freq_hz=center_freq,
                velocity_resolution_kms=velres,
                band=alma_band_from_freq_hz(center_freq),
            )
        )
    for mous_uid in out:
        out[mous_uid].sort(key=lambda template: (template.center_freq_hz, template.spw_id))
    return out


def load_mous_spw_templates_csv(path: str | Path) -> dict[str, list[MousSpwTemplate]]:
    """Load a MOUS/SPW sidecar CSV for distributed WSU projection."""
    csv_path = Path(path)
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return build_templates_from_rows(rows)
