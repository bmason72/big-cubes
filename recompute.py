"""
Config-driven recompute layer for per-MOUS derived quantities.

This module is the "parallel recompute path" described in REFACTOR_PLAN.md:
it reads the existing per-MOUS ecsv's *fixed* inputs (imsize, mfssize,
ntarget, per-milestone nspw and nchan_agg, etc.) and re-derives the
scenario-variable outputs (productsize, datarate, ...) from a
``ScenarioConfig``.  The original ecsv is left untouched; downstream code
(compute_stats, tables, plots) consumes whichever QTable is handed to
it.

Scope (Stages 1-4):
- productsize for BLC + M1/M4/M5 via ``ps = 2 * ntarget * (agg_cube + mfssize * nspw)``;
  pixels_per_beam knob rescales imsize linearly.
- datarate / visrate per milestone (BLC + M1/M4/M5) with per-(array, config,
  milestone) integration-time lookup.
- datavol_{target_tot, cal, total} per milestone, with optional calibrator
  spectral-resolution correction: cal volume scales by 1/R_s and the rate
  columns are time-averaged via M = 1 - f_c + f_c/R_s.  Target volume and
  product/cube sizes are unaffected (imaging is target-only).
- mode="caps": apply per-SPW nchan cap, per-MOUS nchan cap with top-k
  preservation, and cube-size cap to the stored ``nchan_agg_*`` columns
  before downstream recompute; all derived quantities (productsize,
  datarate, visrate, datavol) pick up the reduction automatically.

Memo round-trip contract: with MEMO_CONFIG, every recomputed column
reproduces the stored DB column to within 1e-10 on 99% of rows.  A
handful of BLC rows with ntarget>1 or nonstandard historical
npol/nchan show pre-existing divergences inherited from how the DB was
built; those are documented in test_recompute.py and are not introduced
by this layer.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Dict

import astropy.units as u
import numpy as np
from astropy.table import QTable, Table

from config import MEMO_CONFIG, MitigationCaps, ScenarioConfig
from wsu_projection import load_mous_spw_templates_csv, normalize_mous_uid
from wsu_projection import project_nchan_agg_for_templates

log = logging.getLogger("wsu_pipeline.recompute")


# Constants from wsu_db.calc_datarate / calc_visrate.
_NBYTE_CROSS = 2      # cross-correlation packed to 16 bits
_NBYTE_AUTO = 4       # autocorrelation at 32 bits (the "4*Nant" term)
_NAPC = 1             # number of WVR streams

# Pixels-per-beam baseline baked into the ecsv's stored imsize / mfssize.
# If a scenario sets a different value we rescale geometrically.
_DB_PIXELS_PER_BEAM = 5

# (milestone_key, nchan_agg_col, nspw_col, productsize_col)
_PRODUCTSIZE_COLUMNS = (
    ("BLC", "blc_nchan_agg", "blc_nspw", "blc_productsize"),
    ("M1",  "wsu_nchan_agg_stepped2_initial", "wsu_nspw_initial",
     "wsu_productsize_initial_stepped2"),
    ("M4",  "wsu_nchan_agg_stepped2_ms4", "wsu_nspw_ms4",
     "wsu_productsize_ms4_stepped2"),
    ("M5",  "wsu_nchan_agg_stepped2_goal", "wsu_nspw_goal",
     "wsu_productsize_goal_stepped2"),
)


def _col_value(col, unit: u.Unit) -> np.ndarray:
    """Return ``col`` as a float64 ndarray in ``unit``."""
    if hasattr(col, "to"):
        return np.asarray(col.to(unit).value, dtype=float)
    return np.asarray(col, dtype=float)


def _imsize_mfssize(db: Table, pixels_per_beam: int):
    """Return (imsize_array, mfssize_GB) rescaled to ``pixels_per_beam``.

    Stored imsize / mfssize were computed at _DB_PIXELS_PER_BEAM = 5 px/beam.
    imsize scales linearly with px/beam (cells shrink linearly, FOV
    unchanged), so mfssize (∝ imsize^2) scales quadratically.
    """
    ratio = pixels_per_beam / _DB_PIXELS_PER_BEAM
    imsize = np.asarray(db["imsize"], dtype=float) * ratio
    mfssize_GB = _col_value(db["mfssize"], u.GB) * (ratio ** 2)
    return imsize, mfssize_GB


def recompute_productsize(db: Table, scenario: ScenarioConfig) -> Table:
    """Return a copy of ``db`` with productsize columns recomputed.

    Uses the wsu_db canonical formula:
        productsize = 2 * ntarget * (agg_cube + mfssize * nspw)
        agg_cube    = 4 * imsize^2 * nchan_agg / 1e9    [GB]

    The stored nchan_agg_* columns are used as-is.  In ``mode == "caps"``
    the orchestrator has already mutated them via ``_apply_caps_to_db``,
    so this function's formula is identical regardless of mode.
    """
    out = QTable(db, copy=True)
    imsize, mfssize_GB = _imsize_mfssize(db, scenario.mitigation.pixels_per_beam)
    ntarget = np.asarray(db["ntarget"], dtype=float)

    for key, nchan_col, nspw_col, ps_col in _PRODUCTSIZE_COLUMNS:
        nchan = _col_value(db[nchan_col], u.dimensionless_unscaled) \
            if getattr(db[nchan_col], "unit", None) is not None \
            else np.asarray(db[nchan_col], dtype=float)
        nspw = np.asarray(db[nspw_col], dtype=float)
        agg_cube_GB = 4.0 * imsize ** 2 * nchan / 1.0e9
        ps = 2.0 * ntarget * (agg_cube_GB + mfssize_GB * nspw)
        out[ps_col] = ps * u.GB
        log.debug("recomputed %s (%s rows)", ps_col, len(ps))

    return out


# ---------------------------------------------------------------------------
# Datarate / visrate per milestone (Stage 2)
# ---------------------------------------------------------------------------

# (milestone_key, nchan_agg_col, npol_col, tint_col, datarate_col, visrate_col)
# NB: nchan / npol / tint column names here are only used in
# ``mitigation.mode == "existing"``.  Caps mode (Stage 3) will replace
# the nchan lookup with a recomputed per-MOUS value.
_RATE_COLUMNS = (
    ("BLC", "blc_nchan_agg", "blc_npol", "blc_tint",
     "blc_datarate_typical", "blc_visrate_typical"),
    ("M1",  "wsu_nchan_agg_stepped2_initial", "wsu_npol_initial",
     "wsu_tint_initial",
     "wsu_datarate_initial_stepped2_initial",
     "wsu_visrate_initial_stepped2_initial"),
    ("M4",  "wsu_nchan_agg_stepped2_ms4", "wsu_npol", "wsu_tint",
     "wsu_datarate_ms4_stepped2_typical",
     "wsu_visrate_ms4_stepped2_typical"),
    ("M5",  "wsu_nchan_agg_stepped2_goal", "wsu_npol", "wsu_tint",
     "wsu_datarate_goal_stepped2_typical",
     "wsu_visrate_goal_stepped2_typical"),
)


def _nant_per_row(db: Table, scenario: ScenarioConfig, milestone: str) -> np.ndarray:
    """Nant per MOUS given the scenario's antenna flavor selection.

    Currently driven by the array column alone; when Nant becomes a
    per-MOUS variable (e.g., scenario-driven fractional availability)
    this helper absorbs the logic.
    """
    arr = np.asarray(db["array"])
    out = np.full(arr.shape, np.nan, dtype=float)
    for a in ("12m", "7m"):
        m = arr == a
        if not m.any():
            continue
        out[m] = scenario.antennas.nant(scenario.base, a, milestone)
    return out


def _tint_s_per_row(db: Table, scenario: ScenarioConfig, milestone: str,
                    db_tint_col: str) -> np.ndarray:
    """Integration time per MOUS in seconds.

    If ``scenario.tint.use_db_stored`` is True (MEMO default), reads
    directly from ``db_tint_col`` and applies ``elevation_factor``.
    Otherwise, classifies each MOUS by L80 and looks up via the scenario
    tint rule.
    """
    if scenario.tint.use_db_stored:
        return np.asarray(db[db_tint_col], dtype=float) * scenario.tint.elevation_factor

    arr = np.asarray(db["array"])
    L80 = np.asarray(db["L80"], dtype=float)
    out = np.empty(arr.shape, dtype=float)
    # classify per row — vectorization not worth it, this is ~5k rows
    for i in range(len(arr)):
        cid = scenario.classifier.classify(arr[i], L80[i])
        out[i] = scenario.tint.tint_s(arr[i], cid, milestone)
    return out


# ---------------------------------------------------------------------------
# Shared target-resolution rate + cal-correction helpers
# ---------------------------------------------------------------------------

# Per-milestone target velocity-resolution column (km/s).  Used by the
# calibrator correction to compute R_s = v_cap / velres_target per row.
_VELRES_COLS: Dict[str, str] = {
    "BLC": "blc_velres",
    "M1":  "wsu_velres_stepped2_initial",
    "M4":  "wsu_velres_stepped2",
    "M5":  "wsu_velres_stepped2",
}


# Per-milestone datavol output columns.  The identity
#     datavol_total = datavol_target_tot + datavol_cal
# holds exactly in the base DB; Stage 3 preserves that invariant.
_DATAVOL_COLUMNS = (
    ("BLC",
     "blc_datavol_typical_target_tot",
     "blc_datavol_typical_cal",
     "blc_datavol_typical_total"),
    ("M1",
     "wsu_datavol_initial_stepped2_initial_target_tot",
     "wsu_datavol_initial_stepped2_initial_cal",
     "wsu_datavol_initial_stepped2_initial_total"),
    ("M4",
     "wsu_datavol_ms4_stepped2_typical_target_tot",
     "wsu_datavol_ms4_stepped2_typical_cal",
     "wsu_datavol_ms4_stepped2_typical_total"),
    ("M5",
     "wsu_datavol_goal_stepped2_typical_target_tot",
     "wsu_datavol_goal_stepped2_typical_cal",
     "wsu_datavol_goal_stepped2_typical_total"),
)


# ---------------------------------------------------------------------------
# Mitigation caps (Stage 4)
# ---------------------------------------------------------------------------
# Pure algorithm: given a per-MOUS uniform starting nchan per SPW, apply
# the per-SPW cap, per-MOUS cap with top-k preservation, and cube-size
# cap, in that order.  Returns a per-SPW nchan vector (heterogeneous
# after the per-MOUS preservation step) and the aggregate.
#
# Preservation policy: the top k = ceil(preserve_fraction * nspw) SPWs
# keep their post-per-SPW-cap value (uniform here, since the base DB
# stores per-MOUS velres), and the remaining nspw-k SPWs split the
# leftover budget uniformly via floor division.  If preserving top-k
# already exceeds the MOUS cap, fall back to a uniform clip at
# floor(cap / nspw) and log a warning.


@dataclass(frozen=True)
class _CapResult:
    nchan_per_spw: np.ndarray
    nchan_agg: int


def _apply_caps_to_mous(nchan_per_spw_start, nspw: int,
                        imsize: float, caps: MitigationCaps) -> _CapResult:
    """Apply per-SPW / per-MOUS / cube-size caps to a single MOUS.

    ``nchan_per_spw_start`` is either an integer (uniform SPWs — the
    base DB case, with per-MOUS velres) or a length-``nspw`` array of
    per-SPW nchan values.  The per-MOUS preservation step always
    sorts by nchan descending and preserves the top k, so heterogeneous
    input is handled correctly.

    Returns a per-SPW nchan vector (length nspw) and its sum.
    """
    if nspw <= 0:
        return _CapResult(np.zeros(0, dtype=int), 0)
    if np.isscalar(nchan_per_spw_start):
        if nchan_per_spw_start <= 0:
            return _CapResult(np.zeros(nspw, dtype=int), 0)
        nchan = np.full(nspw, int(nchan_per_spw_start), dtype=int)
    else:
        nchan = np.asarray(nchan_per_spw_start, dtype=int).copy()
        if nchan.shape != (nspw,):
            raise ValueError(
                f"nchan_per_spw_start has shape {nchan.shape}, expected ({nspw},)"
            )

    # 1. Per-SPW cap.
    if caps.nchan_per_spw_cap is not None:
        nchan = np.minimum(nchan, int(caps.nchan_per_spw_cap))

    # 2. Per-MOUS cap with top-k preservation policy.
    # Preserve the k SPWs with the *highest* post-per-SPW-cap nchan —
    # the "highest-spectral-resolution" SPWs in the uniform case this
    # is ill-defined and any k suffice, but for heterogeneous input
    # (future: per-SPW velres, bright-line scenarios) this correctly
    # saves the fine SPWs and reduces the rest.
    if caps.nchan_mous_cap is not None and nchan.sum() > caps.nchan_mous_cap:
        k = max(1, math.ceil(caps.preserve_fraction * nspw))
        k = min(k, nspw)
        order = np.argsort(-nchan, kind="stable")   # descending by nchan
        top = order[:k]
        bot = order[k:]
        preserved = int(nchan[top].sum())
        budget = int(caps.nchan_mous_cap) - preserved
        if budget <= 0 or bot.size == 0:
            # Preservation already exceeds the cap: uniformly clip.
            log.warning("caps: preservation exceeds MOUS cap "
                        "(nspw=%d, k=%d, cap=%d); uniformly clipping",
                        nspw, k, caps.nchan_mous_cap)
            uniform = int(caps.nchan_mous_cap) // nspw
            nchan[:] = uniform
        else:
            nchan[bot] = budget // bot.size

    # 3. Cube-size cap (per SPW): each cube = 4 * imsize^2 * nchan / 1e9 GB.
    if caps.cubesize_cap_GB is not None:
        if imsize > 0 and math.isfinite(imsize):
            cap_bytes = float(caps.cubesize_cap_GB) * 1.0e9
            max_nchan = int(cap_bytes // (4.0 * imsize * imsize))
            if max_nchan < 0:
                max_nchan = 0
            nchan = np.minimum(nchan, max_nchan)

    return _CapResult(nchan_per_spw=nchan, nchan_agg=int(nchan.sum()))


# Per-milestone (nchan_agg column, nspw column) — same keys as _PRODUCTSIZE_COLUMNS.
_MOUS_NCHAN_COLUMNS = {
    "BLC": ("blc_nchan_agg",                  "blc_nspw"),
    "M1":  ("wsu_nchan_agg_stepped2_initial", "wsu_nspw_initial"),
    "M4":  ("wsu_nchan_agg_stepped2_ms4",     "wsu_nspw_ms4"),
    "M5":  ("wsu_nchan_agg_stepped2_goal",    "wsu_nspw_goal"),
}


def _apply_nchan_projection_to_db(db: Table, scenario: ScenarioConfig) -> Table:
    """Return a copy of ``db`` with WSU ``nchan_agg`` columns reprojected.

    The base per-MOUS DB stores one representative WSU spectral resolution
    per MOUS.  ``distributed_binned`` augments that with a MOUS/SPW sidecar
    and recomputes the milestone aggregate channel counts from the current
    within-MOUS SPW resolution mix.

    Only the aggregate ``nchan_agg_*`` columns are overwritten.  Downstream
    derived quantities use those columns directly, so productsize, datarate,
    visrate, datavol, and sysperf pick up the new projection automatically.
    Representative stored per-SPW columns (for example
    ``wsu_nchan_spw_stepped2``) are left unchanged because they no longer
    have a single unambiguous value under heterogeneous SPW projections.
    """
    proj = scenario.nchan_projection
    if proj.mode == "stored":
        return db
    if proj.mode != "distributed_binned":
        raise ValueError(f"unsupported nchan projection mode: {proj.mode}")
    if proj.sidecar_path is None:
        raise ValueError("distributed_binned requires nchan_projection.sidecar_path")
    if scenario.cal_correction.enabled:
        raise ValueError(
            "distributed_binned nchan projection is not supported with "
            "calibrator spectral-resolution correction"
        )
    if scenario.mitigation.mode == "caps":
        raise ValueError(
            "distributed_binned nchan projection is not supported with "
            "mitigation.mode='caps'"
        )
    if proj.mous_id_column not in db.colnames:
        raise KeyError(
            f"MOUS id column {proj.mous_id_column!r} not present in database"
        )

    templates_by_mous = load_mous_spw_templates_csv(proj.sidecar_path)
    out = QTable(db, copy=True)
    mous_values = np.asarray(db[proj.mous_id_column])
    covered = 0

    for ms, (nchan_col, nspw_col) in _MOUS_NCHAN_COLUMNS.items():
        if ms == "BLC":
            continue
        nspw_arr = _col_value(db[nspw_col], u.dimensionless_unscaled) \
            if getattr(db[nspw_col], "unit", None) is not None \
            else np.asarray(db[nspw_col], dtype=float)
        new_agg = _col_value(db[nchan_col], u.dimensionless_unscaled) \
            if getattr(db[nchan_col], "unit", None) is not None \
            else np.asarray(db[nchan_col], dtype=float)
        new_agg = np.asarray(new_agg, dtype=float).copy()

        changed = 0
        for i, mous_value in enumerate(mous_values):
            mous_uid = normalize_mous_uid(mous_value)
            templates = templates_by_mous.get(mous_uid)
            if not templates or not np.isfinite(nspw_arr[i]) or nspw_arr[i] <= 0:
                continue
            new_agg[i] = project_nchan_agg_for_templates(
                mode="distributed_binned",
                milestone=ms,
                templates=templates,
                projected_nspw=float(nspw_arr[i]),
            )
            changed += 1

        if getattr(out[nchan_col], "unit", None) is not None:
            out[nchan_col] = new_agg * u.dimensionless_unscaled
        else:
            out[nchan_col] = new_agg
        covered = max(covered, changed)
        log.debug("projected %s via distributed_binned (%s rows covered)", nchan_col, changed)

    uncovered = len(db) - covered
    if uncovered > 0:
        log.warning(
            "distributed_binned sidecar did not cover %d/%d rows; "
            "uncovered rows kept stored WSU nchan_agg values",
            uncovered, len(db),
        )
    return out


def _apply_caps_to_db(db: Table, scenario: ScenarioConfig) -> Table:
    """Return a copy of ``db`` with nchan_agg_{milestone} re-capped.

    Each milestone's stored ``nchan_agg`` is interpreted as
    ``nchan_per_spw * nspw`` (exact for WSU; floor-approximate for the
    ~8% of BLC rows with heterogeneous historical SPWs).  We derive a
    uniform per-SPW starting nchan and apply the scenario caps per row.

    Downstream recompute functions (productsize, datarate, visrate,
    datavol) consume the mutated ``nchan_agg`` columns, so the cap
    propagates automatically to every derived quantity.
    """
    caps = scenario.mitigation
    if caps.mode != "caps":
        return db
    if (caps.nchan_per_spw_cap is None
            and caps.nchan_mous_cap is None
            and caps.cubesize_cap_GB is None):
        return db  # caps mode with no caps set is a pure no-op

    imsize_scaled, _ = _imsize_mfssize(db, caps.pixels_per_beam)
    out = QTable(db, copy=True)
    n = len(db)
    for ms, (nchan_col, nspw_col) in _MOUS_NCHAN_COLUMNS.items():
        nchan_agg_prev = _col_value(db[nchan_col], u.dimensionless_unscaled) \
            if getattr(db[nchan_col], "unit", None) is not None \
            else np.asarray(db[nchan_col], dtype=float)
        nspw_arr = _col_value(db[nspw_col], u.dimensionless_unscaled) \
            if getattr(db[nspw_col], "unit", None) is not None \
            else np.asarray(db[nspw_col], dtype=float)

        new_agg = np.asarray(nchan_agg_prev, dtype=float).copy()
        for i in range(n):
            ns = int(nspw_arr[i]) if np.isfinite(nspw_arr[i]) else 0
            if ns <= 0 or not np.isfinite(nchan_agg_prev[i]):
                continue
            nchan_start = int(nchan_agg_prev[i] // ns)
            res = _apply_caps_to_mous(nchan_start, ns,
                                      float(imsize_scaled[i]), caps)
            new_agg[i] = float(res.nchan_agg)

        # Preserve the original column's unit flavor (dimensionless vs plain).
        if getattr(out[nchan_col], "unit", None) is not None:
            out[nchan_col] = new_agg * u.dimensionless_unscaled
        else:
            out[nchan_col] = new_agg
        log.debug("capped %s (%s rows changed)", nchan_col,
                  int(np.sum(new_agg != nchan_agg_prev)))
    return out


def _target_rates_per_row(db: Table, scenario: ScenarioConfig, milestone: str,
                          nchan_col: str, npol_col: str, tint_col: str):
    """Return (datarate_GBps, visrate_Gvis_per_hr) at *target* resolution.

    These are the peak rates before the cal-correction multiplier is
    applied; datavol uses these directly (times × rate), and rate
    columns get multiplied by M = 1 - f_c + f_c/R_s downstream.
    """
    nant = _nant_per_row(db, scenario, milestone)
    nbase = nant * (nant - 1.0) / 2.0
    nchan = np.asarray(db[nchan_col], dtype=float)
    npol = np.asarray(db[npol_col], dtype=float)
    tint_s = _tint_s_per_row(db, scenario, milestone, tint_col)
    dr = ((2.0 * _NBYTE_CROSS * _NAPC * nbase + _NBYTE_AUTO * nant)
          * nchan * npol / tint_s / 1.0e9)
    vr = (2.0 * npol * nbase * nchan / 1.0e9) / (tint_s / 3600.0)
    return dr, vr


def _R_per_row(db: Table, scenario: ScenarioConfig, milestone: str) -> np.ndarray:
    """R_s = v_cap / velres_target per row, clipped to ≥ 1.

    Returns an all-ones array when cal correction is disabled or the cap
    callable returns None.  A cap that is coarser-or-equal than the
    target resolution means no reduction (R = 1).
    """
    n = len(db)
    if not scenario.cal_correction.enabled:
        return np.ones(n)
    velres = _col_value(db[_VELRES_COLS[milestone]], u.km / u.s)
    arr = np.asarray(db["array"])
    L80 = np.asarray(db["L80"], dtype=float)
    R = np.ones(n)
    for i in range(n):
        cfg = scenario.classifier.classify(arr[i], L80[i])
        v_cap = scenario.cal_correction.cap_kms(arr[i], cfg, milestone, velres[i])
        if v_cap is None or v_cap <= velres[i] or velres[i] <= 0:
            continue
        R[i] = v_cap / velres[i]
    return R


def _cal_multiplier(db: Table, scenario: ScenarioConfig,
                    milestone: str) -> np.ndarray:
    """Time-weighted average multiplier M = 1 - f_c + f_c / R_s per row."""
    n = len(db)
    if not scenario.cal_correction.enabled:
        return np.ones(n)
    R = _R_per_row(db, scenario, milestone)
    target_time = _col_value(db["target_time_tot"], u.s)
    cal_time = _col_value(db["cal_time"], u.s)
    total = target_time + cal_time
    f_c = np.where(total > 0, cal_time / total, 0.0)
    return 1.0 - f_c + f_c / R


def recompute_datarate(db: Table, scenario: ScenarioConfig) -> Table:
    """Return a copy of ``db`` with datarate columns recomputed.

    Formula (wsu_db.calc_datarate):
        datarate = (2 * Nbyte_cross * Napc * Nbase + Nbyte_auto * Nant)
                   * Nchannels * Npols / Tintegration / 1e9     [GB/s]

    With ``mitigation.mode == "existing"``:
      - Nchannels <- stored nchan_agg_* column for the milestone
      - Npols     <- stored npol column for the milestone (heterogeneous
                     for BLC: {1, 2, 4})
      - Tintegration <- stored tint column (or scenario rule, if
                     use_db_stored=False)
      - Nant      <- scenario.antennas (default: typical = 47/10)

    If ``cal_correction.enabled``, the stored rate is the time-weighted
    average across target + reduced-resolution cal: R_target * M.
    """
    out = QTable(db, copy=True)
    for ms, nchan_col, npol_col, tint_col, dr_col, _ in _RATE_COLUMNS:
        dr, _ = _target_rates_per_row(db, scenario, ms,
                                      nchan_col, npol_col, tint_col)
        M = _cal_multiplier(db, scenario, ms)
        out[dr_col] = (dr * M) * (u.GB / u.s)
        log.debug("recomputed %s", dr_col)
    return out


def recompute_visrate(db: Table, scenario: ScenarioConfig) -> Table:
    """Return a copy of ``db`` with visrate columns recomputed.

    Formula (wsu_db.calc_visrate):
        visrate = (2 * Npols * Nbase * Nchannels / 1e9) / Tintegration_hr
                                                                [Gvis/hr]

    Cal-correction semantics match ``recompute_datarate`` (time-averaged
    via M = 1 - f_c + f_c/R_s).
    """
    try:
        import wsu_db  # noqa: F401 -- registers Gvis
        gvis = u.Unit("Gvis")
    except Exception:
        try:
            gvis = u.def_unit("Gvis")
        except Exception:
            gvis = u.Unit("")

    out = QTable(db, copy=True)
    for ms, nchan_col, npol_col, tint_col, _, vr_col in _RATE_COLUMNS:
        _, vr = _target_rates_per_row(db, scenario, ms,
                                      nchan_col, npol_col, tint_col)
        M = _cal_multiplier(db, scenario, ms)
        out[vr_col] = (vr * M) * gvis / u.hr
        log.debug("recomputed %s", vr_col)
    return out


def recompute_datavol(db: Table, scenario: ScenarioConfig) -> Table:
    """Return a copy of ``db`` with datavol_{target_tot,cal,total} recomputed.

    Formulas (per milestone m):
        target_tot = R_target_m * target_time_tot           [GB]
        cal        = R_target_m * cal_time / R_s            [GB]
        total      = target_tot + cal                        [GB]

    where R_target_m is the target-resolution datarate (same formula as
    ``recompute_datarate`` but without the cal-correction M multiplier)
    and R_s per row comes from the scenario's calibrator-correction cap
    (R_s = 1 when correction is disabled or v_cap >= target velres).

    The base DB sets R_s = 1 for every row, so MEMO_CONFIG reproduces
    the stored columns exactly.  Scenarios that enable cal correction
    shrink the cal contribution without touching target_tot or
    productsize.
    """
    rate_specs = {r[0]: r for r in _RATE_COLUMNS}
    target_time = _col_value(db["target_time_tot"], u.s)
    cal_time = _col_value(db["cal_time"], u.s)

    out = QTable(db, copy=True)
    for ms, tgt_col, cal_col, total_col in _DATAVOL_COLUMNS:
        _, nchan_col, npol_col, tint_col, _, _ = rate_specs[ms]
        dr, _ = _target_rates_per_row(db, scenario, ms,
                                      nchan_col, npol_col, tint_col)
        R = _R_per_row(db, scenario, ms)
        dv_target = dr * target_time
        dv_cal = dr * cal_time / R
        dv_total = dv_target + dv_cal
        out[tgt_col] = dv_target * u.GB
        out[cal_col] = dv_cal * u.GB
        out[total_col] = dv_total * u.GB
        log.debug("recomputed %s / %s / %s", tgt_col, cal_col, total_col)
    return out


_SYSPERF_COLUMNS = (
    # (visrate_col, sysperf_col, mode) where mode ∈ {"allgrid", "aprojonly"}
    # _allgrid: mosaic_aproject=True, wproject=True.
    # _aprojonly: mosaic_aproject=True, wproject=False.
    ("blc_visrate_typical", "blc_sysperf_typical_allgrid", "allgrid"),
    ("wsu_visrate_early_stepped2_typical",
     "wsu_sysperf_early_stepped2_typical_allgrid", "allgrid"),
    ("wsu_visrate_later_2x_stepped2_typical",
     "wsu_sysperf_later_2x_stepped2_typical_allgrid", "allgrid"),
    ("wsu_visrate_later_4x_stepped2_typical",
     "wsu_sysperf_later_4x_stepped2_typical_allgrid", "allgrid"),
    ("blc_visrate_typical", "blc_sysperf_typical_aprojonly", "aprojonly"),
    ("wsu_visrate_early_stepped2_typical",
     "wsu_sysperf_early_stepped2_typical_aprojonly", "aprojonly"),
    ("wsu_visrate_later_2x_stepped2_typical",
     "wsu_sysperf_later_2x_stepped2_typical_aprojonly", "aprojonly"),
    ("wsu_visrate_later_4x_stepped2_typical",
     "wsu_sysperf_later_4x_stepped2_typical_aprojonly", "aprojonly"),
    ("wsu_visrate_initial_stepped2_initial",
     "wsu_sysperf_initial_stepped2_initial_aprojonly", "aprojonly"),
    ("wsu_visrate_ms4_stepped2_typical",
     "wsu_sysperf_ms4_stepped2_typical_aprojonly", "aprojonly"),
    ("wsu_visrate_goal_stepped2_typical",
     "wsu_sysperf_goal_stepped2_typical_aprojonly", "aprojonly"),
)


def _flops_per_vis(db: Table, scenario: ScenarioConfig, mode: str
                   ) -> np.ndarray:
    """Per-row flops/vis coefficient for the given gridding ``mode``.

    ``mode`` ∈ {"allgrid", "aprojonly"}.  Both enable aproject on mosaics;
    "allgrid" additionally enables wproject for Band 1 long baselines
    (and awproject when mosaic-and-long-baseline in Band 1).
    """
    sp = scenario.base.sysperf
    n = len(db)
    flops = np.full(n, sp.flops_per_vis_std, dtype=float)

    mosaic = np.asarray(db["mosaic"]) == "T"
    band = np.asarray(db["band"])
    L80 = _col_value(db["L80"], u.m)
    long_baseline = (band == 1) & (L80 > sp.wproject_L80_threshold_m)

    flops[mosaic] = sp.flops_per_vis_aproj
    if mode == "allgrid":
        flops[long_baseline] = sp.flops_per_vis_wproj
        flops[mosaic & long_baseline] = sp.flops_per_vis_awproj
    return flops


def recompute_sysperf(db: Table, scenario: ScenarioConfig) -> Table:
    """Return a copy of ``db`` with sysperf columns recomputed.

    Formula (wsu_db.calc_sysperf):
        sysperf = visrate_Gvis_hr * 1e9 / 3600
                  * flops_per_vis * k_major_cycles * multiscale_factor
                  / (core_efficiency * parallelization_efficiency)
                  / 1e15                                        [PFLOP/s]

    flops_per_vis is per-row and depends on the gridding mode — see
    ``_flops_per_vis``.  The ``_allgrid`` label enables both aproject
    (mosaics) and wproject (Band 1 long baseline); ``_aprojonly`` enables
    only aproject.  nterms_factor is implicitly 1 (cube imaging).
    """
    sp = scenario.base.sysperf
    coeff = (1e9 / 3600.0) * sp.k_major_cycles * sp.multiscale_factor / (
        sp.core_efficiency * sp.parallelization_efficiency) / 1e15

    out = QTable(db, copy=True)
    flops_cache: Dict[str, np.ndarray] = {}
    for visrate_col, sysperf_col, mode in _SYSPERF_COLUMNS:
        if visrate_col not in db.colnames:
            continue
        if mode not in flops_cache:
            flops_cache[mode] = _flops_per_vis(db, scenario, mode)
        flops = flops_cache[mode]

        col = db[visrate_col]
        visrate = (col.value if hasattr(col, "value")
                   else np.asarray(col, dtype=float))
        sysperf = visrate.astype(float) * flops * coeff
        out[sysperf_col] = sysperf
        log.debug("recomputed %s (mode=%s)", sysperf_col, mode)

    # The stored flops_per_vis_* columns in the DB encode the mode
    # selection used to build it; expose them here so caps or scenario
    # overrides can be compared against the baseline.
    out["flops_per_vis_allgrid"] = flops_cache.get(
        "allgrid", _flops_per_vis(db, scenario, "allgrid"))
    out["flops_per_vis_aprojonly"] = flops_cache.get(
        "aprojonly", _flops_per_vis(db, scenario, "aprojonly"))
    return out


def recompute_db(db: Table, scenario: ScenarioConfig = MEMO_CONFIG) -> Table:
    """Top-level entry point: apply all implemented recomputations.

    Stages compose left-to-right; each function takes a Table and
    returns a new QTable with its target columns overwritten.

    With ``mitigation.mode == "caps"`` the nchan_agg_* columns are
    mutated first (per-SPW / per-MOUS / cube-size caps); downstream
    stages then consume the capped values and propagate the reduction
    to productsize, datarate, visrate, datavol, and sysperf uniformly.
    """
    out = _apply_nchan_projection_to_db(db, scenario)
    out = _apply_caps_to_db(out, scenario)
    out = recompute_productsize(out, scenario)
    out = recompute_datarate(out, scenario)
    out = recompute_visrate(out, scenario)
    out = recompute_datavol(out, scenario)
    out = recompute_sysperf(out, scenario)
    return out
