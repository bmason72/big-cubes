"""
Tests for the config-driven recompute layer.

Stage 1 coverage (productsize only):
- MEMO_CONFIG round-trip matches the stored DB productsize columns.
- Scenario pixels_per_beam != 5 rescales productsize correctly.

Known divergence: stored ``blc_productsize`` is an aggregate of
per-target rows each computed with its own imsize, then summed to the
MOUS level.  Our MOUS-level recompute uses the single ``imsize`` column
as representative for all targets, so multi-target MOUSes diverge
slightly (~2–5% typical, up to ~300% on 5 outliers).  WSU productsize
columns were computed at the MOUS level directly in wsu_db.py and match
to machine precision.  See the comment on
``WSU_size_of_computing_numbers_database.ipynb`` cell 1664 for the
original per-target aggregation decision.
"""

from __future__ import annotations

import os
import sys
import types
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

# Stub ipdb (wsu_db imports it unconditionally; not in venv).
sys.modules.setdefault("ipdb", types.SimpleNamespace(set_trace=lambda: None))

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

DB_PATH = Path(os.environ.get(
    "WSU_FULL_DB",
    HERE / "data" / "wsu_datarates_mit_per_mous_initial_goal_20250423.ecsv",
))


@pytest.mark.skipif(not DB_PATH.exists(),
                    reason=f"full DB not at {DB_PATH}")
class TestProductsizeMemoRoundTrip:
    """With MEMO_CONFIG, recompute must reproduce the DB productsize columns."""

    @pytest.fixture(scope="class")
    def recomputed(self):
        from config import MEMO_CONFIG
        from pipeline import load_database
        from recompute import recompute_db
        db = load_database(str(DB_PATH))
        return db, recompute_db(db, MEMO_CONFIG)

    @pytest.mark.parametrize("col", [
        "wsu_productsize_initial_stepped2",
        "wsu_productsize_ms4_stepped2",
        "wsu_productsize_goal_stepped2",
    ])
    def test_wsu_productsize_exact(self, recomputed, col):
        """WSU productsize columns were derived at MOUS level, so the
        formula reproduction is exact."""
        db, db2 = recomputed
        a = db[col].to("GB").value
        b = db2[col].to("GB").value
        # Use finite subset; NaN-NaN is not a mismatch, NaN-value is.
        fin = np.isfinite(a) & np.isfinite(b)
        rel = np.abs(a[fin] - b[fin]) / np.maximum(np.abs(a[fin]), 1e-9)
        assert rel.max() < 1e-10, (
            f"{col}: max rel diff = {rel.max():.3e} "
            f"(expected exact match)"
        )

    def test_blc_productsize_single_target(self, recomputed):
        """For single-target MOUSes (ntarget == 1) the per-target
        aggregation coincides with the MOUS-level formula, so BLC
        productsize must match exactly."""
        db, db2 = recomputed
        a = db["blc_productsize"].to("GB").value
        b = db2["blc_productsize"].to("GB").value
        ntarget = np.asarray(db["ntarget"])
        m = (ntarget == 1) & np.isfinite(a) & np.isfinite(b)
        rel = np.abs(a[m] - b[m]) / np.maximum(np.abs(a[m]), 1e-9)
        assert rel.max() < 1e-10, (
            f"BLC productsize, ntarget==1: max rel diff = {rel.max():.3e}"
        )

    def test_blc_productsize_multi_target_known_divergence(self, recomputed):
        """For multi-target MOUSes, stored values were built by summing
        per-target productsize where each target had its own imsize,
        so MOUS-level recompute diverges.  We pin the *envelope* here
        rather than per-row agreement — the recompute is
        *formula-consistent*, the DB is *build-history-consistent*,
        and this layer intentionally chooses the former.

        Current envelope (2026-04-21, 1056 multi-target rows):
            median = 1e-16, p75 = 2%, p90 = 7%, p95 = 15%, max = 294%.
        Sum-level agreement is much tighter (see sum-level assertion).
        """
        db, db2 = recomputed
        a = db["blc_productsize"].to("GB").value
        b = db2["blc_productsize"].to("GB").value
        ntarget = np.asarray(db["ntarget"])
        m = (ntarget > 1) & np.isfinite(a) & np.isfinite(b)
        rel = np.abs(a[m] - b[m]) / np.maximum(np.abs(a[m]), 1e-9)
        # Typical row should agree tightly.
        assert np.median(rel) < 1e-3
        assert np.percentile(rel, 90) < 0.12
        # Sum across multi-target MOUSes: systematic bias should be
        # small.  If the recompute formula changes, this catches it.
        sum_ratio = b[m].sum() / a[m].sum()
        assert 0.98 < sum_ratio < 1.05, (
            f"BLC productsize sum ratio (recomputed/stored) = "
            f"{sum_ratio:.4f}; expected ~1.0"
        )


@pytest.mark.skipif(not DB_PATH.exists(),
                    reason=f"full DB not at {DB_PATH}")
class TestPixelsPerBeamScaling:
    """A scenario with pixels_per_beam != 5 must rescale productsize by
    (ratio)^2 on the cube component and (ratio)^2 on the mfs
    component — so the whole productsize scales as (ratio)^2 exactly
    (both terms are ∝ imsize^2)."""

    @pytest.fixture(scope="class")
    def scenario(self):
        from dataclasses import replace
        from config import MEMO_CONFIG
        return replace(MEMO_CONFIG,
                       name="memo_3px",
                       mitigation=replace(MEMO_CONFIG.mitigation,
                                          pixels_per_beam=3))

    def test_m5_scales_quadratically(self, scenario):
        from config import MEMO_CONFIG
        from pipeline import load_database
        from recompute import recompute_db
        db = load_database(str(DB_PATH))
        memo = recompute_db(db, MEMO_CONFIG)
        scaled = recompute_db(db, scenario)
        col = "wsu_productsize_goal_stepped2"
        a = memo[col].to("GB").value
        b = scaled[col].to("GB").value
        m = np.isfinite(a) & (a > 0) & np.isfinite(b)
        ratios = b[m] / a[m]
        expected = (3 / 5) ** 2
        # Every row should hit the expected ratio exactly (same formula,
        # only imsize is rescaled).
        assert np.allclose(ratios, expected, rtol=1e-10), (
            f"mean ratio = {ratios.mean():.6f}, expected {expected:.6f}"
        )


class TestConfigClassifier:
    """Pure-function tests for the Cy13 configuration classifier."""

    def test_every_config_self_classifies(self):
        """Each nominal L80 must classify to its own configuration label."""
        from config import ConfigClassifier
        from data.cycle13_configurations import CONFIG_LABELS_12M, L80_VALUES_12M
        c = ConfigClassifier()
        for label, l80 in zip(CONFIG_LABELS_12M, L80_VALUES_12M):
            assert c.classify("12m", l80) == label

    def test_7m_bypasses_classifier(self):
        from config import ConfigClassifier
        c = ConfigClassifier()
        # 7m ignores L80 -- always '7-m'
        for l80 in (10.0, 30.7, 100.0, 1000.0):
            assert c.classify("7m", l80) == "7-m"

    @pytest.mark.parametrize("l80, expected", [
        (100.0, "C43-1"),     # nearest to 107.1 in log space
        (500.0, "C43-5"),     # log-nearest (rather than C43-4 at 369)
        (1000.0, "C43-6"),    # nearest to 1172.5
        (3000.0, "C43-8"),    # nearest to 3527.3
        (10_000.0, "C43-10"),
    ])
    def test_log_nearest(self, l80, expected):
        from config import ConfigClassifier
        assert ConfigClassifier().classify("12m", l80) == expected

    def test_degenerate_L80_falls_back(self):
        """Non-finite or non-positive L80 must not crash the classifier."""
        import math
        from config import ConfigClassifier
        c = ConfigClassifier()
        assert c.classify("12m", math.nan) == "C43-1"
        assert c.classify("12m", 0.0) == "C43-1"
        assert c.classify("12m", -1.0) == "C43-1"


class TestTintLookup:
    """Default rule + DB-stored + override paths."""

    def test_default_rule_matches_user_directive(self):
        """C43-8+ gets the short tint (~3s); earlier configs the long (~6s)."""
        from config import default_tint_s
        # BLC
        assert default_tint_s("12m", "C43-1", "BLC") == 6.05
        assert default_tint_s("12m", "C43-7", "BLC") == 6.05
        assert default_tint_s("12m", "C43-8", "BLC") == 3.024
        assert default_tint_s("12m", "C43-10", "BLC") == 3.024
        # M1 initial
        assert default_tint_s("12m", "C43-1", "M1") == 6.144
        assert default_tint_s("12m", "C43-8", "M1") == 3.072
        # M4/M5 uniform
        assert default_tint_s("12m", "C43-1", "M4") == 3.072
        assert default_tint_s("12m", "C43-10", "M5") == 3.072
        # 7m
        assert default_tint_s("7m", "7-m", "BLC") == 10.1
        assert default_tint_s("7m", "7-m", "M5") == 9.984

    def test_db_stored_path(self):
        """use_db_stored=True returns the DB tint scaled by elevation_factor."""
        from config import TintLookup
        t = TintLookup(use_db_stored=True, elevation_factor=1.5)
        assert t.tint_s("12m", "C43-3", "BLC", db_tint_s=6.05) == 1.5 * 6.05

    def test_override(self):
        """overrides dict wins over the default rule."""
        from config import TintLookup
        t = TintLookup(use_db_stored=False,
                       overrides={("12m", "C43-1", "BLC"): 42.0})
        assert t.tint_s("12m", "C43-1", "BLC") == 42.0
        # Uncovered keys fall back to default_tint_s.
        assert t.tint_s("12m", "C43-8", "BLC") == 3.024


# ---------------------------------------------------------------------------
# Datarate / visrate memo round-trip (Stage 2)
# ---------------------------------------------------------------------------
# Historical envelope for BLC columns: ~46 of 5192 rows diverge due to
# nonstandard Npol/Nchan used in the archive-era compute path (factors
# of 2/3/4/5 in the ratio).  Bound here by (a) exact match for WSU and
# (b) a permissive p99 on BLC.

@pytest.mark.skipif(not DB_PATH.exists(),
                    reason=f"full DB not at {DB_PATH}")
class TestDatarateMemoRoundTrip:

    @pytest.fixture(scope="class")
    def recomputed(self):
        from config import MEMO_CONFIG
        from pipeline import load_database
        from recompute import recompute_db
        db = load_database(str(DB_PATH))
        return db, recompute_db(db, MEMO_CONFIG)

    @pytest.mark.parametrize("col", [
        "wsu_datarate_initial_stepped2_initial",
        "wsu_datarate_ms4_stepped2_typical",
        "wsu_datarate_goal_stepped2_typical",
    ])
    def test_wsu_exact(self, recomputed, col):
        db, db2 = recomputed
        a = db[col].to("GB/s").value
        b = db2[col].to("GB/s").value
        fin = np.isfinite(a) & np.isfinite(b) & (a > 0)
        rel = np.abs(a[fin] - b[fin]) / a[fin]
        assert rel.max() < 1e-10, f"{col}: max rel diff = {rel.max():.3e}"

    def test_blc_envelope(self, recomputed):
        db, db2 = recomputed
        a = db["blc_datarate_typical"].to("GB/s").value
        b = db2["blc_datarate_typical"].to("GB/s").value
        fin = np.isfinite(a) & np.isfinite(b) & (a > 0)
        rel = np.abs(a[fin] - b[fin]) / a[fin]
        assert np.median(rel) < 1e-10
        # Historical: ~46 outlier rows (0.9%); keep headroom.
        outliers = int(np.sum(rel > 1e-6))
        assert outliers < 80, (
            f"BLC datarate: {outliers} rows > 1e-6 (historical ≈ 46)"
        )


@pytest.mark.skipif(not DB_PATH.exists(),
                    reason=f"full DB not at {DB_PATH}")
class TestVisrateMemoRoundTrip:

    @pytest.fixture(scope="class")
    def recomputed(self):
        from config import MEMO_CONFIG
        from pipeline import load_database
        from recompute import recompute_db
        db = load_database(str(DB_PATH))
        return db, recompute_db(db, MEMO_CONFIG)

    @pytest.mark.parametrize("col", [
        "wsu_visrate_initial_stepped2_initial",
        "wsu_visrate_ms4_stepped2_typical",
        "wsu_visrate_goal_stepped2_typical",
    ])
    def test_wsu_exact(self, recomputed, col):
        db, db2 = recomputed
        a = np.asarray(db[col].value, dtype=float)
        b = np.asarray(db2[col].value, dtype=float)
        fin = np.isfinite(a) & np.isfinite(b) & (a > 0)
        rel = np.abs(a[fin] - b[fin]) / a[fin]
        assert rel.max() < 1e-10, f"{col}: max rel diff = {rel.max():.3e}"

    def test_blc_envelope(self, recomputed):
        db, db2 = recomputed
        a = np.asarray(db["blc_visrate_typical"].value, dtype=float)
        b = np.asarray(db2["blc_visrate_typical"].value, dtype=float)
        fin = np.isfinite(a) & np.isfinite(b) & (a > 0)
        rel = np.abs(a[fin] - b[fin]) / a[fin]
        assert np.median(rel) < 1e-10
        outliers = int(np.sum(rel > 1e-6))
        assert outliers < 80


class TestRecomputeOther:
    """Non-DB-dependent tests."""

    def test_mode_caps_with_no_caps_is_noop(self):
        """mode='caps' with every cap None is a pure no-op — it must
        return the input unchanged (identity check via the _apply_caps_to_db
        helper); downstream recompute then behaves like mode='existing'."""
        from dataclasses import replace
        from astropy.table import QTable
        from config import MEMO_CONFIG, MitigationCaps
        from recompute import _apply_caps_to_db
        sc = replace(MEMO_CONFIG,
                     mitigation=MitigationCaps(mode="caps", pixels_per_beam=5))
        empty = QTable()
        assert _apply_caps_to_db(empty, sc) is empty


# ---------------------------------------------------------------------------
# Datavol memo round-trip + calibrator correction (Stage 3)
# ---------------------------------------------------------------------------
# The DB computes datavol = datarate * {target_time_tot, cal_time, time_tot}
# at target resolution, so:
#   - WSU datavol columns must round-trip exactly with MEMO_CONFIG.
#   - BLC datavol propagates the documented ~46-row datarate divergence
#     (factor-of-N Npol/Nchan mismatches); same envelope applies.

@pytest.mark.skipif(not DB_PATH.exists(),
                    reason=f"full DB not at {DB_PATH}")
class TestDatavolMemoRoundTrip:

    @pytest.fixture(scope="class")
    def recomputed(self):
        from config import MEMO_CONFIG
        from pipeline import load_database
        from recompute import recompute_db
        db = load_database(str(DB_PATH))
        return db, recompute_db(db, MEMO_CONFIG)

    @pytest.mark.parametrize("col", [
        "wsu_datavol_initial_stepped2_initial_target_tot",
        "wsu_datavol_initial_stepped2_initial_cal",
        "wsu_datavol_initial_stepped2_initial_total",
        "wsu_datavol_ms4_stepped2_typical_target_tot",
        "wsu_datavol_ms4_stepped2_typical_cal",
        "wsu_datavol_ms4_stepped2_typical_total",
        "wsu_datavol_goal_stepped2_typical_target_tot",
        "wsu_datavol_goal_stepped2_typical_cal",
        "wsu_datavol_goal_stepped2_typical_total",
    ])
    def test_wsu_exact(self, recomputed, col):
        db, db2 = recomputed
        a = db[col].to("GB").value
        b = db2[col].to("GB").value
        fin = np.isfinite(a) & np.isfinite(b) & (a > 0)
        rel = np.abs(a[fin] - b[fin]) / a[fin]
        assert rel.max() < 1e-10, f"{col}: max rel diff = {rel.max():.3e}"

    @pytest.mark.parametrize("col", [
        "blc_datavol_typical_target_tot",
        "blc_datavol_typical_cal",
        "blc_datavol_typical_total",
    ])
    def test_blc_envelope(self, recomputed, col):
        db, db2 = recomputed
        a = db[col].to("GB").value
        b = db2[col].to("GB").value
        fin = np.isfinite(a) & np.isfinite(b) & (a > 0)
        rel = np.abs(a[fin] - b[fin]) / a[fin]
        assert np.median(rel) < 1e-10
        outliers = int(np.sum(rel > 1e-6))
        assert outliers < 80, (
            f"{col}: {outliers} rows > 1e-6 (historical ≈ 46)"
        )

    def test_identity_target_plus_cal_equals_total(self, recomputed):
        """Recompute must preserve the DB-level identity exactly."""
        _, db2 = recomputed
        for tgt, cal, tot in (
            ("wsu_datavol_ms4_stepped2_typical_target_tot",
             "wsu_datavol_ms4_stepped2_typical_cal",
             "wsu_datavol_ms4_stepped2_typical_total"),
            ("wsu_datavol_goal_stepped2_typical_target_tot",
             "wsu_datavol_goal_stepped2_typical_cal",
             "wsu_datavol_goal_stepped2_typical_total"),
        ):
            t = db2[tgt].to("GB").value
            c = db2[cal].to("GB").value
            s = db2[tot].to("GB").value
            with np.errstate(invalid="ignore", divide="ignore"):
                rel = np.where(s != 0, np.abs((t + c) - s) / np.abs(s), 0)
            assert np.nanmax(rel) < 1e-10


class TestCalibratorCorrectionMath:
    """Hand-calcs for the M = 1 - f_c + f_c/R_s multiplier.

    The config-level ``CalibratorCorrection.cap_kms`` handles None /
    scalar / callable, so these tests pin the contract independently of
    the DB-driven recompute path.
    """

    def test_disabled_returns_none(self):
        from config import CalibratorCorrection
        cc = CalibratorCorrection(enabled=False, v_cap_kms=1.0)
        assert cc.cap_kms("12m", "C43-1", "M5", 0.1) is None

    def test_scalar_cap(self):
        from config import CalibratorCorrection
        cc = CalibratorCorrection(enabled=True, v_cap_kms=1.0)
        assert cc.cap_kms("12m", "C43-1", "M5", 0.1) == 1.0

    def test_callable_cap(self):
        from config import CalibratorCorrection
        # Cap only applies for M5; disable elsewhere.
        fn = lambda arr, cfg, ms, v: 1.0 if ms == "M5" else None
        cc = CalibratorCorrection(enabled=True, v_cap_kms=fn)
        assert cc.cap_kms("12m", "C43-1", "M5", 0.1) == 1.0
        assert cc.cap_kms("12m", "C43-1", "M4", 0.1) is None


@pytest.mark.skipif(not DB_PATH.exists(),
                    reason=f"full DB not at {DB_PATH}")
class TestCalibratorCorrectionApplied:
    """End-to-end: cal correction scales cal by 1/R_s and rates by M,
    and leaves target volume + productsize untouched."""

    @pytest.fixture(scope="class")
    def pair(self):
        from dataclasses import replace
        from config import MEMO_CONFIG, CalibratorCorrection
        from pipeline import load_database
        from recompute import recompute_db
        db = load_database(str(DB_PATH))
        memo = recompute_db(db, MEMO_CONFIG)
        cal_scen = replace(
            MEMO_CONFIG,
            cal_correction=CalibratorCorrection(enabled=True, v_cap_kms=1.0),
        )
        cal = recompute_db(db, cal_scen)
        return db, memo, cal, cal_scen

    def test_target_tot_unchanged(self, pair):
        """Calibrator correction must not touch the target volume."""
        import astropy.units as u
        _, memo, cal, _ = pair
        for col in (
            "wsu_datavol_initial_stepped2_initial_target_tot",
            "wsu_datavol_ms4_stepped2_typical_target_tot",
            "wsu_datavol_goal_stepped2_typical_target_tot",
        ):
            a = memo[col].to(u.GB).value
            b = cal[col].to(u.GB).value
            assert np.max(np.abs(a - b)) == 0.0, col

    def test_productsize_unchanged(self, pair):
        """Imaging products are target-only; cal correction must not touch."""
        import astropy.units as u
        _, memo, cal, _ = pair
        for col in ("blc_productsize", "wsu_productsize_goal_stepped2"):
            a = memo[col].to(u.GB).value
            b = cal[col].to(u.GB).value
            assert np.max(np.abs(a - b)) == 0.0, col

    def test_cal_vol_scales_as_inverse_R(self, pair):
        """cal / memo == 1/R_s per row (M5, v_cap=1 km/s)."""
        import astropy.units as u
        from recompute import _R_per_row
        db, memo, cal, cal_scen = pair
        R = _R_per_row(db, cal_scen, "M5")
        a = memo["wsu_datavol_goal_stepped2_typical_cal"].to(u.GB).value
        b = cal["wsu_datavol_goal_stepped2_typical_cal"].to(u.GB).value
        m = a > 0
        ratio = b[m] / a[m]
        # b/a = 1/R, per row.  Round-off only.
        assert np.max(np.abs(ratio - 1.0 / R[m])) < 1e-10
        # At least some rows must have R>1 for this scenario to be non-trivial.
        assert (R[m] > 1.0).sum() > 100

    def test_rate_scaled_by_M(self, pair):
        """datarate_cal / datarate_memo == M per row."""
        import astropy.units as u
        from recompute import _cal_multiplier
        db, memo, cal, cal_scen = pair
        M = _cal_multiplier(db, cal_scen, "M5")
        a = memo["wsu_datarate_goal_stepped2_typical"].to(u.GB / u.s).value
        b = cal["wsu_datarate_goal_stepped2_typical"].to(u.GB / u.s).value
        m = a > 0
        assert np.max(np.abs(b[m] / a[m] - M[m])) < 1e-10

    def test_v_cap_below_velres_is_noop(self):
        """When v_cap <= velres_target everywhere (cap finer than target),
        cal can't be coarser than target so R = 1 and M = 1 (identity)."""
        from dataclasses import replace
        from config import MEMO_CONFIG, CalibratorCorrection
        from pipeline import load_database
        from recompute import recompute_db
        db = load_database(str(DB_PATH))
        # v_cap=0 is <= every positive velres, so nothing gets reduced.
        noop_cc = CalibratorCorrection(enabled=True, v_cap_kms=0.0)
        memo = recompute_db(db, MEMO_CONFIG)
        noop = recompute_db(db, replace(MEMO_CONFIG, cal_correction=noop_cc))
        for col in (
            "wsu_datavol_goal_stepped2_typical_cal",
            "wsu_datavol_goal_stepped2_typical_total",
            "wsu_datarate_goal_stepped2_typical",
        ):
            a = np.asarray(memo[col].value, dtype=float)
            b = np.asarray(noop[col].value, dtype=float)
            fin = np.isfinite(a) & np.isfinite(b) & (a > 0)
            rel = np.abs(a[fin] - b[fin]) / a[fin]
            assert rel.max() < 1e-10, col


# ---------------------------------------------------------------------------
# Mitigation caps (Stage 4)
# ---------------------------------------------------------------------------
# The cap algorithm is a pure function on a single MOUS; unit-test it
# independently, then verify end-to-end behavior against the DB.

class TestCapAlgorithm:
    """Pure-function tests for ``_apply_caps_to_mous``."""

    def test_no_caps_is_identity(self):
        from config import MitigationCaps
        from recompute import _apply_caps_to_mous
        caps = MitigationCaps(mode="caps")
        res = _apply_caps_to_mous(2904, 8, imsize=1000.0, caps=caps)
        assert res.nchan_agg == 2904 * 8
        assert list(res.nchan_per_spw) == [2904] * 8

    def test_per_spw_cap_clips_uniformly(self):
        from config import MitigationCaps
        from recompute import _apply_caps_to_mous
        caps = MitigationCaps(mode="caps", nchan_per_spw_cap=1000)
        res = _apply_caps_to_mous(2904, 8, imsize=1000.0, caps=caps)
        assert list(res.nchan_per_spw) == [1000] * 8
        assert res.nchan_agg == 8000

    def test_per_spw_cap_above_native_is_noop(self):
        from config import MitigationCaps
        from recompute import _apply_caps_to_mous
        caps = MitigationCaps(mode="caps", nchan_per_spw_cap=10_000)
        res = _apply_caps_to_mous(2904, 8, imsize=1000.0, caps=caps)
        assert res.nchan_agg == 2904 * 8

    def test_mous_cap_top_k_preservation(self):
        """k = ceil(0.25 * 8) = 2 SPWs preserved at 2904; remaining
        6 SPWs share budget = 8000 - 5808 = 2192 → 365 each (floor)."""
        from config import MitigationCaps
        from recompute import _apply_caps_to_mous
        caps = MitigationCaps(mode="caps", nchan_mous_cap=8000,
                              preserve_fraction=0.25)
        res = _apply_caps_to_mous(2904, 8, imsize=1000.0, caps=caps)
        assert list(res.nchan_per_spw) == [2904, 2904] + [365] * 6
        assert res.nchan_agg == 2 * 2904 + 6 * 365   # = 7998 (floor loses 2)

    def test_mous_cap_above_sum_is_noop(self):
        from config import MitigationCaps
        from recompute import _apply_caps_to_mous
        caps = MitigationCaps(mode="caps", nchan_mous_cap=100_000)
        res = _apply_caps_to_mous(2904, 8, imsize=1000.0, caps=caps)
        assert res.nchan_agg == 2904 * 8

    def test_preservation_exceeds_cap_uniform_clip(self):
        """k*nchan_per_spw > cap → degenerate, uniformly clip to cap/nspw."""
        from config import MitigationCaps
        from recompute import _apply_caps_to_mous
        caps = MitigationCaps(mode="caps", nchan_mous_cap=500,
                              preserve_fraction=0.5)
        # nspw=4, k=ceil(0.5*4)=2, preserved=2*2904=5808 > 500 → degenerate.
        res = _apply_caps_to_mous(2904, 4, imsize=1000.0, caps=caps)
        assert list(res.nchan_per_spw) == [125, 125, 125, 125]
        assert res.nchan_agg == 500

    def test_cube_size_cap(self):
        """Cube cap: each cube ≤ 5 GB with imsize=1000 → max nchan=1250."""
        from config import MitigationCaps
        from recompute import _apply_caps_to_mous
        caps = MitigationCaps(mode="caps", cubesize_cap_GB=5.0)
        res = _apply_caps_to_mous(2904, 8, imsize=1000.0, caps=caps)
        assert list(res.nchan_per_spw) == [1250] * 8
        assert res.nchan_agg == 10_000

    def test_combined_caps_apply_in_order(self):
        """per-SPW → per-MOUS (preservation) → cube-size.  Last wins for
        any SPW that exceeds its value after the previous stages."""
        from config import MitigationCaps
        from recompute import _apply_caps_to_mous
        # Per-SPW cap trims 2904 → 2000; per-MOUS cap 6000 with k=2
        # preserves 2*2000=4000, leaves budget 2000 for 6 SPWs → 333 each;
        # cube cap at 3 GB / (4*1000^2) = 750 trims the top two 2000 → 750.
        caps = MitigationCaps(mode="caps",
                              nchan_per_spw_cap=2000,
                              nchan_mous_cap=6000,
                              preserve_fraction=0.25,
                              cubesize_cap_GB=3.0)
        res = _apply_caps_to_mous(2904, 8, imsize=1000.0, caps=caps)
        assert list(res.nchan_per_spw) == [750, 750] + [333] * 6

    def test_nspw_one_preserves_everything(self):
        """nspw=1 always triggers degenerate clip when the mous cap is
        below the preserved value; no 'unpreserved' SPWs to spread to."""
        from config import MitigationCaps
        from recompute import _apply_caps_to_mous
        caps = MitigationCaps(mode="caps", nchan_mous_cap=500)
        res = _apply_caps_to_mous(2904, 1, imsize=1000.0, caps=caps)
        # k=1, preserved=2904, budget=500-2904<0 → degenerate uniform
        assert res.nchan_agg == 500

    def test_zero_start_yields_zero(self):
        from config import MitigationCaps
        from recompute import _apply_caps_to_mous
        caps = MitigationCaps(mode="caps", nchan_per_spw_cap=100)
        res = _apply_caps_to_mous(0, 4, imsize=1000.0, caps=caps)
        assert res.nchan_agg == 0

    def test_heterogeneous_input_preserves_largest_spws(self):
        """Per-MOUS preservation must pick the k *largest* SPWs (highest
        spectral resolution), not the first k by index.  The base DB is
        uniform so this degenerates to any-k, but heterogeneous input
        (future: per-SPW velres, bright-line scenarios) demands the
        top-by-value semantics from REFACTOR_PLAN.md §4.3."""
        from config import MitigationCaps
        from recompute import _apply_caps_to_mous
        # Deliberately put the large SPWs LAST so "first k by index" differs
        # from "top k by value".
        # sum=18500, cap=17600 → fires.  k = ceil(5*0.25) = 2.
        # Top-2 = [9000, 8000]; preserved = 17000; budget = 600 → 200 each.
        start = np.array([500, 500, 500, 8000, 9000], dtype=int)
        caps = MitigationCaps(mode="caps", nchan_mous_cap=17_600,
                              preserve_fraction=0.25)
        res = _apply_caps_to_mous(start, 5, imsize=1000.0, caps=caps)
        # Reduced SPWs are the three low-value ones (indices 0, 1, 2).
        assert list(res.nchan_per_spw) == [200, 200, 200, 8000, 9000]

    def test_per_spw_cap_reshuffles_preservation_order(self):
        """If the per-SPW cap flattens the originally-tallest SPWs so
        they no longer stand out, preservation now picks by the
        *post-per-SPW-cap* distribution — i.e. a stable tie-breaker on
        index is acceptable as long as no capped SPW is preferred over
        an uncapped one."""
        from config import MitigationCaps
        from recompute import _apply_caps_to_mous
        # Before per-SPW cap: [100, 9000, 200, 8000]
        # After cap=500:       [100,  500, 200,  500]
        # k = ceil(4*0.25)=1; top-1 must be one of the 500s, not 200 or 100.
        start = np.array([100, 9000, 200, 8000], dtype=int)
        caps = MitigationCaps(mode="caps",
                              nchan_per_spw_cap=500,
                              nchan_mous_cap=800,
                              preserve_fraction=0.25)
        res = _apply_caps_to_mous(start, 4, imsize=1000.0, caps=caps)
        # One 500 preserved; 3 others share budget 800-500=300 → 100 each.
        # Stable sort → preserves index 1 (first 500).
        assert sorted(res.nchan_per_spw.tolist()) == [100, 100, 100, 500]


class TestDistributedBinnedProjection:
    """Pure-function tests for sidecar-driven WSU aggregate nchan projection."""

    @staticmethod
    def _write_sidecar(path: Path) -> Path:
        path.write_text(
            "\n".join(
                [
                    "mous_uid,spw_id,center_freq_hz,velocity_resolution_kms,band",
                    "uid___A001_X1_X1,1,230000000000.0,0.25,6",
                    "uid___A001_X1_X1,2,231000000000.0,12.0,6",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        return path

    @staticmethod
    def _mini_db() -> Table:
        return Table(
            {
                "mous": ["uid://A001/X1/X1", "uid://A001/X9/X9"],
                "wsu_nchan_agg_stepped2_initial": [999.0, 888.0],
                "wsu_nspw_initial": [6.0, 6.0],
                "wsu_nchan_agg_stepped2_ms4": [999.0, 888.0],
                "wsu_nspw_ms4": [5.5, 5.5],
                "wsu_nchan_agg_stepped2_goal": [999.0, 888.0],
                "wsu_nspw_goal": [16.0, 16.0],
                "blc_nchan_agg": [1.0, 1.0],
                "blc_nspw": [1.0, 1.0],
            }
        )

    def test_sidecar_projection_rewrites_only_covered_rows(self, tmp_path):
        from dataclasses import replace
        from config import MEMO_CONFIG, NchanProjection
        from recompute import _apply_nchan_projection_to_db
        from wsu_projection import load_mous_spw_templates_csv
        from wsu_projection import project_nchan_agg_for_templates

        sidecar = self._write_sidecar(tmp_path / "sidecar.csv")
        scenario = replace(
            MEMO_CONFIG,
            nchan_projection=NchanProjection(
                mode="distributed_binned",
                sidecar_path=str(sidecar),
            ),
        )
        db = self._mini_db()
        out = _apply_nchan_projection_to_db(db, scenario)
        templates = load_mous_spw_templates_csv(sidecar)["uid___A001_X1_X1"]

        expected_m1 = project_nchan_agg_for_templates(
            "distributed_binned", "M1", templates, projected_nspw=6.0
        )
        expected_m4 = project_nchan_agg_for_templates(
            "distributed_binned", "M4", templates, projected_nspw=5.5
        )
        expected_m5 = project_nchan_agg_for_templates(
            "distributed_binned", "M5", templates, projected_nspw=16.0
        )

        assert out["wsu_nchan_agg_stepped2_initial"][0] == pytest.approx(expected_m1)
        assert out["wsu_nchan_agg_stepped2_ms4"][0] == pytest.approx(expected_m4)
        assert out["wsu_nchan_agg_stepped2_goal"][0] == pytest.approx(expected_m5)

        assert out["wsu_nchan_agg_stepped2_initial"][1] == 888.0
        assert out["wsu_nchan_agg_stepped2_ms4"][1] == 888.0
        assert out["wsu_nchan_agg_stepped2_goal"][1] == 888.0

    def test_projection_with_caps_is_rejected(self, tmp_path):
        from dataclasses import replace
        from config import MEMO_CONFIG, MitigationCaps, NchanProjection
        from recompute import recompute_db

        sidecar = self._write_sidecar(tmp_path / "sidecar.csv")
        scenario = replace(
            MEMO_CONFIG,
            mitigation=MitigationCaps(mode="caps", nchan_mous_cap=1000),
            nchan_projection=NchanProjection(
                mode="distributed_binned",
                sidecar_path=str(sidecar),
            ),
        )
        with pytest.raises(ValueError, match="not supported with mitigation.mode='caps'"):
            recompute_db(self._mini_db(), scenario)


@pytest.mark.skipif(not DB_PATH.exists(),
                    reason=f"full DB not at {DB_PATH}")
class TestCapsModeIntegration:
    """mode='caps' applied through the recompute orchestrator."""

    def test_caps_with_no_caps_matches_memo(self):
        """mode='caps' with every cap None must produce identical output
        to mode='existing'.  This is the structural safety net: the mode
        flag alone must not change anything."""
        from dataclasses import replace
        from config import MEMO_CONFIG, MitigationCaps
        from pipeline import load_database
        from recompute import recompute_db
        db = load_database(str(DB_PATH))
        memo = recompute_db(db, MEMO_CONFIG)
        caps = MitigationCaps(mode="caps", pixels_per_beam=5)
        noop = recompute_db(db, replace(MEMO_CONFIG, mitigation=caps))
        for col in (
            "wsu_productsize_goal_stepped2",
            "wsu_nchan_agg_stepped2_goal",
            "wsu_datarate_goal_stepped2_typical",
            "wsu_datavol_goal_stepped2_typical_total",
            "blc_productsize",
        ):
            a = np.asarray(memo[col].value if hasattr(memo[col], "value")
                           else memo[col], dtype=float)
            b = np.asarray(noop[col].value if hasattr(noop[col], "value")
                           else noop[col], dtype=float)
            fin = np.isfinite(a) & np.isfinite(b) & (np.abs(a) > 0)
            rel = np.abs(a[fin] - b[fin]) / np.abs(a[fin])
            assert rel.max() < 1e-12, col

    def test_mous_cap_reduces_nchan_agg(self):
        """An aggressive per-MOUS cap must reduce every row's nchan_agg
        (for WSU, where every row has nchan_agg >> cap)."""
        from dataclasses import replace
        from config import MEMO_CONFIG, MitigationCaps
        from pipeline import load_database
        from recompute import recompute_db
        db = load_database(str(DB_PATH))
        caps = MitigationCaps(mode="caps", nchan_mous_cap=1000)
        out = recompute_db(db, replace(MEMO_CONFIG, mitigation=caps))
        col = "wsu_nchan_agg_stepped2_goal"
        a = np.asarray(db[col].value if hasattr(db[col], "value")
                       else db[col], dtype=float)
        b = np.asarray(out[col].value if hasattr(out[col], "value")
                       else out[col], dtype=float)
        fin = np.isfinite(a) & np.isfinite(b) & (a > 1000)
        # Every WSU M5 row with nchan_agg > cap must have been reduced.
        assert np.all(b[fin] <= 1000 + 1e-9), (
            f"{int((b[fin] > 1000).sum())} rows exceed cap"
        )

    def test_caps_monotonically_reduce(self):
        """No derived quantity may *increase* under any cap scenario
        relative to the memo baseline."""
        from dataclasses import replace
        from config import MEMO_CONFIG, MitigationCaps
        from pipeline import load_database
        from recompute import recompute_db
        db = load_database(str(DB_PATH))
        memo = recompute_db(db, MEMO_CONFIG)
        caps = MitigationCaps(mode="caps",
                              nchan_per_spw_cap=1000,
                              cubesize_cap_GB=5.0)
        out = recompute_db(db, replace(MEMO_CONFIG, mitigation=caps))
        for col in (
            "wsu_nchan_agg_stepped2_goal",
            "wsu_productsize_goal_stepped2",
            "wsu_datarate_goal_stepped2_typical",
            "wsu_datavol_goal_stepped2_typical_total",
        ):
            a = np.asarray(memo[col].value if hasattr(memo[col], "value")
                           else memo[col], dtype=float)
            b = np.asarray(out[col].value if hasattr(out[col], "value")
                           else out[col], dtype=float)
            fin = np.isfinite(a) & np.isfinite(b)
            # Allow 1-bit float slop.
            assert np.all(b[fin] <= a[fin] * (1 + 1e-12) + 1e-9), col

    def test_cube_cap_forces_further_reduction(self):
        """cube-size cap below what per-SPW alone would allow triggers
        an additional reduction in nchan_agg for rows with large imsize."""
        from dataclasses import replace
        from config import MEMO_CONFIG, MitigationCaps
        from pipeline import load_database
        from recompute import recompute_db
        db = load_database(str(DB_PATH))
        # Cap only cube size; no nchan caps.  Any MOUS whose
        # per-SPW cube exceeds 1 GB should see nchan_agg drop.
        caps_cube = MitigationCaps(mode="caps", cubesize_cap_GB=1.0)
        out_cube = recompute_db(db, replace(MEMO_CONFIG, mitigation=caps_cube))
        col = "wsu_nchan_agg_stepped2_goal"
        a = np.asarray(db[col].value if hasattr(db[col], "value")
                       else db[col], dtype=float)
        b = np.asarray(out_cube[col].value if hasattr(out_cube[col], "value")
                       else out_cube[col], dtype=float)
        # At least some rows must have been reduced.
        assert int((b < a).sum()) > 100
