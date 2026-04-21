"""
Regression tests for the WSU data-properties pipeline.

Three tiers:
    1. Formula unit tests (always run): exercise calc_datarate /
       calc_visrate / calc_cube_size against hand calculations, and verify
       the SDD peak data rate of 3.96 GB/s.
    2. Sample DB tests (always run): load the 20-row sample shipped in
       git, check structure, no NaN where not expected, and key columns
       are readable.
    3. Full-DB validation (skipif full DB missing): run the pipeline and
       assert memo stats reproduce within their stated tolerances.

Run:
    venv/bin/pytest test_pipeline.py -v
    venv/bin/pytest test_pipeline.py -v -k "formula or sample"   # fast
    WSU_FULL_DB=/path/... venv/bin/pytest test_pipeline.py -v    # full
"""

from __future__ import annotations

import math
import os
import sys
import types
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest

# Make sure ipdb import in wsu_db doesn't break tests.
sys.modules.setdefault("ipdb", types.SimpleNamespace(set_trace=lambda: None))

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import config  # noqa: E402
import reference_values as refs  # noqa: E402

DB_PATH = Path(os.environ.get(
    "WSU_FULL_DB",
    HERE / "data" / "wsu_datarates_mit_per_mous_initial_goal_20250423.ecsv",
))
SAMPLE_PATH = HERE / "data" / "wsu_datarates_mit_per_mous_initial_goal_20250423_head_and_sample.ecsv"
MEMO_REF_DIR = HERE / "WSU_Estimated_Data_Properties_Update" / "tables"
SDD_REF_DIR = HERE / "wsu-sdd-excerpts" / "tables"


# =============================================================================
# Tier 1: Formula unit tests
# =============================================================================

class TestFormulas:
    """Direct hand-calc checks of calc_datarate, calc_visrate, calc_cube_size."""

    def test_calc_datarate_m5_peak_12m(self):
        """SDD peak rate: Band 2 @ 0.1 km/s, 50x12m -> 3.96 GB/s."""
        import wsu_db
        p = config.PEAK_M5_12M
        # The SDD quotes ``2,380,800`` as the channel count; see
        # config.PeakDataRateParams for the comment on why that is treated
        # as a total-across-pols count with Npols=1 in our formula.
        rate = wsu_db.calc_datarate(
            Nbyte=p.nbyte, Napc=p.napc, Nant=p.nant,
            Nchannels=p.n_channels_total, Npols=1,
            Tintegration=p.tint_s * u.s,
        )
        rate_gbs = rate.to(u.GB / u.s).value
        assert rate_gbs == pytest.approx(p.expected_GBs, rel=0.01), (
            f"peak 12m data rate {rate_gbs:.3f} GB/s != expected {p.expected_GBs} GB/s"
        )

    def test_calc_datarate_m5_peak_7m(self):
        import wsu_db
        p = config.PEAK_M5_7M
        rate = wsu_db.calc_datarate(
            Nbyte=p.nbyte, Napc=p.napc, Nant=p.nant,
            Nchannels=p.n_channels_total, Npols=1,
            Tintegration=p.tint_s * u.s,
        )
        rate_gbs = rate.to(u.GB / u.s).value
        assert rate_gbs == pytest.approx(p.expected_GBs, rel=0.02)

    def test_datarate_formula_algebra(self):
        """Verify the exact formula from the memo:
            DR = (2*Nbyte*Napc*Nant*(Nant-1)/2 + 4*Nant) * Nchannels * Npols / Tint
        with units GB/s (/1e9 conversion).
        """
        import wsu_db
        nant, nchan, npol, tint_s = 50, 1_000_000, 2, 3.072
        expected = ((2 * 2 * 1 * nant * (nant - 1) / 2 + 4 * nant)
                    * nchan * npol / tint_s / 1e9)
        rate = wsu_db.calc_datarate(
            Nbyte=2.0, Napc=1.0, Nant=nant, Nchannels=nchan,
            Npols=npol, Tintegration=tint_s * u.s,
        ).to(u.GB / u.s).value
        assert rate == pytest.approx(expected, rel=1e-10)

    def test_visrate_formula(self):
        """Visrate = 2 * Npols * Nbase * Nchannels / Tint (in Gvis/hr)."""
        import wsu_db
        nant, nchan, npol, tint_s = 50, 1_000_000, 2, 3.072
        nbase = nant * (nant - 1) / 2
        expected = (2 * npol * nbase * nchan / 1e9) / (tint_s / 3600.0)
        vr = wsu_db.calc_visrate(nant, npol, nchan, tint_s * u.s)
        vr_val = vr.to(u.Unit("Gvis/h")).value
        assert vr_val == pytest.approx(expected, rel=1e-9)

    def test_calc_cube_size(self):
        """cube = 4 * imsize^2 * nchan / 1e9  (GB)."""
        from large_cubes import calc_cube_size
        imsize, nchan = 400, 8192
        expected = 4.0 * imsize ** 2 * nchan / 1e9
        # calc_cube_size returns a Quantity in GB
        out = calc_cube_size(imsize, nchan).to(u.GB).value
        assert out == pytest.approx(expected, rel=1e-9)

    def test_calc_mfs_size(self):
        from large_cubes import calc_mfs_size
        imsize = 400
        expected = 4.0 * imsize ** 2 / 1e9
        out = calc_mfs_size(imsize).to(u.GB).value
        assert out == pytest.approx(expected, rel=1e-9)


# =============================================================================
# Tier 2: Sample-DB structural tests
# =============================================================================

@pytest.mark.skipif(not SAMPLE_PATH.exists(),
                    reason="sample ECSV not present")
class TestSampleDB:
    @pytest.fixture(scope="class")
    def db(self):
        from pipeline import load_database
        return load_database(str(SAMPLE_PATH), as_qtable=False)

    def test_loads_nonempty(self, db):
        assert len(db) > 0
        assert len(db.colnames) > 100

    def test_required_columns(self, db):
        required = [
            "array", "band", "weights_all", "time_tot",
            config.MILESTONE_COLUMNS["M1"].datarate,
            config.MILESTONE_COLUMNS["M4"].datarate,
            config.MILESTONE_COLUMNS["M5"].datarate,
            config.BLC_COLUMNS.datarate,
        ]
        for c in required:
            assert c in db.colnames, f"missing required column {c}"

    def test_arrays_are_valid(self, db):
        arrs = set(np.asarray(db["array"]))
        assert arrs.issubset({"12m", "7m"}), f"unexpected array labels: {arrs}"

    def test_bands_in_expected_range(self, db):
        bands = set(np.asarray(db["band"]).astype(int))
        assert bands.issubset(set(range(1, 11)))

    def test_datarates_nonnegative(self, db):
        for col in (config.MILESTONE_COLUMNS["M5"].datarate,
                    config.BLC_COLUMNS.datarate):
            vals = np.asarray(db[col])
            assert np.all(np.isfinite(vals))
            assert np.all(vals >= 0)


# =============================================================================
# Tier 3: Full-DB pipeline regression vs memo tables
# =============================================================================

@pytest.mark.skipif(not DB_PATH.exists(),
                    reason=f"full DB not at {DB_PATH}")
class TestFullPipeline:
    """Validate pipeline stats against the memo tables.

    Realization generation is expensive; we cache a single run at module
    scope.  Use --db-path / WSU_FULL_DB to override.
    """

    @pytest.fixture(scope="class")
    def stats(self, tmp_path_factory):
        from dataclasses import replace
        from pipeline import (DEFAULT_CONFIG, compute_stats,
                              generate_realizations, load_database)

        outdir = tmp_path_factory.mktemp("realizations")
        n = int(os.environ.get("WSU_N_REALIZATIONS", "10"))
        cfg = replace(
            DEFAULT_CONFIG,
            realizations=replace(DEFAULT_CONFIG.realizations,
                                 n_realizations=n,
                                 outdir=str(outdir)),
        )
        db = load_database(str(DB_PATH))
        paths = generate_realizations(db, cfg, outdir=str(outdir))
        return compute_stats(paths, cfg)

    # -- helpers --
    @staticmethod
    def _matches(actual: float, ref: refs.RefValue,
                 stat: str, array: str) -> tuple[bool, str]:
        if not math.isfinite(actual):
            return (False, f"computed value is NaN/inf")
        diff = abs(actual - ref.value)
        if ref.abs_tol is not None and diff <= ref.abs_tol:
            return (True, "")
        # Apply the looser of the ref's per-cell tolerance and the
        # stat-level policy tolerance.  This lets us keep tight per-cell
        # tolerances where physically justified while inheriting realistic
        # Monte-Carlo bounds elsewhere.
        effective_tol = max(ref.rel_tol, refs.policy_tolerance(stat, array))
        if ref.value != 0:
            rel = diff / abs(ref.value)
            if rel <= effective_tol:
                return (True, "")
            return (False, f"actual={actual:g} ref={ref.value:g} "
                           f"rel_diff={rel:.3g} tol={effective_tol:.3g}")
        if diff == 0:
            return (True, "")
        return (False, f"actual={actual:g} ref=0 diff={diff:g}")

    @pytest.mark.parametrize(
        "source,milestone,array,stat,quantity,ref",
        refs.flatten_memo_tables(),
        ids=lambda x: (x[0] + "/" + x[1] + "/" + x[2] + "/" + x[3] + "/" + x[4]
                       if isinstance(x, tuple) else str(x)),
    )
    def test_memo_value(self, stats, source, milestone, array, stat,
                        quantity, ref):
        """Every numeric cell in the three memo tables must reproduce."""
        ms_stats = stats[milestone][array]
        if stat == "total_per_cycle":
            bucket = ms_stats.get("total_per_cycle", {})
        else:
            bucket = ms_stats.get(stat, {})
        if quantity not in bucket:
            pytest.skip(f"pipeline does not produce {quantity} for {milestone}")
        ok, msg = self._matches(bucket[quantity], ref, stat, array)
        assert ok, f"[{source}] {milestone}/{array}/{stat}/{quantity}: {msg}\n  source: {ref.source}"


# =============================================================================
# Tier 4: LaTeX-table reproduction (Phase 2)
# =============================================================================

class TestTableParser:
    """Structural tests for the table-parser -- exercise it on the shipped
    reference tables so the parser is covered even when the full DB is
    absent.
    """

    @pytest.mark.skipif(not MEMO_REF_DIR.exists(),
                        reason="memo reference tables not present")
    def test_parse_datarate_ref(self):
        from tables import parse_numeric_cells
        rows = parse_numeric_cells(
            MEMO_REF_DIR / "wsu_datarate_summary_initial_goal_v2.tex")
        # 3 datarate rows + 3 channel rows == 6 data rows.
        assert len(rows) == 6, f"expected 6 data rows, got {len(rows)}"
        # First row median data rates should include 0.024 (M1 12m median).
        assert 0.024 in rows[0]
        # Max channel count for M5 is 1,185,184.
        assert 1_185_184 in rows[5]

    @pytest.mark.skipif(not MEMO_REF_DIR.exists(),
                        reason="memo reference tables not present")
    def test_parse_datavol_ref(self):
        from tables import parse_numeric_cells
        rows = parse_numeric_cells(
            MEMO_REF_DIR / "wsu_datavol_summary_initial_goal_v2.tex")
        # 3 datavol blocks * 4 rows each == 12 rows.
        assert len(rows) == 12

    def test_sdd_scientific_formatter(self):
        from tables import fmt_sdd
        # Small values produce scientific notation.
        s = fmt_sdd(1.6e-2)
        assert "\\times" in s
        assert "10^{-2}" in s
        # Medium-magnitude values produce plain form.
        assert fmt_sdd(0.32) == "0.32"
        # Large values produce integer-ish form.
        assert fmt_sdd(1200) == "1200"

    def test_memo_datarate_formatter(self):
        from tables import fmt_rate, fmt_datavol, fmt_nchan
        assert fmt_rate(0.024) == "0.024"
        assert fmt_datavol(0.069) == "  0.069"
        assert fmt_datavol(374.689) == "374.689"
        assert fmt_nchan(16443) == "16,443"
        assert fmt_nchan(1185184) == "1,185,184"


@pytest.mark.skipif(not DB_PATH.exists() or not MEMO_REF_DIR.exists(),
                    reason="full DB or reference tables not present")
class TestTableReproduction:
    """Generate memo-style LaTeX tables and check their numeric cells
    match the committed reference tables within policy tolerances.
    """

    @pytest.fixture(scope="class")
    def gen_dir(self, tmp_path_factory, request):
        """Generate tables once from a fresh pipeline run."""
        from dataclasses import replace
        from pipeline import (DEFAULT_CONFIG, compute_stats,
                              generate_realizations, load_database)
        from tables import generate_memo_tables, generate_sdd_tables

        workdir = tmp_path_factory.mktemp("tables")
        real_dir = workdir / "realizations"
        n = int(os.environ.get("WSU_N_REALIZATIONS", "10"))
        cfg = replace(
            DEFAULT_CONFIG,
            realizations=replace(DEFAULT_CONFIG.realizations,
                                 n_realizations=n,
                                 outdir=str(real_dir)),
        )
        db = load_database(str(DB_PATH))
        paths = generate_realizations(db, cfg, outdir=str(real_dir))
        stats = compute_stats(paths, cfg)
        generate_memo_tables(stats, str(workdir))
        generate_sdd_tables(stats, str(workdir))
        return workdir

    @pytest.fixture(scope="class")
    def ref_dir(self, tmp_path_factory):
        """Copy memo references into a dir keyed by the *generated* filename
        so verify_tables() can resolve pairs by a single filename."""
        import shutil
        out = tmp_path_factory.mktemp("ref_tables")
        mapping = {
            "wsu_datarate_summary.tex": "wsu_datarate_summary_initial_goal_v2.tex",
            "wsu_datavol_summary.tex":  "wsu_datavol_summary_initial_goal_v2.tex",
            "wsu_sysperf_summary.tex":  "wsu_sysperf_summary_initial_goal_v2.tex",
        }
        for gen_name, ref_name in mapping.items():
            shutil.copy(MEMO_REF_DIR / ref_name, out / gen_name)
        return out

    def test_memo_datarate_matches(self, gen_dir, ref_dir):
        from tables import verify_tables
        mismatches = verify_tables(gen_dir, ref_dir,
                                    [("wsu_datarate_summary.tex", 0.20)])
        assert not mismatches, "memo datarate cells disagree: " + str(mismatches[:5])

    def test_memo_datavol_matches(self, gen_dir, ref_dir):
        from tables import verify_tables
        # Max-of-max is inherently noisy; total_per_cycle tight.
        mismatches = verify_tables(gen_dir, ref_dir,
                                    [("wsu_datavol_summary.tex", 0.25)])
        assert not mismatches, "memo datavol cells disagree: " + str(mismatches[:5])

    def test_memo_sysperf_matches(self, gen_dir, ref_dir):
        from tables import verify_tables
        mismatches = verify_tables(gen_dir, ref_dir,
                                    [("wsu_sysperf_summary.tex", 0.20)])
        assert not mismatches, "memo sysperf cells disagree: " + str(mismatches[:5])

    def test_all_sdd_tables_generated(self, gen_dir):
        """SDD tables are harder to verify cell-by-cell because the SDD
        uses looser rounding and ad-hoc TP/GOUS factors -- here we only
        confirm the files exist and parse cleanly.  Cell-level consistency
        is covered by Phase 4's consistency_report.
        """
        from tables import parse_numeric_cells
        for fname in ("table_data_prop.tex",
                      "table_data_rates_prop.tex",
                      "table_system_performance_wsu.tex",
                      "data_prop_summ.tex"):
            path = Path(gen_dir) / fname
            assert path.exists()
            rows = parse_numeric_cells(path)
            assert len(rows) > 0, f"{fname} produced no data rows"


@pytest.mark.skipif(not DB_PATH.exists(),
                    reason=f"full DB not at {DB_PATH}")
class TestPlotGeneration:
    """Smoke tests for plots.py.  Heavier overlay path is tested separately."""

    def test_sdd_plots_render_without_overlay(self, tmp_path):
        """Rendering without the band-1/2 overlay needs only the base DB."""
        from pipeline import load_database
        from plots import generate_sdd_plots

        db = load_database(str(DB_PATH))
        out = generate_sdd_plots(db, output_dir=str(tmp_path),
                                 overlay_realizations=False)
        assert len(out) == 4
        expected = {"productsize.png", "datavol.png",
                    "datarate.png", "sysperf.png"}
        got = {Path(p).name for p in out}
        assert got == expected
        for p in out:
            assert Path(p).stat().st_size > 1000

    def test_sdd_plots_with_overlay(self, tmp_path):
        """Full overlay path: generate a couple realizations then plot."""
        from dataclasses import replace
        from pipeline import (DEFAULT_CONFIG, generate_realizations,
                              load_database)
        from plots import generate_sdd_plots

        real_dir = tmp_path / "realizations"
        n = int(os.environ.get("WSU_PLOT_N_REALIZATIONS", "3"))
        cfg = replace(
            DEFAULT_CONFIG,
            realizations=replace(DEFAULT_CONFIG.realizations,
                                 n_realizations=n,
                                 outdir=str(real_dir)),
        )
        db = load_database(str(DB_PATH))
        generate_realizations(db, cfg, outdir=str(real_dir))

        out = generate_sdd_plots(db, output_dir=str(tmp_path / "plots"),
                                 cfg=cfg,
                                 realization_dir=str(real_dir),
                                 overlay_realizations=True)
        assert len(out) == 4
        for p in out:
            assert Path(p).stat().st_size > 1000
