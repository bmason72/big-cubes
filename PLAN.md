# Plan: WSU Data Properties Pipeline

## Context

The repository contains code, notebooks, and documents for projecting ALMA WSU data characteristics (data rates, volumes, cube sizes, system performance) across milestones M1, M4, M5. The current pipeline is notebook-driven, with calculation logic spread across `wsu_db.py` (4600 lines), `large_cubes.py`, and several notebooks. No tests exist. The reference values live in LaTeX tables across three document sets (memo, updated memo, SDD) that partially disagree with each other.

**Goal**: A standalone, tested script that reproduces the memo and SDD tables from the database, identifies discrepancies between documents, and serves as the canonical calculation going forward.

Tasks 1 (pipeline), 2 (error/inconsistency identification), and 3 (regression tests) are treated as one interlocked workstream. Tasks 4 (smarter assumptions) and 5 (GUI) are borne in mind — particularly in keeping assumptions parameterized — but not implemented here.

---

## Phase 1: Config + Reference Values + Stats Engine

### 1a. Create `config.py` — single source of truth for assumptions

Extract from `wsu_db.py` all hardcoded constants into a structured config:

- **Integration times**: 6.144s (M1 12m), 3.072s (M4/M5 12m), 9.984s (7m), BLC values
- **Antenna counts**: typical (47/10), peak (50/12), all (66/16); M1-specific if different
- **Bandwidth per band per milestone**: the dictionaries from `create_initial_wsu_db()` (lines ~4340-4453 of `wsu_db.py`)
- **Channel averaging minima**: `wsu_chanavg_min` and `wsu_chanavg_min_initial` dicts
- **Velocity resolution bins** (stepped2 definition): >10→10, 2-10→2, 0.5-2→0.5, 0.1-0.5→0.1, <0.1→native
- **Imaging parameters**: pixels_per_beam=5, product_factor=2.0
- **System performance parameters**: k=20, flops_per_vis, core_efficiency, parallelization_efficiency
- **SPW bandwidth**: 2.0 GHz; TALON channel: 13.5 kHz

This is the key enabler for Task 4 (experimenting with different assumptions).

### 1b. Create `reference_values.py` — validation targets

Transcribe every numeric cell from these tables:

| Source | File | Content |
|--------|------|---------|
| Memo | `WSU_Estimated_Data_Properties_Update/tables/wsu_datarate_summary_initial_goal_v2.tex` | M1/M4/M5 data rates + channel counts |
| Memo | `WSU_Estimated_Data_Properties_Update/tables/wsu_datavol_summary_initial_goal_v2.tex` | M1/M4/M5 vis volumes + product sizes |
| Memo | `WSU_Estimated_Data_Properties_Update/tables/wsu_sysperf_summary_initial_goal_v2.tex` | M1/M4/M5 data rates + vis rates + sysperf |
| SDD | `wsu-sdd-excerpts/tables/table_data_prop.tex` | M5 vs pre-WSU vis vol + product size |
| SDD | `wsu-sdd-excerpts/tables/table_system_performance_wsu.tex` | M5 vs pre-WSU data rate + vis rate + sysperf |
| SDD | `wsu-sdd-excerpts/tables/data_prop_summ.tex` | Summary: raw vol, product vol, compressed, PFLOP/s |
| SDD | `wsu-sdd-excerpts/tables/table_data_rate_peak.tex` | Peak data rate parameters + result |

Each entry: value, unit, source (document + table + row + column), tolerance.

### 1c. Create `pipeline.py` — stats computation engine

Core function `compute_stats(db_path, config)`:

1. Load database, register `Gvis` unit
2. Generate Band 1/2 realizations via `wsu_db.generate_db_realizations(mydb, n=..., add_initial=True, add_goal=True, add_ms4=True)` — the code for this is in `wsu_db.py` lines 3084-3423 (`add_bands` + `generate_db_realizations`). Historical usage: `n=50` in `WSU_add_band12_to_database.ipynb`, `n=10` in `calculate_initial_WSU_data_properties.ipynb`. Need to determine which was used for the final memo tables.
3. Compute stats across realizations via `wsu_db.calc_wsu_stats_allsamples()` (line 3612) or equivalent logic
4. For each milestone (M1/M4/M5) and BLC:
   - Map milestone to database column names (handling the naming asymmetry: M1 uses `_initial_stepped2_initial`, M4/M5 use `_ms4/_goal_stepped2_typical`)
   - For each array type (12m, 7m, both):
     - Filter rows
     - Compute: median, time-weighted average (using `weights_all` column), maximum
     - Compute totals: sum ÷ 2 for per-cycle (database covers 2 cycles)
   - Quantities: datarate, visrate, nchan_agg, datavol_total, datavol_science, productsize, cubesize, sysperf

Column name mapping (critical detail):

| Milestone | Data rate column | Datavol total column | Product size column |
|-----------|-----------------|---------------------|-------------------|
| M1 | `wsu_datarate_initial_stepped2_initial` | `wsu_datavol_initial_stepped2_initial_total` | `wsu_productsize_initial_stepped2` |
| M4 | `wsu_datarate_ms4_stepped2_typical` | `wsu_datavol_ms4_stepped2_typical_total` | `wsu_productsize_ms4_stepped2` |
| M5 | `wsu_datarate_goal_stepped2_typical` | `wsu_datavol_goal_stepped2_typical_total` | `wsu_productsize_goal_stepped2` |
| BLC | `blc_datarate_typical` | `blc_datavol_typical_total` | `blc_productsize` |

### 1d. Create `test_pipeline.py` — regression tests (pytest)

Three tiers:
1. **Formula unit tests** (always run): verify `calc_datarate`, `calc_visrate`, `calc_cube_size` against hand calculations. Verify peak data rate = 3.96 GB/s for SDD parameters.
2. **Sample tests** (always run, uses 20-row sample in git): verify loading, no NaN, reasonable bounds, correct column access.
3. **Full-DB validation** (conditional, `@pytest.mark.skipif` if DB not present): verify computed stats match every cell in `reference_values.py` memo tables within stated tolerances.

**Deliverable**: Run `pytest test_pipeline.py` and all full-DB tests pass against memo table values.

---

## Phase 2: LaTeX Table Generation

Extend `pipeline.py` with `generate_tables(stats, output_dir)`:

### 2a. Memo-style tables (M1/M4/M5 × 12m/7m/both)

Produce three sideways tables matching existing format:
- `wsu_datarate_summary.tex` — data rates + channel counts
- `wsu_datavol_summary.tex` — vis volumes + product sizes
- `wsu_sysperf_summary.tex` — data rates + vis rates + sysperf

Reuse the LaTeX formatting patterns from `wsu_db.py`'s `make_wsu_stats_table_newstats_*` functions, but read from the new stats dictionary.

### 2b. SDD-style tables (M5 vs pre-WSU)

Produce tables matching `wsu-sdd-excerpts/tables/` format:
- `table_data_prop.tex` — M5 vs pre-WSU vis vol + product size
- `table_data_rates_prop.tex` — M5 data rates only
- `table_system_performance_wsu.tex` — M5 vs pre-WSU rates + vis rates + sysperf
- `data_prop_summ.tex` — top-level summary

**Verification**: Diff generated tables against existing ones.

---

## Phase 3: Distribution Plots

Extend `pipeline.py` with `generate_plots(db, stats, config, output_dir)`:

Call existing `wsu_plots.py` SDD-style functions:
- Visibility data volume CCDF (`plot_datavol_comparison_sdd`)
- Product size CCDF (`plot_productsize_comparison_sdd`)
- Data rate cumulative distribution
- System performance cumulative distribution

Produce both M5-vs-preWSU plots (SDD style) and optionally M1/M4/M5 comparison plots.

---

## Phase 4: Consistency Report

### 4a. Create `consistency_report.py`

Systematically compare:

**Memo vs SDD (M5 values)** — known discrepancies to investigate:

| Quantity | Memo M5 12m | SDD M5 12m | Delta | Notes |
|----------|-------------|------------|-------|-------|
| DR median | 0.121 GB/s | 0.136 GB/s | +12% | SDD "conformance adjustment" for M5 band config |
| DR TWA | 0.478 GB/s | 0.477 GB/s | ~0 | OK |
| Vis rate median | 208.9 Gvis/hr | 230 Gvis/hr | +10% | Same conformance issue? |
| Sysperf TWA | 0.372 PFLOP/s | 0.50 PFLOP/s | +34% | Possibly includes TP and/or GOUS processing |
| Sysperf max | 7.480 PFLOP/s | 10.3 PFLOP/s | +38% | Same — not a simple 3x multiplier |
| Product total | 12.114 PB | 12.0 PB | ~0 | Rounding |

**Important note**: Some SDD vs memo differences may be explained by:
- **TP inclusion**: The SDD `data_prop_summ.tex` explicitly says "7m+12m+TP" for some quantities, while the memo tables show 12m/7m/both without TP. TP adds ~3% per the SDD text.
- **GOUS processing**: The SDD summary table includes a "Prod. volume with GOUS" row (10.7 PB) distinct from the basic product volume (12.1 PB). GOUS scaling adds a ~1.36x factor.
- **Contingency factors**: SDD system performance estimates include "3x overhead/contingency factor" per the table note, but this may be applied differently than a simple multiplier.

These potential explanations should be verified before flagging values as errors.

**Computed vs memo** — pipeline output should match memo values exactly (same database + code).

**Peak data rate self-consistency** — verify SDD formula with adopted parameters yields stated result.

### 4b. Add consistency test cases to `test_pipeline.py`

Document each discrepancy as either:
- An explained difference (with explanation in test docstring)
- An unexplained difference (marked `xfail` pending investigation)

---

## Phase 5: Cleanup

- Add `Makefile` or `pyproject.toml` with: `make pipeline`, `make test`, `make tables`, `make plots`
- Update `CLAUDE.md` with new pipeline usage
- Optionally refactor `wsu_db.py` to import from `config.py` (low priority, careful not to break notebooks)

---

## Files Created / Modified

| File | Action | Purpose |
|------|--------|---------|
| `config.py` | **Create** | Parameterized assumptions (constants, band tables, milestone definitions) |
| `reference_values.py` | **Create** | Validation targets from memo + SDD tables |
| `pipeline.py` | **Create** | Main pipeline: load → compute → tables → plots → report |
| `test_pipeline.py` | **Create** | Regression tests (pytest) |
| `consistency_report.py` | **Create** | Cross-document comparison report |
| `wsu_db.py` | Minimal modify | Only if needed: fix stats function for M1/M4/M5 support |
| `CLAUDE.md` | Update | Add pipeline usage instructions |

---

## Band 1/2 Realizations

The memo statistics include synthetic Band 1/2 observations (bands not present in Cycle 7/8). The code to generate these is already in `wsu_db.py`:

- **`add_bands(mydb, array, band, total_time, ...)`** (line 3084): creates synthetic Band 1/2 rows by randomly sampling from Band 3, scaling frequency-dependent quantities (FOV, resolution, PB, cell), and assigning band-appropriate bandwidths per milestone.
- **`generate_db_realizations(mydb, outDir, n, ...)`** (line 3362): removes 10% 12m / 6% 7m time to make room, generates Band 1+2 rows for both arrays, stacks into N complete realization files in `data/sample_band1_band2/`.
- **`calc_wsu_stats_allsamples(outDir, ...)`** (line 3612): reads all realization files, computes stats averaged across them.

The pipeline will call `generate_db_realizations()` to produce the realization files from the base database, then `calc_wsu_stats_allsamples()` to compute the mean statistics. This should reproduce memo values without any external files. The realization count (`n`) and random seed should be configurable. Historical usage: `WSU_add_band12_to_database.ipynb` used `n=50`; `calculate_initial_WSU_data_properties.ipynb` used `n=10` with `add_initial=True, add_goal=True, add_ms4=True`. Need to determine which was used for the final memo tables.

---

## Future Work Context (Task 4)

The `obsPrepTsiReport.txt` (Section 6) describes how the WSU Observing Tool will handle spectral setup via ATAC Frequency Slices (200 MHz each, 80 in 2x / 160 in 4x), with Channel Averaging Factors (CAFs) and zoom modes. The OT will capture user *science intent* (specific lines, linewidths, priorities) and automatically optimize spectral setups. This is more nuanced than our stepped2 velocity binning — which assigns all projects in a bin to the finest resolution in that bin. For aggregate projections the stepped2 approach is adequate, but Task 4 could benefit from modeling actual spectral line science cases. This document also confirms the ATAC/TPGS technical constraints (13.5 kHz base resolution, 200 MHz FS granularity, zoom factors 2/4/8) that `config.py` should capture.

---

## Verification

After each phase:
- `pytest test_pipeline.py -v` — all tests pass
- `python pipeline.py --db-path data/wsu_datarates_mit_per_mous_initial_goal_20250423.ecsv` — runs without error
- Generated LaTeX tables diff cleanly against reference tables (or differences are documented)
- Generated plots visually match SDD figures
