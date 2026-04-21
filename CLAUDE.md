# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains code, data, and LaTeX documents for projecting data characteristics of ALMA's Wideband Sensitivity Upgrade (WSU) — including data rates, data volumes, cube sizes, and visibility rates — and using those projections for compute and infrastructure cost estimates. The primary input is a database of real ALMA Cycle 7 and 8 observations (`data/wsu_datarates_mit_per_mous_initial_goal_20250423.ecsv`, not in git), from which WSU quantities are extrapolated.

## Running Code

Primary entry points go through `make` (see `make help`):

```
make validation        # reproduce canonical tables + plots in build/validation/
make test              # full pytest suite (~6 min)
make test-fast         # fast tier only (formulas + sample DB + parser; ~1 sec)
make test-keep         # pytest; preserve artifacts in build/test_artifacts/
make list-artifacts    # show what is currently under build/
make clean             # wipe build/ entirely
make clean-test        # remove build/test_artifacts/ only
make clean-validation  # remove build/validation/ only
```

`make validation` writes:
- `build/validation/tables/` — three memo-style + four SDD-style LaTeX tables
- `build/validation/plots/` — SDD-style CCDFs (productsize / datavol / datarate / sysperf)
- `build/validation/realizations/*.ecsv` — N_REAL band-1/2 realization files
- `build/validation/stats.json` — computed medians / TWA / max / totals

Override realization count with `make validation N_REAL=50` (default 10 matches the memo).

The pipeline can also be driven directly (`python pipeline.py --help`) or notebook-style:
```
python sum_m1_prod.py         # legacy product-size summary
jupyter lab                   # exploratory notebooks
```

## Key Terminology

- **BLC** (`blc_*`) — current ALMA Baseline Correlator (12m array) or ACA Correlator (7m array) values
- **WSU** (`wsu_*`) — estimated values for the Wideband Sensitivity Upgrade
- **M1 / M4 / M5** — WSU project milestones (Milestone 1/4/5), roughly equivalent to early/later-2x/later-4x bandwidth stages; the older "IWS/MWS/FWS" terminology is superseded
- **initial / goal** — sub-variants within milestone definitions (see `WSU SDD data volume - Stage Lookup.csv`)
- **MOUS** — Member Observing Unit Set; the primary grouping unit for ALMA observations
- **stepped2** — the preferred velocity-resolution assumption (scenario 2 variant): >10 km/s → 10 km/s, 2–10 → 2, 0.5–2 → 0.5, 0.1–0.5 → 0.1, <0.1 → native resolution

## Core Python Modules

**`wsu_db.py`** — The main calculation library. Contains:
- `create_database()` — builds the per-MOUS WSU parameter database from a cycle 7 archive query table
- `calc_talon_specwidth()` — converts BLC spectral resolution to WSU TALON resolution with per-band minimum channel averaging (`wsu_chanavg_min`)
- Data rate, visibility rate, data volume, and cube/product size computation functions
- Defines a custom astropy unit `gvis` (Gvis/hr) — must be registered in globals before use

**`large_cubes.py`** — Archive query and munging utilities:
- `get_archive_info()` / `read_archive_info()` — TAP query and read from ALMA archive via astroquery
- `munge_archive_info()` — computes imsize, nchan, pb, cell, mosaic flag, and mitigation status from raw archive rows
- `calculate_nchan()` — reverse-engineers number of channels from bandwidth + spectral resolution using BLC correlator rules
- `calc_nchan_max()`, `calc_cube_size()`, `calc_mfs_size()` — key imaging size formulas

**`wsu_plots.py`** — Consolidated plotting code with shared color scheme (`mycolors`) and fiducial Band 2 spectral scan values.

**`sum_m1_prod.py`** — Standalone script summarizing Milestone 1 delivered product sizes by MOUS and ALMA cycle, with optional per-cube size caps.

## Key Equations

```python
cubesize  = 4.0 * imsize**2 * nchan / 1.0e9          # GB
mfssize   = 4.0 * imsize**2 / 1.0e9                   # GB
productsize = 2.0 * (cubesize + mfssize) * nspws

# Output data rate (GB/s)
datarate = (2 * Napc * Nant*(Nant-1)/2 + 4*Nant) * Nchannels * Npols / Tintegration

# Visibility rate (Gvis/hr)
visrate = (2.0 * Npols * Nbase * Nchannels / 1e9) / Tintegration_in_hr
```

WSU integration times: 3.072s (12m array), 9.984s (7m array).  
Typical antenna counts: 47 (12m), 10 (7m).

## Notebook Execution Order

For producing the full WSU database and summary products from scratch:
1. `largeCube_search_v2.ipynb` — archive query and initial data munging
2. `WSU_size_of_computing_numbers_database.ipynb` — builds WSU parameter database
3. `WSU_size_of_computing_numbers_UseCasePlots.ipynb` — use-case plots
4. `WSU_size_of_computing_numbers_CrystalPlots.ipynb` — summary plots

## Data Files

The main database (`data/wsu_datarates_mit_per_mous_initial_goal_20250423.ecsv`) is not in git due to size. A header + 20-row sample is at `data/wsu_datarates_mit_per_mous_initial_goal_20250423_head_and_sample.ecsv`. Read with:
```python
from astropy.table import Table
t = Table.read("data/wsu_datarates_mit_per_mous_initial_goal_20250423.ecsv", format="ascii.ecsv")
```

The README contains a full description of all ~100 column names in the database.

## LaTeX Documents

- `Estimated_WSU_Data_Properties/` — original methodology memo (pre-SDD terminology)
- `WSU_Estimated_Data_Properties_Update/` — updated memo with M1/M4/M5 milestone definitions used in the SDD

Both use `main.tex` as the entry point. Tables are generated separately and `\input{}`-ed from the `tables/` subdirectory.

## Active Tasks (from AGENTS.md)

1. **Core goal**: A script that correctly, reproducibly calculates data rates for M1, M4, M5; produces LaTeX summary tables; and produces distribution plots — matching the style of the WSU SDD data properties appendix
2. Identify errors/inconsistencies in the SDD data properties projections, TSI cost report, and two WSU data properties memos
3. Define regression tests anchored to SDD table entries and current BLC/ACA values
4. Future: smarter mitigation, channel-resolution assumptions, and integration time choices
