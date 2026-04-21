"""
WSU data-properties stats pipeline.

Entry points:
    load_database(path)          -> astropy.table.QTable with Gvis registered
    generate_realizations(db, cfg, outdir) -> list of per-realization paths
    compute_stats(db_or_paths, cfg)         -> nested dict of statistics
    main(argv)                              -> CLI wrapper

The pipeline reproduces the memo tables' M1/M4/M5 statistics by:
  1. Loading the base per-MOUS database (cycle 7/8, bands 3-10).
  2. Generating N realizations that replace some band-3 time with
     synthetic band-1/band-2 rows (via wsu_db.generate_db_realizations).
  3. Computing median / time-weighted-average / max / total per-cycle for
     each (milestone, array, quantity) across all realizations.

Data volume / product size values are converted to TB or PB during stats
collection so downstream code reads directly-comparable scalars.
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import sys
import types
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

# wsu_db.py imports ipdb unconditionally; ipdb is not in the standard venv.
# Stub before any transitive import so the CLI works without test infra.
sys.modules.setdefault("ipdb", types.SimpleNamespace(set_trace=lambda: None))

import astropy.units as u
import numpy as np
from astropy.table import QTable, Table

# Local imports
from config import (BLC_COLUMNS, DEFAULT_CONFIG, MILESTONE_COLUMNS,
                    MilestoneColumns, PipelineConfig)

log = logging.getLogger("wsu_pipeline")


# ---------------------------------------------------------------------------
# Gvis unit registration
# ---------------------------------------------------------------------------
# wsu_db.py defines a custom astropy unit "Gvis" in its module globals.  We
# must do the same before reading any ecsv that references it -- otherwise
# astropy emits UnitsWarning and drops the unit.

_GVIS_REGISTERED = False


def register_gvis_unit() -> u.Unit:
    """Register the custom Gvis unit so stacked tables can read it.

    wsu_db.py registers Gvis at import time; we prefer to piggy-back on
    that so the same unit object lives in both modules' namespaces.
    """
    global _GVIS_REGISTERED
    if not _GVIS_REGISTERED:
        try:
            import wsu_db  # noqa: F401  -- side effect: registers Gvis
        except ImportError:
            try:
                gvis = u.def_unit("Gvis")
                u.add_enabled_units([gvis])
            except Exception:
                pass
        _GVIS_REGISTERED = True
    return u.Unit("Gvis")


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_database(path: str | os.PathLike, as_qtable: bool = True) -> Table:
    """Load a per-MOUS database .ecsv file, registering Gvis first.

    Uses QTable by default so that unit columns become astropy Quantity
    arrays -- required by wsu_db.generate_db_realizations(), which relies
    on Quantity arithmetic.  Pass as_qtable=False for readonly analyses.
    """
    register_gvis_unit()
    log.info("Loading database %s", path)
    cls = QTable if as_qtable else Table
    return cls.read(str(path), format="ascii.ecsv")


# ---------------------------------------------------------------------------
# Realization generation
# ---------------------------------------------------------------------------

def generate_realizations(db: Table, cfg: PipelineConfig,
                          outdir: Optional[str] = None) -> List[str]:
    """Generate N realizations of the database including Band 1/2.

    Thin wrapper around wsu_db.generate_db_realizations().  Calls the
    existing code with add_initial=add_ms4=add_goal=True (required for the
    M1/M4/M5 columns to be populated on the synthesised rows).

    Returns the list of written ecsv paths.  Re-uses existing files if the
    requested count already exists and the caller passes an existing
    outdir -- use ``outdir=None`` to regenerate into the config's default
    directory.
    """
    register_gvis_unit()

    # Deferred import: wsu_db is large and pulls in matplotlib etc.
    import wsu_db

    out = outdir or cfg.realizations.outdir
    os.makedirs(out, exist_ok=True)

    existing = sorted(Path(out).glob(f"{cfg.realizations.filename}_*.ecsv"))
    if len(existing) >= cfg.realizations.n_realizations:
        log.info("Reusing %d existing realizations in %s",
                 len(existing), out)
        return [str(p) for p in existing[: cfg.realizations.n_realizations]]

    # Seed for reproducibility; wsu_db uses np.random.default_rng() which
    # reads from a module-level BitGenerator, so we reseed numpy's default.
    np.random.seed(cfg.realizations.seed)

    log.info("Generating %d realizations into %s",
             cfg.realizations.n_realizations, out)
    wsu_db.generate_db_realizations(
        db,
        outDir=out,
        filename=cfg.realizations.filename,
        frac_12m=cfg.realizations.frac_12m_replaced,
        frac_7m=cfg.realizations.frac_7m_replaced,
        n=cfg.realizations.n_realizations,
        add_initial=True,
        add_goal=True,
        add_ms4=True,
    )
    written = sorted(Path(out).glob(f"{cfg.realizations.filename}_*.ecsv"))
    return [str(p) for p in written]


# ---------------------------------------------------------------------------
# Stats computation
# ---------------------------------------------------------------------------

# Quantities computed for each milestone / array.  Units listed are what we
# *report* -- internal conversions happen below.
QUANTITY_SPEC: Dict[str, Dict[str, str]] = {
    # name        : { "column": <attr of MilestoneColumns>, "report_unit": <str> }
    "datarate":        {"column": "datarate",         "report_unit": "GB/s"},
    "visrate":         {"column": "visrate",          "report_unit": "Gvis/hr"},
    "nchan_agg":       {"column": "nchan_agg",        "report_unit": ""},
    "sysperf":         {"column": "sysperf",          "report_unit": "PFLOP/s"},
    "datavol_total":   {"column": "datavol_total",    "report_unit": "TB"},
    "datavol_science": {"column": "datavol_target_tot","report_unit": "TB"},
    "productsize":     {"column": "productsize",      "report_unit": "TB"},
}


_ARRAY_MASKS = ("12m", "7m", "both")


def _array_mask(db: Table, array: str) -> np.ndarray:
    if array == "both":
        return (db["array"] == "12m") | (db["array"] == "7m")
    return db["array"] == array


def _as_value(col, unit: str):
    """Return column values converted to ``unit``.

    ``unit`` is a human-readable report label.  Columns that astropy can't
    interpret with a unit (e.g. sysperf, stored as raw PFLOP/s floats; or
    dimensionless channel counts) are returned as-is.
    """
    col_unit = getattr(col, "unit", None)
    raw = np.asarray(getattr(col, "value", col)).astype(float)

    if unit == "" or unit is None or col_unit is None:
        return raw

    target = u.Unit("Gvis/h") if unit == "Gvis/hr" else None
    try:
        if target is None:
            target = u.Unit(unit)
    except Exception:
        # Target not representable in astropy (e.g. PFLOP/s). Assume the
        # column is stored natively in the target unit already.
        return raw

    try:
        return col.to(target).value
    except Exception:
        return raw


def _compute_per_sample_stats(db: Table, col_name: str, unit: str
                              ) -> Dict[str, Dict[str, float]]:
    """Compute median/twa/max/total for each array mask on one realization.

    Returns: { array: { "median", "twa", "max", "total" } }
    """
    results: Dict[str, Dict[str, float]] = {}
    weights_all = np.asarray(db["weights_all"], dtype=float)
    values_full = _as_value(db[col_name], unit)

    for arr in _ARRAY_MASKS:
        mask = _array_mask(db, arr)
        values = values_full[mask]
        wts = weights_all[mask]
        # guard against empty / all-nan
        finite = np.isfinite(values)
        v_f = values[finite]
        w_f = wts[finite]
        stats = {}
        if len(v_f) == 0:
            stats = {"median": np.nan, "twa": np.nan,
                     "max": np.nan, "total": np.nan}
        else:
            stats["median"] = float(np.nanmedian(v_f))
            # Weighted average using the per-array weights (unnormalized over
            # the full database, so effectively a time-weighted average
            # restricted to the subset of rows).  Matches wsu_db code.
            w_sum = float(np.nansum(w_f))
            if w_sum > 0:
                stats["twa"] = float(np.nansum(v_f * w_f) / w_sum)
            else:
                stats["twa"] = float(np.nanmean(v_f))
            stats["max"] = float(np.nanmax(v_f))
            stats["total"] = float(np.nansum(v_f))
        results[arr] = stats
    return results


def compute_stats(paths: Iterable[str], cfg: PipelineConfig = DEFAULT_CONFIG
                  ) -> Dict:
    """Compute milestone stats by averaging over realization files.

    Output:
        stats[milestone][array][stat][quantity] = float
    where
        milestone : "M1" | "M4" | "M5" | "BLC"
        array     : "12m" | "7m" | "both"
        stat      : "median" | "twa" | "max" | "total" | "total_per_cycle"
    """
    register_gvis_unit()

    # milestone -> { column_attr: (actual_col, report_unit) }
    ms_map: Dict[str, MilestoneColumns] = {
        **{k: v for k, v in MILESTONE_COLUMNS.items()},
        "BLC": BLC_COLUMNS,
    }

    # Accumulator: per-sample arrays of values.
    # accum[ms][array][stat][quantity] = list of floats, one per realization.
    accum: Dict[str, Dict[str, Dict[str, Dict[str, List[float]]]]] = {
        ms: {a: {s: {q: [] for q in QUANTITY_SPEC}
                 for s in ("median", "twa", "max", "total")}
             for a in _ARRAY_MASKS}
        for ms in ms_map
    }

    paths = list(paths)
    if not paths:
        raise ValueError("no realization paths supplied to compute_stats")

    for i, path in enumerate(paths):
        log.info("Computing stats from %s (%d/%d)", path, i + 1, len(paths))
        db = Table.read(str(path), format="ascii.ecsv")
        for ms, cols in ms_map.items():
            for qname, qspec in QUANTITY_SPEC.items():
                col_attr = qspec["column"]
                report_unit = qspec["report_unit"]
                col_name = getattr(cols, col_attr)
                if col_name not in db.colnames:
                    continue
                per_arr = _compute_per_sample_stats(db, col_name, report_unit)
                for arr, stats in per_arr.items():
                    for s, v in stats.items():
                        accum[ms][arr][s][qname].append(v)

    n_real = len(paths)

    # Average the per-sample stats across realizations.
    final: Dict[str, Dict[str, Dict[str, Dict[str, float]]]] = {}
    for ms in accum:
        final[ms] = {}
        for arr in accum[ms]:
            final[ms][arr] = {}
            for stat in accum[ms][arr]:
                final[ms][arr][stat] = {}
                for q, arr_vals in accum[ms][arr][stat].items():
                    if not arr_vals:
                        continue
                    arr_vals_np = np.array(arr_vals, dtype=float)
                    final[ms][arr][stat][q] = float(np.nanmean(arr_vals_np))
            # "total_per_cycle": sum / 2 (database covers two ALMA cycles)
            # converted to PB for datavol_* and productsize.
            final[ms][arr]["total_per_cycle"] = {}
            for q in ("datavol_total", "datavol_science", "productsize"):
                totals = accum[ms][arr]["total"].get(q)
                if not totals:
                    continue
                # totals is in TB (report_unit); /2 -> per cycle; /1000 -> PB
                final[ms][arr]["total_per_cycle"][q] = \
                    float(np.nanmean(totals)) / 2.0 / 1000.0

    final["_metadata"] = {
        "n_realizations": n_real,
        "realization_paths": [str(p) for p in paths],
    }
    return final


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-path",
                        default="data/wsu_datarates_mit_per_mous_initial_goal_20250423.ecsv",
                        help="Path to per-MOUS database")
    parser.add_argument("--outdir",
                        default=None,
                        help="Directory for realization files (default from config)")
    parser.add_argument("--n-realizations", type=int, default=None,
                        help="Override number of realizations")
    parser.add_argument("--skip-realizations", action="store_true",
                        help="Reuse existing realization files if present")
    parser.add_argument("--stats-out",
                        default=None,
                        help="Optional path to write a stats JSON")
    parser.add_argument("--generate-tables",
                        default=None,
                        metavar="DIR",
                        help="Write memo-style and SDD-style LaTeX tables to DIR")
    parser.add_argument("--generate-plots",
                        default=None,
                        metavar="DIR",
                        help="Write SDD-style CCDF plots (png) to DIR")
    parser.add_argument("--no-plot-overlay", action="store_true",
                        help="Skip band-1/2 realization overlay on plots")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    cfg = DEFAULT_CONFIG
    if args.n_realizations is not None:
        from dataclasses import replace
        cfg = replace(cfg, realizations=replace(cfg.realizations,
                                                n_realizations=args.n_realizations))

    if args.skip_realizations:
        outdir = args.outdir or cfg.realizations.outdir
        paths = sorted(Path(outdir).glob(f"{cfg.realizations.filename}_*.ecsv"))
        paths = [str(p) for p in paths[: cfg.realizations.n_realizations]]
        if not paths:
            print(f"No realizations found in {outdir}", file=sys.stderr)
            return 2
    else:
        db = load_database(args.db_path)
        paths = generate_realizations(db, cfg, outdir=args.outdir)

    stats = compute_stats(paths, cfg)

    if args.stats_out:
        import json
        # drop non-serializable metadata arrays, keep scalars
        ser = {k: v for k, v in stats.items() if k != "_metadata"}
        ser["_metadata"] = stats["_metadata"]
        with open(args.stats_out, "w") as fh:
            json.dump(ser, fh, indent=2)

    if args.generate_tables:
        from tables import generate_memo_tables, generate_sdd_tables
        memo = generate_memo_tables(stats, args.generate_tables)
        sdd = generate_sdd_tables(stats, args.generate_tables)
        print(f"Wrote {len(memo) + len(sdd)} tables to {args.generate_tables}")

    if args.generate_plots:
        from plots import generate_sdd_plots
        db_for_plots = load_database(args.db_path)
        written = generate_sdd_plots(
            db_for_plots,
            output_dir=args.generate_plots,
            cfg=cfg,
            realization_dir=args.outdir,
            overlay_realizations=not args.no_plot_overlay,
        )
        print(f"Wrote {len(written)} plots to {args.generate_plots}")

    # Quick summary
    for ms in ("M1", "M4", "M5", "BLC"):
        if ms not in stats:
            continue
        for arr in _ARRAY_MASKS:
            twa = stats[ms][arr].get("twa", {})
            dr = twa.get("datarate")
            vr = twa.get("visrate")
            if dr is not None and vr is not None:
                print(f"{ms} {arr}: TWA datarate={dr:.4f} GB/s, "
                      f"visrate={vr:.2f} Gvis/hr")

    return 0


if __name__ == "__main__":
    sys.exit(main())
