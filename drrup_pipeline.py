#!/usr/bin/env python3
"""Top-level driver for the main WSU data-properties pipeline."""

from __future__ import annotations

import argparse
import logging
from dataclasses import replace
from pathlib import Path

from config import DEFAULT_CONFIG
from pipeline import compute_stats, find_realizations, generate_realizations, load_database
from plots import generate_sdd_plots
from run_layout import next_run_dir, write_config_json
from tables import generate_memo_tables, generate_sdd_tables, write_all_tables_document


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--db-path",
        default="data/wsu_datarates_mit_per_mous_initial_goal_20250423.ecsv",
        help="Path to per-MOUS database",
    )
    parser.add_argument("--run-dir", default=None, help="Root output directory")
    parser.add_argument(
        "--validation-run",
        action="store_true",
        help="Use out/valid_drrupPl_<date><suffix>/ for default outputs",
    )
    parser.add_argument("--n-realizations", type=int, default=None)
    parser.add_argument(
        "--skip-realizations",
        action="store_true",
        help="Reuse an existing realization directory instead of generating new ones",
    )
    parser.add_argument(
        "--realization-dir",
        default=None,
        help="Existing or target realization directory (default: <run-dir>/realizations)",
    )
    parser.add_argument(
        "--no-plot-overlay",
        action="store_true",
        help="Skip band-1/2 realization overlay on plots",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    cfg = DEFAULT_CONFIG
    if args.n_realizations is not None:
        cfg = replace(
            DEFAULT_CONFIG,
            realizations=replace(DEFAULT_CONFIG.realizations, n_realizations=args.n_realizations),
        )

    run_dir = Path(args.run_dir) if args.run_dir else next_run_dir(
        "drrupPl", validation=args.validation_run
    )
    real_dir = Path(args.realization_dir) if args.realization_dir else run_dir / "realizations"
    tables_dir = run_dir / "tables"
    plots_dir = run_dir / "plots"
    stats_path = run_dir / "stats.json"
    run_dir.mkdir(parents=True, exist_ok=True)
    real_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    if args.skip_realizations:
        paths = find_realizations(cfg, real_dir)
    else:
        db = load_database(args.db_path)
        paths = generate_realizations(db, cfg, outdir=str(real_dir))

    stats = compute_stats(paths, cfg)

    import json

    ser = {k: v for k, v in stats.items() if k != "_metadata"}
    ser["_metadata"] = stats["_metadata"]
    stats_path.write_text(json.dumps(ser, indent=2) + "\n", encoding="utf-8")

    memo = generate_memo_tables(stats, str(tables_dir))
    sdd = generate_sdd_tables(stats, str(tables_dir))
    write_all_tables_document(memo + sdd, tables_dir / "all_tables.tex")

    db_for_plots = load_database(args.db_path)
    generate_sdd_plots(
        db_for_plots,
        output_dir=str(plots_dir),
        cfg=cfg,
        realization_dir=str(real_dir),
        overlay_realizations=not args.no_plot_overlay,
    )

    write_config_json(
        run_dir / "config.json",
        {
            "tool": "drrup_pipeline",
            "db_path": args.db_path,
            "run_dir": run_dir,
            "realization_dir": real_dir,
            "tables_dir": tables_dir,
            "plots_dir": plots_dir,
            "stats_path": stats_path,
            "skip_realizations": args.skip_realizations,
            "validation_run": args.validation_run,
            "effective_config": cfg,
            "argv": argv,
        },
    )

    print(f"Wrote outputs to {run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
