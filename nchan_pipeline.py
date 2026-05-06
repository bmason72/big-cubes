#!/usr/bin/env python3
"""Top-level driver for the ryanCy11 SPW/nchan estimation workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

import ryan_cy11_wsu_projection
import spw_setup_summary
import spw_tabulate
from config import DEFAULT_CONFIG
from run_layout import next_run_dir, write_config_json


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        default="ryanCy11",
        help="Root directory containing one MOUS directory per entry",
    )
    parser.add_argument("--run-dir", default=None, help="Root output directory")
    parser.add_argument(
        "--validation-run",
        action="store_true",
        help="Use out/valid_nchanEst_<date><suffix>/ for default outputs",
    )
    parser.add_argument("--metadata-csv", default=None)
    parser.add_argument("--max-mous", type=int, default=None)
    parser.add_argument("--coverage-fraction", type=float, default=0.85)
    parser.add_argument(
        "--bandwidth-cut-mode",
        choices=("fixed", "percentile"),
        default="fixed",
    )
    parser.add_argument("--bandwidth-low-cut-mhz", type=float, default=90.0)
    parser.add_argument("--bandwidth-high-cut-mhz", type=float, default=1200.0)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    run_dir = Path(args.run_dir) if args.run_dir else next_run_dir(
        "nchanEst", validation=args.validation_run
    )
    tab_dir = run_dir / "tabulation"
    setup_dir = run_dir / "setup_summary"
    wsu_dir = run_dir / "wsu_projection"
    for directory in (run_dir, tab_dir, setup_dir, wsu_dir):
        directory.mkdir(parents=True, exist_ok=True)

    tab_args = [args.input_dir, str(tab_dir)]
    if args.metadata_csv:
        tab_args.extend(["--metadata-csv", args.metadata_csv])
    if args.max_mous is not None:
        tab_args.extend(["--max-mous", str(args.max_mous)])
    spw_tabulate.main(tab_args)

    setup_args = [
        str(tab_dir),
        str(setup_dir),
        "--coverage-fraction", str(args.coverage_fraction),
        "--bandwidth-cut-mode", args.bandwidth_cut_mode,
        "--bandwidth-low-cut-mhz", str(args.bandwidth_low_cut_mhz),
        "--bandwidth-high-cut-mhz", str(args.bandwidth_high_cut_mhz),
    ]
    spw_setup_summary.main(setup_args)

    ryan_cy11_wsu_projection.main([str(tab_dir), str(wsu_dir)])

    write_config_json(
        run_dir / "config.json",
        {
            "tool": "nchan_pipeline",
            "input_dir": args.input_dir,
            "run_dir": run_dir,
            "tabulation_dir": tab_dir,
            "setup_summary_dir": setup_dir,
            "wsu_projection_dir": wsu_dir,
            "validation_run": args.validation_run,
            "tabulation_args": tab_args,
            "setup_summary_args": setup_args,
            "wsu_projection_args": [str(tab_dir), str(wsu_dir)],
            "effective_config": DEFAULT_CONFIG,
            "argv": argv,
        },
    )

    print(f"Wrote outputs to {run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
