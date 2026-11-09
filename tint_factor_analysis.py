#!/usr/bin/env python3
"""Visibility-volume factor for "6.144 s except long baselines" at M4/M5.

The PDR baseline assumes all 12m WSU observations at M4/M5 use 3.072 s
integrations.  Current practice is ~6 s except in long-baseline (LB)
configurations; this script computes the visibility-volume factor implied by
keeping that practice (12m non-LB EBs at 6.144 s = 0.5x volume; 12m LB and
all 7m unchanged), from two independent samples:

1. The Cy7/8 per-MOUS database (one LB cycle + one no-LB cycle), weighting
   each MOUS by its projected per-milestone WSU visibility volume
   (``wsu_datavol_<ms>_stepped2_typical_total``) and classifying LB by
   ``L80 > threshold``.
2. The Cy11 EB sample, weighting each EB by its projected
   distributed_binned volume from ``eb_projection_factors.csv`` and
   classifying LB from the listobs antenna table (cy11_geometry).

Outputs a per-cycle and combined factor table (CSV + markdown).
"""

from __future__ import annotations

import argparse
import csv
import math
import warnings
from collections import defaultdict
from pathlib import Path

LB_L80_THRESHOLD_M = 3000.0  # matches blc_tint_breakpt_m (current 3s/6s practice split)

DB_MILESTONE_COLUMNS = {
    "M1": "wsu_datavol_initial_stepped2_initial_total",
    "M4": "wsu_datavol_ms4_stepped2_typical_total",
    "M5": "wsu_datavol_goal_stepped2_typical_total",
}


def db_tint_factors(db_path: Path) -> list[dict[str, object]]:
    from astropy.table import Table

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        db = Table.read(db_path, format="ascii.ecsv")

    rows: list[dict[str, object]] = []
    cycles = sorted(set(str(c) for c in db["cycle_info"]))
    milestone_cols = {
        ms: col for ms, col in DB_MILESTONE_COLUMNS.items() if col in db.colnames
    }
    for ms, col in milestone_cols.items():
        per_cycle_factors = {}
        for cycle in cycles + ["all"]:
            sel = (
                [str(c) == cycle for c in db["cycle_info"]]
                if cycle != "all"
                else [True] * len(db)
            )
            total = 0.0
            weighted = 0.0
            lb_vol = 0.0
            n = n_lb = 0
            for i, keep in enumerate(sel):
                if not keep:
                    continue
                vol = float(db[col][i])
                if not math.isfinite(vol) or vol <= 0:
                    continue
                is_12m = str(db["array"][i]) == "12m"
                l80 = float(db["L80"][i])
                is_lb = is_12m and math.isfinite(l80) and l80 > LB_L80_THRESHOLD_M
                # M1 already assumes 6.144 s for all 12m EBs: factor 1
                if ms == "M1":
                    w = 1.0
                else:
                    w = 1.0 if (is_lb or not is_12m) else 0.5
                total += vol
                weighted += vol * w
                n += 1
                if is_lb:
                    n_lb += 1
                    lb_vol += vol
            factor = weighted / total if total > 0 else math.nan
            per_cycle_factors[cycle] = factor
            rows.append(
                {
                    "sample": f"cy78db_{cycle}",
                    "milestone": ms,
                    "n": n,
                    "n_lb": n_lb,
                    "lb_volume_fraction": lb_vol / total if total > 0 else math.nan,
                    "tint_volume_factor": factor,
                }
            )
        cycle_only = [per_cycle_factors[c] for c in cycles]
        rows.append(
            {
                "sample": "cy78db_cycle_mean",
                "milestone": ms,
                "n": "",
                "n_lb": "",
                "lb_volume_fraction": "",
                "tint_volume_factor": sum(cycle_only) / len(cycle_only),
            }
        )
    return rows


def cy11_tint_factors(
    projection_csv: Path,
    tabulation_dir: Path,
    search_roots: list[Path],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    from cy11_geometry import parse_listobs_geometry
    from ryan_cy11_wsu_projection import resolve_local_path

    proj_rows = list(csv.DictReader(projection_csv.open()))
    geometry_cache: dict[str, object] = {}
    eb_geometry_rows: list[dict[str, object]] = []

    # Per-MOUS aggregation: LB classified at MOUS level (max L80 over its
    # EBs), volumes summed over EBs, so counts are MOUS counts as in the
    # Cy7/8 DB rows.
    mous_l80: dict[str, float] = defaultdict(float)
    mous_array: dict[str, str] = {}
    mous_vols: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    eb_records: list[tuple[str, dict[str, str], object]] = []
    for row in proj_rows:
        path = resolve_local_path(row["listobs_path"], search_roots)
        key = str(path)
        if key not in geometry_cache:
            geometry_cache[key] = parse_listobs_geometry(path)
        geom = geometry_cache[key]
        mous_uid = row["mous_uid"]
        if math.isfinite(geom.l80_m):
            mous_l80[mous_uid] = max(mous_l80[mous_uid], geom.l80_m)
        mous_array[mous_uid] = row["array_type"]
        for ms in ("M1", "M4", "M5"):
            vol = float(row[f"{ms.lower()}_distributed_binned_volume_gb"])
            if math.isfinite(vol) and vol > 0:
                mous_vols[mous_uid][ms] += vol
        eb_records.append((mous_uid, row, geom))

    rows: list[dict[str, object]] = []
    for ms in ("M1", "M4", "M5"):
        total = weighted = lb_vol = 0.0
        n = n_lb = 0
        for mous_uid, vols in mous_vols.items():
            vol = vols.get(ms, 0.0)
            if vol <= 0:
                continue
            is_12m = mous_array[mous_uid] == "12m"
            is_lb = is_12m and mous_l80[mous_uid] > LB_L80_THRESHOLD_M
            if ms == "M1":
                w = 1.0
            else:
                w = 1.0 if (is_lb or not is_12m) else 0.5
            total += vol
            weighted += vol * w
            n += 1
            if is_lb:
                n_lb += 1
                lb_vol += vol
        factor = weighted / total if total > 0 else math.nan
        rows.append(
            {
                "sample": "cy11",
                "milestone": ms,
                "n": n,
                "n_lb": n_lb,
                "lb_volume_fraction": lb_vol / total if total > 0 else math.nan,
                "tint_volume_factor": factor,
            }
        )

    for mous_uid, row, geom in eb_records:
        is_12m = row["array_type"] == "12m"
        eb_geometry_rows.append(
            {
                "mous_uid": mous_uid,
                "eb_uid": row["eb_uid"],
                "array_type": row["array_type"],
                "l80_m": geom.l80_m,
                "lmax_m": geom.lmax_m,
                "n_ant": geom.n_ant,
                "is_long_baseline": is_12m and mous_l80[mous_uid] > LB_L80_THRESHOLD_M,
                "m4_distributed_binned_volume_gb": float(
                    row["m4_distributed_binned_volume_gb"]
                ),
            }
        )
    return rows, eb_geometry_rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--db-path",
        default="data/wsu_datarates_mit_per_mous_initial_goal_20250423.ecsv",
    )
    parser.add_argument(
        "--projection-csv",
        default="out/valid_nchanEst_10jun26a/wsu_projection/eb_projection_factors.csv",
    )
    parser.add_argument(
        "--tabulation-dir",
        default="out/valid_nchanEst_10jun26a/tabulation",
    )
    parser.add_argument("--output-dir", default="out/tint_factor")
    args = parser.parse_args(argv)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = db_tint_factors(Path(args.db_path))
    cy11_rows, eb_geometry_rows = cy11_tint_factors(
        Path(args.projection_csv),
        Path(args.tabulation_dir),
        [Path.cwd(), Path(args.tabulation_dir).parent, Path(__file__).resolve().parent],
    )
    rows.extend(cy11_rows)

    fields = ["sample", "milestone", "n", "n_lb", "lb_volume_fraction", "tint_volume_factor"]
    with (out_dir / "tint_factors.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    geo_fields = [
        "mous_uid", "eb_uid", "array_type", "l80_m", "lmax_m", "n_ant",
        "is_long_baseline", "m4_distributed_binned_volume_gb",
    ]
    with (out_dir / "cy11_eb_geometry.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=geo_fields)
        writer.writeheader()
        writer.writerows(eb_geometry_rows)

    lines = [
        "# Tint Factor Analysis (6.144 s except long baselines)",
        "",
        f"- LB classification: 12m array with L80 > {LB_L80_THRESHOLD_M:.0f} m",
        "- Factor = (vis volume with 12m non-LB at 6.144 s) / (vis volume all-12m at 3.072 s)",
        "- M1 factor is 1.0 by construction (M1 already assumes 6.144 s).",
        "",
        "| Sample | Milestone | N | N(LB) | LB vol frac | Tint factor |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        lbf = row["lb_volume_fraction"]
        lbf_s = f"{lbf:.3f}" if isinstance(lbf, float) else ""
        lines.append(
            f"| {row['sample']} | {row['milestone']} | {row['n']} | {row['n_lb']} "
            f"| {lbf_s} | {float(row['tint_volume_factor']):.4f} |"
        )
    (out_dir / "tint_factors.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print((out_dir / "tint_factors.md").read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
