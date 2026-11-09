#!/usr/bin/env python3
"""Cross-check of the 100 GB cube-cap factor against the Cy7/8 database.

Computes, per cycle, the M1 product-volume factor from capping each
projected WSU M1 cube at 100 GB (cap, not discard), using the DB's stored
``imsize``, ``wsu_nchan_spw_stepped2_initial``, ``wsu_nspw_initial`` and
``ntarget``.  No line fractions exist for Cy7/8 (no cont.dat), so this
checks the cap-only factor and its long-baseline-cycle sensitivity.
"""

from __future__ import annotations

import argparse
import csv
import math
import warnings
from pathlib import Path

# The DB imsize is at 5 pix/beam; the refined scenario images at 3.5 pix/beam,
# so a 3.5-pixel cube capped at 100 GB is equivalent to capping the 5-pixel
# cube at 100/(3.5/5)^2 ~= 204 GB.  Capping the stored 5-pixel cubes at that
# threshold reproduces the "3.5 pix then 100 GB cap" scenario.
CUBE_CAP_GB = 100.0 / (3.5 / 5.0) ** 2


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--db-path",
        default="data/wsu_datarates_mit_per_mous_initial_goal_20250423.ecsv",
    )
    parser.add_argument("--output-dir", default="out/psi_cap")
    parser.add_argument("--lb-threshold-m", type=float, default=3000.0)
    args = parser.parse_args(argv)

    from astropy.table import Table

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        db = Table.read(args.db_path, format="ascii.ecsv")

    out_rows = []
    cycles = sorted(set(str(c) for c in db["cycle_info"]))
    subsets = [("all", lambda i: True)]
    subsets += [(c, lambda i, c=c: str(db["cycle_info"][i]) == c) for c in cycles]
    subsets += [
        ("long_baseline", lambda i: str(db["array"][i]) == "12m"
            and float(db["L80"][i]) > args.lb_threshold_m),
        ("non_long_baseline", lambda i: not (str(db["array"][i]) == "12m"
            and float(db["L80"][i]) > args.lb_threshold_m)),
    ]

    imsize = db["imsize"]
    nchan_spw = db["wsu_nchan_spw_stepped2_initial"]
    nspw = db["wsu_nspw_initial"]
    ntarget = db["ntarget"]

    for label, selector in subsets:
        raw = capped = 0.0
        n = 0
        for i in range(len(db)):
            if not selector(i):
                continue
            try:
                im = float(imsize[i])
                nc = float(nchan_spw[i])
                ns = float(nspw[i])
                nt = float(ntarget[i])
            except (TypeError, ValueError):
                continue
            if not all(map(math.isfinite, (im, nc, ns, nt))) or min(im, nc, ns, nt) <= 0:
                continue
            cube_gb = 4.0 * im * im * nc / 1.0e9
            mfs_gb = 4.0 * im * im / 1.0e9
            raw += 2.0 * (cube_gb + mfs_gb) * ns * nt
            capped += 2.0 * (min(cube_gb, CUBE_CAP_GB) + mfs_gb) * ns * nt
            n += 1
        out_rows.append(
            {
                "subset": label,
                "n_mous": n,
                "raw_prodvol_tb": raw / 1.0e3,
                "f_cap_only": capped / raw if raw > 0 else math.nan,
            }
        )

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "cy78_cap_crosscheck.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)

    lines = [
        "# Cy7/8 DB 100 GB cap cross-check (M1, stored stepped2-initial projection)",
        "",
        "| Subset | N MOUS | raw prod vol [TB] | f(cap only) |",
        "| --- | --- | --- | --- |",
    ]
    for row in out_rows:
        lines.append(
            f"| {row['subset']} | {row['n_mous']} | {row['raw_prodvol_tb']:.1f} "
            f"| {row['f_cap_only']:.4f} |"
        )
    (out_dir / "cy78_cap_crosscheck.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print((out_dir / "cy78_cap_crosscheck.md").read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
