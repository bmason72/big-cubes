#!/usr/bin/env python3
"""Cy11 projected M1 cube sizes, prioritized-subset-imaging, and 100 GB cap.

Pipeline:
1. Per MOUS, take the representative-EB science SPW templates from the
   tuning-aware tabulation (``mous_tuning_spw_templates.csv``).
2. Synthesized beam per SPW from the listobs antenna geometry
   (theta = 0.574 lambda / L80) and primary beam / cell / imsize using the
   same conventions as ``large_cubes.munge_archive_info`` (5 pix/beam;
   single-field imsize = 1.539 * PB / cell; mosaic imsize =
   (mosaic FOV + 0.7 PB) / cell with hex-packed pointings at PB/2 spacing).
3. Project each current SPW to WSU M1 2-GHz SPWs (stepped2-binned
   resolution, M1 chanavg minima); each current template represents an equal
   share of the projected whole-SPW count (ceil convention).
4. Prioritized subset imaging (PSI): retained cube volume =
   line_fraction * cubesize (no floor); MFS images always kept.
5. 100 GB per-cube cap applied *after* PSI (cap, not discard).

Outputs per-MOUS and per-cube tables, aggregate product-volume factors
(PSI-only, cap-only, PSI+cap; split by long-baseline), CCDF plots, and a
markdown summary.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

from cy11_geometry import parse_listobs_geometry
from ryan_cy11_wsu_projection import resolve_local_path
from wsu_projection import (
    alma_band_from_freq_hz,
    normalize_mous_uid,
    projected_nchan_per_spw,
    projected_spw_equivalents,
    stepped2_floor_velocity,
)
from config import DEFAULT_CONFIG

C_M_S = 299792458.0
BEAM_L80_COEFF = 0.574          # theta = 0.574 lambda / L80   (ALMA Technical Handbook)
PIX_PER_BEAM = 3.5             # refined-scenario sampling: cubes are sized at 3.5 pix/beam so
                              # the 100 GB cap acts AFTER the 5->3.5 pixel reduction (the correct
                              # 3.5pix -> PSI -> cap order).  The 5->3.5 reduction itself is applied
                              # separately in refined_rollout_table.py (PIX_FACTOR); the PSI/cap
                              # factor produced here is a ratio (capped/uncapped at this sampling),
                              # which is pixel-scale invariant except through the cap nonlinearity,
                              # so it is not double counted.  Use 5.0 to reproduce pre-2026-06-25.
SF_IMSIZE_SCALE = 1.1 * (1.12 / 1.22) * math.sqrt(-math.log(0.2) / math.log(2.0))  # = 1.539
MOSAIC_PB_PAD = 0.70            # as in munge_archive_info
CUBE_CAP_GB = 100.0
MIN_IMSIZE_PIX = 64.0


def primary_beam_arcsec(array_type: str, freq_ghz: float) -> float:
    if array_type == "7m":
        return 33.3 * 300.0 / freq_ghz
    return 19.4 * 300.0 / freq_ghz


def synthesized_beam_arcsec(l80_m: float, freq_ghz: float) -> float:
    lam_m = C_M_S / (freq_ghz * 1.0e9)
    theta_rad = BEAM_L80_COEFF * lam_m / l80_m
    return math.degrees(theta_rad) * 3600.0


def imsize_pixels(array_type: str, freq_ghz: float, l80_m: float, n_pointings: int) -> tuple[float, float, float]:
    """Return (imsize, cell_arcsec, beam_arcsec) per DRRUP-style conventions."""
    beam = synthesized_beam_arcsec(l80_m, freq_ghz)
    cell = beam / PIX_PER_BEAM
    pb = primary_beam_arcsec(array_type, freq_ghz)
    if n_pointings <= 1:
        imsize = SF_IMSIZE_SCALE * pb / cell
    else:
        # hex-packed pointing centers at PB/2 spacing; equivalent circular
        # diameter of n pointing centers, plus the munge-style PB padding
        fov = 1.13 * (pb / 2.0) * math.sqrt(float(n_pointings))
        imsize = (fov + MOSAIC_PB_PAD * pb) / cell
    imsize = max(round(imsize, -1), MIN_IMSIZE_PIX)
    return imsize, cell, beam


def load_pipeline_stats(mous_dir: Path) -> tuple[int, list[int]]:
    """Return (n_target, list of n_pointings per target) from pipeline_stats."""
    stats_files = sorted(mous_dir.glob("pipeline_stats_*.json"))
    if not stats_files:
        return 1, [1]
    try:
        data = json.loads(stats_files[0].read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return 1, [1]
    key = next((k for k in data if k != "header"), None)
    if key is None:
        return 1, [1]
    top = data[key]
    targets = top.get("TARGET", {})
    pointings = []
    for tgt in targets.values():
        np_val = tgt.get("n_pointings", {}).get("value", 1)
        pointings.append(int(np_val) if np_val else 1)
    if not pointings:
        n_target = int(top.get("n_target", {}).get("value", 1) or 1)
        pointings = [1] * max(n_target, 1)
    return len(pointings), pointings


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tabulation-dir", default="out/valid_nchanEst_10jun26a/tabulation")
    parser.add_argument("--setup-summary-dir", default="out/valid_nchanEst_10jun26a/setup_summary")
    parser.add_argument("--input-dir", default="ryanCy11")
    parser.add_argument("--output-dir", default="out/psi_cap")
    parser.add_argument("--lb-threshold-m", type=float, default=3000.0)
    args = parser.parse_args(argv)

    tab_dir = Path(args.tabulation_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    input_dir = Path(args.input_dir)
    search_roots = [Path.cwd(), tab_dir.parent, Path(__file__).resolve().parent]

    # --- inputs -----------------------------------------------------------
    template_rows = list(csv.DictReader((tab_dir / "mous_tuning_spw_templates.csv").open()))
    by_mous: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in template_rows:
        by_mous[row["mous_uid"]].append(row)

    linefrac: dict[tuple[str, int], float] = {}
    lf_path = Path(args.setup_summary_dir) / "mous_spw_line_fraction.csv"
    lf_bucket_sums: dict[str, list[float]] = defaultdict(list)
    for row in csv.DictReader(lf_path.open()):
        try:
            value = float(row["line_fraction_percent"]) / 100.0
        except ValueError:
            continue
        if math.isfinite(value):
            value = min(max(value, 0.0), 1.0)
            linefrac[(row["mous_uid"], int(row["spw_id"]))] = value
            lf_bucket_sums[row["resolution_bucket"]].append(value)
    lf_bucket_mean = {
        bucket: sum(vals) / len(vals) for bucket, vals in lf_bucket_sums.items() if vals
    }
    lf_global_mean = (
        sum(v for vals in lf_bucket_sums.values() for v in vals)
        / max(sum(len(vals) for vals in lf_bucket_sums.values()), 1)
    )

    # --- per-MOUS cube computation ---------------------------------------
    geometry_cache: dict[str, object] = {}
    mous_rows: list[dict[str, object]] = []
    cube_rows: list[dict[str, object]] = []
    n_lf_fallback = 0
    n_geom_missing = 0

    for mous_uid, rows in sorted(by_mous.items()):
        listobs_path = rows[0]["listobs_path"]
        try:
            resolved = resolve_local_path(listobs_path, search_roots)
        except FileNotFoundError:
            n_geom_missing += 1
            continue
        key = str(resolved)
        if key not in geometry_cache:
            geometry_cache[key] = parse_listobs_geometry(resolved)
        geom = geometry_cache[key]
        if not math.isfinite(geom.l80_m) or geom.l80_m <= 0:
            n_geom_missing += 1
            continue
        array_type = geom.array_type
        is_lb = array_type == "12m" and geom.l80_m > args.lb_threshold_m

        mous_dir = input_dir / mous_uid
        n_target, pointings = load_pipeline_stats(mous_dir)

        band = alma_band_from_freq_hz(float(rows[0]["center_freq_hz"]))
        n_tunings = int(rows[0].get("n_tuning_groups_in_mous") or 1)
        total_bw_hz = float(rows[0].get("mous_total_science_bw_hz") or 0.0)
        milestone_bw_hz = DEFAULT_CONFIG.milestones["M1"].bandwidth_for(band) * 1.0e9
        projected_tunings = max(1, math.ceil(total_bw_hz / milestone_bw_hz)) if total_bw_hz > 0 else 1
        nspw_wsu = projected_tunings * projected_spw_equivalents("M1", band)
        spw_weight = nspw_wsu / len(rows)  # WSU SPWs represented by each current SPW

        vol_raw = vol_psi = vol_psi_cap = vol_cap_only = 0.0
        mfs_vol = 0.0
        max_cube = 0.0
        for n_point in pointings:
            for row in rows:
                freq_ghz = float(row["center_freq_hz"]) / 1.0e9
                velres = float(row["velocity_resolution_kms"])
                imsize, cell, beam = imsize_pixels(array_type, freq_ghz, geom.l80_m, n_point)
                nchan_wsu = projected_nchan_per_spw(
                    requested_vel_kms=stepped2_floor_velocity(velres),
                    original_vel_kms=velres,
                    center_freq_hz=float(row["center_freq_hz"]),
                    band=band,
                    milestone="M1",
                )
                cube_gb = 4.0 * imsize**2 * nchan_wsu / 1.0e9
                mfs_gb = 4.0 * imsize**2 / 1.0e9
                lf = linefrac.get((mous_uid, int(row["spw_id"])))
                if lf is None:
                    bucket = "CONT"  # fallback key resolved below
                    lf = lf_bucket_mean.get(bucket, lf_global_mean)
                    n_lf_fallback += 1
                cube_psi_gb = lf * cube_gb
                cube_psi_cap_gb = min(cube_psi_gb, CUBE_CAP_GB)
                cube_cap_only_gb = min(cube_gb, CUBE_CAP_GB)

                vol_raw += spw_weight * cube_gb
                vol_psi += spw_weight * cube_psi_gb
                vol_psi_cap += spw_weight * cube_psi_cap_gb
                vol_cap_only += spw_weight * cube_cap_only_gb
                mfs_vol += spw_weight * mfs_gb
                max_cube = max(max_cube, cube_gb)
                cube_rows.append(
                    {
                        "mous_uid": mous_uid,
                        "spw_id": row["spw_id"],
                        "array_type": array_type,
                        "is_long_baseline": is_lb,
                        "n_pointings": n_point,
                        "freq_ghz": round(freq_ghz, 3),
                        "l80_m": round(geom.l80_m, 1),
                        "beam_arcsec": round(beam, 4),
                        "imsize_pix": imsize,
                        "nchan_wsu_m1": nchan_wsu,
                        "line_fraction": round(lf, 4),
                        "wsu_spw_weight": round(spw_weight, 3),
                        "cube_gb": round(cube_gb, 3),
                        "cube_psi_gb": round(cube_psi_gb, 3),
                        "cube_psi_cap_gb": round(cube_psi_cap_gb, 3),
                    }
                )

        mous_rows.append(
            {
                "mous_uid": mous_uid,
                "array_type": array_type,
                "is_long_baseline": is_lb,
                "l80_m": round(geom.l80_m, 1),
                "n_target": n_target,
                "n_current_spw": len(rows),
                "nspw_wsu_m1": nspw_wsu,
                "max_cube_gb": round(max_cube, 2),
                "prodvol_raw_gb": 2.0 * (vol_raw + mfs_vol),
                "prodvol_psi_gb": 2.0 * (vol_psi + mfs_vol),
                "prodvol_psi_cap_gb": 2.0 * (vol_psi_cap + mfs_vol),
                "prodvol_cap_only_gb": 2.0 * (vol_cap_only + mfs_vol),
            }
        )

    # --- aggregate factors -------------------------------------------------
    def factor_rows(selector, label):
        sel = [m for m in mous_rows if selector(m)]
        raw = sum(m["prodvol_raw_gb"] for m in sel)
        psi = sum(m["prodvol_psi_gb"] for m in sel)
        psi_cap = sum(m["prodvol_psi_cap_gb"] for m in sel)
        cap_only = sum(m["prodvol_cap_only_gb"] for m in sel)
        return {
            "subset": label,
            "n_mous": len(sel),
            "raw_prodvol_tb": raw / 1.0e3,
            "f_psi": psi / raw if raw > 0 else math.nan,
            "f_cap_only": cap_only / raw if raw > 0 else math.nan,
            "f_psi_cap": psi_cap / raw if raw > 0 else math.nan,
        }

    factors = [
        factor_rows(lambda m: True, "all"),
        factor_rows(lambda m: m["is_long_baseline"], "long_baseline"),
        factor_rows(lambda m: not m["is_long_baseline"], "non_long_baseline"),
        factor_rows(lambda m: m["array_type"] == "12m", "12m"),
        factor_rows(lambda m: m["array_type"] == "7m", "7m"),
    ]
    f_m1 = factors[0]["f_psi_cap"]
    f_m1_psi = factors[0]["f_psi"]
    f_m4 = (1.0 + f_m1_psi) / 2.0  # PSI midpoint at M4; no cap at M4

    # --- outputs ------------------------------------------------------------
    def write_csv(path: Path, rows, fields):
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)

    write_csv(out_dir / "mous_psi_cap.csv", mous_rows, list(mous_rows[0].keys()))
    write_csv(out_dir / "cube_psi_cap.csv", cube_rows, list(cube_rows[0].keys()))
    write_csv(out_dir / "psi_cap_factors.csv", factors, list(factors[0].keys()))

    # --- plots ---------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        sizes_raw = np.array([c["cube_gb"] for c in cube_rows])
        sizes_psi = np.array([c["cube_psi_gb"] for c in cube_rows])
        sizes_cap = np.array([c["cube_psi_cap_gb"] for c in cube_rows])
        weights = np.array([c["wsu_spw_weight"] for c in cube_rows])

        fig, ax = plt.subplots(figsize=(7, 5))
        for vals, label, color in (
            (sizes_raw, "M1 cubes (unmitigated)", "C0"),
            (sizes_psi, "after PSI (line-fraction)", "C1"),
            (sizes_cap, "after PSI + 100 GB cap", "C2"),
        ):
            order = np.argsort(vals)
            v = vals[order]
            w = weights[order]
            ccdf = 1.0 - np.cumsum(w) / np.sum(w)
            ax.plot(v, ccdf, label=label, color=color)
        ax.axvline(CUBE_CAP_GB, color="k", ls=":", lw=1, label="100 GB cap")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("cube size [GB]")
        ax.set_ylabel("fraction of WSU SPW-cubes larger than size")
        ax.set_title("Cy11-projected M1 cube sizes (distributed_binned)")
        ax.legend(loc="lower left", fontsize=9)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(out_dir / "cube_size_ccdf.png", dpi=150)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(7, 5))
        lf_vals = np.array([c["line_fraction"] for c in cube_rows])
        ax.hist(lf_vals, bins=40, color="C0", alpha=0.8)
        ax.set_xlabel("line fraction (1 - line-free fraction)")
        ax.set_ylabel("EB-SPW count")
        ax.set_title("cont.dat line fractions used as PSI retention")
        fig.tight_layout()
        fig.savefig(out_dir / "line_fraction_hist.png", dpi=150)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(7, 5))
        l80s = np.array([m["l80_m"] for m in mous_rows if m["array_type"] == "12m"])
        ax.hist(np.log10(l80s), bins=40, color="C0", alpha=0.8)
        ax.axvline(math.log10(args.lb_threshold_m), color="k", ls=":", label="LB threshold (3 km)")
        ax.set_xlabel("log10 L80 [m] (12m MOUS)")
        ax.set_ylabel("MOUS count")
        ax.legend()
        fig.tight_layout()
        fig.savefig(out_dir / "l80_hist.png", dpi=150)
        plt.close(fig)
        plots_ok = True
    except ImportError:
        plots_ok = False

    # --- markdown summary -----------------------------------------------------
    lines = [
        "# Cy11 Prioritized Subset Imaging + 100 GB Cap Analysis (M1)",
        "",
        "Conventions: theta = 0.574 lambda/L80 (listobs antenna table); 5 pix/beam;",
        "imsize per munge_archive_info (SF: 1.539 PB/cell; mosaic: hex-packed pointings",
        "at PB/2 + 0.7 PB pad); WSU M1 SPWs at 2 GHz, ceil SPW count, stepped2-binned",
        "resolution with M1 chanavg minima; PSI retention = cont.dat line fraction",
        "(no floor); MFS always kept; cap = min(cube, 100 GB) after PSI.",
        "",
        f"- MOUS used: {len(mous_rows)} (geometry missing/invalid: {n_geom_missing})",
        f"- line-fraction fallbacks (no cont.dat match): {n_lf_fallback} "
        f"of {len(cube_rows)} cube rows (sample-mean fallback {lf_global_mean:.3f})",
        "",
        "| Subset | N MOUS | raw prod vol [TB] | f(PSI) | f(cap only) | f(PSI+cap) |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for f in factors:
        lines.append(
            f"| {f['subset']} | {f['n_mous']} | {f['raw_prodvol_tb']:.1f} "
            f"| {f['f_psi']:.4f} | {f['f_cap_only']:.4f} | {f['f_psi_cap']:.4f} |"
        )
    lines.extend(
        [
            "",
            f"**Adopted M1 product factor (PSI then cap): {f_m1:.4f}**",
            f"**Adopted M4 product factor (PSI midpoint (1+f_PSI)/2, no cap): {f_m4:.4f}**",
            "**M5: 1.0 (no PSI, no cap).**",
            "",
            "Note: PSI retention uses each current SPW's cont.dat line fraction applied",
            "to its projected 2 GHz WSU descendants. For narrow current SPWs this is",
            "conservative (lines fill a smaller fraction of a 2 GHz window than of the",
            "narrow window they were targeted with).",
        ]
    )
    if plots_ok:
        lines.append("")
        lines.append("Plots: cube_size_ccdf.png, line_fraction_hist.png, l80_hist.png")
    (out_dir / "psi_cap_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print((out_dir / "psi_cap_summary.md").read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
