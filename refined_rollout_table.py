#!/usr/bin/env python3
"""Build the refined rollout projection table from measured factors.

Reads the measured factor analyses (tint_factor_analysis,
cy11_psi_cap_analysis, nchan summaries) and applies them as multiplicative
factors to the PDR rollout baseline, emitting a LaTeX table and a markdown
derivation trail.

Factor chain (on the PDR baseline):
- visibility side: tint factor (M4/M5; two-cycle Cy7/8 mean) x empirical
  Nchan factor (per milestone, aggregate-channel ratio
  distributed_binned / memo_uniform_binned).
- product side: (3.5/5)^2 pixel factor x empirical Nchan factor.
- mitigation: PSI+100GB-cap factor at M1 (two-cycle average: Cy11 all =
  LB cycle, Cy11 non-LB subset = no-LB-cycle proxy); PSI midpoint
  (1+f_PSI)/2 at M4 (no cap); none at M5.
- system performance and M4/M5 cores scale with the visibility factor;
  M1 cores are bottoms-up against the updated current product volume.

Updated current (2026) anchors from F. Stoehr: raw 157.7 TB/yr, mitigated
compressed products 109 TB/yr (sum matches the measured 267 TB/yr total).
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

# --- PDR rollout baseline (memo Table "rollout-pdr") -----------------------
PDR = {
    # quantity: (2026, M1, M4, M5, unit)
    "raw_vol":        (0.099, 0.613, 2.06, 5.00, "PB/yr"),
    "raw_rate":       (0.010, 0.0660, 0.221, 0.540, "GB/s"),
    "peak_rate":      (None, 0.218, 1.74, 3.48, "GB/s"),
    "max_rate":       (None, 0.309, 2.24, 4.43, "GB/s"),
    "prod_vol":       (0.61, 4.50, 5.07, 12.1, "PB/yr"),
    "prod_compr":     (0.40, 2.92, 3.30, 7.87, "PB/yr"),
    "mitig_compr":    (0.070, 0.585, 3.30, 7.87, "PB/yr"),
    "mitig_gous":     (None, 0.585, 4.48, 10.7, "PB/yr"),
    "total_vol":      (0.17, 1.20, 6.55, 15.7, "PB/yr"),
    "sysperf":        (0.0065, 0.00990, 0.263, 0.570, "PFLOP/s"),
    "cores":          (384, 3080, 7160, 15500, "cores"),
}

ROW_LABELS = {
    "raw_vol": "Raw data vol (12m+7m+TP)",
    "raw_rate": "Time averaged raw data rate (12m+7m+TP)",
    "peak_rate": "Expected peak data rate (not limiting)",
    "max_rate": "Max data rate (12m+7m+TP, Sys-Req-610)",
    "prod_vol": "Prod Vol (12m+7m+TP)",
    "prod_compr": "Prod vol (compressed)",
    "mitig_compr": "Mitigated prod vol (compressed)",
    "mitig_gous": "Mitig.\\ prod vol with GOUS",
    "total_vol": "Total data volume",
    "sysperf": "Estimated system performance",
    "cores": "Cores",
}

PIX_FACTOR = (3.5 / 5.0) ** 2

# Updated current-state anchors (F. Stoehr, 2026)
STOEHR_RAW_PB = 0.1577
STOEHR_MITIG_COMPR_PB = 0.109
CURRENT_CORES = 384


def read_tint_factor(path: Path) -> dict[str, float]:
    rows = list(csv.DictReader(path.open()))
    out = {}
    for row in rows:
        if row["sample"] == "cy78db_cycle_mean":
            out[row["milestone"]] = float(row["tint_volume_factor"])
    return out


def read_nchan_factors(path: Path) -> dict[str, float]:
    rows = list(csv.DictReader(path.open()))
    by = {(r["milestone"], r["method"]): float(r["projected_total_nchan_agg"]) for r in rows}
    return {
        ms: by[(ms, "distributed_binned")] / by[(ms, "memo_uniform_binned")]
        for ms in ("M1", "M4", "M5")
    }


def read_calcap_factors(path: Path) -> dict[str, float]:
    """Measured 1 km/s calibrator-resolution-cap factor (capped/uncapped
    projected visibility volume over the science+bandpass+phase+check
    buckets), distributed_binned scenario."""
    rows = list(csv.DictReader(path.open()))
    out = {}
    for row in rows:
        if row["method"] == "distributed_binned":
            out[row["milestone"]] = (
                float(row["projected_selected_total_volume_gb"])
                / float(row["projected_selected_total_uncapped_volume_gb"])
            )
    return out


def read_psi_factors(path: Path) -> dict[str, float]:
    rows = {r["subset"]: r for r in csv.DictReader(path.open())}
    f_all_psicap = float(rows["all"]["f_psi_cap"])
    f_nlb_psicap = float(rows["non_long_baseline"]["f_psi_cap"])
    f_all_psi = float(rows["all"]["f_psi"])
    f_nlb_psi = float(rows["non_long_baseline"]["f_psi"])
    # two-cycle average: Cy11 (LB cycle) + non-LB subset (no-LB-cycle proxy)
    f_m1 = (f_all_psicap + f_nlb_psicap) / 2.0
    f_psi_avg = (f_all_psi + f_nlb_psi) / 2.0
    return {
        "m1_psi_cap": f_m1,
        "psi_avg": f_psi_avg,
        "m4_psi": (1.0 + f_psi_avg) / 2.0,
        "cy11_lb_cycle_psi_cap": f_all_psicap,
        "no_lb_proxy_psi_cap": f_nlb_psicap,
    }


def fmt(value, sig=3):
    if value is None:
        return "--"
    if value == 0:
        return "0"
    if abs(value) >= 100:
        return f"{value:.0f}"
    return f"{value:.{sig}g}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tint-csv", default="out/tint_factor_10jun26a/tint_factors.csv")
    parser.add_argument(
        "--nchan-csv",
        default="out/valid_nchanEst_10jun26a/wsu_projection/summary_overall_nchan.csv",
    )
    parser.add_argument("--psi-csv", default="out/psi_cap_10jun26a/psi_cap_factors.csv")
    parser.add_argument(
        "--calcap-csv",
        default="out/calcap_10jun26a/summary_overall_intent_aware.csv",
    )
    parser.add_argument("--output-dir", default="out/refined_rollout")
    args = parser.parse_args(argv)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tint = read_tint_factor(Path(args.tint_csv))
    nchan = read_nchan_factors(Path(args.nchan_csv))
    psi = read_psi_factors(Path(args.psi_csv))
    calcap = read_calcap_factors(Path(args.calcap_csv))

    vis = {ms: tint[ms] * nchan[ms] * calcap[ms] for ms in ("M1", "M4", "M5")}
    prod = {ms: PIX_FACTOR * nchan[ms] for ms in ("M1", "M4", "M5")}
    mitig = {"M1": psi["m1_psi_cap"], "M4": psi["m4_psi"], "M5": 1.0}

    cols = {"2026": 0, "M1": 1, "M4": 2, "M5": 3}
    refined: dict[str, dict[str, float | None]] = {q: {} for q in PDR}

    gous_ratio = {
        ms: (PDR["mitig_gous"][cols[ms]] / PDR["mitig_compr"][cols[ms]])
        for ms in ("M1", "M4", "M5")
    }

    for ms in ("M1", "M4", "M5"):
        i = cols[ms]
        refined["raw_vol"][ms] = PDR["raw_vol"][i] * vis[ms]
        refined["raw_rate"][ms] = PDR["raw_rate"][i] * vis[ms]
        refined["peak_rate"][ms] = PDR["peak_rate"][i]
        refined["max_rate"][ms] = PDR["max_rate"][i]
        refined["prod_vol"][ms] = PDR["prod_vol"][i] * prod[ms]
        refined["prod_compr"][ms] = PDR["prod_compr"][i] * prod[ms]
        refined["mitig_compr"][ms] = refined["prod_compr"][ms] * mitig[ms]
        refined["mitig_gous"][ms] = refined["mitig_compr"][ms] * gous_ratio[ms]
        refined["total_vol"][ms] = refined["raw_vol"][ms] + refined["mitig_gous"][ms]
        if ms == "M1":
            refined["sysperf"][ms] = None  # bottoms-up milestone
            refined["cores"][ms] = CURRENT_CORES * (
                refined["mitig_compr"][ms] / STOEHR_MITIG_COMPR_PB
            )
        else:
            refined["sysperf"][ms] = PDR["sysperf"][i] * vis[ms]
            refined["cores"][ms] = PDR["cores"][i] * vis[ms]

    current = {
        "raw_vol": STOEHR_RAW_PB,
        # PDR's own rate/volume convention (effective observing time):
        # 0.010 GB/s per 0.099 PB/yr in the 2026 column
        "raw_rate": STOEHR_RAW_PB * (PDR["raw_rate"][0] / PDR["raw_vol"][0]),
        "peak_rate": None,
        "max_rate": None,
        "prod_vol": PDR["prod_vol"][0],
        "prod_compr": PDR["prod_compr"][0],
        "mitig_compr": STOEHR_MITIG_COMPR_PB,
        "mitig_gous": None,
        "total_vol": STOEHR_RAW_PB + STOEHR_MITIG_COMPR_PB,
        "sysperf": PDR["sysperf"][0],
        "cores": CURRENT_CORES,
    }

    # ---- LaTeX table -------------------------------------------------------
    lines = [
        r"\begin{table}[p]",
        r"\centering",
        r"\caption{Refined projection: PDR baseline scaled by the measured "
        r"integration-time, beam-sampling, empirical-Nchan, and PSI/cap factors. "
        r"The 2026 raw and mitigated-product anchors are the updated measured "
        r"values (F.\ Stoehr); 2026 unmitigated product volumes retain the PDR "
        r"estimates.}",
        r"\label{tab:rollout-refined}",
        r"\small",
        r"\setlength{\tabcolsep}{4pt}",
        r"\begin{tabular}{>{\raggedright\arraybackslash}p{2.9in}ccccc}",
        r"\toprule",
        r"Quantity & 2026 & M1 & M4 & M5 & Units \\",
        r"\midrule",
    ]
    for q in PDR:
        unit = PDR[q][4]
        vals = [current[q]] + [refined[q][ms] for ms in ("M1", "M4", "M5")]
        cells = " & ".join(fmt(v) for v in vals)
        lines.append(f"{ROW_LABELS[q]} & {cells} & {unit} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    (out_dir / "rollout_refined.tex").write_text("\n".join(lines), encoding="utf-8")

    # ---- factor table (memo Table 'mitig' replacement) ----------------------
    flines = [
        r"\begin{table}[h]",
        r"\centering",
        r"\small",
        r"\setlength{\tabcolsep}{4pt}",
        r"\begin{tabular}{>{\raggedright\arraybackslash}p{2.6in}cccll}",
        r"\toprule",
        r"Refinement / mitigation & M1 & M4 & M5 & Impacts vis? & Impacts prod? \\",
        r"\midrule",
        rf"Tint = 6\,s except long baseline (measured) & 1.00 & {tint['M4']:.3f} & {tint['M5']:.3f} & yes & no \\",
        rf"5 pix/beam $\rightarrow$ 3.5 pix/beam & {PIX_FACTOR:.3f} & {PIX_FACTOR:.3f} & {PIX_FACTOR:.3f} & no & yes \\",
        rf"Empirical Nchan projection (measured) & {nchan['M1']:.3f} & {nchan['M4']:.3f} & {nchan['M5']:.3f} & yes & yes \\",
        rf"PSI + 100 GB cap (measured; replaces assumed mitigation) & {mitig['M1']:.3f} & {mitig['M4']:.3f} & 1.00 & no & yes \\",
        rf"Cap calibrators at 1 km/s (measured) & {calcap['M1']:.3f} & {calcap['M4']:.3f} & {calcap['M5']:.3f} & yes & no \\",
        r"\midrule",
        rf"Net visibility factor & {vis['M1']:.3f} & {vis['M4']:.3f} & {vis['M5']:.3f} & & \\",
        rf"Net product factor (before mitigation) & {prod['M1']:.3f} & {prod['M4']:.3f} & {prod['M5']:.3f} & & \\",
        r"\bottomrule",
        r"\end{tabular}",
        r"\caption{Measured refinement and mitigation factors. "
        r"Tint: two-cycle Cy7/8 mean (LB fraction by projected volume); "
        r"Nchan: aggregate-channel ratio distributed\_binned / memo\_uniform\_binned (Cy11); "
        r"PSI+cap at M1: two-cycle average of the Cy11 line-fraction + 100 GB cap factor; "
        r"M4 uses the PSI midpoint $(1+f_{\rm PSI})/2$ without the cap; "
        r"calibrator cap: capped/uncapped projected visibility volume from the "
        r"intent-aware Cy11 projection.}",
        r"\label{tab:mitig-measured}",
        r"\end{table}",
        "",
    ]
    (out_dir / "factors_measured.tex").write_text("\n".join(flines), encoding="utf-8")

    # ---- markdown derivation trail ------------------------------------------
    md = [
        "# Refined Rollout Table Derivation",
        "",
        "## Measured input factors",
        "",
        f"- Tint (6.144 s except LB; two-cycle Cy7/8 mean): M1 1.0, M4 {tint['M4']:.4f}, M5 {tint['M5']:.4f}",
        f"- Calibrator 1 km/s cap (intent-aware Cy11, distributed_binned): "
        f"M1 {calcap['M1']:.4f}, M4 {calcap['M4']:.4f}, M5 {calcap['M5']:.4f}",
        f"- Empirical Nchan (distributed_binned / memo_uniform_binned, aggregate channels): "
        f"M1 {nchan['M1']:.4f}, M4 {nchan['M4']:.4f}, M5 {nchan['M5']:.4f}",
        f"- Pixel sampling (3.5/5)^2 = {PIX_FACTOR:.4f}",
        f"- PSI+cap M1: Cy11 LB-cycle {psi['cy11_lb_cycle_psi_cap']:.4f}, "
        f"no-LB-cycle proxy {psi['no_lb_proxy_psi_cap']:.4f}, two-cycle mean {psi['m1_psi_cap']:.4f}",
        f"- PSI M4 midpoint: (1 + {psi['psi_avg']:.4f})/2 = {psi['m4_psi']:.4f}; M5 1.0",
        "",
        "## Net chain factors (vs PDR baseline)",
        "",
        f"- Visibility: M1 {vis['M1']:.4f}, M4 {vis['M4']:.4f}, M5 {vis['M5']:.4f}",
        f"- Product (pre-mitigation): M1 {prod['M1']:.4f}, M4 {prod['M4']:.4f}, M5 {prod['M5']:.4f}",
        f"- Mitigation: M1 {mitig['M1']:.4f}, M4 {mitig['M4']:.4f}, M5 1.0",
        "",
        "## 2026 anchors (F. Stoehr)",
        "",
        f"- Raw {STOEHR_RAW_PB} PB/yr; mitigated compressed products {STOEHR_MITIG_COMPR_PB} PB/yr "
        f"(sum {STOEHR_RAW_PB + STOEHR_MITIG_COMPR_PB:.3f} PB/yr = measured 267 TB/yr total)",
        f"- M1 cores bottoms-up: {CURRENT_CORES} x (M1 mitigated compressed / {STOEHR_MITIG_COMPR_PB}) "
        f"= {refined['cores']['M1']:.0f}",
        "",
        "## Refined table",
        "",
        "| Quantity | 2026 | M1 | M4 | M5 | Units |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for q in PDR:
        vals = [current[q]] + [refined[q][ms] for ms in ("M1", "M4", "M5")]
        md.append(
            f"| {ROW_LABELS[q].replace(chr(92), '')} | "
            + " | ".join(fmt(v) for v in vals)
            + f" | {PDR[q][4]} |"
        )
    (out_dir / "rollout_refined.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print((out_dir / "rollout_refined.md").read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
