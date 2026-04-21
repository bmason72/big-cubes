"""
SDD-style distribution plots for the WSU data-properties pipeline.

Thin wrapper around the SDD plot functions in ``wsu_plots.py``.  The four
figures produced match the SDD data-properties appendix:

    productsize.png  -- MOUS product-size CCDF, Pre-WSU vs Milestone 5
    datavol.png      -- Visibility data-volume CCDF
    datarate.png     -- Data-rate CCDF
    sysperf.png      -- System-performance CCDF

Each plot optionally overlays a shaded band showing the min/max envelope
across band-1/band-2 Monte-Carlo realizations (the same
``hist_cumulative`` structure the original notebooks feed into these
functions).
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from astropy.table import Table

from config import DEFAULT_CONFIG, PipelineConfig
from pipeline import load_database, register_gvis_unit

log = logging.getLogger("wsu_pipeline.plots")


# M5/goal quantity list for calc_wsu_stats_allsamples -- the SDD plots
# reference these column names as band1_band2_estimate keys.
# Sysperf flavor -- only "aprojonly" exists on the M5/goal columns; the
# allgrid variant of wsu_sysperf_goal_stepped2_* is not in the base DB.
_SYSPERF_LABEL = "aprojonly"

_SDD_PLOT_QUANTITIES: List[str] = [
    "wsu_productsize_goal_stepped2",
    "wsu_datavol_goal_stepped2_typical_total",
    "wsu_datarate_goal_stepped2_typical",
    f"wsu_sysperf_goal_stepped2_typical_{_SYSPERF_LABEL}",
]


def _compute_realization_aggregate(
    outdir: str,
    filename: str,
    quantity_list: Optional[List[str]] = None,
) -> Dict:
    """Run wsu_db.calculate_dist and return its hist_cumulative dict.

    ``calculate_dist`` reads every realization in ``outdir``, builds
    cumulative histograms per quantity, and stores the min/max/median
    envelope across realizations.  We pull the hist_cumulative dict
    that the SDD plot functions expect as ``band1_band2_estimate``.
    """
    import matplotlib
    matplotlib.use("Agg")  # calculate_dist calls plt.hist under the hood
    import wsu_db

    qlist = quantity_list or _SDD_PLOT_QUANTITIES
    n_files = len(list(Path(outdir).glob(f"{filename}_*.ecsv")))
    log.info("Aggregating %d realizations in %s", n_files, outdir)
    agg = wsu_db.calculate_dist(
        outDir=outdir,
        filename=filename,
        quantity_list=qlist,
    )
    return agg.get("hist_cumulative", {})


def generate_sdd_plots(
    db: Table,
    output_dir: str,
    cfg: PipelineConfig = DEFAULT_CONFIG,
    realization_dir: Optional[str] = None,
    overlay_realizations: bool = True,
) -> List[str]:
    """Produce the four SDD-style CCDF plots.

    Parameters
    ----------
    db
        Base per-MOUS database (the one loaded by ``pipeline.load_database``).
    output_dir
        Directory to write .png files into; created if absent.
    cfg
        Pipeline config; controls realization dir / filename defaults.
    realization_dir
        Override the realization directory.  Defaults to
        ``cfg.realizations.outdir``.  Only consulted when
        ``overlay_realizations`` is True.
    overlay_realizations
        If True, compute the band-1/2 envelope from realizations and
        overlay it as a shaded band.  Set False for a quick single-pass
        plot.
    """
    register_gvis_unit()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import wsu_plots

    os.makedirs(output_dir, exist_ok=True)

    estimate: Optional[Dict] = None
    if overlay_realizations:
        outdir = realization_dir or cfg.realizations.outdir
        if Path(outdir).exists():
            estimate = _compute_realization_aggregate(
                outdir, cfg.realizations.filename
            ) or None
        else:
            log.warning("Realization dir %s not found; skipping overlay", outdir)

    written: List[str] = []

    specs = [
        ("productsize.png", wsu_plots.plot_productsize_comparison_sdd,
         {"plot_title": "Product Size"}),
        ("datavol.png", wsu_plots.plot_datavol_comparison_sdd,
         {"plot_title": "Visibility Data Volume", "datatype": "total"}),
        ("datarate.png", wsu_plots.plot_datarate_comparison_sdd,
         {"plot_title": "Data Rate"}),
        ("sysperf.png", wsu_plots.plot_soc_result_cumulative_sdd,
         {"plot_title": "System Performance", "label": _SYSPERF_LABEL}),
    ]

    for fname, func, kwargs in specs:
        figpath = str(Path(output_dir) / fname)
        plt.figure(figsize=(7, 5))
        # sysperf uses pltname= instead of figname=.
        save_kw = "pltname" if func is wsu_plots.plot_soc_result_cumulative_sdd else "figname"
        try:
            func(db, band1_band2_estimate=estimate, **{save_kw: figpath}, **kwargs)
            plt.close("all")
            written.append(figpath)
            log.info("Wrote %s", figpath)
        except Exception as exc:
            plt.close("all")
            log.error("Failed to generate %s: %s", fname, exc)

    return written


def generate_plots_cli(
    db_path: str,
    output_dir: str,
    cfg: PipelineConfig = DEFAULT_CONFIG,
    realization_dir: Optional[str] = None,
    overlay_realizations: bool = True,
) -> List[str]:
    """Convenience entry point: load DB from path, then call generate_sdd_plots."""
    db = load_database(db_path)
    return generate_sdd_plots(
        db,
        output_dir=output_dir,
        cfg=cfg,
        realization_dir=realization_dir,
        overlay_realizations=overlay_realizations,
    )
