"""
LaTeX table generators for the WSU data-properties pipeline.

Two families:

    * ``generate_memo_tables(stats, out_dir)``
      Produces the three sideways tables in the memo format:
        - wsu_datarate_summary.tex
        - wsu_datavol_summary.tex
        - wsu_sysperf_summary.tex

    * ``generate_sdd_tables(stats, out_dir)``
      Produces the four SDD-style tables (M5 vs pre-WSU):
        - table_data_prop.tex
        - table_data_rates_prop.tex
        - table_system_performance_wsu.tex
        - data_prop_summ.tex

Both families accept the dict returned by ``pipeline.compute_stats``.

For verification use ``verify_tables`` which parses the numeric cells out of
a generated .tex file and compares them against the corresponding reference
table in the repo, applying the same tolerance policy used by the test suite.
"""

from __future__ import annotations

import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# ---------------------------------------------------------------------------
# Number formatting helpers
# ---------------------------------------------------------------------------

def fmt_rate(x: float) -> str:
    """Data-rate formatter: ``0.024`` (3 decimal places, unpadded)."""
    return f"{x:.3f}"


def fmt_datavol(x: float) -> str:
    """Data-volume formatter: width-7 right-aligned, 3 decimals.

    Yields strings like ``"  0.069"`` / ``" 13.249"`` / ``"374.689"`` /
    ``"999.172"`` to match the memo's visual alignment.
    """
    return f"{x:>7.3f}"


def fmt_visrate(x: float) -> str:
    """Vis-rate formatter: min-width 5, 1 decimal -> ``" 41.7"`` / ``"3002.8"``."""
    return f"{x:>5.1f}"


def fmt_sysperf(x: float) -> str:
    """Sysperf formatter: min-width 7, 3 decimals, trailing space.

    Yields ``"  0.012 "`` / ``"  7.480 "`` to match the memo cells, which
    include an intentional space before the closing ``$``.
    """
    return f"{x:>7.3f} "


def fmt_nchan(x: float) -> str:
    """Integer with thousands separator: ``16,443`` / ``1,185,184``."""
    return f"{x:,.0f}"


# SDD-style: mixed scientific/plain notation depending on magnitude.
def fmt_sdd(x: float, sig_decimals: int = 2) -> str:
    """Format x for SDD-style tables.

    * abs(x) < 0.1 : scientific notation ``"X.X \\\\times 10^{-n}"``
    * otherwise    : plain with 2-3 sig figs

    Returns a LaTeX-ready math string (without surrounding ``$``).
    """
    if x == 0 or not math.isfinite(x):
        return f"{x:g}"
    absx = abs(x)
    if absx < 0.1:
        # scientific
        exp = int(math.floor(math.log10(absx)))
        mantissa = x / (10 ** exp)
        return f"{mantissa:.1f} \\times 10^{{{exp}}}"
    if absx >= 100:
        # plain integer-ish
        if absx >= 1000:
            return f"{x:.0f}"
        return f"{x:.0f}"
    if absx >= 10:
        return f"{x:.1f}"
    if absx >= 1:
        return f"{x:.2f}"
    # 0.1 <= x < 1
    return f"{x:.2f}"


# ---------------------------------------------------------------------------
# Cell assembly
# ---------------------------------------------------------------------------

def _row_cells_dollar(values: Sequence[float], fmt) -> str:
    """Join values as ``"$v$"`` cells separated by ``"& "`` for memo tables."""
    return "& ".join(f"${fmt(v)}$" for v in values)


def _row_cells_plain(values: Sequence[float], fmt) -> str:
    """Join values without dollars (for nchan row)."""
    return "& ".join(fmt(v) for v in values)


def _row_cells_sdd(values: Sequence[float]) -> str:
    return "   & ".join(f"${fmt_sdd(v)}$" for v in values)


# ---------------------------------------------------------------------------
# Memo tables
# ---------------------------------------------------------------------------

_MEMO_HEADER = r"""\definecolor{myband}{RGB}{255,235,205}
\definecolor{myinitial}{RGB}{181,110,110}
\definecolor{myms4}{RGB}{115, 112, 138}
\definecolor{mygoal}{RGB}{152,168,214}
\definecolor{my2x}{RGB}{115, 112, 138}
\definecolor{my4x}{RGB}{251, 206, 177}

"""

_MEMO_COLSPEC = (
    r"\begin{tblr}{colspec={|[2pt]Q[l]|Q[l]|[1pt]c|c|c|c|c|c|c|c|c|[2pt]}," "\n"
    r"width=\textwidth," "\n"
    r"cells = {font=\scriptsize}}" "\n"
)

_MEMO_TOP = (
    r"\hline[2pt]" "\n"
    r"\SetCell[c=2,r=2]{c} & & \SetCell[c=3]{c,bg=myinitial} Milestone 1 & & & "
    r"\SetCell[c=3]{c,bg=myms4} Milestone 4 & & & \SetCell[c=3]{c,bg=mygoal} "
    r"Milestone 5 & & & \\ \hline[1pt]" "\n"
    r"& & 12m & 7m & both & 12m & 7m & both & 12m & 7m & both  \\ \hline[1pt]" "\n"
)

_MS_ORDER = ("M1", "M4", "M5")
_ARR_ORDER = ("12m", "7m", "both")


def _values_flat(stats: Dict, stat: str, quantity: str,
                 fallback: float = math.nan) -> List[float]:
    """Return a length-9 list: [M1-12m, M1-7m, M1-both, M4-..., M5-...]."""
    out: List[float] = []
    for ms in _MS_ORDER:
        for arr in _ARR_ORDER:
            bucket = stats[ms][arr].get(stat, {})
            out.append(bucket.get(quantity, fallback))
    return out


def _total_per_cycle_flat(stats: Dict, quantity: str) -> List[float]:
    out: List[float] = []
    for ms in _MS_ORDER:
        for arr in _ARR_ORDER:
            bucket = stats[ms][arr].get("total_per_cycle", {})
            out.append(bucket.get(quantity, math.nan))
    return out


def write_datarate_summary(stats: Dict, path: str | os.PathLike) -> str:
    """Write the memo-style data-rate + channel-count table."""
    dr_med = _values_flat(stats, "median", "datarate")
    dr_twa = _values_flat(stats, "twa", "datarate")
    dr_max = _values_flat(stats, "max", "datarate")
    ch_med = _values_flat(stats, "median", "nchan_agg")
    ch_twa = _values_flat(stats, "twa", "nchan_agg")
    ch_max = _values_flat(stats, "max", "nchan_agg")

    body = (
        f"Data Rate & {{Median (GB/s)}}& {_row_cells_dollar(dr_med, fmt_rate)}\\\\ \n"
        f" & {{Time Weighted \\\\ Average (GB/s)}}& "
        f"{_row_cells_dollar(dr_twa, fmt_rate)}\\\\ \n"
        f" & Maximum (GB/s)& {_row_cells_dollar(dr_max, fmt_rate)}\\\\ \n"
        "\\hline \n"
        f"{{Number of \\\\ of Channels}} & {{Median }}& "
        f"{_row_cells_plain(ch_med, fmt_nchan)}\\\\ \n"
        f" & {{Time Weighted \\\\ Average }}& "
        f"{_row_cells_plain(ch_twa, fmt_nchan)}\\\\ \n"
        f" & Maximum & {_row_cells_plain(ch_max, fmt_nchan)}\\\\ \n"
    )

    content = (
        _MEMO_HEADER
        + r"\begin{sidewaystable}" + "\n"
        + r"\centering" + "\n"
        + r"\caption{Overview of Data Rate Properties for  WSU "
          r"\label{tab:overview_datarates}}" + "\n"
        + _MEMO_COLSPEC
        + _MEMO_TOP
        + body
        + "\n" + r"\hline[2pt]" + "\n"
        + r"\end{tblr}" + "\n"
        + r"\end{sidewaystable}" + "\n"
    )
    Path(path).write_text(content)
    return str(path)


def write_datavol_summary(stats: Dict, path: str | os.PathLike) -> str:
    """Write memo-style data-volume + product-size table.

    The memo groups data volume into three blocks (total vis, science vis,
    product size); each block has median/TWA/max rows plus a final
    per-cycle total row (``\\cline{2-14}`` under the "max" row).
    """
    rows: List[str] = []
    for label, key in [
        ("Visibility Data \\\\ Volume (Total)",  "datavol_total"),
        ("Visibility Data \\\\ Volume (Science)", "datavol_science"),
        ("Product Size \\\\ (Total)",            "productsize"),
    ]:
        med = _values_flat(stats, "median", key)
        twa = _values_flat(stats, "twa",    key)
        mx  = _values_flat(stats, "max",    key)
        per_cycle = _total_per_cycle_flat(stats, key)

        # Median / TWA / Max / per-cycle total
        block = (
            f"{{{label} }} & {{Median (TB)}}& "
            f"{_row_cells_dollar(med, fmt_datavol)}\\\\ \n"
            f" & {{Time Weighted \\\\ Average (TB)}} & "
            f"{_row_cells_dollar(twa, fmt_datavol)}\\\\ \n"
            f" & Maximum (TB)& {_row_cells_dollar(mx, fmt_datavol)}\\\\ \n"
            "\\cline{2-14} \n"
            f"  & {{{{ {{\\bf Total per cycle (PB)}}}}}}& "
            f"{_row_cells_dollar(per_cycle, fmt_datavol)}\\\\ \n"
            "\\hline \n"
        )
        rows.append(block)

    content = (
        _MEMO_HEADER
        + r"\begin{sidewaystable}" + "\n"
        + r"\centering" + "\n"
        + r"\caption{Overview of Data Volume Properties for WSU "
          r"\label{tab:overview_datavol}}" + "\n"
        + _MEMO_COLSPEC
        + _MEMO_TOP
        + rows[0] + rows[1] + rows[2]
        + "\n" + r"\hline[2pt]" + "\n"
        + r"\end{tblr}" + "\n"
        + r"\end{sidewaystable}" + "\n"
    )
    Path(path).write_text(content)
    return str(path)


def write_sysperf_summary(stats: Dict, path: str | os.PathLike) -> str:
    """Write memo-style data-rate + vis-rate + sysperf table."""
    dr_med = _values_flat(stats, "median", "datarate")
    dr_twa = _values_flat(stats, "twa",    "datarate")
    dr_max = _values_flat(stats, "max",    "datarate")
    vr_med = _values_flat(stats, "median", "visrate")
    vr_twa = _values_flat(stats, "twa",    "visrate")
    vr_max = _values_flat(stats, "max",    "visrate")
    sp_med = _values_flat(stats, "median", "sysperf")
    sp_twa = _values_flat(stats, "twa",    "sysperf")
    sp_max = _values_flat(stats, "max",    "sysperf")

    body = (
        f"Data Rate & {{Median (GB/s)}}& {_row_cells_dollar(dr_med, fmt_rate)}\\\\ \n"
        f" & {{Time Weighted \\\\ Average (GB/s)}}& "
        f"{_row_cells_dollar(dr_twa, fmt_rate)}\\\\ \n"
        f" & Maximum (GB/s)& {_row_cells_dollar(dr_max, fmt_rate)}\\\\ \n"
        "\\hline \n"
        f"Visibility Rate & {{Median (Gvis/hr)}}& "
        f"{_row_cells_dollar(vr_med, fmt_visrate)}\\\\ \n"
        f" & {{Time Weighted \\\\ Average (Gvis/hr)}}& "
        f"{_row_cells_dollar(vr_twa, fmt_visrate)}\\\\ \n"
        f" & Maximum (Gvis/hr)& {_row_cells_dollar(vr_max, fmt_visrate)}\\\\ \n"
        "\\hline \n"
        f"{{System \\\\ Performance}} & {{Median (PFLOP/s)}}& "
        f"{_row_cells_dollar(sp_med, fmt_sysperf)}\\\\ \n"
        f" & {{Time Weighted \\\\ Average (PFLOP/s)}}& "
        f"{_row_cells_dollar(sp_twa, fmt_sysperf)}\\\\ \n"
        f" & Maximum (PFLOP/s)& {_row_cells_dollar(sp_max, fmt_sysperf)}\\\\ \n"
    )

    content = (
        _MEMO_HEADER
        + r"\begin{sidewaystable}" + "\n"
        + r"\centering" + "\n"
        + r"\caption{Overview of System Performance Related Quantities for "
          r"WSU \label{tab:overview_sysperf}}" + "\n"
        + _MEMO_COLSPEC
        + _MEMO_TOP
        + body
        + "\n" + r"\hline[2pt]" + "\n"
        + r"\end{tblr}" + "\n"
        + r"\end{sidewaystable}" + "\n"
    )
    Path(path).write_text(content)
    return str(path)


def generate_memo_tables(stats: Dict, out_dir: str | os.PathLike) -> List[str]:
    """Write all three memo tables into ``out_dir``.  Returns list of paths."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    return [
        write_datarate_summary(stats, out / "wsu_datarate_summary.tex"),
        write_datavol_summary(stats,  out / "wsu_datavol_summary.tex"),
        write_sysperf_summary(stats,  out / "wsu_sysperf_summary.tex"),
    ]


# ---------------------------------------------------------------------------
# SDD-style tables (M5 vs pre-WSU)
# ---------------------------------------------------------------------------

def write_sdd_data_prop(stats: Dict, path: str | os.PathLike) -> str:
    """Write SDD-style visibility-volume + product-size table.

    Columns: WSU M5 12m, WSU M5 7m, pre-WSU (BLC) 12m, pre-WSU (BLC) 7m.
    """
    # Pull the four columns per stat/quantity.
    def row(stat: str, quantity: str) -> List[float]:
        return [
            stats["M5"]["12m"][stat].get(quantity, math.nan),
            stats["M5"]["7m"][stat].get(quantity, math.nan),
            stats["BLC"]["12m"][stat].get(quantity, math.nan),
            stats["BLC"]["7m"][stat].get(quantity, math.nan),
        ]

    def row_per_cycle(quantity: str) -> List[float]:
        return [
            stats["M5"]["12m"]["total_per_cycle"].get(quantity, math.nan),
            stats["M5"]["7m"]["total_per_cycle"].get(quantity, math.nan),
            stats["BLC"]["12m"]["total_per_cycle"].get(quantity, math.nan),
            stats["BLC"]["7m"]["total_per_cycle"].get(quantity, math.nan),
        ]

    def fmt_row(values: Sequence[float]) -> str:
        return "   & ".join(f"${fmt_sdd(v)}$" for v in values)

    lines = [
        r"\begin{table}[hbt!]",
        r"\centering",
        r"\begin{talltblr}[  ",
        r"  caption={Estimated  MOUS-level interferometric data volume properties.},",
        r"  label={tab:dataprop},",
        r"  remark{Note}={\small From \cite{drru_mean}, with WSU  updated to conform to the  configuration at  Milestone 5.}",
        r"  ]{",
        r"  colspec={Q[l,m]Q[l,m]*{4}{X[c,m]}},",
        r"  width=\textwidth,",
        r"  hlines,",
        r"  vlines,",
        r"  vline{1-2} = {1-2}{0pt},",
        r"  hline{1} = {1-2}{0pt},",
        r"  cells = {font=\small},",
        r"  cell{1-2}{3-4} = {my2x, font=\small\bfseries},",
        r"  cell{1-2}{5-6} = {my4x, font=\small\bfseries},",
        r"  row{6,10} = {bg=highlight2, font=\small\bfseries},",
        r"  }",
        r"\SetCell[c=2,r=2]{c} & & \SetCell[c=2]{c}{WSU at Milestone 5} & &  \SetCell[c=2]{c}{pre-WSU} &  \\",
        r"& & 12 m & 7 m & 12 m & 7 m \\",
        r"\SetCell[r=4]{bg=lightgray}{Visibility Data \\",
        r"Volume (Total)} & {Median (TB)}                     & "
        + fmt_row(row("median", "datavol_total")) + r"     \\",
        r" & {Observing Time\\Weighted Average (TB)} & "
        + fmt_row(row("twa", "datavol_total")) + r"   \\",
        r" & Maximum (TB)                           & "
        + fmt_row(row("max", "datavol_total")) + r"      \\",
        r" & Total per cycle (PB)                   & "
        + fmt_row(row_per_cycle("datavol_total")) + r"      \\",
        r"\SetCell[r=4]{bg=lightgray}{Product Size\\(Total)}        & {Median (TB)}                     & "
        + fmt_row(row("median", "productsize")) + r"       \\",
        r" & {Observing Time\\Weighted Average (TB)} & "
        + fmt_row(row("twa", "productsize")) + r"      \\",
        r" & Maximum (TB)                           & "
        + fmt_row(row("max", "productsize")) + r"      \\",
        r" & Total per cycle (PB)                   & "
        + fmt_row(row_per_cycle("productsize")) + r"      \\",
        r"\end{talltblr}",
        r"\end{table}",
        "",
    ]
    content = "\n".join(lines)
    Path(path).write_text(content)
    return str(path)


def write_sdd_data_rates_prop(stats: Dict, path: str | os.PathLike) -> str:
    """Write SDD-style M5-only data-rate table (12m / 7m / both)."""
    def row(stat: str) -> List[float]:
        return [
            stats["M5"]["12m"][stat]["datarate"],
            stats["M5"]["7m"][stat]["datarate"],
            stats["M5"]["both"][stat]["datarate"],
        ]

    def fmt_row(values):
        return " & ".join(f"${fmt_sdd(v)}$" for v in values)

    lines = [
        r"\begin{center}",
        r"\begin{talltblr}",
        r"  [caption={Estimated data rate properties for the WSU at Milestone 5\TblrNote{a}. Note that the total power antennas will increase the averages rates by $\sim 3\%$ \cite{drru_mean}.},",
        r"   note{a}={As per \cite{drru_mean}, these quantities assume a realistic number of antennas: 47$\times$12-m and 10$\times$7-m.},",
        r"   label={tab:overview_datarates}]",
        r"   {",
        r"    colspec={Q[l,m,3cm] *{3}{X[c,m]}},",
        r"    width=\textwidth,",
        r"    cells = {font=\small},",
        r"    hlines,",
        r"    vlines,",
        r"    cell{1-2}{2-4} = {my4x, font=\small\bfseries},",
        r"    cell{1}{1} = {lightgray, font=\small\bfseries}",
        r"  }",
        r"\SetCell[r=2]{c}{Data Rate\\ (\GBs)} & \SetCell[c=3]{c}{WSU at Milestone 5} \\",
        r"& 12 m & 7 m & both \\",
        r"{Median} & " + fmt_row(row("median")) + r" \\ ",
        r"{Observing Time \\ Weighted Average} & " + fmt_row(row("twa")) + r"\\",
        r"Maximum & " + fmt_row(row("max")) + r" \\ ",
        r"\end{talltblr}",
        r"\end{center}",
        "",
    ]
    content = "\n".join(lines)
    Path(path).write_text(content)
    return str(path)


def write_sdd_system_performance(stats: Dict, path: str | os.PathLike) -> str:
    """Write SDD-style rate + vis-rate + sysperf (M5 vs pre-WSU)."""
    def row(stat: str, quantity: str) -> List[float]:
        return [
            stats["M5"]["12m"][stat].get(quantity, math.nan),
            stats["M5"]["7m"][stat].get(quantity, math.nan),
            stats["BLC"]["12m"][stat].get(quantity, math.nan),
            stats["BLC"]["7m"][stat].get(quantity, math.nan),
        ]

    def fmt_row(values):
        return "   & ".join(f"${fmt_sdd(v)}$" for v in values)

    lines = [
        r"\begin{table}[hbt!]",
        r"\centering",
        r"\begin{talltblr}[%",
        r"  caption={Estimated MOUS-level interferometric data rate and required compute system performance.},",
        r"  label={tab:computeSizeTable1},",
        r"  remark{Note}={\small From \cite{dpSizeOfCompute} and \cite{drru_mean}, with WSU updated to conform with Milestone 5. Note that there is a typo in the units given in the middle 3 rows of Table 2 of \cite{dpSizeOfCompute}, which should be MVis/hr. FLOP/s estimates assume 4000 hr/year of science observing and include the $3\times$ overhead/contingency factor,  as discussed in \cite{computeCostTsi}.}",
        r"]{",
        r"  colspec={Q[l,m]Q[l,m]*{4}{X[c,m]}},",
        r"  width=\textwidth,",
        r"  hlines,",
        r"  vlines,",
        r"  vline{1-2} = {1-2}{0pt},",
        r"  hline{1}            = {1-2}{0pt},",
        r"  hline{4,5,7,8,10,11} = {2-6}{0pt},",
        r"  cells = {font=\small},",
        r"  cell{1-2}{3-4} = {my2x, font=\small\bfseries},",
        r"  cell{1-2}{5-6} = {my4x, font=\small\bfseries},",
        r"}",
        r"\SetCell[c=2,r=2]{c}      & & \SetCell[c=2]{c}{WSU at Milestone 5} & & \SetCell[c=2]{c}{pre\u2011WSU} & \\",
        r"                         & & 12 m      & 7 m      & 12 m      & 7 m      \\",
        r"\SetCell[r=3]{bg=lightgray}{Data Rate}",
        r"  & {Median (GB/s)}                   & " + fmt_row(row("median", "datarate")) + r"   \\",
        r"  & {Time Weighted\\Average (GB/s)}   & " + fmt_row(row("twa",    "datarate")) + r"   \\",
        r"  & Maximum (GB/s)                    & " + fmt_row(row("max",    "datarate")) + r"   \\",
        r"\SetCell[r=3]{bg=lightgray}{Visibility Rate}",
        r"  & {Median (Gvis/hr)}                & " + fmt_row(row("median", "visrate"))  + r"     \\",
        r"  & {Time Weighted\\Average (Gvis/hr)}& " + fmt_row(row("twa",    "visrate"))  + r"     \\",
        r"  & Maximum (Gvis/hr)                 & " + fmt_row(row("max",    "visrate"))  + r"    \\",
        r"\SetCell[r=3]{bg=lightgray}{System\\Performance}",
        r"  & {Median (PFLOP/s)}                & " + fmt_row(row("median", "sysperf"))  + r"   \\",
        r"  & {Time Weighted\\Average (PFLOP/s)}& " + fmt_row(row("twa",    "sysperf"))  + r"   \\",
        r"  & Maximum (PFLOP/s)                 & " + fmt_row(row("max",    "sysperf"))  + r"   \\",
        r"\end{talltblr}",
        r"\end{table}",
        "",
    ]
    content = "\n".join(lines)
    Path(path).write_text(content)
    return str(path)


# Additional factors needed for the summary table; user-configurable.
@dataclass
class SummaryFactors:
    """Ad-hoc multipliers applied to raw pipeline outputs to reproduce the
    SDD summary table ``data_prop_summ.tex``.

    These capture quantities that the base pipeline does not compute:
      - ``tp_uplift`` : fraction added to rates / vols for total-power
        antennas (the SDD text quotes ~3%).
      - ``compression_factor`` : product-volume-compressed / product-volume.
      - ``gous_factor`` : prod-volume-with-gous / prod-volume.
      - ``contingency_factor`` : PFLOP/s contingency (3x per SDD note).
    """
    tp_uplift: float = 0.03
    compression_factor: float = 0.66    # 7.9/12.1 from SDD quoted numbers
    gous_factor: float = 0.884           # 10.7/12.1 ~ compressed+GOUS routing
    contingency_factor: float = 3.0


def _sum_over_arrays(stats: Dict, ms: str, stat: str, quantity: str) -> float:
    """Sum of 12m + 7m (approx TP contribution applied separately)."""
    a = stats[ms]["12m"][stat].get(quantity, math.nan)
    b = stats[ms]["7m"][stat].get(quantity, math.nan)
    return a + b


def write_sdd_summary(stats: Dict, path: str | os.PathLike,
                       factors: Optional[SummaryFactors] = None) -> str:
    """Write SDD-style top-level summary table (``data_prop_summ.tex``).

    The summary table mixes computed pipeline values with ad-hoc factors
    (TP inclusion, compression, GOUS routing, FLOP contingency).  The
    factors are parameterised via :class:`SummaryFactors` so downstream
    consumers can re-run with alternative assumptions.
    """
    f = factors or SummaryFactors()

    # Raw visibility volume (total_per_cycle).
    m5_raw = _sum_over_arrays(stats, "M5", "total_per_cycle", "datavol_total") * (1 + f.tp_uplift)
    blc_raw = _sum_over_arrays(stats, "BLC", "total_per_cycle", "datavol_total") * (1 + f.tp_uplift)

    # Time-averaged data rate (TWA from 12m, since "both" already pools).
    m5_dr = stats["M5"]["both"]["twa"]["datarate"] * (1 + f.tp_uplift)
    blc_dr = stats["BLC"]["both"]["twa"]["datarate"] * (1 + f.tp_uplift)

    # Product volume (total_per_cycle).
    m5_prod = _sum_over_arrays(stats, "M5", "total_per_cycle", "productsize") * (1 + f.tp_uplift)
    blc_prod = _sum_over_arrays(stats, "BLC", "total_per_cycle", "productsize") * (1 + f.tp_uplift)

    # Derived products.
    m5_prod_compressed = m5_prod * f.compression_factor
    blc_prod_compressed = blc_prod * f.compression_factor
    m5_prod_gous = m5_prod * f.gous_factor
    m5_total = m5_raw + m5_prod
    blc_total = blc_raw + blc_prod

    # Sysperf with contingency.
    m5_sp = stats["M5"]["both"]["twa"]["sysperf"] * f.contingency_factor
    blc_sp = stats["BLC"]["both"]["twa"]["sysperf"] * f.contingency_factor

    def c(x: float) -> str:
        return f"${fmt_sdd(x)}$"

    lines = [
        r"\begin{table}[hbt!]",
        r"\centering",
        r"\begin{talltblr}",
        r"  [caption={Summary of key WSU data properties compared to current values.},",
        r"   label={tab:data_summary},",
        r"   remark{Notes}={\small Except as indicated all quantities are from \cite{drru_mean} with ``Later WSU'' adjusted  for conformance to Milestone 5, and with other corrections as described in the text above. {\it Exceptions}: \textsuperscript{a}  Pre-WSU raw data volume and aggregate mitigated product volumes were measured in the archive \cite{computeCostTsi}; \textsuperscript{b} From \cite{dpSizeOfCompute} with WSU updated to Milestone 5, assuming 4000 hrs/year of science observing, with $3\times$ contingency factor \cite{computeCostTsi}.}]",
        r"   {",
        r"    colspec={Q[l,m,4.5cm] X[c,m] X[c,m] X[c,m]},",
        r"    width=\textwidth,",
        r"    cells = {font=\small},",
        r"    hlines,",
        r"    vlines,",
        r"    cell{1}{2} = {my2x, font=\small\bfseries},",
        r"    cell{1}{3} = {my4x, font=\small\bfseries},",
        r"    cell{1}{1} = {lightgray, font=\small\bfseries},",
        r"    cell{1}{4} = {lightgray, font=\small\bfseries}",
        r"  }",
        r"Quantity & WSU at Milestone 5 & Pre-WSU Value & Units \\ ",
        f"{{Raw data volume (7m+12m+TP)\\textsuperscript{{a}}}} & {c(m5_raw)} &  {c(blc_raw)} & PB/cycle  \\\\",
        f"{{Time-averaged data rate (7m+12m+TP)\\textsuperscript{{b}}}} & {c(m5_dr)} & {c(blc_dr)}  & GB/s  \\\\",
        f"{{Product volume (7m+12m+TP)}} & {c(m5_prod)} &  {c(blc_prod)} & PB/cycle  \\\\",
        f"{{Product volume (compressed)}} & {c(m5_prod_compressed)} & {c(blc_prod_compressed)} & PB/cycle   \\\\",
        r"{Mitigated product volume\textsuperscript{a}} & --- & $7.0 \times 10^{-2}$  & PB/cycle  \\",
        f"{{Prod. volume with GOUS}} & {c(m5_prod_gous)} & --- & PB/cycle \\\\",
        f"{{Total data volume}} & {c(m5_total)} & {c(blc_total)} & PB/cycle \\\\",
        f"{{Estimated Processing System Performance\\textsuperscript{{b}}}} & {c(m5_sp)} & {c(blc_sp)} & PFLOP/s \\\\",
        r"\end{talltblr}",
        r"\end{table}",
        "",
    ]
    content = "\n".join(lines)
    Path(path).write_text(content)
    return str(path)


def generate_sdd_tables(stats: Dict, out_dir: str | os.PathLike,
                         factors: Optional[SummaryFactors] = None) -> List[str]:
    """Write all four SDD-style tables into ``out_dir``."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    return [
        write_sdd_data_prop(stats,          out / "table_data_prop.tex"),
        write_sdd_data_rates_prop(stats,    out / "table_data_rates_prop.tex"),
        write_sdd_system_performance(stats, out / "table_system_performance_wsu.tex"),
        write_sdd_summary(stats,            out / "data_prop_summ.tex", factors),
    ]


# ---------------------------------------------------------------------------
# Verification: parse numeric cells from a LaTeX table
# ---------------------------------------------------------------------------

# Recognise either a plain number ``123.45`` or scientific ``1.6 \times
# 10^{-2}``.  Whitespace inside the ``\\times`` form is flexible.
_SCIENTIFIC_RE = re.compile(
    r"(?P<mantissa>-?\d+(?:\.\d+)?)\s*\\times\s*10\^\{(?P<exp>-?\d+)\}"
)
# Two alternatives: thousands-separated (must have at least one ``,``) OR
# a plain number.  Requiring ``+`` on the comma group stops a value like
# ``1168.000`` from matching as ``116`` via the thousands-separator path.
_PLAIN_RE = re.compile(r"-?\d{1,3}(?:,\d{3})+(?:\.\d+)?|-?\d+(?:\.\d+)?")


def _parse_cell(cell: str) -> Optional[float]:
    """Extract a float from one LaTeX cell, skipping non-numeric cells."""
    s = cell.strip()
    # Strip leading/trailing ``$``.
    s = s.strip("$ \t\n\\")
    # Short-circuit on obviously-non-numeric cells.
    if not s or s in {"---", "N/A", ""}:
        return None
    # Try scientific first.
    m = _SCIENTIFIC_RE.search(s)
    if m:
        return float(m.group("mantissa")) * (10 ** int(m.group("exp")))
    # Plain number.
    m = _PLAIN_RE.search(s)
    if m:
        try:
            return float(m.group(0).replace(",", ""))
        except ValueError:
            return None
    return None


# Tokens that, when present anywhere in a row, indicate the row is
# structural (header / rule / preamble) and its numeric content is not a
# data cell -- e.g. ``\\SetCell[c=3,r=2]``, ``\\cline{2-14}``, the
# color/alignment bits.  Data rows may still reference the colored cells
# by name, but never contain these macros.
_STRUCTURAL_TOKENS = (
    "\\SetCell", "\\cline", "\\hline",
    "\\caption", "\\begin", "\\end", "\\definecolor",
    "\\TblrNote", "\\SetCell[", "\\multicolumn",
)


def parse_numeric_cells(tex_path: str | os.PathLike) -> List[List[Optional[float]]]:
    """Parse a LaTeX table file and return ``[row [cell, cell, ...]]``.

    Strips commented lines and preamble/header rows, finds rows terminated
    by ``\\\\``, and splits each row on unescaped ``&``.  Rows containing
    any of the tokens in ``_STRUCTURAL_TOKENS`` are skipped so that digits
    hiding inside ``\\SetCell[c=3,r=2]`` or ``{RGB}{255,235,205}`` are not
    mistaken for data cells.
    """
    content = Path(tex_path).read_text()
    lines = [l for l in content.splitlines()
             if not l.lstrip().startswith("%")]
    stripped = "\n".join(lines)

    # Split rows on ``\\`` only when it sits at end-of-line (optionally
    # followed by whitespace and a newline).  An inline ``\\`` inside a
    # cell -- e.g. ``{Time Weighted \\ Average (GB/s)}`` -- must not
    # cleave the row in half.
    rows: List[List[Optional[float]]] = []
    for raw in re.split(r"\\\\\s*(?:\n|$)", stripped, flags=re.MULTILINE):
        # Each chunk may contain preamble (``\hline[1pt]``, ``\cline{}``,
        # ``\SetCell ...``) on preceding lines plus the data row on the
        # last line.  Strip to the last newline-separated piece so the
        # structural-token filter doesn't drop a valid data row that
        # happens to share a chunk with a header line.
        chunk_lines = [ln for ln in raw.splitlines() if ln.strip()]
        if not chunk_lines:
            continue
        raw = chunk_lines[-1].strip()
        if "&" not in raw:
            continue
        if any(tok in raw for tok in _STRUCTURAL_TOKENS):
            continue
        cells = raw.split("&")
        numeric = [_parse_cell(c) for c in cells]
        if any(v is not None for v in numeric):
            rows.append(numeric)
    return rows


@dataclass
class CellMismatch:
    table: str
    row_idx: int
    col_idx: int
    ref_value: float
    gen_value: Optional[float]
    rel_diff: Optional[float]
    tolerance: float


def verify_tables(generated_dir: str | os.PathLike,
                   reference_dir: str | os.PathLike,
                   pairs: Iterable[Tuple[str, float]],
                   abs_floor: float = 0.001,
                   ) -> List[CellMismatch]:
    """Compare generated vs reference .tex files cell by cell.

    ``pairs`` is an iterable of ``(filename, tolerance)`` -- tolerance is
    the relative tolerance applied to every numeric cell.  ``abs_floor``
    is an absolute-difference floor below which cells are considered
    equal regardless of the relative difference (protects against the
    25% relative error that appears between ``0.004`` and ``0.005``
    simply because the memo rounds tiny values to 3 decimals).
    """
    gen = Path(generated_dir)
    ref = Path(reference_dir)
    mismatches: List[CellMismatch] = []

    for filename, tol in pairs:
        gen_cells = parse_numeric_cells(gen / filename)
        ref_cells = parse_numeric_cells(ref / filename)

        for r_idx in range(max(len(gen_cells), len(ref_cells))):
            g_row = gen_cells[r_idx] if r_idx < len(gen_cells) else []
            r_row = ref_cells[r_idx] if r_idx < len(ref_cells) else []
            for c_idx in range(max(len(g_row), len(r_row))):
                g_val = g_row[c_idx] if c_idx < len(g_row) else None
                r_val = r_row[c_idx] if c_idx < len(r_row) else None
                if r_val is None:
                    continue  # non-numeric ref cell
                if g_val is None:
                    mismatches.append(CellMismatch(
                        filename, r_idx, c_idx, r_val, None, None, tol))
                    continue
                if abs(g_val - r_val) <= abs_floor:
                    continue
                if r_val == 0:
                    rel = 0.0 if g_val == 0 else float("inf")
                else:
                    rel = abs(g_val - r_val) / abs(r_val)
                if rel > tol:
                    mismatches.append(CellMismatch(
                        filename, r_idx, c_idx, r_val, g_val, rel, tol))
    return mismatches
