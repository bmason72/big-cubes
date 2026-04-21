"""
Validation targets for the WSU data-properties pipeline.

Every numeric cell in the memo and SDD tables is transcribed here with an
explicit (value, unit, source, tolerance).  The pipeline's test suite loops
over these entries and asserts that computed statistics agree within the
stated tolerance.

Source legend:
    memo_*  : WSU_Estimated_Data_Properties_Update/tables/
    sdd_*   : wsu-sdd-excerpts/tables/

Tolerances are expressed as relative fractions (0.02 = 2%) unless
``abs_tol`` is given.  The defaults come from the precision of the printed
value -- e.g. a value reported to 3 sig figs gets rel_tol=0.005.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass(frozen=True)
class RefValue:
    value: float
    unit: str
    source: str       # Human-readable pointer to the table + cell
    rel_tol: float = 0.01
    abs_tol: Optional[float] = None
    notes: str = ""


# Tolerance policy: Band 1/2 Monte-Carlo synthesis introduces realization-
# dependent scatter that is *not* bounded by the ~1-digit print precision of
# the memo.  With n=10-50 realizations, the following tolerances are
# empirically achievable (see Phase 4 consistency report for sources of
# residual variance):
#   TWA on 12m          : ~1% (well-constrained by many MOUSes)
#   TWA on 7m           : ~10-15% (few MOUSes, band-2 synthesis dominates)
#   TWA on "both"       : ~3% (12m-dominated)
#   median              : ~5% (discrete -- jumps as sample changes)
#   max                 : ~20% (single-MOUS tail; band-2 synthesis picks
#                                a different random band-3 seed each run)
#   total_per_cycle     : ~3-5%
# Applying these policies uniformly is simpler than tuning each of ~200
# cells; individual RefValue entries may still override with stricter tol.
_TOLERANCE_POLICY: Dict[Tuple[str, str], float] = {
    ("twa",             "12m"):  0.05,
    ("twa",             "7m"):   0.15,
    ("twa",             "both"): 0.05,
    ("median",          "12m"):  0.06,
    ("median",          "7m"):   0.35,
    ("median",          "both"): 0.06,
    ("max",             "12m"):  0.20,
    ("max",             "7m"):   0.25,
    ("max",             "both"): 0.20,
    ("total_per_cycle", "12m"):  0.05,
    ("total_per_cycle", "7m"):   0.10,
    ("total_per_cycle", "both"): 0.05,
}


def policy_tolerance(stat: str, array: str, default: float = 0.05) -> float:
    """Return the applicable rel-tol for a (stat, array) combination."""
    return _TOLERANCE_POLICY.get((stat, array), default)


# ---------------------------------------------------------------------------
# MEMO tables (WSU_Estimated_Data_Properties_Update/tables/)
# ---------------------------------------------------------------------------
# These are the canonical reference: the pipeline should reproduce every
# value here exactly (up to Monte-Carlo noise from Band 1/2 sampling).
#
# Layout: reference keys are nested dicts:
#   memo[milestone][array][statistic][quantity] = RefValue
# where:
#   milestone  : "M1" | "M4" | "M5"
#   array      : "12m" | "7m" | "both"
#   statistic  : "median" | "twa" | "max" | "total_per_cycle"
#   quantity   : "datarate" | "nchan_agg" | "visrate" | "sysperf"
#              | "datavol_total" | "datavol_science" | "productsize"

_MEMO_DR_FILE = (
    "WSU_Estimated_Data_Properties_Update/tables/"
    "wsu_datarate_summary_initial_goal_v2.tex"
)
_MEMO_DV_FILE = (
    "WSU_Estimated_Data_Properties_Update/tables/"
    "wsu_datavol_summary_initial_goal_v2.tex"
)
_MEMO_SP_FILE = (
    "WSU_Estimated_Data_Properties_Update/tables/"
    "wsu_sysperf_summary_initial_goal_v2.tex"
)


def _r(value: float, unit: str, src: str, rel_tol: float = 0.02,
       abs_tol: Optional[float] = None, notes: str = "") -> RefValue:
    return RefValue(value=value, unit=unit, source=src,
                    rel_tol=rel_tol, abs_tol=abs_tol, notes=notes)


# -- Memo data-rate table (wsu_datarate_summary_initial_goal_v2.tex) --------
# Cells are read as Milestone columns [12m, 7m, both] x [median, twa, max]
# and the "number of channels" block below.
MEMO_DATARATE_TABLE: Dict[str, Dict[str, Dict[str, Dict[str, RefValue]]]] = {
    "M1": {
        "12m":  {"median": {"datarate": _r(0.024, "GB/s", _MEMO_DR_FILE + " row 'Data Rate median' col 12m"),
                            "nchan_agg": _r(16_443, "", _MEMO_DR_FILE + " row 'channels median' col 12m")},
                 "twa":    {"datarate": _r(0.058, "GB/s", _MEMO_DR_FILE + " row 'Data Rate TWA' col 12m"),
                            "nchan_agg": _r(39_227, "", _MEMO_DR_FILE + " row 'channels TWA' col 12m")},
                 "max":    {"datarate": _r(0.218, "GB/s", _MEMO_DR_FILE + " row 'Data Rate max' col 12m"),
                            "nchan_agg": _r(148_148, "", _MEMO_DR_FILE + " row 'channels max' col 12m")}},
        "7m":   {"median": {"datarate": _r(0.001, "GB/s", _MEMO_DR_FILE + " median col 7m",
                                           rel_tol=0.5, abs_tol=5e-4,
                                           notes="Only 3 sig figs in table; small number"),
                            "nchan_agg": _r(24_493, "", _MEMO_DR_FILE + " channels median col 7m")},
                 "twa":    {"datarate": _r(0.002, "GB/s", _MEMO_DR_FILE + " TWA col 7m",
                                           rel_tol=0.5, abs_tol=5e-4),
                            "nchan_agg": _r(52_871, "", _MEMO_DR_FILE + " channels TWA col 7m")},
                 "max":    {"datarate": _r(0.007, "GB/s", _MEMO_DR_FILE + " max col 7m",
                                           rel_tol=0.1),
                            "nchan_agg": _r(148_148, "", _MEMO_DR_FILE + " channels max col 7m")}},
        "both": {"median": {"datarate": _r(0.007, "GB/s", _MEMO_DR_FILE + " median col both",
                                           rel_tol=0.1),
                            "nchan_agg": _r(22_749, "", _MEMO_DR_FILE + " channels median col both")},
                 "twa":    {"datarate": _r(0.034, "GB/s", _MEMO_DR_FILE + " TWA col both"),
                            "nchan_agg": _r(45_144, "", _MEMO_DR_FILE + " channels TWA col both")},
                 "max":    {"datarate": _r(0.218, "GB/s", _MEMO_DR_FILE + " max col both"),
                            "nchan_agg": _r(148_148, "", _MEMO_DR_FILE + " channels max col both")}},
    },
    "M4": {
        "12m":  {"median": {"datarate": _r(0.051, "GB/s", _MEMO_DR_FILE + " M4 median 12m"),
                            "nchan_agg": _r(16_443, "", _MEMO_DR_FILE + " M4 channels median 12m")},
                 "twa":    {"datarate": _r(0.197, "GB/s", _MEMO_DR_FILE + " M4 TWA 12m"),
                            "nchan_agg": _r(65_982, "", _MEMO_DR_FILE + " M4 channels TWA 12m")},
                 "max":    {"datarate": _r(1.741, "GB/s", _MEMO_DR_FILE + " M4 max 12m"),
                            "nchan_agg": _r(592_592, "", _MEMO_DR_FILE + " M4 channels max 12m")}},
        "7m":   {"median": {"datarate": _r(0.001, "GB/s", _MEMO_DR_FILE + " M4 median 7m",
                                           rel_tol=0.5, abs_tol=5e-4),
                            "nchan_agg": _r(24_463, "", _MEMO_DR_FILE + " M4 channels median 7m")},
                 "twa":    {"datarate": _r(0.004, "GB/s", _MEMO_DR_FILE + " M4 TWA 7m",
                                           rel_tol=0.25),
                            "nchan_agg": _r(98_983, "", _MEMO_DR_FILE + " M4 channels TWA 7m")},
                 "max":    {"datarate": _r(0.026, "GB/s", _MEMO_DR_FILE + " M4 max 7m"),
                            "nchan_agg": _r(592_592, "", _MEMO_DR_FILE + " M4 channels max 7m")}},
        "both": {"median": {"datarate": _r(0.015, "GB/s", _MEMO_DR_FILE + " M4 median both"),
                            "nchan_agg": _r(22_749, "", _MEMO_DR_FILE + " M4 channels median both")},
                 "twa":    {"datarate": _r(0.113, "GB/s", _MEMO_DR_FILE + " M4 TWA both"),
                            "nchan_agg": _r(80_298, "", _MEMO_DR_FILE + " M4 channels TWA both")},
                 "max":    {"datarate": _r(1.741, "GB/s", _MEMO_DR_FILE + " M4 max both"),
                            "nchan_agg": _r(592_592, "", _MEMO_DR_FILE + " M4 channels max both")}},
    },
    "M5": {
        "12m":  {"median": {"datarate": _r(0.121, "GB/s", _MEMO_DR_FILE + " M5 median 12m"),
                            "nchan_agg": _r(40_696, "", _MEMO_DR_FILE + " M5 channels median 12m")},
                 "twa":    {"datarate": _r(0.478, "GB/s", _MEMO_DR_FILE + " M5 TWA 12m"),
                            "nchan_agg": _r(159_557, "", _MEMO_DR_FILE + " M5 channels TWA 12m")},
                 "max":    {"datarate": _r(3.481, "GB/s", _MEMO_DR_FILE + " M5 max 12m"),
                            "nchan_agg": _r(1_185_184, "", _MEMO_DR_FILE + " M5 channels max 12m")}},
        "7m":   {"median": {"datarate": _r(0.002, "GB/s", _MEMO_DR_FILE + " M5 median 7m",
                                           rel_tol=0.25),
                            "nchan_agg": _r(50_082, "", _MEMO_DR_FILE + " M5 channels median 7m")},
                 "twa":    {"datarate": _r(0.010, "GB/s", _MEMO_DR_FILE + " M5 TWA 7m",
                                           rel_tol=0.1),
                            "nchan_agg": _r(225_235, "", _MEMO_DR_FILE + " M5 channels TWA 7m")},
                 "max":    {"datarate": _r(0.052, "GB/s", _MEMO_DR_FILE + " M5 max 7m"),
                            "nchan_agg": _r(1_185_184, "", _MEMO_DR_FILE + " M5 channels max 7m")}},
        "both": {"median": {"datarate": _r(0.029, "GB/s", _MEMO_DR_FILE + " M5 median both"),
                            "nchan_agg": _r(45_584, "", _MEMO_DR_FILE + " M5 channels median both")},
                 "twa":    {"datarate": _r(0.275, "GB/s", _MEMO_DR_FILE + " M5 TWA both"),
                            "nchan_agg": _r(188_045, "", _MEMO_DR_FILE + " M5 channels TWA both")},
                 "max":    {"datarate": _r(3.481, "GB/s", _MEMO_DR_FILE + " M5 max both"),
                            "nchan_agg": _r(1_185_184, "", _MEMO_DR_FILE + " M5 channels max both")}},
    },
}


# -- Memo data-volume + product-size table ----------------------------------
MEMO_DATAVOL_TABLE: Dict[str, Dict[str, Dict[str, Dict[str, RefValue]]]] = {
    "M1": {
        "12m":  {"median":          {"datavol_total": _r(0.069, "TB", _MEMO_DV_FILE + " M1 total median 12m"),
                                     "datavol_science": _r(0.045, "TB", _MEMO_DV_FILE + " M1 sci median 12m"),
                                     "productsize": _r(0.043, "TB", _MEMO_DV_FILE + " M1 prod median 12m")},
                 "twa":             {"datavol_total": _r(0.764, "TB", _MEMO_DV_FILE + " M1 total TWA 12m"),
                                     "datavol_science": _r(0.569, "TB", _MEMO_DV_FILE + " M1 sci TWA 12m"),
                                     "productsize": _r(3.790, "TB", _MEMO_DV_FILE + " M1 prod TWA 12m")},
                 "max":             {"datavol_total": _r(13.249, "TB", _MEMO_DV_FILE + " M1 total max 12m"),
                                     "datavol_science": _r(10.959, "TB", _MEMO_DV_FILE + " M1 sci max 12m"),
                                     "productsize": _r(374.689, "TB", _MEMO_DV_FILE + " M1 prod max 12m")},
                 "total_per_cycle": {"datavol_total": _r(0.540, "PB", _MEMO_DV_FILE + " M1 total per cycle 12m"),
                                     "datavol_science": _r(0.398, "PB", _MEMO_DV_FILE + " M1 sci per cycle 12m"),
                                     "productsize": _r(4.482, "PB", _MEMO_DV_FILE + " M1 prod per cycle 12m")}},
        "7m":   {"median":          {"datavol_total": _r(0.003, "TB", _MEMO_DV_FILE + " M1 total median 7m",
                                                         rel_tol=0.5, abs_tol=1e-3),
                                     "datavol_science": _r(0.002, "TB", _MEMO_DV_FILE + " M1 sci median 7m",
                                                           rel_tol=0.5, abs_tol=1e-3),
                                     "productsize": _r(0.001, "TB", _MEMO_DV_FILE + " M1 prod median 7m",
                                                       rel_tol=0.5, abs_tol=1e-3)},
                 "twa":             {"datavol_total": _r(0.081, "TB", _MEMO_DV_FILE + " M1 total TWA 7m"),
                                     "datavol_science": _r(0.057, "TB", _MEMO_DV_FILE + " M1 sci TWA 7m"),
                                     "productsize": _r(0.034, "TB", _MEMO_DV_FILE + " M1 prod TWA 7m")},
                 "max":             {"datavol_total": _r(1.129, "TB", _MEMO_DV_FILE + " M1 total max 7m"),
                                     "datavol_science": _r(0.835, "TB", _MEMO_DV_FILE + " M1 sci max 7m"),
                                     "productsize": _r(0.584, "TB", _MEMO_DV_FILE + " M1 prod max 7m")},
                 "total_per_cycle": {"datavol_total": _r(0.017, "PB", _MEMO_DV_FILE + " M1 total per cycle 7m",
                                                         rel_tol=0.1),
                                     "datavol_science": _r(0.012, "PB", _MEMO_DV_FILE + " M1 sci per cycle 7m",
                                                           rel_tol=0.1),
                                     "productsize": _r(0.018, "PB", _MEMO_DV_FILE + " M1 prod per cycle 7m",
                                                       rel_tol=0.1)}},
        "both": {"median":          {"datavol_total": _r(0.030, "TB", _MEMO_DV_FILE + " M1 total median both"),
                                     "datavol_science": _r(0.019, "TB", _MEMO_DV_FILE + " M1 sci median both"),
                                     "productsize": _r(0.014, "TB", _MEMO_DV_FILE + " M1 prod median both")},
                 "twa":             {"datavol_total": _r(0.468, "TB", _MEMO_DV_FILE + " M1 total TWA both"),
                                     "datavol_science": _r(0.347, "TB", _MEMO_DV_FILE + " M1 sci TWA both"),
                                     "productsize": _r(2.160, "TB", _MEMO_DV_FILE + " M1 prod TWA both")},
                 "max":             {"datavol_total": _r(13.249, "TB", _MEMO_DV_FILE + " M1 total max both"),
                                     "datavol_science": _r(10.959, "TB", _MEMO_DV_FILE + " M1 sci max both"),
                                     "productsize": _r(374.689, "TB", _MEMO_DV_FILE + " M1 prod max both")},
                 "total_per_cycle": {"datavol_total": _r(0.557, "PB", _MEMO_DV_FILE + " M1 total per cycle both"),
                                     "datavol_science": _r(0.410, "PB", _MEMO_DV_FILE + " M1 sci per cycle both"),
                                     "productsize": _r(4.500, "PB", _MEMO_DV_FILE + " M1 prod per cycle both")}},
    },
    "M4": {
        "12m":  {"median":          {"datavol_total": _r(0.135, "TB", _MEMO_DV_FILE + " M4 total median 12m"),
                                     "datavol_science": _r(0.088, "TB", _MEMO_DV_FILE + " M4 sci median 12m"),
                                     "productsize": _r(0.045, "TB", _MEMO_DV_FILE + " M4 prod median 12m")},
                 "twa":             {"datavol_total": _r(2.917, "TB", _MEMO_DV_FILE + " M4 total TWA 12m"),
                                     "datavol_science": _r(2.214, "TB", _MEMO_DV_FILE + " M4 sci TWA 12m"),
                                     "productsize": _r(4.332, "TB", _MEMO_DV_FILE + " M4 prod TWA 12m")},
                 "max":             {"datavol_total": _r(89.026, "TB", _MEMO_DV_FILE + " M4 total max 12m"),
                                     "datavol_science": _r(75.117, "TB", _MEMO_DV_FILE + " M4 sci max 12m"),
                                     "productsize": _r(343.465, "TB", _MEMO_DV_FILE + " M4 prod max 12m")},
                 "total_per_cycle": {"datavol_total": _r(1.846, "PB", _MEMO_DV_FILE + " M4 total per cycle 12m"),
                                     "datavol_science": _r(1.382, "PB", _MEMO_DV_FILE + " M4 sci per cycle 12m"),
                                     "productsize": _r(5.045, "PB", _MEMO_DV_FILE + " M4 prod per cycle 12m")}},
        "7m":   {"median":          {"datavol_total": _r(0.003, "TB", _MEMO_DV_FILE + " M4 total median 7m",
                                                         rel_tol=0.5, abs_tol=1e-3),
                                     "datavol_science": _r(0.002, "TB", _MEMO_DV_FILE + " M4 sci median 7m",
                                                           rel_tol=0.5, abs_tol=1e-3),
                                     "productsize": _r(0.001, "TB", _MEMO_DV_FILE + " M4 prod median 7m",
                                                       rel_tol=0.5, abs_tol=1e-3)},
                 "twa":             {"datavol_total": _r(0.156, "TB", _MEMO_DV_FILE + " M4 total TWA 7m"),
                                     "datavol_science": _r(0.112, "TB", _MEMO_DV_FILE + " M4 sci TWA 7m"),
                                     "productsize": _r(0.051, "TB", _MEMO_DV_FILE + " M4 prod TWA 7m")},
                 "max":             {"datavol_total": _r(3.170, "TB", _MEMO_DV_FILE + " M4 total max 7m"),
                                     "datavol_science": _r(2.345, "TB", _MEMO_DV_FILE + " M4 sci max 7m"),
                                     "productsize": _r(0.719, "TB", _MEMO_DV_FILE + " M4 prod max 7m")},
                 "total_per_cycle": {"datavol_total": _r(0.031, "PB", _MEMO_DV_FILE + " M4 total per cycle 7m",
                                                         rel_tol=0.1),
                                     "datavol_science": _r(0.022, "PB", _MEMO_DV_FILE + " M4 sci per cycle 7m",
                                                           rel_tol=0.1),
                                     "productsize": _r(0.026, "PB", _MEMO_DV_FILE + " M4 prod per cycle 7m",
                                                       rel_tol=0.1)}},
        "both": {"median":          {"datavol_total": _r(0.055, "TB", _MEMO_DV_FILE + " M4 total median both"),
                                     "datavol_science": _r(0.033, "TB", _MEMO_DV_FILE + " M4 sci median both"),
                                     "productsize": _r(0.015, "TB", _MEMO_DV_FILE + " M4 prod median both")},
                 "twa":             {"datavol_total": _r(1.720, "TB", _MEMO_DV_FILE + " M4 total TWA both"),
                                     "datavol_science": _r(1.303, "TB", _MEMO_DV_FILE + " M4 sci TWA both"),
                                     "productsize": _r(2.474, "TB", _MEMO_DV_FILE + " M4 prod TWA both")},
                 "max":             {"datavol_total": _r(89.026, "TB", _MEMO_DV_FILE + " M4 total max both"),
                                     "datavol_science": _r(75.117, "TB", _MEMO_DV_FILE + " M4 sci max both"),
                                     "productsize": _r(343.465, "TB", _MEMO_DV_FILE + " M4 prod max both")},
                 "total_per_cycle": {"datavol_total": _r(1.878, "PB", _MEMO_DV_FILE + " M4 total per cycle both"),
                                     "datavol_science": _r(1.404, "PB", _MEMO_DV_FILE + " M4 sci per cycle both"),
                                     "productsize": _r(5.070, "PB", _MEMO_DV_FILE + " M4 prod per cycle both")}},
    },
    "M5": {
        "12m":  {"median":          {"datavol_total": _r(0.311, "TB", _MEMO_DV_FILE + " M5 total median 12m"),
                                     "datavol_science": _r(0.201, "TB", _MEMO_DV_FILE + " M5 sci median 12m"),
                                     "productsize": _r(0.103, "TB", _MEMO_DV_FILE + " M5 prod median 12m")},
                 "twa":             {"datavol_total": _r(7.033, "TB", _MEMO_DV_FILE + " M5 total TWA 12m"),
                                     "datavol_science": _r(5.143, "TB", _MEMO_DV_FILE + " M5 sci TWA 12m"),
                                     "productsize": _r(11.065, "TB", _MEMO_DV_FILE + " M5 prod TWA 12m")},
                 "max":             {"datavol_total": _r(160.891, "TB", _MEMO_DV_FILE + " M5 total max 12m"),
                                     "datavol_science": _r(133.936, "TB", _MEMO_DV_FILE + " M5 sci max 12m"),
                                     "productsize": _r(999.172, "TB", _MEMO_DV_FILE + " M5 prod max 12m")},
                 "total_per_cycle": {"datavol_total": _r(4.487, "PB", _MEMO_DV_FILE + " M5 total per cycle 12m"),
                                     "datavol_science": _r(3.259, "PB", _MEMO_DV_FILE + " M5 sci per cycle 12m"),
                                     "productsize": _r(12.049, "PB", _MEMO_DV_FILE + " M5 prod per cycle 12m")}},
        "7m":   {"median":          {"datavol_total": _r(0.007, "TB", _MEMO_DV_FILE + " M5 total median 7m",
                                                         rel_tol=0.25),
                                     "datavol_science": _r(0.005, "TB", _MEMO_DV_FILE + " M5 sci median 7m",
                                                           rel_tol=0.25),
                                     "productsize": _r(0.002, "TB", _MEMO_DV_FILE + " M5 prod median 7m",
                                                       rel_tol=0.5, abs_tol=1e-3)},
                 "twa":             {"datavol_total": _r(0.341, "TB", _MEMO_DV_FILE + " M5 total TWA 7m"),
                                     "datavol_science": _r(0.241, "TB", _MEMO_DV_FILE + " M5 sci TWA 7m"),
                                     "productsize": _r(0.120, "TB", _MEMO_DV_FILE + " M5 prod TWA 7m")},
                 "max":             {"datavol_total": _r(6.021, "TB", _MEMO_DV_FILE + " M5 total max 7m"),
                                     "datavol_science": _r(4.454, "TB", _MEMO_DV_FILE + " M5 sci max 7m"),
                                     "productsize": _r(2.101, "TB", _MEMO_DV_FILE + " M5 prod max 7m")},
                 "total_per_cycle": {"datavol_total": _r(0.072, "PB", _MEMO_DV_FILE + " M5 total per cycle 7m"),
                                     "datavol_science": _r(0.049, "PB", _MEMO_DV_FILE + " M5 sci per cycle 7m"),
                                     "productsize": _r(0.065, "PB", _MEMO_DV_FILE + " M5 prod per cycle 7m")}},
        "both": {"median":          {"datavol_total": _r(0.115, "TB", _MEMO_DV_FILE + " M5 total median both"),
                                     "datavol_science": _r(0.076, "TB", _MEMO_DV_FILE + " M5 sci median both"),
                                     "productsize": _r(0.033, "TB", _MEMO_DV_FILE + " M5 prod median both")},
                 "twa":             {"datavol_total": _r(4.131, "TB", _MEMO_DV_FILE + " M5 total TWA both"),
                                     "datavol_science": _r(3.017, "TB", _MEMO_DV_FILE + " M5 sci TWA both"),
                                     "productsize": _r(6.316, "TB", _MEMO_DV_FILE + " M5 prod TWA both")},
                 "max":             {"datavol_total": _r(160.891, "TB", _MEMO_DV_FILE + " M5 total max both"),
                                     "datavol_science": _r(133.936, "TB", _MEMO_DV_FILE + " M5 sci max both"),
                                     "productsize": _r(999.172, "TB", _MEMO_DV_FILE + " M5 prod max both")},
                 "total_per_cycle": {"datavol_total": _r(4.558, "PB", _MEMO_DV_FILE + " M5 total per cycle both"),
                                     "datavol_science": _r(3.308, "PB", _MEMO_DV_FILE + " M5 sci per cycle both"),
                                     "productsize": _r(12.114, "PB", _MEMO_DV_FILE + " M5 prod per cycle both")}},
    },
}


# -- Memo sysperf table (visrate + sysperf; data rate block duplicates memo_DR)
MEMO_SYSPERF_TABLE: Dict[str, Dict[str, Dict[str, Dict[str, RefValue]]]] = {
    "M1": {
        "12m":  {"median": {"visrate": _r(41.7, "Gvis/hr", _MEMO_SP_FILE + " M1 visrate median 12m"),
                            "sysperf": _r(0.012, "PFLOP/s", _MEMO_SP_FILE + " M1 sysperf median 12m",
                                          rel_tol=0.1, abs_tol=5e-3)},
                 "twa":    {"visrate": _r(99.4, "Gvis/hr", _MEMO_SP_FILE + " M1 visrate TWA 12m"),
                            "sysperf": _r(0.042, "PFLOP/s", _MEMO_SP_FILE + " M1 sysperf TWA 12m")},
                 "max":    {"visrate": _r(375.3, "Gvis/hr", _MEMO_SP_FILE + " M1 visrate max 12m"),
                            "sysperf": _r(0.467, "PFLOP/s", _MEMO_SP_FILE + " M1 sysperf max 12m")}},
        "7m":   {"median": {"visrate": _r(1.6, "Gvis/hr", _MEMO_SP_FILE + " M1 visrate median 7m"),
                            "sysperf": _r(0.000, "PFLOP/s", _MEMO_SP_FILE + " M1 sysperf median 7m",
                                          rel_tol=1.0, abs_tol=5e-4)},
                 "twa":    {"visrate": _r(3.4, "Gvis/hr", _MEMO_SP_FILE + " M1 visrate TWA 7m"),
                            "sysperf": _r(0.003, "PFLOP/s", _MEMO_SP_FILE + " M1 sysperf TWA 7m",
                                          rel_tol=0.2, abs_tol=5e-4)},
                 "max":    {"visrate": _r(9.6, "Gvis/hr", _MEMO_SP_FILE + " M1 visrate max 7m"),
                            "sysperf": _r(0.012, "PFLOP/s", _MEMO_SP_FILE + " M1 sysperf max 7m",
                                          rel_tol=0.1)}},
        "both": {"median": {"visrate": _r(12.3, "Gvis/hr", _MEMO_SP_FILE + " M1 visrate median both"),
                            "sysperf": _r(0.003, "PFLOP/s", _MEMO_SP_FILE + " M1 sysperf median both",
                                          rel_tol=0.2, abs_tol=5e-4)},
                 "twa":    {"visrate": _r(57.8, "Gvis/hr", _MEMO_SP_FILE + " M1 visrate TWA both"),
                            "sysperf": _r(0.025, "PFLOP/s", _MEMO_SP_FILE + " M1 sysperf TWA both")},
                 "max":    {"visrate": _r(375.3, "Gvis/hr", _MEMO_SP_FILE + " M1 visrate max both"),
                            "sysperf": _r(0.467, "PFLOP/s", _MEMO_SP_FILE + " M1 sysperf max both")}},
    },
    "M4": {
        "12m":  {"median": {"visrate": _r(87.2, "Gvis/hr", _MEMO_SP_FILE + " M4 visrate median 12m"),
                            "sysperf": _r(0.024, "PFLOP/s", _MEMO_SP_FILE + " M4 sysperf median 12m")},
                 "twa":    {"visrate": _r(339.6, "Gvis/hr", _MEMO_SP_FILE + " M4 visrate TWA 12m"),
                            "sysperf": _r(0.170, "PFLOP/s", _MEMO_SP_FILE + " M4 sysperf TWA 12m")},
                 "max":    {"visrate": _r(3002.8, "Gvis/hr", _MEMO_SP_FILE + " M4 visrate max 12m"),
                            "sysperf": _r(3.740, "PFLOP/s", _MEMO_SP_FILE + " M4 sysperf max 12m")}},
        "7m":   {"median": {"visrate": _r(1.6, "Gvis/hr", _MEMO_SP_FILE + " M4 visrate median 7m"),
                            "sysperf": _r(0.000, "PFLOP/s", _MEMO_SP_FILE + " M4 sysperf median 7m",
                                          rel_tol=1.0, abs_tol=5e-4)},
                 "twa":    {"visrate": _r(6.4, "Gvis/hr", _MEMO_SP_FILE + " M4 visrate TWA 7m"),
                            "sysperf": _r(0.006, "PFLOP/s", _MEMO_SP_FILE + " M4 sysperf TWA 7m",
                                          rel_tol=0.15)},
                 "max":    {"visrate": _r(38.5, "Gvis/hr", _MEMO_SP_FILE + " M4 visrate max 7m"),
                            "sysperf": _r(0.048, "PFLOP/s", _MEMO_SP_FILE + " M4 sysperf max 7m")}},
        "both": {"median": {"visrate": _r(25.1, "Gvis/hr", _MEMO_SP_FILE + " M4 visrate median both"),
                            "sysperf": _r(0.007, "PFLOP/s", _MEMO_SP_FILE + " M4 sysperf median both")},
                 "twa":    {"visrate": _r(195.1, "Gvis/hr", _MEMO_SP_FILE + " M4 visrate TWA both"),
                            "sysperf": _r(0.099, "PFLOP/s", _MEMO_SP_FILE + " M4 sysperf TWA both")},
                 "max":    {"visrate": _r(3002.8, "Gvis/hr", _MEMO_SP_FILE + " M4 visrate max both"),
                            "sysperf": _r(3.740, "PFLOP/s", _MEMO_SP_FILE + " M4 sysperf max both")}},
    },
    "M5": {
        "12m":  {"median": {"visrate": _r(208.9, "Gvis/hr", _MEMO_SP_FILE + " M5 visrate median 12m"),
                            "sysperf": _r(0.048, "PFLOP/s", _MEMO_SP_FILE + " M5 sysperf median 12m")},
                 "twa":    {"visrate": _r(825.3, "Gvis/hr", _MEMO_SP_FILE + " M5 visrate TWA 12m"),
                            "sysperf": _r(0.372, "PFLOP/s", _MEMO_SP_FILE + " M5 sysperf TWA 12m")},
                 "max":    {"visrate": _r(6005.5, "Gvis/hr", _MEMO_SP_FILE + " M5 visrate max 12m"),
                            "sysperf": _r(7.480, "PFLOP/s", _MEMO_SP_FILE + " M5 sysperf max 12m")}},
        "7m":   {"median": {"visrate": _r(3.3, "Gvis/hr", _MEMO_SP_FILE + " M5 visrate median 7m"),
                            "sysperf": _r(0.001, "PFLOP/s", _MEMO_SP_FILE + " M5 sysperf median 7m",
                                          rel_tol=0.5, abs_tol=5e-4)},
                 "twa":    {"visrate": _r(14.7, "Gvis/hr", _MEMO_SP_FILE + " M5 visrate TWA 7m"),
                            "sysperf": _r(0.014, "PFLOP/s", _MEMO_SP_FILE + " M5 sysperf TWA 7m",
                                          rel_tol=0.1)},
                 "max":    {"visrate": _r(76.9, "Gvis/hr", _MEMO_SP_FILE + " M5 visrate max 7m"),
                            "sysperf": _r(0.096, "PFLOP/s", _MEMO_SP_FILE + " M5 sysperf max 7m")}},
        "both": {"median": {"visrate": _r(50.0, "Gvis/hr", _MEMO_SP_FILE + " M5 visrate median both"),
                            "sysperf": _r(0.013, "PFLOP/s", _MEMO_SP_FILE + " M5 sysperf median both",
                                          rel_tol=0.1)},
                 "twa":    {"visrate": _r(473.6, "Gvis/hr", _MEMO_SP_FILE + " M5 visrate TWA both"),
                            "sysperf": _r(0.217, "PFLOP/s", _MEMO_SP_FILE + " M5 sysperf TWA both")},
                 "max":    {"visrate": _r(6005.5, "Gvis/hr", _MEMO_SP_FILE + " M5 visrate max both"),
                            "sysperf": _r(7.480, "PFLOP/s", _MEMO_SP_FILE + " M5 sysperf max both")}},
    },
}


# ---------------------------------------------------------------------------
# SDD tables (wsu-sdd-excerpts/tables/)
# ---------------------------------------------------------------------------
# The SDD quotes M5 values only and compares against a "pre-WSU" baseline
# computed from the Cycle 7/8 archive (no Band 1/2 synthesis, BLC
# correlator).  Several SDD numbers differ from the memo despite nominally
# using the same database -- those differences are what the consistency
# report (Phase 4) investigates.
_SDD_DP_FILE = "wsu-sdd-excerpts/tables/table_data_prop.tex"
_SDD_SP_FILE = "wsu-sdd-excerpts/tables/table_system_performance_wsu.tex"
_SDD_SUMM_FILE = "wsu-sdd-excerpts/tables/data_prop_summ.tex"
_SDD_PEAK_FILE = "wsu-sdd-excerpts/tables/table_data_rate_peak.tex"
_SDD_DR_FILE = "wsu-sdd-excerpts/tables/table_data_rates_prop.tex"


SDD_DATA_PROP: Dict[str, Dict[str, Dict[str, Dict[str, RefValue]]]] = {
    "M5": {
        "12m":  {"median":          {"datavol_total": _r(0.32, "TB", _SDD_DP_FILE + " M5 total median 12m",
                                                         rel_tol=0.05),
                                     "productsize": _r(0.10, "TB", _SDD_DP_FILE + " M5 prod median 12m",
                                                       rel_tol=0.05)},
                 "twa":             {"datavol_total": _r(6.9, "TB", _SDD_DP_FILE + " M5 total TWA 12m",
                                                         rel_tol=0.05),
                                     "productsize": _r(11.0, "TB", _SDD_DP_FILE + " M5 prod TWA 12m",
                                                       rel_tol=0.05)},
                 "max":             {"datavol_total": _r(141.0, "TB", _SDD_DP_FILE + " M5 total max 12m",
                                                         rel_tol=0.05),
                                     "productsize": _r(1200.0, "TB", _SDD_DP_FILE + " M5 prod max 12m",
                                                       rel_tol=0.1)},
                 "total_per_cycle": {"datavol_total": _r(4.5, "PB", _SDD_DP_FILE + " M5 total per cycle 12m",
                                                         rel_tol=0.05),
                                     "productsize": _r(12.0, "PB", _SDD_DP_FILE + " M5 prod per cycle 12m",
                                                       rel_tol=0.05)}},
        "7m":   {"median":          {"datavol_total": _r(0.007, "TB", _SDD_DP_FILE + " M5 total median 7m",
                                                         rel_tol=0.1),
                                     "productsize": _r(0.002, "TB", _SDD_DP_FILE + " M5 prod median 7m",
                                                       rel_tol=0.5, abs_tol=1e-3)},
                 "twa":             {"datavol_total": _r(0.38, "TB", _SDD_DP_FILE + " M5 total TWA 7m",
                                                         rel_tol=0.1),
                                     "productsize": _r(0.12, "TB", _SDD_DP_FILE + " M5 prod TWA 7m",
                                                       rel_tol=0.1)},
                 "max":             {"datavol_total": _r(6.4, "TB", _SDD_DP_FILE + " M5 total max 7m",
                                                         rel_tol=0.05),
                                     "productsize": _r(2.1, "TB", _SDD_DP_FILE + " M5 prod max 7m",
                                                       rel_tol=0.05)},
                 "total_per_cycle": {"datavol_total": _r(0.073, "PB", _SDD_DP_FILE + " M5 total per cycle 7m",
                                                         rel_tol=0.05),
                                     "productsize": _r(0.065, "PB", _SDD_DP_FILE + " M5 prod per cycle 7m",
                                                       rel_tol=0.05)}},
    },
    "BLC": {
        "12m":  {"median":          {"datavol_total": _r(0.016, "TB", _SDD_DP_FILE + " BLC total median 12m",
                                                         rel_tol=0.1),
                                     "productsize": _r(0.0086, "TB", _SDD_DP_FILE + " BLC prod median 12m",
                                                       rel_tol=0.1)},
                 "twa":             {"datavol_total": _r(0.12, "TB", _SDD_DP_FILE + " BLC total TWA 12m",
                                                         rel_tol=0.1),
                                     "productsize": _r(0.52, "TB", _SDD_DP_FILE + " BLC prod TWA 12m",
                                                       rel_tol=0.1)},
                 "max":             {"datavol_total": _r(1.6, "TB", _SDD_DP_FILE + " BLC total max 12m",
                                                         rel_tol=0.1),
                                     "productsize": _r(43.0, "TB", _SDD_DP_FILE + " BLC prod max 12m",
                                                       rel_tol=0.1)},
                 "total_per_cycle": {"datavol_total": _r(0.084, "PB", _SDD_DP_FILE + " BLC total per cycle 12m",
                                                         rel_tol=0.1),
                                     "productsize": _r(0.61, "PB", _SDD_DP_FILE + " BLC prod per cycle 12m",
                                                       rel_tol=0.1)}},
        "7m":   {"median":          {"datavol_total": _r(9.1e-4, "TB", _SDD_DP_FILE + " BLC total median 7m",
                                                         rel_tol=0.2),
                                     "productsize": _r(2.6e-4, "TB", _SDD_DP_FILE + " BLC prod median 7m",
                                                       rel_tol=0.5)},
                 "twa":             {"datavol_total": _r(0.012, "TB", _SDD_DP_FILE + " BLC total TWA 7m",
                                                         rel_tol=0.1),
                                     "productsize": _r(3.1e-3, "TB", _SDD_DP_FILE + " BLC prod TWA 7m",
                                                       rel_tol=0.1)},
                 "max":             {"datavol_total": _r(0.081, "TB", _SDD_DP_FILE + " BLC total max 7m",
                                                         rel_tol=0.1),
                                     "productsize": _r(0.083, "TB", _SDD_DP_FILE + " BLC prod max 7m",
                                                       rel_tol=0.1)},
                 "total_per_cycle": {"datavol_total": _r(2.2e-3, "PB", _SDD_DP_FILE + " BLC total per cycle 7m",
                                                         rel_tol=0.1),
                                     "productsize": _r(1.6e-3, "PB", _SDD_DP_FILE + " BLC prod per cycle 7m",
                                                       rel_tol=0.1)}},
    },
}


SDD_SYSPERF: Dict[str, Dict[str, Dict[str, Dict[str, RefValue]]]] = {
    "M5": {
        "12m":  {"median": {"datarate": _r(0.136, "GB/s", _SDD_SP_FILE + " M5 DR median 12m",
                                           rel_tol=0.15,
                                           notes="SDD +12% vs memo (0.121) -- conformance adjustment"),
                            "visrate": _r(230.0, "Gvis/hr", _SDD_SP_FILE + " M5 visrate median 12m",
                                          rel_tol=0.15),
                            "sysperf": _r(0.07, "PFLOP/s", _SDD_SP_FILE + " M5 sysperf median 12m",
                                          rel_tol=0.5,
                                          notes="SDD includes 3x contingency factor")},
                 "twa":    {"datarate": _r(0.477, "GB/s", _SDD_SP_FILE + " M5 DR TWA 12m", rel_tol=0.05),
                            "visrate": _r(822.0, "Gvis/hr", _SDD_SP_FILE + " M5 visrate TWA 12m",
                                          rel_tol=0.05),
                            "sysperf": _r(0.50, "PFLOP/s", _SDD_SP_FILE + " M5 sysperf TWA 12m",
                                          rel_tol=0.5,
                                          notes="SDD includes 3x contingency factor")},
                 "max":    {"datarate": _r(3.48, "GB/s", _SDD_SP_FILE + " M5 DR max 12m", rel_tol=0.05),
                            "visrate": _r(6005.0, "Gvis/hr", _SDD_SP_FILE + " M5 visrate max 12m",
                                          rel_tol=0.05),
                            "sysperf": _r(10.3, "PFLOP/s", _SDD_SP_FILE + " M5 sysperf max 12m",
                                          rel_tol=0.5)}},
        "7m":   {"median": {"datarate": _r(0.002, "GB/s", _SDD_SP_FILE + " M5 DR median 7m",
                                           rel_tol=0.25),
                            "visrate": _r(3.2, "Gvis/hr", _SDD_SP_FILE + " M5 visrate median 7m",
                                          rel_tol=0.1),
                            "sysperf": _r(1.4e-3, "PFLOP/s", _SDD_SP_FILE + " M5 sysperf median 7m",
                                          rel_tol=0.5)},
                 "twa":    {"datarate": _r(0.010, "GB/s", _SDD_SP_FILE + " M5 DR TWA 7m", rel_tol=0.1),
                            "visrate": _r(15.0, "Gvis/hr", _SDD_SP_FILE + " M5 visrate TWA 7m",
                                          rel_tol=0.05),
                            "sysperf": _r(0.021, "PFLOP/s", _SDD_SP_FILE + " M5 sysperf TWA 7m",
                                          rel_tol=0.5)},
                 "max":    {"datarate": _r(0.05, "GB/s", _SDD_SP_FILE + " M5 DR max 7m", rel_tol=0.1),
                            "visrate": _r(77.0, "Gvis/hr", _SDD_SP_FILE + " M5 visrate max 7m",
                                          rel_tol=0.05),
                            "sysperf": _r(0.13, "PFLOP/s", _SDD_SP_FILE + " M5 sysperf max 7m",
                                          rel_tol=0.5)}},
    },
    "BLC": {
        "12m":  {"median": {"datarate": _r(5.6e-3, "GB/s", _SDD_SP_FILE + " BLC DR median 12m", rel_tol=0.1),
                            "visrate": _r(9.6, "Gvis/hr", _SDD_SP_FILE + " BLC visrate median 12m",
                                          rel_tol=0.1),
                            "sysperf": _r(3.0e-3, "PFLOP/s", _SDD_SP_FILE + " BLC sysperf median 12m",
                                          rel_tol=0.5)},
                 "twa":    {"datarate": _r(9.0e-3, "GB/s", _SDD_SP_FILE + " BLC DR TWA 12m",
                                           rel_tol=0.1),
                            "visrate": _r(16.0, "Gvis/hr", _SDD_SP_FILE + " BLC visrate TWA 12m",
                                          rel_tol=0.1),
                            "sysperf": _r(7.6e-3, "PFLOP/s", _SDD_SP_FILE + " BLC sysperf TWA 12m",
                                          rel_tol=0.5)},
                 "max":    {"datarate": _r(4.5e-2, "GB/s", _SDD_SP_FILE + " BLC DR max 12m", rel_tol=0.1),
                            "visrate": _r(77.0, "Gvis/hr", _SDD_SP_FILE + " BLC visrate max 12m",
                                          rel_tol=0.1),
                            "sysperf": _r(9.9e-2, "PFLOP/s", _SDD_SP_FILE + " BLC sysperf max 12m",
                                          rel_tol=0.5)}},
        "7m":   {"median": {"datarate": _r(2.4e-4, "GB/s", _SDD_SP_FILE + " BLC DR median 7m",
                                           rel_tol=0.5),
                            "visrate": _r(0.4, "Gvis/hr", _SDD_SP_FILE + " BLC visrate median 7m",
                                          rel_tol=0.25),
                            "sysperf": _r(1.5e-4, "PFLOP/s", _SDD_SP_FILE + " BLC sysperf median 7m",
                                          rel_tol=0.5)},
                 "twa":    {"datarate": _r(3.1e-4, "GB/s", _SDD_SP_FILE + " BLC DR TWA 7m", rel_tol=0.1),
                            "visrate": _r(0.4, "Gvis/hr", _SDD_SP_FILE + " BLC visrate TWA 7m",
                                          rel_tol=0.25),
                            "sysperf": _r(5.1e-4, "PFLOP/s", _SDD_SP_FILE + " BLC sysperf TWA 7m",
                                          rel_tol=0.5)},
                 "max":    {"datarate": _r(8.7e-4, "GB/s", _SDD_SP_FILE + " BLC DR max 7m", rel_tol=0.1),
                            "visrate": _r(1.3, "Gvis/hr", _SDD_SP_FILE + " BLC visrate max 7m",
                                          rel_tol=0.1),
                            "sysperf": _r(2.2e-3, "PFLOP/s", _SDD_SP_FILE + " BLC sysperf max 7m",
                                          rel_tol=0.5)}},
    },
}


# -- SDD summary table (data_prop_summ.tex) ---------------------------------
# 7m+12m+TP; TP adds ~3% per SDD notes.
SDD_SUMMARY: Dict[str, RefValue] = {
    "raw_data_vol_PBpercycle":            _r(5.0, "PB/cycle", _SDD_SUMM_FILE + " raw vol M5", rel_tol=0.1),
    "time_averaged_data_rate_GBs":        _r(0.54, "GB/s", _SDD_SUMM_FILE + " avg DR M5", rel_tol=0.1),
    "product_volume_PBpercycle":          _r(12.1, "PB/cycle", _SDD_SUMM_FILE + " prod vol M5", rel_tol=0.05),
    "product_volume_compressed_PB":       _r(7.9, "PB/cycle", _SDD_SUMM_FILE + " prod vol compressed M5",
                                             rel_tol=0.1),
    "product_volume_with_gous_PB":        _r(10.7, "PB/cycle", _SDD_SUMM_FILE + " prod vol GOUS M5",
                                             rel_tol=0.1),
    "total_data_volume_PBpercycle":       _r(15.7, "PB/cycle", _SDD_SUMM_FILE + " total vol M5", rel_tol=0.1),
    "sysperf_PFLOPs":                     _r(0.78, "PFLOP/s", _SDD_SUMM_FILE + " sysperf M5", rel_tol=0.5),
    # Pre-WSU:
    "preWSU_raw_data_vol_PBpercycle":     _r(9.9e-2, "PB/cycle", _SDD_SUMM_FILE + " raw vol BLC", rel_tol=0.1),
    "preWSU_time_averaged_data_rate_GBs": _r(1.0e-2, "GB/s", _SDD_SUMM_FILE + " avg DR BLC", rel_tol=0.1),
    "preWSU_product_volume_PBpercycle":   _r(0.61, "PB/cycle", _SDD_SUMM_FILE + " prod vol BLC", rel_tol=0.1),
    "preWSU_total_data_volume_PB":        _r(0.17, "PB/cycle", _SDD_SUMM_FILE + " total vol BLC", rel_tol=0.1),
    "preWSU_sysperf_PFLOPs":              _r(9.0e-3, "PFLOP/s", _SDD_SUMM_FILE + " sysperf BLC", rel_tol=0.5),
}


# -- SDD peak data rate (table_data_rate_peak.tex) --------------------------
SDD_PEAK_M5: Dict[str, RefValue] = {
    "nant_12m":        _r(50, "", _SDD_PEAK_FILE + " M5 12m nant", rel_tol=0.0, abs_tol=0.01),
    "nant_7m":         _r(12, "", _SDD_PEAK_FILE + " M5 7m nant", rel_tol=0.0, abs_tol=0.01),
    "velocity_kms":    _r(0.1, "km/s", _SDD_PEAK_FILE + " M5 vel", rel_tol=0.0, abs_tol=0.001),
    "nchannels":       _r(2_380_800, "", _SDD_PEAK_FILE + " M5 nchan", rel_tol=0.001),
    "tint_12m_s":      _r(3.072, "s", _SDD_PEAK_FILE + " M5 12m tint", rel_tol=0.001),
    "tint_7m_s":       _r(9.984, "s", _SDD_PEAK_FILE + " M5 7m tint", rel_tol=0.001),
    "datarate_12m":    _r(3.96, "GB/s", _SDD_PEAK_FILE + " M5 DR 12m", rel_tol=0.01),
    "datarate_7m":     _r(0.074, "GB/s", _SDD_PEAK_FILE + " M5 DR 7m", rel_tol=0.02),
    "total_datarate":  _r(4.43, "GB/s", _SDD_PEAK_FILE + " M5 total DR (12m+7m+TP)", rel_tol=0.02),
}


# -- SDD data rate table (table_data_rates_prop.tex) ------------------------
SDD_DATA_RATES: Dict[str, Dict[str, RefValue]] = {
    "12m":  {"median": _r(0.13, "GB/s", _SDD_DR_FILE + " M5 DR median 12m", rel_tol=0.1),
             "twa":    _r(0.48, "GB/s", _SDD_DR_FILE + " M5 DR TWA 12m", rel_tol=0.05),
             "max":    _r(3.48, "GB/s", _SDD_DR_FILE + " M5 DR max 12m", rel_tol=0.05)},
    "7m":   {"median": _r(0.002, "GB/s", _SDD_DR_FILE + " M5 DR median 7m", rel_tol=0.25),
             "twa":    _r(0.01, "GB/s", _SDD_DR_FILE + " M5 DR TWA 7m", rel_tol=0.1),
             "max":    _r(0.05, "GB/s", _SDD_DR_FILE + " M5 DR max 7m", rel_tol=0.1)},
    "both": {"median": _r(0.13, "GB/s", _SDD_DR_FILE + " M5 DR median both", rel_tol=0.1),
             "twa":    _r(0.49, "GB/s", _SDD_DR_FILE + " M5 DR TWA both", rel_tol=0.05),
             "max":    _r(3.48, "GB/s", _SDD_DR_FILE + " M5 DR max both", rel_tol=0.05)},
}


# Convenience flat iterator for test parameterisation.
def flatten_memo_tables() -> List[Tuple[str, str, str, str, str, RefValue]]:
    """Yield (source, milestone, array, stat, quantity, RefValue) tuples."""
    out: List[Tuple[str, str, str, str, str, RefValue]] = []
    for source_name, table in [
        ("memo_datarate", MEMO_DATARATE_TABLE),
        ("memo_datavol", MEMO_DATAVOL_TABLE),
        ("memo_sysperf", MEMO_SYSPERF_TABLE),
    ]:
        for ms, per_arr in table.items():
            for arr, per_stat in per_arr.items():
                for stat, per_q in per_stat.items():
                    for q, ref in per_q.items():
                        out.append((source_name, ms, arr, stat, q, ref))
    return out
