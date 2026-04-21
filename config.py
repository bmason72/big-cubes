"""
Single source of truth for WSU data-properties projections.

All constants and per-band / per-milestone parameters live here so the
pipeline can be re-run with alternate assumptions without touching the
calculation code. Values are transcribed from wsu_db.py (as of 2026-04-20).

Sources for each block are called out in comments so future updates can be
traced back to the code path they came from.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import astropy.units as u


# ---------------------------------------------------------------------------
# TALON / ATAC base parameters
# ---------------------------------------------------------------------------
# From wsu_db.py:calc_talon_specwidth (line ~56) and obsPrepTsiReport.txt.
TALON_CHANNEL_KHZ = 13.5  # base channel size out of the correlator
SPW_BANDWIDTH_GHZ = 2.0   # per-spw bandwidth

# Minimum channel-averaging factor per band.  Two flavours: the generic
# stepped2 minimum used at M4/M5 (wsu_chanavg_min) and the tighter minimum
# that applies to the Milestone-1 (initial) correlator restrictions.
# Extracted from wsu_db.py lines 17-39.
WSU_CHANAVG_MIN: Dict[int, float] = {
    1: 1.0, 2: 2.0, 3: 3.0, 4: 5.0, 5: 6.0,
    6: 8.0, 7: 10.0, 8: 15.0, 9: 20.0, 10: 32.0,
}

WSU_CHANAVG_MIN_INITIAL: Dict[int, float] = {
    1: 4.0, 2: 8.0, 3: 8.0, 4: 5.0, 5: 6.0,
    6: 8.0, 7: 10.0, 8: 15.0, 9: 20.0, 10: 32.0,
}


# ---------------------------------------------------------------------------
# Antennas and integration times
# ---------------------------------------------------------------------------
# Typical / peak / all(=12m+7m+TP) counts are used for different scenarios.
# wsu_db.py lines 230-242 and 4327-4336.
@dataclass(frozen=True)
class ArrayParams:
    name: str
    nant_typical: int
    nant_array: int   # peak configuration (50 for 12m, 12 for 7m)
    nant_all: int     # 12m+7m+TP for 12m; 7m+TP for 7m
    tint_wsu_s: float   # M4/M5 WSU integration time (seconds)
    tint_initial_s: float  # M1 integration time
    tint_blc_7m_s: float = 10.1  # only used for 7m
    tint_blc_12m_short_s: float = 6.05
    tint_blc_12m_long_s: float = 3.024
    blc_tint_breakpt_m: float = 3000.0   # L80 split between short/long


ARRAY_12M = ArrayParams(
    name="12m",
    nant_typical=47,
    nant_array=50,
    nant_all=66,
    tint_wsu_s=3.072,
    tint_initial_s=6.144,
)

ARRAY_7M = ArrayParams(
    name="7m",
    nant_typical=10,
    nant_array=12,
    nant_all=16,
    tint_wsu_s=9.984,
    tint_initial_s=9.984,
)

ARRAYS: Dict[str, ArrayParams] = {"12m": ARRAY_12M, "7m": ARRAY_7M}


# ---------------------------------------------------------------------------
# Milestone definitions (bandwidth per band, npol, nant/tint override)
# ---------------------------------------------------------------------------
# Bandwidth-per-band dictionaries from wsu_db.py create_initial_wsu_db()
# (lines ~4340-4454).  Values apply to 12m and 7m unless noted.
#
# Bands 1 and 2 are *not* in the base database; they are synthesised by
# add_bands() (wsu_db.py line 3084) and use different per-milestone
# bandwidths from the band 3-10 dictionaries below, so we keep a separate
# lookup for them.

@dataclass(frozen=True)
class MilestoneBands:
    """Bandwidth (GHz) per band for a given milestone.  Missing bands raise."""

    name: str
    bandwidth_ghz: Dict[int, float]
    npol: int = 2
    uses_initial_chanavg_min: bool = False  # True for M1 only
    # Integration-time / antenna overrides relative to ArrayParams.
    # M1 uses longer 6.144s at 12m + 47 ants; M4/M5 use the WSU default.
    use_initial_tint: bool = False

    def bandwidth_for(self, band: int) -> float:
        if band not in self.bandwidth_ghz:
            raise KeyError(f"band {band} missing from milestone {self.name}")
        return self.bandwidth_ghz[band]


# Milestone 1 (initial).
M1 = MilestoneBands(
    name="M1",
    bandwidth_ghz={
        1: 8.0,    # band 1  (from add_bands)
        2: 16.0,   # band 2  (from add_bands)
        3: 16.0,   # "band 2 (high)" / band 3 proper
        4: 8.0, 5: 8.0, 6: 11.0, 7: 8.0, 8: 8.0,
        9: 16.0, 10: 16.0,
    },
    npol=2,
    uses_initial_chanavg_min=True,
    use_initial_tint=True,
)

# Milestone 4 (mid-stage; 2x bandwidth for most bands vs BLC).
M4 = MilestoneBands(
    name="M4",
    bandwidth_ghz={
        1: 8.0, 2: 16.0,
        3: 16.0, 4: 8.0, 5: 8.0, 6: 11.0, 7: 8.0, 8: 8.0,
        9: 16.0, 10: 16.0,
    },
)

# Milestone 5 (goal; bands 2, 6, 7, 8 at 32 GHz).
M5 = MilestoneBands(
    name="M5",
    bandwidth_ghz={
        1: 8.0,   # stays at 8
        2: 32.0,  # band 2 low freq side
        3: 32.0,  # upper part of band 2 / band 3
        4: 8.0, 5: 8.0,
        6: 32.0, 7: 32.0, 8: 32.0,
        9: 16.0, 10: 16.0,
    },
)

MILESTONES: Dict[str, MilestoneBands] = {"M1": M1, "M4": M4, "M5": M5}


# ---------------------------------------------------------------------------
# Database column name mapping
# ---------------------------------------------------------------------------
# The base DB uses asymmetric column suffixes for M1 (uses "_initial" for
# both antenna count and channel min) vs M4/M5 (uses "_typical" antenna
# count).  Keep the mapping in one place so downstream code doesn't need to
# know about the naming quirk.
@dataclass(frozen=True)
class MilestoneColumns:
    datarate: str
    visrate: str
    datavol_total: str
    datavol_target_tot: str  # science-only data volume
    productsize: str
    sysperf: str
    nchan_agg: str


MILESTONE_COLUMNS: Dict[str, MilestoneColumns] = {
    "M1": MilestoneColumns(
        datarate="wsu_datarate_initial_stepped2_initial",
        visrate="wsu_visrate_initial_stepped2_initial",
        datavol_total="wsu_datavol_initial_stepped2_initial_total",
        datavol_target_tot="wsu_datavol_initial_stepped2_initial_target_tot",
        productsize="wsu_productsize_initial_stepped2",
        sysperf="wsu_sysperf_initial_stepped2_initial_aprojonly",
        nchan_agg="wsu_nchan_agg_stepped2_initial",
    ),
    "M4": MilestoneColumns(
        datarate="wsu_datarate_ms4_stepped2_typical",
        visrate="wsu_visrate_ms4_stepped2_typical",
        datavol_total="wsu_datavol_ms4_stepped2_typical_total",
        datavol_target_tot="wsu_datavol_ms4_stepped2_typical_target_tot",
        productsize="wsu_productsize_ms4_stepped2",
        sysperf="wsu_sysperf_ms4_stepped2_typical_aprojonly",
        nchan_agg="wsu_nchan_agg_stepped2_ms4",
    ),
    "M5": MilestoneColumns(
        datarate="wsu_datarate_goal_stepped2_typical",
        visrate="wsu_visrate_goal_stepped2_typical",
        datavol_total="wsu_datavol_goal_stepped2_typical_total",
        datavol_target_tot="wsu_datavol_goal_stepped2_typical_target_tot",
        productsize="wsu_productsize_goal_stepped2",
        sysperf="wsu_sysperf_goal_stepped2_typical_aprojonly",
        nchan_agg="wsu_nchan_agg_stepped2_goal",
    ),
}

BLC_COLUMNS = MilestoneColumns(
    datarate="blc_datarate_typical",
    visrate="blc_visrate_typical",
    datavol_total="blc_datavol_typical_total",
    datavol_target_tot="blc_datavol_typical_target_tot",
    productsize="blc_productsize",
    sysperf="blc_sysperf_typical_aprojonly",
    nchan_agg="blc_nchan_agg",
)


# ---------------------------------------------------------------------------
# Velocity-resolution binning (stepped2)
# ---------------------------------------------------------------------------
# The "stepped2" binning maps each MOUS's requested velocity resolution to
# the floor of the bin it falls in: >10 -> 10; 2-10 -> 2; 0.5-2 -> 0.5;
# 0.1-0.5 -> 0.1; <0.1 -> native.  Stored as (upper-bound, value-assigned).
STEPPED2_BINS_KMS: List[Tuple[float, float]] = [
    (float("inf"), 10.0),
    (10.0, 2.0),
    (2.0, 0.5),
    (0.5, 0.1),
    (0.1, 0.0),  # 0.0 = use native (finest) resolution
]


# ---------------------------------------------------------------------------
# Imaging parameters (cube / product size formulas)
# ---------------------------------------------------------------------------
# From wsu_db.py calc_cube_size / calc_mfs_size; values are float32 cubes
# plus a factor-of-2 to include calibrators + weblog etc.
PIXELS_PER_BEAM = 5
CUBE_BYTES_PER_VOXEL = 4.0  # float32
PRODUCT_FACTOR = 2.0        # final product vs raw cube


# ---------------------------------------------------------------------------
# System-performance parameters
# ---------------------------------------------------------------------------
# All from wsu_db.py calc_sysperf (line ~3947).  Defaults mirror the memo.
@dataclass(frozen=True)
class SysperfParams:
    k_major_cycles: int = 20
    multiscale_factor: float = 1.2
    core_efficiency: float = 0.05
    parallelization_efficiency: float = 0.8
    flops_per_vis_std: float = 1280.8
    flops_per_vis_aproj: float = 7472.8
    flops_per_vis_wproj: float = 21768.4
    flops_per_vis_awproj: float = 39704.8
    # Band-1 long-baseline wproject threshold (m).
    wproject_L80_threshold_m: float = 6200.0


SYSPERF = SysperfParams()


# ---------------------------------------------------------------------------
# Peak data rate check parameters (SDD table_data_rate_peak.tex)
# ---------------------------------------------------------------------------
# These reproduce the SDD Milestone-5 peak-rate entry for unit tests.
@dataclass(frozen=True)
class PeakDataRateParams:
    nbyte: float
    napc: float
    nant: int
    velocity_kms: float
    # SDD table_data_rate_peak.tex labels this "Total channels per
    # polarization = 2,380,800".  That label is misleading: the number
    # reproduces the quoted 3.96 GB/s only when used as the *total* channel
    # count (i.e. already summed across both pols).  We therefore store it
    # as n_channels_total and pass Npols=1 to calc_datarate.
    n_channels_total: int
    tint_s: float
    expected_GBs: float


PEAK_M5_12M = PeakDataRateParams(
    nbyte=2.0, napc=1.0, nant=50, velocity_kms=0.1,
    n_channels_total=2_380_800, tint_s=3.072,
    expected_GBs=3.96,
)
PEAK_M5_7M = PeakDataRateParams(
    nbyte=2.0, napc=1.0, nant=12, velocity_kms=0.1,
    n_channels_total=2_380_800, tint_s=9.984,
    expected_GBs=0.074,
)


# ---------------------------------------------------------------------------
# Band 1/2 Monte-Carlo realizations
# ---------------------------------------------------------------------------
# generate_db_realizations() replaces some band 3/6/7 time with synthetic
# band 1/2 observations.  Defaults mirror the notebook
# calculate_initial_WSU_data_properties.ipynb which produced the memo tables.
@dataclass(frozen=True)
class RealizationParams:
    n_realizations: int = 10
    frac_12m_replaced: float = 0.10   # 10% of 12m time -> band1+band2
    frac_7m_replaced: float = 0.06    # 6% of 7m time  -> band1+band2
    outdir: str = "data/sample_band1_band2"
    filename: str = "wsu_datarates_initial_goal"
    seed: int = 20250423


REALIZATIONS = RealizationParams()


# ---------------------------------------------------------------------------
# Full bundled config
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class PipelineConfig:
    talon_channel_kHz: float = TALON_CHANNEL_KHZ
    spw_bandwidth_GHz: float = SPW_BANDWIDTH_GHZ
    chanavg_min: Dict[int, float] = field(default_factory=lambda: dict(WSU_CHANAVG_MIN))
    chanavg_min_initial: Dict[int, float] = field(
        default_factory=lambda: dict(WSU_CHANAVG_MIN_INITIAL)
    )
    arrays: Dict[str, ArrayParams] = field(default_factory=lambda: dict(ARRAYS))
    milestones: Dict[str, MilestoneBands] = field(default_factory=lambda: dict(MILESTONES))
    milestone_columns: Dict[str, MilestoneColumns] = field(
        default_factory=lambda: dict(MILESTONE_COLUMNS)
    )
    blc_columns: MilestoneColumns = field(default_factory=lambda: BLC_COLUMNS)
    stepped2_bins_kms: List[Tuple[float, float]] = field(
        default_factory=lambda: list(STEPPED2_BINS_KMS)
    )
    pixels_per_beam: int = PIXELS_PER_BEAM
    cube_bytes_per_voxel: float = CUBE_BYTES_PER_VOXEL
    product_factor: float = PRODUCT_FACTOR
    sysperf: SysperfParams = field(default_factory=lambda: SYSPERF)
    realizations: RealizationParams = field(default_factory=lambda: REALIZATIONS)


DEFAULT_CONFIG = PipelineConfig()
