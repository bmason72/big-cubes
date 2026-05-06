"""
ALMA array configuration parameters from the Cycle-13 Technical Handbook,
Table 7.3.

Used by the scenario recompute path to classify 12m MOUSes by their L80
(80th-percentile projected baseline, meters) into a configuration label,
which in turn drives the integration-time lookup.

Baselines are projected for a source transiting at zenith.  Scenarios
may apply an elevation factor on top of the tint lookup to approximate
lower-elevation targets.
"""

from __future__ import annotations

from typing import Dict, List


# Each entry: (config_label, min_baseline_m, L05_m, L80_m, max_baseline_m)
# Source: Cy13 Technical Handbook §7.3.
CY13_CONFIGS: List[tuple] = [
    ("7-m",    8.7,   9.1,    30.7,   45.0),
    ("C43-1", 14.6,  21.4,   107.1,  160.7),
    ("C43-2", 14.6,  27.0,   143.8,  313.7),
    ("C43-3", 14.6,  37.6,   235.4,  500.2),
    ("C43-4", 14.6,  54.1,   369.2,  783.5),
    ("C43-5", 14.6,  90.9,   623.8, 1397.9),
    ("C43-6", 14.6, 148.6,  1172.5, 2516.9),
    ("C43-7", 64.0, 235.2,  1673.1, 3637.8),
    ("C43-8", 110.4, 427.3, 3527.3, 8547.7),
    ("C43-9", 367.6, 746.9, 6482.6, 13894.2),
    ("C43-10", 244.0, 1228.1, 8685.9, 16194.0),
]


# Convenience dict: label -> L80
L80_BY_CONFIG: Dict[str, float] = {label: l80 for label, *_, l80, _ in CY13_CONFIGS}


# 12m-only labels (excludes 7-m); used when classifying 12m MOUSes.
CONFIG_LABELS_12M: List[str] = [c[0] for c in CY13_CONFIGS if c[0] != "7-m"]

# 12m L80 values aligned with CONFIG_LABELS_12M.
L80_VALUES_12M: List[float] = [L80_BY_CONFIG[c] for c in CONFIG_LABELS_12M]


# "Long-baseline" label set per user directive (C-8 or longer configurations
# use the short integration time).  Short-baseline is the complement.
LONG_BASELINE_CONFIGS_12M: set = {"C43-8", "C43-9", "C43-10"}
SHORT_BASELINE_CONFIGS_12M: set = set(CONFIG_LABELS_12M) - LONG_BASELINE_CONFIGS_12M
