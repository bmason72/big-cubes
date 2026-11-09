"""Antenna-geometry extraction for the ryanCy11 listobs sample.

Parses the ``Antennas:`` table of a listobs file and computes the projected
baseline-length distribution (notably L80, the 80th-percentile baseline
length) per EB.  This is the geometry input for long-baseline classification
and synthesized-beam estimates in the Cy11-based WSU refinements; it can be
cross-checked against the pipeline_stats JSON L80 values where present.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

# "  0    DA41  A004      12.0 m   -067.45.15.9  -22.53.28.0   52.6602  -704.4161  21.7753  ..."
ANTENNA_ROW_RE = re.compile(
    r"^\s+\d+\s+(?P<name>\S+)\s+(?P<station>\S+)\s+(?P<diam>[\d.]+)\s*m\s+"
    r"\S+\s+\S+\s+"
    r"(?P<east>-?[\d.]+)\s+(?P<north>-?[\d.]+)\s+(?P<elev>-?[\d.]+)\s+"
)

ANTENNAS_HEADER_RE = re.compile(r"^Antennas:\s*(?P<count>\d+):")


@dataclass(frozen=True)
class EbGeometry:
    listobs_path: str
    n_ant: int
    n_ant_12m: int
    n_ant_7m: int
    l80_m: float
    lmax_m: float

    @property
    def array_type(self) -> str:
        return "12m" if self.n_ant_12m >= self.n_ant_7m else "7m"


def baseline_percentile(positions: list[tuple[float, float]], percentile: float) -> float:
    lengths = []
    for i in range(len(positions)):
        e1, n1 = positions[i]
        for j in range(i + 1, len(positions)):
            e2, n2 = positions[j]
            lengths.append(math.hypot(e1 - e2, n1 - n2))
    if not lengths:
        return math.nan
    lengths.sort()
    if len(lengths) == 1:
        return lengths[0]
    # linear-interpolated percentile, matching numpy's default
    rank = (percentile / 100.0) * (len(lengths) - 1)
    lo = int(math.floor(rank))
    hi = int(math.ceil(rank))
    frac = rank - lo
    return lengths[lo] * (1.0 - frac) + lengths[hi] * frac


def parse_listobs_geometry(path: Path) -> EbGeometry:
    """Extract antenna positions from a listobs file and compute L80/Lmax."""
    positions: list[tuple[float, float]] = []
    diams: list[float] = []
    in_antennas = False
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if ANTENNAS_HEADER_RE.match(raw):
            in_antennas = True
            continue
        if not in_antennas:
            continue
        match = ANTENNA_ROW_RE.match(raw)
        if match is None:
            if positions:
                break  # end of antenna table
            continue
        positions.append((float(match.group("east")), float(match.group("north"))))
        diams.append(float(match.group("diam")))

    n_12m = sum(1 for d in diams if d >= 10.0)
    n_7m = sum(1 for d in diams if d < 10.0)
    return EbGeometry(
        listobs_path=str(path),
        n_ant=len(positions),
        n_ant_12m=n_12m,
        n_ant_7m=n_7m,
        l80_m=baseline_percentile(positions, 80.0),
        lmax_m=baseline_percentile(positions, 100.0),
    )
