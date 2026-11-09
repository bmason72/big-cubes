#!/usr/bin/env python3
"""Tabulate Cycle-11 science integration time by array, band, and configuration."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any

from data.cycle13_configurations import CONFIG_LABELS_12M


ARRAY_CONFIGS = {
    "12m": CONFIG_LABELS_12M,
    "7m": ["7-m"],
}
BANDS = [str(band) for band in range(1, 11)]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def calculate_stats(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for row in rows:
        key = (row["array"], row["band"], row["estimated_config"])
        grouped[key].append(float(row["integration_time_s"]))

    result: list[dict[str, Any]] = []
    for array in ("12m", "7m"):
        for band in BANDS:
            for config in ARRAY_CONFIGS[array]:
                values = grouped.get((array, band, config), [])
                if not values:
                    continue
                result.append(
                    {
                        "array": array,
                        "band": band,
                        "estimated_config": config,
                        "n_mous": len(values),
                        "mean_integration_time_s": statistics.fmean(values),
                        "stddev_integration_time_s": statistics.pstdev(values),
                    }
                )
    return result


def write_long_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames = [
        "array",
        "band",
        "estimated_config",
        "n_mous",
        "mean_integration_time_s",
        "stddev_integration_time_s",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(array: str, stats_rows: list[dict[str, Any]]) -> list[str]:
    configs = ARRAY_CONFIGS[array]
    lookup = {
        (str(row["band"]), str(row["estimated_config"])): row
        for row in stats_rows
        if row["array"] == array
    }
    lines = [
        f"## {array}",
        "",
        "| Band | " + " | ".join(configs) + " |",
        "|---:|" + "|".join("---:" for _ in configs) + "|",
    ]
    for band in BANDS:
        cells: list[str] = []
        for config in configs:
            row = lookup.get((band, config))
            if row is None:
                cells.append("--")
            else:
                cells.append(
                    f"{float(row['mean_integration_time_s']):.3f} ± "
                    f"{float(row['stddev_integration_time_s']):.3f} "
                    f"(n={int(row['n_mous'])})"
                )
        lines.append(f"| {band} | " + " | ".join(cells) + " |")
    return lines


def write_markdown(path: Path, stats_rows: list[dict[str, Any]]) -> None:
    lines = [
        "# Cycle-11 Science Integration Time by Band and Configuration",
        "",
        "Cells give mean ± population standard deviation in seconds, followed by the number of MOUS-band rows.",
        "Configuration is estimated from the median L80 across each MOUS's EBs.",
        "",
    ]
    lines.extend(markdown_table("12m", stats_rows))
    lines.append("")
    lines.extend(markdown_table("7m", stats_rows))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_csv",
        type=Path,
        nargs="?",
        default=Path("ryanCy11out/integration_times/cy11_mous_science_setups.csv"),
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        nargs="?",
        default=Path("ryanCy11out/integration_times"),
    )
    args = parser.parse_args(argv)

    stats_rows = calculate_stats(read_rows(args.input_csv))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_long_csv(args.output_dir / "integration_time_by_band_config.csv", stats_rows)
    write_markdown(args.output_dir / "integration_time_by_band_config.md", stats_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
