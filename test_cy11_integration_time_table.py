from __future__ import annotations

import math

from cy11_integration_time_table import calculate_stats, markdown_table


def test_calculate_stats_groups_mous_band_rows():
    rows = [
        {"array": "12m", "band": "6", "estimated_config": "C43-5", "integration_time_s": "2.0"},
        {"array": "12m", "band": "6", "estimated_config": "C43-5", "integration_time_s": "4.0"},
        {"array": "7m", "band": "3", "estimated_config": "7-m", "integration_time_s": "6.0"},
    ]

    stats = calculate_stats(rows)

    twelve = next(row for row in stats if row["array"] == "12m")
    seven = next(row for row in stats if row["array"] == "7m")
    assert twelve["n_mous"] == 2
    assert twelve["mean_integration_time_s"] == 3.0
    assert twelve["stddev_integration_time_s"] == 1.0
    assert seven["n_mous"] == 1
    assert seven["stddev_integration_time_s"] == 0.0


def test_markdown_table_has_band_rows_and_empty_cells():
    stats = [
        {
            "array": "7m",
            "band": "3",
            "estimated_config": "7-m",
            "n_mous": 2,
            "mean_integration_time_s": 6.05,
            "stddev_integration_time_s": 0.0,
        }
    ]

    text = "\n".join(markdown_table("7m", stats))

    assert "| 3 | 6.050 ± 0.000 (n=2) |" in text
    assert "| 2 | -- |" in text
