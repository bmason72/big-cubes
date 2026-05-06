from __future__ import annotations

import csv
import json
import math
import sys

from ryan_cy11_wsu_projection import (
    METHODS,
    current_eb_rate_gbps,
    parse_listobs_runtime,
    projected_channel_pol_per_spw,
    main,
)
from wsu_projection import (
    build_templates_from_rows,
    projected_spw_equivalents,
    stepped2_floor_velocity,
)


def test_stepped2_floor_velocity_matches_memo_bins():
    assert stepped2_floor_velocity(12.0) == 10.0
    assert stepped2_floor_velocity(3.2) == 2.0
    assert stepped2_floor_velocity(0.7) == 0.5
    assert stepped2_floor_velocity(0.2) == 0.1
    assert math.isclose(stepped2_floor_velocity(0.05), 0.05)


def test_projected_spw_equivalents_match_current_rules():
    assert projected_spw_equivalents("M1", 6) == 6
    assert math.isclose(projected_spw_equivalents("M4", 6), 5.5)
    assert math.isclose(projected_spw_equivalents("M5", 6), 16.0)


def test_distribution_and_exact_modes_reduce_channels_relative_to_memo():
    rows = [
        {
            "mous_uid": "m1",
            "eb_uid": "e1",
            "spw_id": "1",
            "center_freq_hz": str(230e9),
            "velocity_resolution_kms": "0.25",
            "corr_products": "XX YY",
        },
        {
            "mous_uid": "m1",
            "eb_uid": "e1",
            "spw_id": "2",
            "center_freq_hz": str(231e9),
            "velocity_resolution_kms": "0.30",
            "corr_products": "XX YY",
        },
        {
            "mous_uid": "m1",
            "eb_uid": "e1",
            "spw_id": "3",
            "center_freq_hz": str(232e9),
            "velocity_resolution_kms": "2.50",
            "corr_products": "XX YY",
        },
        {
            "mous_uid": "m1",
            "eb_uid": "e1",
            "spw_id": "4",
            "center_freq_hz": str(233e9),
            "velocity_resolution_kms": "12.0",
            "corr_products": "XX YY",
        },
    ]
    templates = build_templates_from_rows(rows)["m1"]
    corr_counts = [2, 2, 2, 2]

    memo = projected_channel_pol_per_spw("memo_uniform_binned", "M4", templates, corr_counts)
    uniform_exact = projected_channel_pol_per_spw("memo_uniform_exact", "M4", templates, corr_counts)
    distributed_binned = projected_channel_pol_per_spw("distributed_binned", "M4", templates, corr_counts)
    distributed_exact = projected_channel_pol_per_spw("distributed_exact", "M4", templates, corr_counts)

    assert memo > uniform_exact
    assert memo > distributed_binned
    assert distributed_binned > distributed_exact


def test_current_eb_rate_uses_actual_antennas_and_intervals(tmp_path):
    listobs = tmp_path / "eb.txt"
    listobs.write_text(
        "\n".join(
            [
                "Observation: ALMA",
                "Antennas: 10:",
                "  Date        Timerange (UTC)          Scan  FldId FieldName             nRows     SpwIds   Average Interval(s)    ScanIntent",
                "  01-Jan-2025/00:00:00.0 - 00:01:00.0     1      2 target                 100  [1,2]  [10.0, 10.0] [OBSERVE_TARGET#ON_SOURCE]",
                "           (nRows = Total number of rows per scan)",
            ]
        ),
        encoding="utf-8",
    )
    runtime = parse_listobs_runtime(listobs)
    rows = [
        {"spw_id": "1", "nchan": "128", "corr_products": "XX YY"},
        {"spw_id": "2", "nchan": "240", "corr_products": "XX YY"},
    ]
    rate = current_eb_rate_gbps(rows, runtime)
    nbase = 10 * 9 / 2
    cross_auto = 2.0 * 2.0 * 1.0 * nbase + 4.0 * 10
    expected = cross_auto * ((128 * 2) / 10.0 + (240 * 2) / 10.0) / 1.0e9
    assert math.isclose(rate, expected)


def test_main_writes_projection_outputs(tmp_path):
    tab = tmp_path / "tab"
    tab.mkdir()
    listobs = tmp_path / "eb.txt"
    listobs.write_text(
        "\n".join(
            [
                "Observation: ALMA",
                "Antennas: 47:",
                "  Date        Timerange (UTC)          Scan  FldId FieldName             nRows     SpwIds   Average Interval(s)    ScanIntent",
                "  01-Jan-2025/00:00:00.0 - 00:01:00.0     1      2 target                 100  [1,2,3,4]  [6.0, 6.0, 6.0, 6.0] [OBSERVE_TARGET#ON_SOURCE]",
                "           (nRows = Total number of rows per scan)",
            ]
        ),
        encoding="utf-8",
    )
    rows = [
        {
            "mous_uid": "m1",
            "eb_uid": "e1",
            "listobs_path": str(listobs),
            "spw_id": "1",
            "center_freq_hz": str(100e9),
            "nchan": "128",
            "velocity_resolution_kms": "0.25",
            "corr_products": "XX YY",
        },
        {
            "mous_uid": "m1",
            "eb_uid": "e1",
            "listobs_path": str(listobs),
            "spw_id": "2",
            "center_freq_hz": str(101e9),
            "nchan": "128",
            "velocity_resolution_kms": "0.30",
            "corr_products": "XX YY",
        },
        {
            "mous_uid": "m1",
            "eb_uid": "e1",
            "listobs_path": str(listobs),
            "spw_id": "3",
            "center_freq_hz": str(102e9),
            "nchan": "128",
            "velocity_resolution_kms": "2.50",
            "corr_products": "XX YY",
        },
        {
            "mous_uid": "m1",
            "eb_uid": "e1",
            "listobs_path": str(listobs),
            "spw_id": "4",
            "center_freq_hz": str(103e9),
            "nchan": "128",
            "velocity_resolution_kms": "12.0",
            "corr_products": "XX YY",
        },
    ]
    with (tab / "spw_level_table.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "mous_uid",
                "eb_uid",
                "listobs_path",
                "spw_id",
                "center_freq_hz",
                "nchan",
                "velocity_resolution_kms",
                "corr_products",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    out = tmp_path / "out"
    argv = sys.argv
    try:
        sys.argv = ["ryan_cy11_wsu_projection.py", str(tab), str(out)]
        rc = main()
    finally:
        sys.argv = argv

    assert rc == 0
    assert (out / "eb_projection_factors.csv").exists()
    assert (out / "mous_spw_templates.csv").exists()
    assert (out / "summary_overall.csv").exists()
    assert (out / "summary_by_array.csv").exists()
    assert (out / "summary.md").exists()

    summary_text = (out / "summary.md").read_text(encoding="utf-8")
    assert "memo_uniform_binned" in summary_text
    assert "distributed_exact" in summary_text

    metadata = json.loads((out / "run_metadata.json").read_text(encoding="utf-8"))
    assert metadata["methods"] == list(METHODS)
    assert metadata["mous_spw_sidecar"].endswith("mous_spw_templates.csv")
