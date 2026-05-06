from __future__ import annotations

import csv
import json
import math
import sys

from spw_setup_summary import (
    EbRuntimeMetadata,
    build_eb_assignments,
    build_mous_assignments,
    choose_gap_cut,
    compute_eb_data_rate_gbps,
    dedupe_eb_spw_rows,
    derive_cut_metadata,
    main,
    select_top_rows_by_coverage,
    summarize_setups,
)


def test_choose_gap_cut_moves_boundary_off_spike():
    values = [0.0585938] * 6 + [0.9375] * 4 + [1.875] * 10 + [2.0] * 8
    cut = choose_gap_cut(values, 2.0 / 3.0)

    assert math.isclose(cut["raw_quantile"], 1.875)
    assert math.isclose(cut["cut"], (1.875 + 2.0) / 2.0)
    assert cut["lower_value"] == 1.875
    assert cut["upper_value"] == 2.0


def test_dedupe_eb_spw_rows_uses_eb_spw_not_target_spw():
    rows = [
        {"mous_uid": "m1", "eb_uid": "e1", "spw_id": "17", "total_bw_hz": "1", "velocity_resolution_kms": "1"},
        {"mous_uid": "m1", "eb_uid": "e1", "spw_id": "17", "total_bw_hz": "1", "velocity_resolution_kms": "1"},
        {"mous_uid": "m1", "eb_uid": "e1", "spw_id": "19", "total_bw_hz": "1", "velocity_resolution_kms": "1"},
        {"mous_uid": "m1", "eb_uid": "e2", "spw_id": "17", "total_bw_hz": "1", "velocity_resolution_kms": "1"},
    ]

    deduped = dedupe_eb_spw_rows(rows)

    assert len(deduped) == 3


def test_derive_cut_metadata_defaults_to_fixed_bandwidth_cuts():
    rows = [
        {"total_bw_hz": str(58_593_800.0), "velocity_resolution_kms": "0.1"},
        {"total_bw_hz": str(117_187_500.0), "velocity_resolution_kms": "1.0"},
        {"total_bw_hz": str(2_000_000_000.0), "velocity_resolution_kms": "10.0"},
    ]

    metadata = derive_cut_metadata(
        rows,
        bandwidth_cut_mode="fixed",
        bandwidth_low_cut_mhz=90.0,
        bandwidth_high_cut_mhz=1200.0,
    )

    assert metadata["bandwidth"]["mode"] == "fixed"
    assert math.isclose(metadata["bandwidth"]["low_cut"], 0.09)
    assert math.isclose(metadata["bandwidth"]["high_cut"], 1.2)
    assert math.isnan(metadata["bandwidth"]["q33"])


def test_derive_cut_metadata_keeps_percentile_bandwidth_mode_available():
    rows = [
        {"total_bw_hz": str(58_593_800.0), "velocity_resolution_kms": "0.1"},
        {"total_bw_hz": str(937_500_000.0), "velocity_resolution_kms": "1.0"},
        {"total_bw_hz": str(1_875_000_000.0), "velocity_resolution_kms": "10.0"},
        {"total_bw_hz": str(2_000_000_000.0), "velocity_resolution_kms": "20.0"},
    ]

    metadata = derive_cut_metadata(
        rows,
        bandwidth_cut_mode="percentile",
        bandwidth_low_cut_mhz=90.0,
        bandwidth_high_cut_mhz=1200.0,
    )

    assert metadata["bandwidth"]["mode"] == "percentile"
    assert metadata["bandwidth"]["low_cut"] > 0.09
    assert metadata["bandwidth"]["high_cut"] > 1.2


def test_compute_eb_data_rate_uses_fiducial_antennas_and_corr_count():
    rows = [
        {"spw_id": "1", "nchan": "128", "corr_products": "XX  YY"},
        {"spw_id": "2", "nchan": "480", "corr_products": "XX  YY"},
    ]
    metadata = EbRuntimeMetadata(
        listobs_path="x",
        actual_antenna_count=47,
        fiducial_nant=43,
        array_type="12m",
        science_interval_median_s=6.0,
        science_interval_by_spw_s={1: 6.0, 2: 6.0},
    )

    rate = compute_eb_data_rate_gbps(rows, metadata)
    nbase = 43 * 42 / 2
    expected = nbase * ((128 * 8 / 6.0) + (480 * 8 / 6.0)) / 1.0e9

    assert math.isclose(rate, expected)


def test_setup_summary_uses_eb_coverage_and_mean_data_rate():
    mous_assignments = [
        {"mous_uid": "m1", "n_ebs": 2, "n_spws": 2, "setup_signature": "[WC,WC]", "setup_tokens": "WC WC", "line_fraction_percent": 20.0},
        {"mous_uid": "m2", "n_ebs": 1, "n_spws": 1, "setup_signature": "[NF]", "setup_tokens": "NF", "line_fraction_percent": 10.0},
    ]
    eb_assignments = [
        {"mous_uid": "m1", "eb_uid": "e1", "setup_signature": "[WC,WC]", "actual_n_spws": 2, "data_rate_gbps": 0.30},
        {"mous_uid": "m1", "eb_uid": "e2", "setup_signature": "[WC,WC]", "actual_n_spws": 2, "data_rate_gbps": 0.50},
        {"mous_uid": "m2", "eb_uid": "e3", "setup_signature": "[NF]", "actual_n_spws": 1, "data_rate_gbps": 0.05},
    ]

    summaries = summarize_setups(mous_assignments, eb_assignments)
    top = summaries[0]

    assert top["setup_signature"] == "[WC,WC]"
    assert top["eb_count"] == 2
    assert top["eb_spw_count"] == 4
    assert math.isclose(top["eb_spw_percent"], 80.0)
    assert math.isclose(top["mean_eb_data_rate_gbps"], 0.40)

    selected = select_top_rows_by_coverage(summaries, 0.66)
    assert len(selected) == 1
    assert selected[0]["setup_signature"] == "[WC,WC]"


def test_main_writes_summary_plots_and_eb_assignments(tmp_path):
    tabulation_dir = tmp_path / "tab"
    tabulation_dir.mkdir()

    listobs_dir = tmp_path / "listobs"
    listobs_dir.mkdir()
    listobs1 = listobs_dir / "eb1.txt"
    listobs1.write_text(
        "\n".join(
            [
                "Observation: ALMA",
                "Antennas: 47:",
                "  Date        Timerange (UTC)          Scan  FldId FieldName             nRows     SpwIds   Average Interval(s)    ScanIntent",
                "  01-Jan-2025/00:00:00.0 - 00:01:00.0     1      2 target                 100  [1,2]  [6.0, 6.0] [OBSERVE_TARGET#ON_SOURCE]",
                "           (nRows = Total number of rows per scan)",
            ]
        ),
        encoding="utf-8",
    )
    listobs2 = listobs_dir / "eb2.txt"
    listobs2.write_text(
        "\n".join(
            [
                "Observation: ALMA",
                "Antennas: 9:",
                "  Date        Timerange (UTC)          Scan  FldId FieldName             nRows     SpwIds   Average Interval(s)    ScanIntent",
                "  01-Jan-2025/00:00:00.0 - 00:01:00.0     1      2 target                 100  [3]  [10.0] [OBSERVE_TARGET#ON_SOURCE]",
                "           (nRows = Total number of rows per scan)",
            ]
        ),
        encoding="utf-8",
    )

    spw_level_rows = [
        {
            "mous_uid": "m1",
            "eb_uid": "e1",
            "listobs_path": str(listobs1),
            "spw_id": "1",
            "center_freq_hz": "100e9",
            "nchan": "128",
            "total_bw_hz": str(2.0e9),
            "velocity_resolution_kms": "20.0",
            "corr_products": "XX  YY",
        },
        {
            "mous_uid": "m1",
            "eb_uid": "e1",
            "listobs_path": str(listobs1),
            "spw_id": "2",
            "center_freq_hz": "110e9",
            "nchan": "128",
            "total_bw_hz": str(2.0e9),
            "velocity_resolution_kms": "20.0",
            "corr_products": "XX  YY",
        },
        {
            "mous_uid": "m1",
            "eb_uid": "e1",
            "listobs_path": str(listobs1),
            "spw_id": "2",
            "center_freq_hz": "110e9",
            "nchan": "128",
            "total_bw_hz": str(2.0e9),
            "velocity_resolution_kms": "20.0",
            "corr_products": "XX  YY",
        },
        {
            "mous_uid": "m2",
            "eb_uid": "e2",
            "listobs_path": str(listobs2),
            "spw_id": "3",
            "center_freq_hz": "90e9",
            "nchan": "240",
            "total_bw_hz": str(0.1e9),
            "velocity_resolution_kms": "0.1",
            "corr_products": "XX  YY",
        },
    ]
    with (tabulation_dir / "spw_level_table.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "mous_uid",
                "eb_uid",
                "listobs_path",
                "spw_id",
                "center_freq_hz",
                "nchan",
                "total_bw_hz",
                "velocity_resolution_kms",
                "corr_products",
            ],
        )
        writer.writeheader()
        writer.writerows(spw_level_rows)

    mous_spw_rows = [
        {
            "mous_uid": "m1",
            "spw_id": "1",
            "center_freq_hz": "100e9",
            "total_bw_hz": str(2.0e9),
            "velocity_resolution_kms": "20.0",
            "mean_summed_line_fraction_of_spw": "0.20",
            "n_ebs_parsed": "1",
        },
        {
            "mous_uid": "m1",
            "spw_id": "2",
            "center_freq_hz": "110e9",
            "total_bw_hz": str(2.0e9),
            "velocity_resolution_kms": "20.0",
            "mean_summed_line_fraction_of_spw": "0.10",
            "n_ebs_parsed": "1",
        },
        {
            "mous_uid": "m2",
            "spw_id": "3",
            "center_freq_hz": "90e9",
            "total_bw_hz": str(0.1e9),
            "velocity_resolution_kms": "0.1",
            "mean_summed_line_fraction_of_spw": "0.05",
            "n_ebs_parsed": "1",
        },
    ]
    with (tabulation_dir / "mous_spw_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "mous_uid",
                "spw_id",
                "center_freq_hz",
                "total_bw_hz",
                "velocity_resolution_kms",
                "mean_summed_line_fraction_of_spw",
                "n_ebs_parsed",
            ],
        )
        writer.writeheader()
        writer.writerows(mous_spw_rows)

    output_dir = tmp_path / "out"
    argv = sys.argv
    try:
        sys.argv = ["spw_setup_summary.py", str(tabulation_dir), str(output_dir), "--coverage-fraction", "0.5"]
        rc = main()
    finally:
        sys.argv = argv

    assert rc == 0
    assert (output_dir / "setup_summary.csv").exists()
    assert (output_dir / "setup_summary.md").exists()
    assert (output_dir / "eb_setup_assignments.csv").exists()
    assert (output_dir / "bandwidth_hist_eb_spw.png").exists()
    assert (output_dir / "resolution_hist_eb_spw.png").exists()

    summary_text = (output_dir / "setup_summary.md").read_text(encoding="utf-8")
    assert "Bandwidth cuts (GHz, fixed): low=0.090000, high=1.200000" in summary_text
    assert "Top 1 Setup Signatures" in summary_text
    assert "covering 50% of EBs target" in summary_text
    assert "EB-SPWs (66.67%)" in summary_text
    assert "mean EB data rate" in summary_text

    cut_metadata = json.loads((output_dir / "cut_metadata.json").read_text(encoding="utf-8"))
    assert cut_metadata["n_unique_eb_spw"] == 3
    assert cut_metadata["selected_setup_count"] == 1
    assert cut_metadata["bandwidth"]["mode"] == "fixed"
