from __future__ import annotations

import csv
import json
import math
import sys

from ryan_cy11_wsu_projection import (
    METHODS,
    current_eb_rate_gbps,
    parse_listobs_runtime,
    projected_channel_pol_agg_for_templates,
    projected_channel_pol_per_spw,
    main,
)
from wsu_projection import (
    MousSpwTemplate,
    build_templates_from_rows,
    projected_spw_equivalents,
    stepped2_floor_velocity,
)


def write_minimal_tuning_inputs(tab, rows, total_scan_time_s=60.0, tuning_total_bw_hz=None):
    tuning_total_bw_hz = tuning_total_bw_hz or float(len(rows) * 2.0e9)
    current_nchan = sum(int(row["nchan"]) for row in rows)

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

    with (tab / "eb_tuning_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "mous_uid",
                "eb_uid",
                "listobs_path",
                "tuning_group_id",
                "tuning_group_index",
                "tuning_signature",
                "target_names",
                "n_spw_in_tuning",
                "tuning_total_bw_hz",
                "tuning_center_freq_hz",
                "n_scans_in_tuning",
                "total_scan_time_s",
                "n_tuning_groups_in_eb",
                "is_spectral_scan",
                "eb_total_science_bw_hz",
                "eb_current_nchan_agg",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "mous_uid": rows[0]["mous_uid"],
                "eb_uid": rows[0]["eb_uid"],
                "listobs_path": rows[0]["listobs_path"],
                "tuning_group_id": "TG01",
                "tuning_group_index": 1,
                "tuning_signature": ",".join(row["spw_id"] for row in rows),
                "target_names": "target",
                "n_spw_in_tuning": len(rows),
                "tuning_total_bw_hz": tuning_total_bw_hz,
                "tuning_center_freq_hz": sum(float(row["center_freq_hz"]) for row in rows) / len(rows),
                "n_scans_in_tuning": 1,
                "total_scan_time_s": total_scan_time_s,
                "n_tuning_groups_in_eb": 1,
                "is_spectral_scan": False,
                "eb_total_science_bw_hz": tuning_total_bw_hz,
                "eb_current_nchan_agg": current_nchan,
            }
        )

    with (tab / "eb_tuning_spw_table.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "mous_uid",
                "eb_uid",
                "listobs_path",
                "tuning_group_id",
                "tuning_group_index",
                "tuning_signature",
                "target_names",
                "n_spw_in_tuning",
                "tuning_total_bw_hz",
                "spw_id",
                "spw_name",
                "nchan",
                "center_freq_hz",
                "freq_min_hz",
                "freq_max_hz",
                "total_bw_hz",
                "abs_chan_width_hz",
                "velocity_resolution_kms",
                "corr_products",
                "median_interval_s",
                "n_scans_in_tuning",
                "total_scan_time_s",
                "n_tuning_groups_in_eb",
                "is_spectral_scan",
                "eb_total_science_bw_hz",
                "eb_current_nchan_agg",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "mous_uid": row["mous_uid"],
                    "eb_uid": row["eb_uid"],
                    "listobs_path": row["listobs_path"],
                    "tuning_group_id": "TG01",
                    "tuning_group_index": 1,
                    "tuning_signature": ",".join(r["spw_id"] for r in rows),
                    "target_names": "target",
                    "n_spw_in_tuning": len(rows),
                    "tuning_total_bw_hz": tuning_total_bw_hz,
                    "spw_id": row["spw_id"],
                    "spw_name": f"spw{row['spw_id']}",
                    "nchan": row["nchan"],
                    "center_freq_hz": row["center_freq_hz"],
                    "freq_min_hz": row["center_freq_hz"],
                    "freq_max_hz": row["center_freq_hz"],
                    "total_bw_hz": 2.0e9,
                    "abs_chan_width_hz": 1.0e5,
                    "velocity_resolution_kms": row["velocity_resolution_kms"],
                    "corr_products": row["corr_products"],
                    "median_interval_s": 6.0,
                    "n_scans_in_tuning": 1,
                    "total_scan_time_s": total_scan_time_s,
                    "n_tuning_groups_in_eb": 1,
                    "is_spectral_scan": False,
                    "eb_total_science_bw_hz": tuning_total_bw_hz,
                    "eb_current_nchan_agg": current_nchan,
                }
            )

    with (tab / "mous_tuning_spw_templates.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "mous_uid",
                "representative_eb_uid",
                "tuning_group_id",
                "tuning_group_index",
                "tuning_total_bw_hz",
                "spw_id",
                "center_freq_hz",
                "velocity_resolution_kms",
                "n_tuning_groups_in_mous",
                "is_spectral_scan_mous",
                "mous_total_science_bw_hz",
                "mous_current_nchan_agg",
                "eb_tuning_consistency_flag",
                "eb_tuning_consistency_notes",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "mous_uid": row["mous_uid"],
                    "representative_eb_uid": row["eb_uid"],
                    "tuning_group_id": "TG01",
                    "tuning_group_index": 1,
                    "tuning_total_bw_hz": tuning_total_bw_hz,
                    "spw_id": row["spw_id"],
                    "center_freq_hz": row["center_freq_hz"],
                    "velocity_resolution_kms": row["velocity_resolution_kms"],
                    "n_tuning_groups_in_mous": 1,
                    "is_spectral_scan_mous": False,
                    "mous_total_science_bw_hz": tuning_total_bw_hz,
                    "mous_current_nchan_agg": current_nchan,
                    "eb_tuning_consistency_flag": "consistent",
                    "eb_tuning_consistency_notes": "",
                }
            )


def test_stepped2_floor_velocity_matches_memo_bins():
    assert stepped2_floor_velocity(12.0) == 10.0
    assert stepped2_floor_velocity(3.2) == 2.0
    assert stepped2_floor_velocity(0.7) == 0.5
    assert stepped2_floor_velocity(0.2) == 0.1
    assert math.isclose(stepped2_floor_velocity(0.05), 0.05)


def test_projected_spw_equivalents_use_ceil_at_all_milestones():
    assert projected_spw_equivalents("M1", 6) == 6
    assert projected_spw_equivalents("M4", 6) == 6  # ceil(11/2), not 5.5
    assert math.isclose(projected_spw_equivalents("M5", 6), 16.0)
    assert projected_spw_equivalents("M4", 3) == 8
    assert projected_spw_equivalents("M5", 3) == 16


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


def test_projected_channel_pol_agg_cap_reduces_high_resolution_projection():
    templates = [
        MousSpwTemplate("m1", 1, 100e9, 0.25, 3),
        MousSpwTemplate("m1", 2, 101e9, 0.30, 3),
    ]
    corr_counts = [2, 2]
    uncapped = projected_channel_pol_agg_for_templates(
        "distributed_exact", "M4", templates, corr_counts, projected_nspw=1.0
    )
    capped = projected_channel_pol_agg_for_templates(
        "distributed_exact", "M4", templates, corr_counts, projected_nspw=1.0, velocity_cap_kms=1.0
    )
    assert capped < uncapped


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
    write_minimal_tuning_inputs(tab, rows)

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
    assert (out / "mous_tuning_spw_templates.csv").exists()
    assert (out / "summary_overall.csv").exists()
    assert (out / "summary_by_array.csv").exists()
    assert (out / "summary_overall_nchan.csv").exists()
    assert (out / "summary_overall_nchan.tex").exists()
    assert (out / "summary.md").exists()

    summary_text = (out / "summary.md").read_text(encoding="utf-8")
    assert "memo_uniform_binned" in summary_text
    assert "distributed_exact" in summary_text

    metadata = json.loads((out / "run_metadata.json").read_text(encoding="utf-8"))
    assert metadata["methods"] == list(METHODS)
    assert metadata["mous_spw_sidecar"].endswith("mous_tuning_spw_templates.csv")


def test_intent_aware_projection_writes_cap_outputs_and_reduces_selected_volume(tmp_path):
    tab = tmp_path / "tab"
    tab.mkdir()
    target_listobs = tmp_path / "sessionsession_1__uid___A002_X1_X2_targets.ms__listobs.txt"
    full_listobs = tmp_path / "sessionsession_1__uid___A002_X1_X2.ms__listobs.txt"
    target_listobs.write_text(
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
    full_listobs.write_text(
        "\n".join(
            [
                "Observation: ALMA",
                "Antennas: 47:",
                "  Date        Timerange (UTC)          Scan  FldId FieldName             nRows     SpwIds   Average Interval(s)    ScanIntent",
                "  01-Jan-2025/00:00:00.0 - 00:01:00.0     1      2 target                 100  [1,2]  [6.0, 6.0] [OBSERVE_TARGET#ON_SOURCE]",
                "  01-Jan-2025/00:01:00.0 - 00:01:30.0     2      3 bpcal                  100  [1,2]  [6.0, 6.0] [CALIBRATE_BANDPASS#ON_SOURCE]",
                "  01-Jan-2025/00:01:30.0 - 00:02:00.0     3      4 phcal                  100  [1,2]  [6.0, 6.0] [CALIBRATE_PHASE#ON_SOURCE]",
                "  01-Jan-2025/00:02:00.0 - 00:02:15.0     4      5 check                  100  [1,2]  [6.0, 6.0] [OBSERVE_CHECK_SOURCE#ON_SOURCE]",
                "           (nRows = Total number of rows per scan)",
                "Spectral Windows:",
                "SpwID  Name                                  #Chans Frame Ch0(MHz) ChanWid(kHz) TotBW(kHz) CtrFreq(MHz) BBC Num  Corrs",
                "0      none                                     1 TOPO 0.0 1.0 1.0 0.0 0 XX",
                "1      BB_1#FULL_RES                          128 TOPO 100000.0 83.4 10675.2 100500.0 1 XX YY",
                "2      BB_2#FULL_RES                          128 TOPO 101000.0 83.4 10675.2 101500.0 2 XX YY",
                "Sources:",
            ]
        ),
        encoding="utf-8",
    )
    rows = [
        {
            "mous_uid": "m1",
            "eb_uid": "e1",
            "listobs_path": str(target_listobs),
            "spw_id": "1",
            "center_freq_hz": str(100.5e9),
            "nchan": "128",
            "velocity_resolution_kms": "0.25",
            "corr_products": "XX YY",
        },
        {
            "mous_uid": "m1",
            "eb_uid": "e1",
            "listobs_path": str(target_listobs),
            "spw_id": "2",
            "center_freq_hz": str(101.5e9),
            "nchan": "128",
            "velocity_resolution_kms": "0.25",
            "corr_products": "XX YY",
        },
    ]
    write_minimal_tuning_inputs(tab, rows)

    out = tmp_path / "out"
    rc = main(
        [
            str(tab),
            str(out),
            "--intent-aware",
            "--apply-calibrator-cap",
        ]
    )
    assert rc == 0
    assert (out / "eb_intent_projection_factors.csv").exists()
    assert (out / "summary_overall_intent_aware.csv").exists()
    assert (out / "intent_aware_scan_details.csv").exists()

    with (out / "summary_overall_intent_aware.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    row = next(item for item in rows if item["milestone"] == "M4" and item["method"] == "distributed_exact")
    assert float(row["projected_selected_total_volume_gb"]) < float(row["projected_selected_total_uncapped_volume_gb"])
    assert float(row["cap_delta_total_volume_gb"]) < 0.0

    summary_text = (out / "summary.md").read_text(encoding="utf-8")
    assert "Intent-Aware Projection" in summary_text


def test_intent_aware_without_cap_skips_out_of_band_calibrator_spws(tmp_path):
    tab = tmp_path / "tab"
    tab.mkdir()
    target_listobs = tmp_path / "sessionsession_1__uid___A002_X1_X2_targets.ms__listobs.txt"
    full_listobs = tmp_path / "sessionsession_1__uid___A002_X1_X2.ms__listobs.txt"
    target_listobs.write_text(
        "\n".join(
            [
                "Observation: ALMA",
                "Antennas: 47:",
                "  Date        Timerange (UTC)          Scan  FldId FieldName             nRows     SpwIds   Average Interval(s)    ScanIntent",
                "  01-Jan-2025/00:00:00.0 - 00:01:00.0     1      2 target                 100  [1]  [6.0] [OBSERVE_TARGET#ON_SOURCE]",
                "           (nRows = Total number of rows per scan)",
            ]
        ),
        encoding="utf-8",
    )
    full_listobs.write_text(
        "\n".join(
            [
                "Observation: ALMA",
                "Antennas: 47:",
                "  Date        Timerange (UTC)          Scan  FldId FieldName             nRows     SpwIds   Average Interval(s)    ScanIntent",
                "  01-Jan-2025/00:00:00.0 - 00:01:00.0     1      2 target                 100  [1]  [6.0] [OBSERVE_TARGET#ON_SOURCE]",
                "  01-Jan-2025/00:01:00.0 - 00:01:30.0     2      3 bpcal                  100  [2]  [6.0] [CALIBRATE_BANDPASS#ON_SOURCE]",
                "           (nRows = Total number of rows per scan)",
                "Spectral Windows:",
                "SpwID  Name                                  #Chans Frame Ch0(MHz) ChanWid(kHz) TotBW(kHz) CtrFreq(MHz) BBC Num  Corrs",
                "1      BB_1#FULL_RES                          128 TOPO 100000.0 83.4 10675.2 100500.0 1 XX YY",
                "2      BAD#FULL_RES                           128 TOPO 118300.0 83.4 10675.2 118345.0 2 XX YY",
                "Sources:",
            ]
        ),
        encoding="utf-8",
    )
    rows = [
        {
            "mous_uid": "m1",
            "eb_uid": "e1",
            "listobs_path": str(target_listobs),
            "spw_id": "1",
            "center_freq_hz": str(100.5e9),
            "nchan": "128",
            "velocity_resolution_kms": "0.25",
            "corr_products": "XX YY",
        },
    ]
    write_minimal_tuning_inputs(tab, rows)

    out = tmp_path / "out"
    rc = main(
        [
            str(tab),
            str(out),
            "--intent-aware",
        ]
    )
    assert rc == 0
    assert (out / "summary_overall_intent_aware.csv").exists()


def test_write_csv_ignores_extra_row_fields(tmp_path):
    from ryan_cy11_wsu_projection import write_csv

    path = tmp_path / "out.csv"
    write_csv(path, [{"a": "1", "b": "2", "extra": "x"}], ["a", "b"])
    rows = list(csv.DictReader(path.open()))
    assert rows == [{"a": "1", "b": "2"}]


def test_write_latex_factor_table_escapes_underscores_once(tmp_path):
    from ryan_cy11_wsu_projection import MILESTONES, write_latex_factor_table

    rows = [
        {
            "milestone": milestone,
            "method": method,
            "sample_total_factor": 2.0,
            "mean_eb_factor": 1.5,
            "median_eb_factor": 1.0,
        }
        for milestone in MILESTONES
        for method in METHODS
    ]
    path = tmp_path / "table.tex"
    write_latex_factor_table(path, rows, "caption", "tab:label")
    text = path.read_text(encoding="utf-8")
    assert r"memo\_uniform\_binned" in text
    assert r"\\_" not in text


def test_memo_uniform_band_from_nearest_template_when_mean_in_interband_gap():
    from ryan_cy11_wsu_projection import projected_channel_pol_agg_for_templates
    from wsu_projection import project_nchan_agg_for_templates

    # band 3 (112 GHz) + band 4 (125 GHz) SPWs: mean = 118.5 GHz falls in the
    # 116-125 GHz inter-band gap and must not raise
    templates = [
        MousSpwTemplate("m", 1, 112.0e9, 5.0, 3),
        MousSpwTemplate("m", 2, 125.0e9, 5.0, 4),
    ]
    value = projected_channel_pol_agg_for_templates(
        "memo_uniform_binned", "M1", templates, [2, 2], projected_nspw=1.0
    )
    assert math.isfinite(value) and value > 0
    value2 = project_nchan_agg_for_templates(
        "memo_uniform_binned", "M1", templates, projected_nspw=1.0
    )
    assert math.isfinite(value2) and value2 > 0
