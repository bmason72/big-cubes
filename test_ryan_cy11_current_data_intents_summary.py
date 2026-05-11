from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from ryan_cy11_current_data_intents_summary import classify_bucket, main, matched_full_listobs_path


def test_classify_bucket_uses_anchor_intent_rules():
    assert classify_bucket(["OBSERVE_TARGET#ON_SOURCE", "CALIBRATE_WVR#ON_SOURCE"]) == "science"
    assert classify_bucket(["CALIBRATE_BANDPASS#ON_SOURCE", "CALIBRATE_FLUX#ON_SOURCE", "CALIBRATE_WVR#ON_SOURCE"]) == "bandpass"
    assert classify_bucket(["CALIBRATE_PHASE#ON_SOURCE", "CALIBRATE_WVR#ON_SOURCE"]) == "phase"
    assert classify_bucket(["OBSERVE_CHECK_SOURCE#ON_SOURCE", "CALIBRATE_WVR#ON_SOURCE"]) == "check_source"
    assert classify_bucket(["CALIBRATE_DIFFGAIN#REFERENCE", "CALIBRATE_WVR#REFERENCE"]) == "diffgain_reference"
    assert classify_bucket(["CALIBRATE_ATMOSPHERE#OFF_SOURCE", "CALIBRATE_WVR#OFF_SOURCE"]) == "atmosphere"


def test_matched_full_listobs_path_replaces_target_suffixes(tmp_path):
    listobs_dir = tmp_path / "listobs"
    listobs_dir.mkdir()
    full = listobs_dir / "sessionsession_1__uid___A002_X1_X2.ms__listobs.txt"
    full.write_text("", encoding="utf-8")
    target = listobs_dir / "sessionsession_1__uid___A002_X1_X2_targets.ms__listobs.txt"
    target.write_text("", encoding="utf-8")

    resolved = matched_full_listobs_path(str(target), [tmp_path])
    assert resolved == full


def test_main_writes_current_summary_outputs(tmp_path):
    tab = tmp_path / "tab"
    tab.mkdir()
    listobs_dir = tmp_path / "sample" / "mous1" / "listobs"
    listobs_dir.mkdir(parents=True)

    target_listobs = listobs_dir / "sessionsession_1__uid___A002_X123_X456_targets.ms__listobs.txt"
    target_listobs.write_text("placeholder\n", encoding="utf-8")

    full_listobs = listobs_dir / "sessionsession_1__uid___A002_X123_X456.ms__listobs.txt"
    full_listobs.write_text(
        "\n".join(
            [
                "Observation: ALMA",
                "Antennas: 47:",
                "  Date        Timerange (UTC)          Scan  FldId FieldName             nRows     SpwIds   Average Interval(s)    ScanIntent",
                "  01-Jan-2025/00:00:00.0 - 00:01:00.0     1      0 bpcal                  100  [1,2,3]  [6.0, 6.0, 6.0] [CALIBRATE_BANDPASS#ON_SOURCE,CALIBRATE_FLUX#ON_SOURCE,CALIBRATE_WVR#ON_SOURCE]",
                "              00:01:05.0 - 00:02:05.0     2      1 target                 100  [3]  [6.0] [OBSERVE_TARGET#ON_SOURCE]",
                "              00:02:10.0 - 00:03:10.0     3      2 phasecal               100  [3]  [6.0] [CALIBRATE_PHASE#ON_SOURCE,CALIBRATE_WVR#ON_SOURCE]",
                "           (nRows = Total number of rows per scan)",
                "Spectral Windows:  (3 unique spectral windows and 1 unique polarization setups)",
                "  SpwID  Name                                      #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) CtrFreq(MHz) BBC Num  Corrs",
                "  1      X#BB_1#SQLD                                  1   TOPO  100000.000   2000000.000   2000000.0 100000.0000        1  XX  YY",
                "  2      X#BB_1#SW-01#CH_AVG                          1   TOPO  100000.000   1875000.000   1875000.0 100000.0000        1  XX  YY",
                "  3      X#BB_1#SW-01#FULL_RES                      240   TOPO  100000.000      3906.250    937500.0 100468.7500        1  XX  YY",
                "Sources: 1",
            ]
        ),
        encoding="utf-8",
    )

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
        writer.writerow(
            {
                "mous_uid": "m1",
                "eb_uid": "e1",
                "listobs_path": str(target_listobs),
                "spw_id": "3",
                "center_freq_hz": str(100468.75e6),
                "nchan": "240",
                "velocity_resolution_kms": "0.011658",
                "corr_products": "XX YY",
            }
        )

    out = tmp_path / "out"
    rc = main([str(tab), str(out), "--skip-plots"])

    assert rc == 0
    assert (out / "scan_level_table.csv").exists()
    assert (out / "bucket_by_eb.csv").exists()
    assert (out / "bucket_by_mous.csv").exists()
    assert (out / "bucket_fraction_summary.csv").exists()
    assert (out / "summary.md").exists()
    assert (out / "sanity_checks.csv").exists()

    with (out / "bucket_by_eb.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    bandpass_row = next(row for row in rows if row["bucket"] == "bandpass")
    science_row = next(row for row in rows if row["bucket"] == "science")
    phase_row = next(row for row in rows if row["bucket"] == "phase")
    assert float(bandpass_row["volume_gb"]) > 0.0
    assert float(science_row["volume_gb"]) > 0.0
    assert float(phase_row["volume_gb"]) > 0.0

    with (out / "scan_level_table.csv").open(newline="", encoding="utf-8") as handle:
        scan_rows = list(csv.DictReader(handle))
    bandpass_scan = next(row for row in scan_rows if row["bucket"] == "bandpass")
    assert bandpass_scan["excluded_sqld_count"] == "1"
    assert bandpass_scan["excluded_ch_avg_count"] == "1"
    assert bandpass_scan["retained_spw_count"] == "1"

    with (out / "sanity_checks.csv").open(newline="", encoding="utf-8") as handle:
        sanity_row = next(csv.DictReader(handle))
    assert math.isclose(float(sanity_row["time_fraction_sum"]), 1.0)
    assert math.isclose(float(sanity_row["volume_fraction_sum"]), 1.0)

    metadata = json.loads((out / "run_metadata.json").read_text(encoding="utf-8"))
    assert metadata["n_sample_eb"] == 1
    assert metadata["plots"]["plots_written"] is False
