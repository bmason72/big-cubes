from __future__ import annotations

import math

from spw_tabulate import (
    SpwInfo,
    buffer_intervals,
    complement_intervals,
    extract_spw_level_rows,
    is_target_listobs_path,
    interval_total_width,
    merge_intervals,
    parse_contdat,
    parse_listobs,
    process_mous_dir,
    summarize_contdat_entry,
)


def test_parse_listobs_extracts_science_spws(tmp_path):
    path = tmp_path / "sessionsession_1__uid___A002_X123_X456_targets.ms__listobs.txt"
    path.write_text(
        "\n".join(
            [
                "  Date        Timerange (UTC)          Scan  FldId FieldName             nRows     SpwIds   Average Interval(s)    ScanIntent",
                "  04-Jan-2025/14:07:55.2 - 14:17:31.9     6      2 gp1547+2203             428640  [5,9,16,21]  [6.05, 6.05, 6.05, 6.05] [OBSERVE_TARGET#ON_SOURCE]",
                "           (nRows = Total number of rows per scan)",
                "Fields: 1",
                "  ID   Code Name                RA               Decl           Epoch   SrcId      nRows",
                "  2    none gp1547+2203         15:47:48.990000 +22.03.03.20000 ICRS    2        1737120",
                "Spectral Windows:  (4 unique spectral windows and 1 unique polarization setups)",
                "  SpwID  Name                                       #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) CtrFreq(MHz) BBC Num  Corrs",
                "  5      X#BB_1#SW-01#FULL_RES    128   TOPO   99038.315    -15625.000   2000000.0  98046.1280        1  XX  YY",
                "  9      X#BB_3#SW-01#FULL_RES    128   TOPO  109053.940     15625.000   2000000.0 110046.1280        3  XX  YY",
                "  16     X#BB_2#SW-01#FULL_RES    128   TOPO  100913.315    -15625.000   2000000.0  99921.1280        2  XX  YY",
                "  21     X#BB_4#SW-01#FULL_RES    240   TOPO  111302.132      3906.250    937500.0 111768.9290        4  XX  YY",
                "Sources: 1",
            ]
        ),
        encoding="utf-8",
    )
    parsed = parse_listobs(path, "uid___A001_X1_X2")
    rows = extract_spw_level_rows(parsed)

    assert len(rows) == 4
    assert {row["spw_id"] for row in rows} == {5, 9, 16, 21}
    assert all(row["science_target_flag"] for row in rows)


def test_parse_contdat_extracts_multiple_ranges(tmp_path):
    path = tmp_path / "cont.dat"
    path.write_text(
        "\n".join(
            [
                "Field: gp1547+2203",
                "",
                "SpectralWindow: 5 X#BB_1#SW-01#FULL_RES",
                "97.1750867149~97.9562523699GHz LSRK",
                "98.1906020664~98.8780278428GHz LSRK",
                "",
                "SpectralWindow: 9 X#BB_3#SW-01#FULL_RES",
                "Flags: ALLCONT",
                "109.1736455316~110.8765892210GHz LSRK",
            ]
        ),
        encoding="utf-8",
    )
    entries, issues = parse_contdat(path, "uid___A001_X1_X2")

    assert not issues
    assert len(entries) == 2
    assert entries[0].target_name == "gp1547+2203"
    assert entries[0].spw_id == 5
    assert len(entries[0].continuum_ranges_hz) == 2
    assert entries[1].flags == ["ALLCONT"]


def test_line_metrics_use_continuum_complement_and_buffer():
    spw_row = {
        "freq_min_hz": 100.0,
        "freq_max_hz": 200.0,
        "total_bw_hz": 100.0,
        "abs_chan_width_hz": 1.0,
    }
    entry = type("Entry", (), {
        "mous_uid": "uid___A001_X1_X2",
        "path": "cont.dat",
        "target_name": "target",
        "spw_id": 1,
        "flags": [],
        "continuum_ranges_hz": [(120.0, 130.0), (150.0, 160.0)],
    })()

    issues = []
    row = summarize_contdat_entry(entry, spw_row, issues)

    assert math.isclose(row["line_range_width_hz"], 100.0)
    assert math.isclose(row["summed_line_width_hz"], 80.0)
    assert math.isclose(row["buffered_summed_line_width_hz"], 100.0)
    assert row["contdat_parse_warning"] == ""


def test_interval_helpers_merge_and_complement():
    merged = merge_intervals([(1.0, 2.0), (1.5, 3.0), (5.0, 6.0)])
    assert merged == [(1.0, 3.0), (5.0, 6.0)]

    complement = complement_intervals([(2.0, 4.0), (6.0, 7.0)], 1.0, 8.0)
    assert complement == [(1.0, 2.0), (4.0, 6.0), (7.0, 8.0)]
    assert math.isclose(interval_total_width(complement), 4.0)

    buffered = buffer_intervals([(4.0, 6.0)], 2.0, 1.0, 8.0)
    assert buffered == [(2.0, 8.0)]


def test_spwinfo_frequency_edges_and_velocity():
    spw = SpwInfo(
        spw_id=21,
        spw_name="X",
        nchan=240,
        frame="TOPO",
        ch0_mhz=111302.132,
        chan_width_khz=3906.250,
        total_bw_khz=937500.0,
        center_mhz=111768.9290,
        corrs="XX YY",
    )

    assert spw.freq_max_hz > spw.freq_min_hz
    assert spw.abs_chan_width_hz == 3906.25 * 1e3
    assert spw.velocity_resolution_kms > 0


def test_is_target_listobs_path_only_accepts_targets_ms():
    from pathlib import Path

    assert is_target_listobs_path(Path("sessionsession_1__uid___A002_X1_X2_targets.ms__listobs.txt"))
    assert not is_target_listobs_path(Path("sessionsession_1__uid___A002_X1_X2_targets_line.ms__listobs.txt"))
    assert not is_target_listobs_path(Path("sessionsession_1__uid___A002_X1_X2.ms__listobs.txt"))


def test_process_mous_dir_ignores_full_and_targets_line_listobs(tmp_path):
    mous_dir = tmp_path / "uid___A001_X1_X2"
    listobs_dir = mous_dir / "listobs"
    listobs_dir.mkdir(parents=True)

    target_text = "\n".join(
        [
            "  Date        Timerange (UTC)          Scan  FldId FieldName             nRows     SpwIds   Average Interval(s)    ScanIntent",
            "  04-Jan-2025/14:07:55.2 - 14:17:31.9     6      2 science                100  [5]  [6.05] [OBSERVE_TARGET#ON_SOURCE]",
            "           (nRows = Total number of rows per scan)",
            "Fields: 1",
            "  ID   Code Name                RA               Decl           Epoch   SrcId      nRows",
            "  2    none science             15:47:48.990000 +22.03.03.20000 ICRS    2        1737120",
            "Spectral Windows:  (1 unique spectral windows and 1 unique polarization setups)",
            "  SpwID  Name                                       #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) CtrFreq(MHz) BBC Num  Corrs",
            "  5      X#BB_1#SW-01#FULL_RES    128   TOPO   99038.315    -15625.000   2000000.0  98046.1280        1  XX  YY",
            "Sources: 1",
        ]
    )
    full_text = "\n".join(
        [
            "  Date        Timerange (UTC)          Scan  FldId FieldName             nRows     SpwIds   Average Interval(s)    ScanIntent",
            "  04-Jan-2025/14:07:55.2 - 14:17:31.9     6      2 fullscience            100  [9]  [6.05] [OBSERVE_TARGET#ON_SOURCE]",
            "           (nRows = Total number of rows per scan)",
            "Fields: 1",
            "  ID   Code Name                RA               Decl           Epoch   SrcId      nRows",
            "  2    none fullscience         15:47:48.990000 +22.03.03.20000 ICRS    2        1737120",
            "Spectral Windows:  (1 unique spectral windows and 1 unique polarization setups)",
            "  SpwID  Name                                       #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) CtrFreq(MHz) BBC Num  Corrs",
            "  9      X#BB_3#SW-01#FULL_RES    128   TOPO  109053.940     15625.000   2000000.0 110046.1280        3  XX  YY",
            "Sources: 1",
        ]
    )
    line_text = full_text.replace("fullscience", "linetarget").replace("[9]", "[16]").replace("  9      ", "  16     ")

    (listobs_dir / "sessionsession_1__uid___A002_X1_X2_targets.ms__listobs.txt").write_text(target_text, encoding="utf-8")
    (listobs_dir / "sessionsession_1__uid___A002_X1_X2.ms__listobs.txt").write_text(full_text, encoding="utf-8")
    (listobs_dir / "sessionsession_1__uid___A002_X1_X2_targets_line.ms__listobs.txt").write_text(line_text, encoding="utf-8")

    (mous_dir / "cont.dat").write_text("", encoding="utf-8")

    issues = []
    spw_rows, cont_rows = process_mous_dir(mous_dir, issues)

    assert len(spw_rows) == 1
    assert spw_rows[0]["spw_id"] == 5
    assert spw_rows[0]["target_name"] == "science"
    assert cont_rows == []
