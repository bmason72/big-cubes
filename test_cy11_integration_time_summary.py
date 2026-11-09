from __future__ import annotations

import csv
import json
from pathlib import Path

from cy11_integration_time_summary import build_compact_rows, main


def listobs_text(*, full: bool = False) -> str:
    spw_ids = "[0,5]" if full else "[5]"
    intervals = "[0.016, 6.05]" if full else "[6.05]"
    return "\n".join(
        [
            "  Date        Timerange (UTC)          Scan  FldId FieldName             nRows     SpwIds   Average Interval(s)    ScanIntent",
            f"  04-Jan-2025/14:07:55.2 - 14:17:31.9     6      2 target                 100  {spw_ids}  {intervals} [OBSERVE_TARGET#ON_SOURCE]",
            "           (nRows = Total number of rows per scan)",
            "Fields: 1",
            "  ID   Code Name                RA               Decl           Epoch   SrcId      nRows",
            "  2    none target              15:47:48.990000 +22.03.03.20000 ICRS    2        100",
            "Spectral Windows:  (2 unique spectral windows and 1 unique polarization setups)",
            "  SpwID  Name                                       #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) CtrFreq(MHz) BBC Num  Corrs",
            "  0      X#ALMA_RB_03#BB_1#SQLD                         1   TOPO   98000.000   2000000.000   2000000.0  98000.0000        1  XX  YY",
            "  5      X#ALMA_RB_03#BB_1#SW-01#FULL_RES             128   TOPO   99000.000    -15625.000   2000000.0  98000.0000        1  XX  YY",
            "Sources: 1",
            "Antennas: 3:",
            "  ID   Name  Station   Diam.    Long.         Lat.                Offset from array center (m)                ITRF Geocentric coordinates (m)",
            "   0    DA41  A001      12.0 m   -067.45.15.9  -22.53.28.0    0.0       0.0       0.0  1 2 3",
            "   1    DA42  A002      12.0 m   -067.45.15.9  -22.53.28.0  100.0       0.0       0.0  1 2 3",
            "   2    DA43  A003      12.0 m   -067.45.15.9  -22.53.28.0    0.0     100.0       0.0  1 2 3",
        ]
    )


def write_ppr(path: Path) -> None:
    path.write_text(
        "<SciPipeRequest><ProjectSummary><ProposalCode>2024.1.00001.S</ProposalCode></ProjectSummary>"
        "<ProcedureTitle>hifa_cal</ProcedureTitle><SchedBlockRef entityId='uid://A001/X1/X2'/>"
        "<SBTitle>Undefined</SBTitle></SciPipeRequest>",
        encoding="utf-8",
    )


def test_main_uses_full_ms_fallback_and_excludes_sqld(tmp_path):
    source = tmp_path / "source"
    mous = source / "uid___A001_X1_X1"
    listobs = mous / "listobs"
    listobs.mkdir(parents=True)
    path = listobs / "sessionsession_1__uid___A002_X1_X2.ms__listobs.txt"
    path.write_text(listobs_text(full=True), encoding="utf-8")
    write_ppr(mous / "member.pprequest.xml")
    (mous / "pipeline_stats_uid___A001_X1_X1.json").write_text(
        json.dumps({"uid___A001_X1_X1": {"EB": {"uid___A002_X1_X2.ms": {"L80": {"value": 235.4}}}}}),
        encoding="utf-8",
    )
    with (source / "listobs_extraction_index.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["mous_uid", "status"])
        writer.writeheader()
        writer.writerow({"mous_uid": "uid___A001_X1_X1", "status": "OK"})

    output = tmp_path / "output"
    assert main([str(source), str(output)]) == 0

    with (output / "mous_integration_times.csv").open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    assert row["project_code"] == "2024.1.00001.S"
    assert row["sb_uid"] == "uid://A001/X1/X2"
    assert row["target_names"] == "target"
    assert row["bands"] == "3"
    assert row["input_kinds"] == "full_ms_fallback"
    assert float(row["integration_time_mode_s"]) == 6.05
    assert "128ch@2GHz" in row["spectral_signature"]
    assert "1ch" not in row["spectral_signature"]
    assert row["estimated_configurations"] == "C43-3"

    with (output / "cy11_mous_science_setups.csv").open(newline="", encoding="utf-8") as handle:
        compact = next(csv.DictReader(handle))
    assert list(compact) == [
        "mous_uid",
        "integration_time_s",
        "array",
        "band",
        "spectral_signature",
        "l80_m",
        "estimated_config",
    ]
    assert compact["mous_uid"] == "uid___A001_X1_X1"
    assert float(compact["integration_time_s"]) == 6.05
    assert compact["spectral_signature"] == "TG01[128ch@2GHz]"


def test_compact_rows_split_band_to_band_mous():
    rows = []
    for band, spw_id, nchan, bandwidth in [("4", 1, 128, 2e9), ("9", 2, 1920, 1.875e9)]:
        rows.append(
            {
                "mous_uid": "m1",
                "eb_uid": "e1",
                "array": "12m",
                "band": band,
                "tuning_group": "TG01",
                "spw_id": spw_id,
                "nchan": nchan,
                "total_bw_hz": bandwidth,
                "median_interval_s": 6.14,
                "l80_m": 235.4,
            }
        )

    compact = build_compact_rows(rows)

    assert len(compact) == 2
    assert [row["band"] for row in compact] == ["4", "9"]
    assert {row["integration_time_s"] for row in compact} == {6.14}
    assert compact[0]["spectral_signature"] == "TG01[128ch@2GHz]"
    assert compact[1]["spectral_signature"] == "TG01[1920ch@1.875GHz]"
