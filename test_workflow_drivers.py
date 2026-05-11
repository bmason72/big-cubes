from __future__ import annotations

import json
from pathlib import Path

from run_layout import next_run_dir, write_config_json


HERE = Path(__file__).resolve().parent
SAMPLE_DB = HERE / "data" / "wsu_datarates_mit_per_mous_initial_goal_20250423_head_and_sample.ecsv"


def test_next_run_dir_advances_suffix(tmp_path):
    first = next_run_dir("drrupPl", out_root=tmp_path)
    first.mkdir(parents=True)
    second = next_run_dir("drrupPl", out_root=tmp_path)
    assert first.name.endswith("a")
    assert second.name.endswith("b")


def test_write_config_json_stringifies_nonfinite_floats(tmp_path):
    target = tmp_path / "config.json"
    write_config_json(target, {"limit": float("inf"), "value": 1.5})
    raw = target.read_text(encoding="utf-8")
    assert '"Infinity"' in raw
    assert ": Infinity" not in raw


def test_pipeline_skip_realizations_cli_writes_stats(tmp_path, monkeypatch):
    import pipeline

    real_dir = tmp_path / "realizations"
    real_dir.mkdir()

    seen = {}

    def fake_find_realizations(cfg, outdir):
        seen["outdir"] = outdir
        seen["n_realizations"] = cfg.realizations.n_realizations
        return ["a.ecsv", "b.ecsv"]

    def fake_compute_stats(paths, cfg):
        seen["paths"] = paths
        return {"_metadata": {"n_paths": len(paths)}}

    monkeypatch.setattr(pipeline, "find_realizations", fake_find_realizations)
    monkeypatch.setattr(pipeline, "compute_stats", fake_compute_stats)
    monkeypatch.setattr(
        pipeline,
        "generate_realizations",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("generate_realizations should not run")),
    )

    stats_out = tmp_path / "stats.json"
    rc = pipeline.main(
        [
            "--skip-realizations",
            "--outdir", str(real_dir),
            "--n-realizations", "2",
            "--stats-out", str(stats_out),
        ]
    )
    assert rc == 0
    assert stats_out.exists()
    assert seen["outdir"] == str(real_dir)
    assert seen["n_realizations"] == 2
    assert seen["paths"] == ["a.ecsv", "b.ecsv"]


def test_drrup_pipeline_driver_creates_root_outputs(tmp_path, monkeypatch):
    import drrup_pipeline

    run_dir = tmp_path / "run"

    monkeypatch.setattr(drrup_pipeline, "find_realizations", lambda cfg, outdir: ["r1.ecsv"])
    monkeypatch.setattr(drrup_pipeline, "compute_stats", lambda paths, cfg: {"_metadata": {"n": len(paths)}})
    monkeypatch.setattr(drrup_pipeline, "load_database", lambda path: object())

    def fake_generate_memo_tables(stats, out_dir):
        out = Path(out_dir)
        path = out / "memo.tex"
        path.write_text("memo", encoding="utf-8")
        return [str(path)]

    def fake_generate_sdd_tables(stats, out_dir):
        out = Path(out_dir)
        path = out / "sdd.tex"
        path.write_text("sdd", encoding="utf-8")
        return [str(path)]

    def fake_write_all_tables_document(paths, out_path):
        Path(out_path).write_text("all", encoding="utf-8")
        return str(out_path)

    def fake_generate_sdd_plots(db, output_dir, cfg, realization_dir, overlay_realizations):
        out = Path(output_dir)
        plot = out / "plot.png"
        plot.write_bytes(b"png")
        return [str(plot)]

    monkeypatch.setattr(drrup_pipeline, "generate_memo_tables", fake_generate_memo_tables)
    monkeypatch.setattr(drrup_pipeline, "generate_sdd_tables", fake_generate_sdd_tables)
    monkeypatch.setattr(drrup_pipeline, "write_all_tables_document", fake_write_all_tables_document)
    monkeypatch.setattr(drrup_pipeline, "generate_sdd_plots", fake_generate_sdd_plots)

    rc = drrup_pipeline.main(
        [
            "--run-dir", str(run_dir),
            "--skip-realizations",
            "--realization-dir", str(tmp_path / "realizations"),
        ]
    )
    assert rc == 0
    assert (run_dir / "stats.json").exists()
    assert (run_dir / "tables" / "all_tables.tex").exists()
    assert (run_dir / "plots" / "plot.png").exists()
    cfg = json.loads((run_dir / "config.json").read_text(encoding="utf-8"))
    assert cfg["tool"] == "drrup_pipeline"


def test_nchan_pipeline_driver_creates_root_outputs(tmp_path, monkeypatch):
    import nchan_pipeline

    run_dir = tmp_path / "run"
    calls = []
    monkeypatch.setattr(nchan_pipeline.spw_tabulate, "main", lambda argv=None: calls.append(("tab", argv)) or 0)
    monkeypatch.setattr(nchan_pipeline.spw_setup_summary, "main", lambda argv=None: calls.append(("setup", argv)) or 0)
    monkeypatch.setattr(nchan_pipeline.ryan_cy11_current_data_intents_summary, "main", lambda argv=None: calls.append(("current", argv)) or 0)
    monkeypatch.setattr(nchan_pipeline.ryan_cy11_wsu_projection, "main", lambda argv=None: calls.append(("wsu", argv)) or 0)

    rc = nchan_pipeline.main(["--run-dir", str(run_dir), "--input-dir", "ryanCy11"])
    assert rc == 0
    assert (run_dir / "tabulation").exists()
    assert (run_dir / "setup_summary").exists()
    assert (run_dir / "current_data_properties").exists()
    assert (run_dir / "wsu_projection").exists()
    assert (run_dir / "config.json").exists()
    assert [name for name, _ in calls] == ["tab", "setup", "current", "wsu"]


def test_nchan_pipeline_passes_intent_aware_projection_flags(tmp_path, monkeypatch):
    import nchan_pipeline

    run_dir = tmp_path / "run"
    calls = []
    monkeypatch.setattr(nchan_pipeline.spw_tabulate, "main", lambda argv=None: 0)
    monkeypatch.setattr(nchan_pipeline.spw_setup_summary, "main", lambda argv=None: 0)
    monkeypatch.setattr(nchan_pipeline.ryan_cy11_current_data_intents_summary, "main", lambda argv=None: 0)
    monkeypatch.setattr(nchan_pipeline.ryan_cy11_wsu_projection, "main", lambda argv=None: calls.append(argv) or 0)

    rc = nchan_pipeline.main(
        [
            "--run-dir", str(run_dir),
            "--input-dir", "ryanCy11",
            "--intent-aware",
            "--apply-calibrator-cap",
            "--calibrator-cap-kms", "0.8",
        ]
    )
    assert rc == 0
    assert calls == [[
        str(run_dir / "tabulation"),
        str(run_dir / "wsu_projection"),
        "--intent-aware",
        "--apply-calibrator-cap",
        "--calibrator-cap-kms",
        "0.8",
    ]]


def test_recompute_cli_writes_db_and_config(tmp_path):
    import recompute

    run_dir = tmp_path / "recompute_run"
    rc = recompute.main(
        [
            "--db-path", str(SAMPLE_DB),
            "--scenario", "VALIDATION_TWEAKED_CONFIG",
            "--run-dir", str(run_dir),
        ]
    )
    assert rc == 0
    assert (run_dir / "recomputed_validation_tweaked.ecsv").exists()
    assert (run_dir / "config.json").exists()


def test_recompute_realization_dir_uses_existing_files_without_default_count(tmp_path):
    import recompute
    from config import DEFAULT_CONFIG

    in_dir = tmp_path / "realizations"
    out_dir = tmp_path / "recomputed"
    in_dir.mkdir()

    for idx in range(2):
        target = in_dir / f"{DEFAULT_CONFIG.realizations.filename}_{idx:04d}.ecsv"
        target.write_text(SAMPLE_DB.read_text(encoding="utf-8"), encoding="utf-8")

    written = recompute.recompute_realization_dir(
        str(in_dir),
        str(out_dir),
        recompute._load_named_scenario("VALIDATION_TWEAKED_CONFIG"),
    )

    assert len(written) == 2
    assert all(Path(path).exists() for path in written)
