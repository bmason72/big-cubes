# Refactor Plan: Config-Driven Scenario Recomputation

Author context: drafted 2026-04-21 in a Claude Code session that completed
pipeline Phases 1–3 (see `PLAN.md`).  Handoff target: Codex, implementing
this afternoon.  A running context / memory dump is included in
§[Context Handoff](#context-handoff) at the end so nothing is lost in
the switch.

---

## 1. Goal

Make `config.py` *load-bearing* so that editing it (or passing an
alternate config object) re-derives the projected quantities — data
rates, visibility rates, data volumes, product sizes, system
performance — from the per-MOUS ecsv **without** touching the existing
`wsu_db.py` notebook-era computation.  The current ecsv + existing
BLC/M1/M4/M5 columns remain the canonical validation target; a
*parallel* recompute layer drives scenarios.

**Success criterion:** running the new recompute path with a
`MEMO_CONFIG` that encodes the memo's assumptions must reproduce each
derived column already in the ecsv to within a documented numerical
tolerance.  Deviations are bugs, not features.

---

## 2. Architecture

```
                         ┌──────────────────────────┐
 scenario: ScenarioCfg → │ recompute.py             │
                         │                          │
 data/...ecsv           → │ fixed per-MOUS inputs     │
                         │  (L80, band/freq,        │
                         │   velres per SPW, nspw,   │
                         │   target_tos, cal_tos,   │
                         │   pb/imsize/fov, ...)    │
                         │                          │
                         │ → derived columns        │
                         │   (datarate, visrate,    │
                         │    datavol, productsize, │
                         │    sysperf, nchan_agg)   │
                         └──────────────────────────┘
                                       │
                                       ▼
                         pipeline.compute_stats()
                              tables.py
                              plots.py
```

- **Keep `wsu_db.py` untouched** except if we discover a bug.
- **Keep `pipeline.compute_stats` unchanged** on the read-columns path.
- **New module `recompute.py`** produces a *new* QTable with the same
  schema as the existing ecsv but with scenario-driven values in the
  derived columns.  Downstream code (`compute_stats`, `tables`,
  `plots`) is unchanged — it just reads from the recomputed table.
- Scenarios are Python dataclass instances (sub-configs of
  `PipelineConfig`), not YAML.  Rationale: type-checked, discoverable
  in an IDE, diffable in git, and lets the user plug in callables
  (e.g. for the configuration→tint heuristic) without inventing a DSL.

---

## 3. What each scenario knob does

| Knob | Section | Affects |
|------|---------|---------|
| 12m configuration from L80 | §4.1 | integration time, baseline count? (see note) |
| Integration-time table per (array, configuration, milestone) | §4.1 | datarate, visrate, datavol, sysperf |
| Elevation factor on tint (default 1.0) | §4.1 | same as above |
| Antenna counts per (array, milestone) | §4.2 | datarate, visrate, sysperf |
| Channel-averaging minimums per band | §4.3 | nchan_agg, downstream |
| Stepped2 velocity bins | §4.3 | nchan_agg, downstream |
| Per-SPW nchan cap | §4.3 | nchan_agg, cubesize, productsize |
| Per-MOUS nchan cap + preserve-fraction X | §4.3 | same |
| Cube-size cap | §4.4 | productsize |
| Pixels per beam (default 5) | §4.4 | imsize, cubesize, mfssize, productsize |
| Mitigation mode: "existing" vs "caps" | §4.4 | which path decides nchan_product, imsize |
| Calibrator spectral-resolution cap `v_cal_cap(array, config, target_velres)` | §4.5 | datarate, visrate, datavol (NOT productsize) |
| Sysperf flops / efficiency knobs | §4.6 | sysperf |

---

## 4. Derived quantity recompute — formulas

Notation: each MOUS row has fixed inputs (L80, band, freq, nspw,
`target_tos`, `cal_tos`, and a per-SPW requested velocity-resolution
vector `velres_req[s]`).  Per milestone `m ∈ {M1, M4, M5}` and BLC we
produce the full set of derived columns.

### 4.1 Integration time

**Configuration classification from L80** (12m only).  Lookup table
from the ALMA Cycle-13 Technical Handbook §7 (user supplied):

| Configuration | L80 (m) | Max baseline (m) |
|---------------|---------|------------------|
| 7-m           | 30.7    | 45.0             |
| C43-1         | 107.1   | 160.7            |
| C43-2         | 143.8   | 313.7            |
| C43-3         | 235.4   | 500.2            |
| C43-4         | 369.2   | 783.5            |
| C43-5         | 623.8   | 1397.9           |
| C43-6         | 1172.5  | 2516.9           |
| C43-7         | 1673.1  | 3637.8           |
| C43-8         | 3527.3  | 8547.7           |
| C43-9         | 6482.6  | 13894.2          |
| C43-10        | 8685.9  | 16194.0          |

Classification rule: 12m MOUS gets the configuration whose L80 is
*nearest in log space* to the MOUS's L80.  7m MOUSes bypass this —
they get the single "7-m" configuration.  The classifier should be a
pure function `classify_config(array, L80) -> str` so it can be
swapped out for fancier heuristics later.

**Integration time lookup:**

```
tint_i = elevation_factor * tint_table[(array_i, config_i, milestone)]
```

`tint_table` is a dict in config.py; the implementer should seed it
with the current memo defaults (12m: 3.072 s all configs at M4/M5,
6.144 s at M1; 7m: 9.984 s all milestones) and leave slots for the
user to override per configuration.

`elevation_factor` is a single scalar in config; default `1.0`.

### 4.2 Antenna counts

Already a table in `config.py`
(`ArrayParams.nant_typical / nant_array / nant_all`).  Scenario
selects one flavor per milestone.  Default: memo uses `nant_typical`
everywhere.

Number of baselines: `Nbase = Nant * (Nant - 1) / 2`.

### 4.3 Channels (nchan_agg) + caps

Per SPW s in MOUS i:

1. **Target velocity resolution** starts at `velres_req[s]`.  Apply
   the stepped2 binning (config-controlled) to get `velres_binned[s]`.
   Stepped2 default: `>10 → 10`, `2–10 → 2`, `0.5–2 → 0.5`, `0.1–0.5
   → 0.1`, `<0.1 → native`.
2. Convert velocity resolution to channel width in kHz via
   `chanwidth_kHz = (velres_binned / c) * freq_GHz * 1e6`, then enforce
   the per-band floor:
   `chanwidth_kHz = max(chanwidth_kHz, talon_channel_kHz * chanavg_min[band, milestone])`.
3. `nchan_raw[s] = floor(spw_bandwidth_GHz * 1e6 / chanwidth_kHz)`.
4. **Per-SPW cap**: `nchan[s] = min(nchan_raw[s], nchan_per_spw_cap)`
   (if the cap is set; None = no cap).  This is a post-hoc ceiling —
   it effectively forces additional channel averaging in that SPW.
5. **Per-MOUS cap** `N_mous_cap`, with preservation policy
   `preserve_fraction = 1 - X` (user default `X = 0.75`, so 25 % of
   SPWs are preserved at their post-per-SPW-cap resolution):
   - Sort SPWs by `nchan[s]` descending.
   - Preserve the top `k = ceil(preserve_fraction * nspw)` SPWs as-is.
   - Sum of preserved nchans = `P`.
   - Remaining budget: `B = N_mous_cap - P`.
   - If `B <= 0`: preservation would already exceed the cap;
     implementer choice — either fail loudly (config is internally
     inconsistent) or fall back to uniform averaging of everything
     to fit.  Recommend: log a warning, uniformly clip.
   - Otherwise distribute `B` across the remaining `nspw - k` SPWs
     uniformly: `nchan[s ∈ unpreserved] = floor(B / (nspw - k))`.

6. `nchan_agg_i = sum_s nchan[s]`.

**Note on per-SPW velres vector:** the current ecsv stores one
`velres` per MOUS, not per SPW.  In practice WSU SPWs have uniform
requested resolution per MOUS (the memo assumes it).  If we later
admit heterogeneous SPWs, this formula already generalizes.  For now
treat `velres_req[s] = velres_req[MOUS]` for all s.

### 4.4 Imaging: cubesize, mfssize, productsize

**Pixels per beam** is configurable; default 5.  imsize already in
DB was computed at 5 px/beam.  To reflect a config change:

```
imsize_scenario = imsize_db * (pixels_per_beam_scenario / 5)
```

(Integer-rounded; watch for off-by-one on small MOUSes.)

**Cube size (per SPW):**

```
cubesize[s] = 4 * imsize^2 * nchan[s] / 1e9    # GB
mfssize     = 4 * imsize^2 / 1e9                # GB
```

**Product size:** `productsize = 2 * sum_s (cubesize[s] + mfssize) = 2 * nspw * (cubesize_avg + mfssize)`
but sum-form is exact when SPWs differ.  Factor of 2 is the
calibrators/weblog allowance.

**Mitigation modes:**

- `mode == "existing"`:  use the DB's `mitigatedprodsize`,
  `initialprodsize`, and per-milestone `mitigated/unmitigated`
  columns directly.  Scenario caps are ignored for the imaging
  path.
- `mode == "caps"`:  ignore DB mitigation entirely; apply only
  the scenario caps:
  - `cubesize_cap`: if `cubesize[s] > cubesize_cap`, reduce `nchan[s]`
    for that SPW to fit.  (This is where cube cap interacts with
    MOUS nchan cap — see note below.)
  - Per-SPW and per-MOUS nchan caps (§4.3) already applied upstream.

**Interaction between cube-size cap and nchan caps:** apply caps in
this order so there are no surprises: per-SPW nchan → per-MOUS nchan
preservation-policy → cube-size cap → recompute `nchan_agg`.  If the
cube cap triggers further reductions, those propagate back to
`nchan[s]` and thus to data rate / data volume.  Document that the
cube cap is the *last* constraint applied so the nchan caps don't
mask it.

### 4.5 Calibrator spectral-resolution correction

**Rationale:** in real observations, calibrators don't need the same
spectral resolution as the target.  If the target is at very fine
resolution (e.g. < 1 km/s) we can (and should) record calibrator
visibilities at coarser resolution, shrinking the cal visibility
stream.  This reduces visibility volumes and data rates but not
product sizes (product comes from imaging the target, not the cal).

Using approach (b) from the clarification thread — a post-hoc
multiplier — because the base DB has `datavol_*_cal` and
`datavol_*_target_tot` columns but (per the user's read) does *not*
carry a distinct cal-visibility-resolution axis we could drive
directly.  If implementer finds those columns are in fact broken out
per SPW res, prefer approach (a) and note the change.

**Config surface:**

```python
cal_velres_cap_kms: Callable[(array, config_id, target_velres_kms)] -> float | None
```

Default: single-scalar cap `1.0 km/s` for all (array, config),
irrespective of target resolution.  The callable returns the cap
velocity resolution `v_cap` to use for calibrator recording, or
`None` to disable correction.

**Math** (per MOUS i, per SPW s, but collapsible to per-MOUS if
velres_req is uniform across SPWs):

```
if v_cap is None or velres_req[s] >= v_cap:
    R_s = 1.0                          # no reduction
else:
    R_s = v_cap / velres_req[s]        # e.g. v_cap=1 km/s, target=0.1 → R=10
```

`f_c_i = cal_tos_i / total_tos_i` (or `cal_vol_i / total_vol_i`; user
derivation uses the fraction of visibilities).  Verify in
implementation which convention matches the DB; they should agree
modulo tiny rounding.

Volume multiplier for MOUS i (averaged over SPWs in the uniform
case): `M_i = 1 - f_c_i + f_c_i / R_i`.

Apply to:
- `datarate_total_i  *= M_i`   ← user confirmed: yes, rate is affected
- `visrate_total_i   *= M_i`   ← same
- `datavol_total_i   *= M_i`   ← same (science + (reduced cal))
- `datavol_cal_i     *= 1/R_i` (if we choose to expose separately)

Do **not** apply to:
- `productsize`, `cubesize`, `mfssize`   ← imaging uses target only

Downstream of visibility (sysperf): the flops/vis coefficient stays
the same, but the vis count drops → sysperf goes down proportionally.

### 4.6 Sysperf

```
sysperf_i = (k_major * multiscale_factor * visrate_i * flops_per_vis) /
            (core_efficiency * parallelization_efficiency)
```

`flops_per_vis` is the imaging mode coefficient (std / aproj / wproj /
awproj).  Existing `calc_sysperf` logic in `wsu_db.py` does
band-specific mode selection (Band 1 long baselines → wproject).
Plan: reimplement `calc_sysperf` as a pure function in
`recompute.py`, taking the same inputs plus a `SysperfParams` from
config.  Keep the wproject L80 threshold configurable (currently
6200 m).

---

## 5. Config surface (extension of `config.py`)

New dataclasses (sketch):

```python
@dataclass(frozen=True)
class ConfigClassifier:
    """L80 -> Cycle-13 configuration label."""
    configs: Dict[str, float]   # {"7-m": 30.7, "C43-1": 107.1, ...}
    def classify(self, array: str, L80_m: float) -> str: ...

@dataclass(frozen=True)
class TintLookup:
    """(array, config_id, milestone) -> tint in seconds."""
    table: Dict[Tuple[str, str, str], float]
    elevation_factor: float = 1.0
    def tint(self, array, config_id, milestone) -> float: ...

@dataclass(frozen=True)
class MitigationCaps:
    mode: str = "existing"          # or "caps"
    nchan_per_spw_cap: Optional[int] = None
    nchan_mous_cap: Optional[int] = None
    preserve_fraction: float = 0.25    # top 25% of SPWs preserved (X=0.75)
    cubesize_cap_GB: Optional[float] = None
    pixels_per_beam: int = 5

@dataclass(frozen=True)
class CalibratorCorrection:
    enabled: bool = False
    v_cap_kms: Union[float, Callable] = 1.0   # callable: (array, config, v_target) -> v_cap
    @property
    def cap_fn(self) -> Callable: ...

@dataclass(frozen=True)
class ScenarioConfig:
    name: str
    base: PipelineConfig          # existing config wrapper
    configs: ConfigClassifier
    tint: TintLookup
    mitigation: MitigationCaps
    cal_correction: CalibratorCorrection
    # Everything else inherited from `base`.
```

Ship with two canonical scenarios:

- `MEMO_CONFIG` — reproduces the memo/current ecsv (validation anchor).
- `MEMO_3PX_CONFIG` — memo defaults but `pixels_per_beam=3`, to give
  the user a first meaningful exploration knob out of the box.

---

## 6. File layout

New:
- `recompute.py` — core recompute engine.
- `scenarios.py` — catalog of `ScenarioConfig` instances; imported by
  CLI.  Keep each scenario as a module-level constant.
- `data/cycle13_configurations.py` — the L80 table as Python data (also
  consulted by `ConfigClassifier`).
- `test_recompute.py` — unit tests (config classifier, tint lookup,
  cap algorithm, cal correction math) + integration test (memo
  round-trip).

Extend:
- `config.py` — add the new dataclasses from §5.  Keep the existing
  `PipelineConfig` untouched; compose via `ScenarioConfig`.
- `pipeline.py` — accept `--scenario NAME`; if set, call
  `recompute.recompute_db(db, scenario)` before `compute_stats`.
- `tables.py` / `plots.py` — no code changes needed; they consume the
  stats dict.  Only the output paths change (per scenario).
- `test_pipeline.py` — add memo-scenario round-trip test; downgrade
  existing reference-value tests to `@pytest.mark.parametrize` on
  the memo scenario so non-memo scenarios don't hit them.
- `Makefile` — add `make validation SCENARIO=name` target.  Default
  scenario remains `memo`.

---

## 7. Implementation order (2–4 hour target)

1. **Config extension** (§5).  Pure plumbing; no behavior change.
2. **Configuration classifier** (§4.1).  Unit test against the
   Cycle-13 table.
3. **Recompute: datarate, visrate, datavol** (§4.1, 4.2, 4.3 with
   mode="existing" and `cal_correction.enabled=False`).  At this
   point the memo-scenario round-trip test should pass.
4. **Mitigation caps** (§4.3, 4.4, mode="caps").  Unit tests on the
   preservation-policy algorithm with hand-built SPW vectors.
5. **Cube-size cap + pixels_per_beam** (§4.4).
6. **Calibrator correction** (§4.5).  Unit tests on the `M_i` formula
   with hand calcs.
7. **Sysperf recompute** (§4.6).  Verify against existing column.
8. **Band 1/2 realizations**: wrap `wsu_db.generate_db_realizations`
   with an optional recompute pass after synthesis.  See §8 — this
   may be deferrable.
9. **CLI + Makefile target**.
10. **Smoke tests** for `MEMO_3PX_CONFIG` (pixels_per_beam=3) and
    one caps-mode scenario.

---

## 8. Band 1/2 realizations — minimal path

`add_bands()` synthesises Band 1/2 rows from Band 3 with rescaled
frequency-dependent quantities and milestone-specific bandwidths.
The *structural* assumptions (bandwidths, chanavg mins) are
scenario-independent per user §6 response.  So: Band 1/2 synthesis
can run **before** the recompute pass, producing synthetic rows
with the same *input* columns; the recompute pass then derives the
scenario-specific output columns for all rows uniformly.

Acceptable caveat: if a scenario changes `nant`, `tint`, caps, or
calibrator correction, the Band 1/2 rows reflect those changes
automatically because recompute runs after synthesis.  User is fine
skipping Band 1/2 refinement for now.

---

## 9. Validation

### 9.1 Memo round-trip (the acceptance test)

Pseudocode:

```python
db = load_database(DB_PATH)
db_recomputed = recompute_db(db, MEMO_CONFIG)

for col in DERIVED_COLUMNS:
    assert_allclose(db_recomputed[col], db[col], rtol=1e-3, atol=1e-6)
```

`DERIVED_COLUMNS` = every `blc_*`, `wsu_*_initial_*`, `wsu_*_ms4_*`,
`wsu_*_goal_*` column.  A discrepancy means either (a) `MEMO_CONFIG`
doesn't actually encode the memo assumptions or (b) `wsu_db.py` has
a bug.  Tighten `rtol` as confidence builds.

### 9.2 Existing reference-value tests

The current `TestTableReproduction` class asserts generated tables
match memo reference .tex cells.  After the refactor these run with
`scenario=memo` only.  Non-memo scenarios are only smoke-tested.

### 9.3 Unit tests to add

- `test_classify_config` — every row of the Cycle-13 table classifies
  to itself.
- `test_classify_config_interpolation` — L80 = 500 m → C43-4 (log
  nearest).
- `test_mous_cap_preserves_topk` — given a fake SPW vector, confirm
  the preservation policy.
- `test_cal_correction_multiplier` — hand calcs at `f_c=0.3, R=5`.
- `test_cal_correction_nontrigger` — v_target > v_cap → M = 1.

---

## 10. Open questions deferred to implementation

1. **`_target_tot` column semantics.**  The DB has
   `*_target_tot` and `*_cal` columns.  Before wiring (§4.5) confirm
   `target_tot + cal == total` and `target_tot` is target-only.  If
   the column actually carries a cal-resolution-reducible component
   already, approach (a) becomes viable (use the columns directly
   instead of the multiplier).
2. **Integration-time defaults per configuration.**  User has not yet
   specified per-configuration `tint` values; the plan seeds them to
   the current memo values (3.072 / 6.144 / 9.984 s) and leaves the
   override slots explicit.  User will fill in when exploring.
3. **"both"-array aggregation bug.**  Currently
   `pipeline._compute_per_sample_stats` averages rather than sums 12m
   + 7m for the "both" row of datarate/visrate/sysperf.  Flagged
   separately in memory; **not** in scope for this refactor, but the
   new recompute path should at least not *worsen* the behavior.
   Decide after the refactor lands.
4. **Sysperf imaging-mode selection.**  The memo logic picks
   std/aproj/wproj based on band + L80.  Port faithfully for the
   memo round-trip; make the thresholds configurable afterwards.

---

## 11. Risks

| Risk | Mitigation |
|------|------------|
| Memo round-trip fails by >1 %. | Inspect each formula against `wsu_db.py`; most likely culprits are stepped2 binning edge cases, chanavg_min rounding, and integer nchan truncation. |
| Cap algorithm mishandles `nspw=1` or `nspw<k`. | Unit tests with nspw ∈ {1, 2, 4, 8}. |
| Cal correction silently doubles when chained with existing `_cal` vs `_target_tot` math. | Explicit toggle; default OFF.  Document in module docstring that it is a *replacement* for per-column cal handling, not an addition. |
| Per-SPW velres vector: DB only has per-MOUS. | Ship code that accepts per-SPW input but defaults to broadcasting the MOUS-level value. |

---

## Context handoff

*Everything below is a verbatim dump of the state I've built up this
session.  The next tool (Codex) starts cold; this is the read-in.*

### User profile
- ALMA WSU domain expert.  Cares about correctness, reproducibility,
  and discoverable validation artifacts.  Values concise answers.
- Expects to iterate on scenarios; the *entire point* of this
  refactor is to let him explore config-driven variants without
  rewriting calculation code.  Canonical memo results must remain
  reproducible.

### Repo state (2026-04-21)
- `PLAN.md` tracks the five-phase pipeline plan; Phases 1–3 done
  (config skeleton, pipeline+tables, SDD plots).
- `make validation` reproduces canonical outputs into
  `build/validation/`.  `make test-keep` preserves test artifacts in
  `build/test_artifacts/`.  `make clean[-{test,validation}]` cleans.
- Full test suite (`make test`) = 237 tests, ~6 min; fast tier
  (`make test-fast`) = 15 tests, <1 s.  All green at time of
  drafting.
- Current ecsv uses column suffixes `_initial_stepped2_initial` for
  M1, `_ms4_stepped2_typical` for M4, `_goal_stepped2_typical` for
  M5.  BLC uses `_typical`.  See `config.MILESTONE_COLUMNS`.

### Known quirks
- `wsu_db.py` imports `ipdb` unconditionally and it isn't in the
  venv; `pipeline.py` now stubs it at import time.  Don't undo that.
- The custom `Gvis` astropy unit must be registered before reading
  the ecsv; `pipeline.register_gvis_unit()` handles this by
  importing `wsu_db` (which registers at module load).
- Base DB must be loaded as `QTable` (not `Table`) or
  `generate_db_realizations` fails on unitless arithmetic.
- The SDD sysperf plot uses `_aprojonly` sysperf columns; the
  `_allgrid` variant of `wsu_sysperf_goal_stepped2_*` does not exist
  on the base DB.  See `plots._SYSPERF_LABEL`.
- Band 1/2 realization aggregation for plots comes from
  `wsu_db.calculate_dist()`, *not* `calc_wsu_stats_allsamples` (that
  one returns a different-shaped dict; spent 20 min on this).

### Known open issue (deferred)
- `pipeline._compute_per_sample_stats` computes the "both" row by
  masking 12m ∪ 7m and taking a time-weighted average across all
  rows.  User (2026-04-21): this should be a *sum* for instantaneous
  quantities (datarate, visrate, sysperf).  Not in scope for this
  refactor.  See memory `project_both_array_datarate_bug.md`.

### Python env
- Venv: `/Users/briantest/Code/alma/venv/` (Python 3.13.5).
- Key packages: astropy 7.2, numpy 2.4, matplotlib 3.10, pytest 9.
- Full DB: `data/wsu_datarates_mit_per_mous_initial_goal_20250423.ecsv`
  (~5192 rows, ungitted).
- Sample DB (20 rows, in git):
  `data/wsu_datarates_mit_per_mous_initial_goal_20250423_head_and_sample.ecsv`.

### Source map (lines approximate)
- `wsu_db.py` (4600 lines) — notebook-era compute.  Key functions:
  `calc_datarate` (1246), `calc_visrate` (1295), `calc_cube_size`
  (near imaging section), `calc_sysperf` (3947), `calc_talon_specwidth`
  (42), `add_bands` (3084), `generate_db_realizations` (3362),
  `calculate_dist` (3427), `calc_wsu_stats_allsamples` (3612),
  `create_initial_wsu_db` (4228 — M1/4/5 bandwidth dicts),
  `calc_sysperf_average_over_subset` (4561).
- `config.py` — declarative only; not yet wired into computation.
- `pipeline.py` — `compute_stats` reads pre-computed ecsv columns;
  `generate_realizations` wraps `generate_db_realizations`.
- `tables.py` — memo + SDD LaTeX generation; tolerant cell-compare.
- `plots.py` — four SDD CCDFs.

### Memory cross-references (stored at
`/Users/briantest/.claude/projects/-Users-briantest-Code-alma-big-cubes/memory/`)
- `user_role.md`
- `project_pipeline_plan.md`
- `project_both_array_datarate_bug.md`
- `feedback_config_not_wired.md` — *this plan addresses it*
- `reference_python_env.md`

### Tone guidance for the next tool
- Keep commits small and scoped per implementation-order item.
- Run the memo round-trip test after every step; it catches most
  regressions immediately.
- When a formula diverges from `wsu_db.py`, **port the current
  behavior first** and file a follow-up; don't silently "fix" things
  on the way through, as the memo references anchor to the current
  behavior.
