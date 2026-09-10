# Changelog

All notable changes to WrinkleFE are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project aims to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

In addition to the standard categories, a **Numerical results** category
calls out any change that shifts predictions (failure-criterion fixes,
CZM-law changes, default-parameter changes, field-composition changes).
Those are the entries to scan when upgrading an engineering analysis
tool — and a results JSON's `provenance` block lets you detect which
version produced a given file.

## [Unreleased]

### Fixed
- Solver — **FE local-frame stress and strain were transformed with the
  wrong matrix and in the wrong order** (found by audit).
  `StaticSolver.recover_element_results` built one transform,
  `T_total = T_wrinkle @ T_ply` from `stress_transformation_3d`, and applied
  it to *both* the stress and the engineering-strain vector. Two distinct
  errors:
  1. **Wrong matrix for strain.** Engineering strain rotates with
     `T_eps = R T_sigma R^-1` (`R = diag(1,1,1,2,2,2)`), not `T_sigma`.
     Reusing the stress matrix mis-scaled every shear component of
     `strain_local` by a factor of two, and mixed that error into the
     normal components wherever the rotation couples normal to shear.
  2. **Wrong composition order.** The element stiffness is built as
     `C_bar = R_y(R_z(C)) = T_s(phi)⁻¹ T_s(th)⁻¹ C T_e(th) T_e(phi)`, so
     `sigma_g = C_bar eps_g` implies
     `sigma_local = T_ply @ T_wrinkle @ sigma_global` — the **ply**
     transform on the left, because the wrinkle rotation is the outer one
     and is therefore the first undone. The code had the two factors
     reversed, so `stress_local` was not the material-frame stress of the
     `stress_global` it had just computed.

  The order error vanishes identically when either angle is zero
  (`T(0) = I`), so **0 deg plies and pristine elements were always exact**
  — which is why the UD validation ledger never drifted and still shows
  zero drift after the fix. It bites on off-axis plies of a wrinkled
  laminate: measured on a `[0/45/-45/90]s` IM7/8552 coupon at A = 0.3 mm,
  `sigma_local` was off by up to **43.9 MPa (8.0 %) on the 90 deg plies**
  and ~1 % on the ±45 deg plies, against a transverse strength of order
  80 MPa. Because the failure criteria read `stress_local`, this reached
  the matrix-mode failure indices.

  Guarded going forward by an invariant that does not restate the transform
  algebra: in the material frame the *unrotated* stiffness must hold,
  `sigma_local == C_material @ eps_local`. That residual was 7.1e-3 of peak
  stress before the fix and is 6.9e-16 after. A second test records the
  physical reason the order is not arbitrary — for a 90 deg ply the wrinkle
  rotation is about that ply's own fibre axis, so material-frame `sigma_11`
  must be independent of the wrinkle angle, which only the correct order
  satisfies.

### Added
- Analysis — **the FE path now assembles the thermal initial-strain load
  vector** (issue #273, Stage 2 — *Fixes #273*). Stage 1 made
  `AnalysisConfig.delta_T` reachable on the CLT path and deliberately
  **rejected** it on the FE path, because the element formulation had no
  place to put it. It does now:
  `Hex8Element.thermal_force_vector()` integrates `∫ Bᵀ C̄ ε_th dV`,
  `GlobalAssembler.assemble_thermal_force()` scatters it into the global
  right-hand side, and stress recovery subtracts the thermal strain so
  `σ = C̄ (B u − ε_th)`. The CTE vector `[α₁, α₂, α₃, 0, 0, 0]` is rotated
  into global axes with the **inverse strain** transformation, applied in
  the same ply-z then wrinkle-y sequence the stiffness uses — because the
  wrinkle rotates the fibre frame, the CTE mismatch concentrates in
  exactly the elements where the failure criteria are evaluated, which is
  the reason to carry the term into the FE path rather than leaving it at
  laminate level. `Hex8IElement` partitions the load into `(f_u, f_a)` and
  condenses it consistently with its stiffness, and its recovered
  incompatible modes pick up `+K_aa⁻¹ f_a^th`. Internal-force assembly
  subtracts the thermal load (the element internal force is
  `K_e u_e − f^th_e`), so the CZM and progressive-damage Newton residuals
  are correct by construction. The Stage 1 refusals — config validation,
  the `run(analytical_only=False)` mirror guard, the Streamlit run
  handler, and the CLI's exit-2 — are all gone; `--delta-T` no longer
  requires `--analytical-only`. New `examples/15_cure_residual_stress.py`.
  **Verified**, not assumed: a flat `[0/90]s` IM7/8552 laminate solved
  with a statically determinate restraint reproduces the closed-form CLT
  ply stresses to **0.3 %** (issue #273 acceptance criterion 3), and free
  thermal expansion of a single element develops zero stress under
  combined ply and wrinkle rotation — the check that catches a
  stress-instead-of-strain CTE rotation, a mismatched stiffness inside the
  load integral, or a stress recovery that forgets `ε_th`.
- Docs — **the validation figures are now linked from the validation
  ledger** (issue #278). `docs/internal/VALIDATION.md` gains an
  *Interlaminar (CZM) validation evidence — Phase 7* section embedding all
  seven NASA-TM benchmark comparisons (DCB, ENF, 4PB, MMB 25/50/75, mixity
  synthesis) beside the test that regenerates each one, plus the two legacy
  knockdown snapshots marked explicitly as archival (no script in the tree
  regenerates them). Previously the directory was the repository's largest
  content while no document referenced any of it.

### Numerical results
- FE `stress_local` / `strain_local` — and therefore every failure index,
  failure mode and retention factor derived from them — change for
  laminates with **off-axis plies and a non-zero wrinkle angle**, per the
  transform fix above. UD `[0]_n` results and all pristine elements are
  **bit-identical**, and `python scripts/validate.py` shows zero ledger
  drift. On a `[0/45/-45/90]s` coupon at A = 0.3 mm the LaRC05 maximum
  failure index moves 0.8506 -> 0.8459 (-0.55 %) and the retention factor
  0.84075 -> 0.84089: small at the maximum, because that maximum sits in
  the fibre-dominated 0 deg plies which were never affected, but up to
  8 % per element in the 90 deg plies where matrix modes are assessed.
- FE stresses, failure indices and retention factors change when
  `delta_T != 0` — previously such a run was refused outright, so nothing
  that ran before produces a different number now (`delta_T == 0` is
  bit-identical, pinned by a test). Measured on a wrinkled `[0/90]s`
  IM7/8552 coupon, a ΔT = −155 cool-down adds **+35.3 MPa** of matrix
  (σ₂) tension, and adds the *same* amount under tension and compression
  — a load-independent residual, as it should be. The effect on failure is
  therefore **signed**: the LaRC05 maximum index goes 0.497 → 0.754 in
  tension (residual adds to the mechanical matrix tension) and
  0.348 → 0.214 in compression (residual relieves it). Cure residual
  stress is not a blanket penalty, and a model that reported one would be
  wrong for the compression cases this package mostly analyses.
- Two deliberate asymmetries in what carries the temperature: the pristine
  **retention baseline** is solved at the same `delta_T`, so retention
  factors compare like with like; the **measured global modulus** is
  solved at `delta_T = 0`, because a residual load adds a
  strain-*independent* offset to the reaction force and reporting that as
  a stiffness change would be a wrong number. (The legacy local-σ₁₁
  modulus proxy is derived from the single thermally-loaded field, so it
  inherits the thermal state on both sides of its ratio.)

### Changed
- Tests — **generated validation figures no longer land in a git-tracked
  path** (issue #278). The seven Phase-7 tests resolve their output through
  `tests/integration/_figure_output.validation_figure_path()`, which
  defaults to the git-ignored `figures/_generated/` and honours
  `WRINKLEFE_FIGURE_DIR`; refreshing the committed evidence is now the
  explicit `WRINKLEFE_FIGURE_DIR=figures pytest tests/integration -m
  integration`. Running the suite no longer dirties the working tree.

### Removed
- `figures/fig_page4.png` and `figures/fig_page5.png` (issue #278) — 16.2 MB
  of the directory's 18 MB, referenced by nothing. They were not PNGs: no
  PNG signature, exactly 8,400,000 bytes each (2000x1400 raw RGB), so they
  could not have rendered in GitHub, Markdown or Sphinx even if linked. The
  directory drops from 18 MB to 1.5 MB; the blobs remain in git history.
- Analysis — **thermal / cure-residual loading is reachable from
  `AnalysisConfig`** (issue #273, Stage 1 — the FE initial-strain term is
  Stage 2). The CLT machinery (`LoadState.delta_T`,
  `Laminate.thermal_resultants`, the thermal branch of
  `midplane_strains`) had been implemented and tested for a long time but
  dead-ended one layer below the pipeline, which built a hardcoded
  `LoadState(Nx=applied_strain * 1000.0)`. There is now an
  `AnalysisConfig.delta_T` (default `0.0`), a CLI `analyze --delta-T`
  flag, and an expert-mode *Temperature change from cure ΔT [°C]* input
  in the Streamlit app, so a user can finally ask "what does my knockdown
  look like including cure-induced residual stress?".
  **Sign convention, stated in the field docstring, the CLI help, the app
  tooltip, the run summary, the README and the units page:** `delta_T` is
  the temperature change *from the stress-free (cure) state*
  (`T_service − T_stress_free`), so a **cure cool-down is negative** — a
  177 °C cure taken to 22 °C service is `delta_T = -155`.
  **Scope boundary (lifted later in this same unreleased cycle by Stage
  2, above).** Stage 1 covered the analytical/CLT path only: with no FE
  thermal initial-strain load vector, a non-zero `delta_T` with
  `analytical_only=False` was **rejected at construction and again at
  `run()` time** rather than silently dropped — an FE result that quietly
  omits a load this large is a wrong number, not a missing feature. On
  the analytical path a CLT first-ply-failure report is now produced when
  `delta_T != 0` (the closed-form knockdown carries no temperature term,
  so without it the feature would have no output at all); `delta_T == 0`
  leaves that path exactly as it was, `failure_report is None` included.
  Moisture (`delta_C`, `beta1/2/3`) is deliberately **deferred**: nothing
  in the CLT solve consumes `delta_C`, so exposing it would create the
  silent no-op this change removes.
- CI — **Python 3.13 is now a required test cell, and 3.14 runs as a
  contained experiment** (issue #279 — Fixes #279). The matrix stopped
  at 3.12, so a package whose users increasingly run a 2026 distro
  default had zero signal on the two newest interpreters. 3.13 joins the
  required matrix on both Ubuntu and macOS and gains its trove
  classifier; the release workflow's wheel smoke test picks it up too, so
  a published wheel is exercised on every version the classifiers claim.
  3.14 runs in a separate non-blocking `test-experimental` job — Ubuntu
  only, `continue-on-error`, installed without the `vtk` extra, because
  VTK wheels are the expected gap on a just-released CPython and pyvista
  is lazily imported behind `pytest.importorskip`. It deliberately gets
  **no** trove classifier while its cell is allowed to fail; the matrix
  cell and the classifier get promoted together once VTK ships wheels.
  `requires-python` is unchanged: 3.10 stays the floor.
- Docs — a **Python support policy** in `CONTRIBUTING.md` (issue #279):
  every CPython in bugfix/security status is a required matrix cell and a
  classifier, the newest release may run as a non-blocking experiment
  until the scientific stack's wheels stabilise, no version gets a
  classifier while its cell can fail, and the whole thing is reviewed
  each October.
- Release automation — **tagging is now the whole release procedure**
  (issue #264 — Fixes #264). Pushing a `vX.Y.Z` tag runs
  `.github/workflows/release.yml`: build → version-check and wheel-test →
  publish → GitHub Release.
  - The **version check compares the tag against the built wheel's
    metadata**, not against a grep of `pyproject.toml` — it installs
    `dist/*.whl` into a throwaway venv and reads the version back through
    `importlib.metadata`. That proves what the artifact about to be
    uploaded declares, closing the issue #21 class of skew (a published
    artifact disagreeing with the repository) permanently.
  - The **smoke test exercises the installed wheel**, not the source
    tree, on Python 3.10/3.11/3.12: the checkout supplies only the tests,
    the package comes from a clean-venv wheel install, and a step asserts
    `wrinklefe.__file__` resolves inside site-packages before the tests
    run.
  - **PyPI publishing uses trusted publishing** (OIDC, `id-token: write`,
    a `pypi` environment) — there is no API token in repository secrets.
    The one-time trusted-publisher registration and the Zenodo archiving
    toggle are maintainer actions on those sites; their exact fields are
    documented under "Release procedure" in `CONTRIBUTING.md`.
- CI — a **`build` job on every pull request** running `python -m build`
  and `twine check dist/*` (issue #264), so packaging breakage (a bad
  classifier, a dropped `py.typed`, a MANIFEST.in typo) fails a PR instead
  of surfacing by hand at upload time.
- CI — a **`citation` job validating `CITATION.cff`** against the CFF
  1.2.0 schema with `cffconvert --validate` (issue #284). A malformed
  citation file fails silently today: GitHub's "Cite this repository"
  button just stops offering BibTeX.
- Tests — `tests/test_version.py` now also locks **`CITATION.cff`'s
  `version:` to `pyproject.toml`** (issue #284), so a version bump that
  forgets the citation file fails CI rather than advertising a version
  that was never released.
- Docs — a **"Release procedure" section in `CONTRIBUTING.md`** (issues
  #264, #284): the bump/changelog/tag sequence, the one-time PyPI
  trusted-publisher fields and Zenodo toggle, where a minted DOI gets
  inserted, and how to rehearse an upload against TestPyPI.
- Docs — **interpretation guidance, a units/conventions reference and a
  worked tutorial** (issue #378 — Fixes #378). Three new pages, wired into
  the Sphinx toctree after *Getting started*:
  - **Interpreting results** (`docs/interpreting_results.md`) — what each
    headline number is and, just as importantly, is not:
    `analytical_knockdown` (and why `analytical_onset_knockdown` is the
    tension first-load-drop), the three stiffness numbers
    (`analytical_modulus_knockdown`, the local-σ₁₁ `modulus_retention`, and
    the coupon-level `modulus_retention_global` that should be preferred),
    the FE first-ply-failure `retention_factors` and why they are degenerate
    for pristine UD compression, `progressive_knockdown` as the ultimate,
    the CZM outputs with `czm_converged` as the precondition for quoting any
    of them, and the safe-side, forward-verified semantics of the goal-seek
    acceptance limit (issue #280). Closes with the severity bands
    **transcribed from `_SEVERITY_BANDS` in `wrinklefe.io.export`**, which
    remains authoritative, carrying the export's existing non-binding MRB
    language verbatim rather than any new wording.
  - **Units & conventions** (`docs/units_conventions.md`) — one sourced
    table (mm / MPa; toughness in N/mm, numerically kJ/m²; strain and `Vf`
    as fractions; angles in degrees but misalignment stored on
    `AnalysisResults` in radians), the x/y/z and material-frame coordinate
    conventions, the sign conventions (compression-negative
    `applied_strain`, unsigned strength allowables, `delta_T` from the
    stress-free state), and the one place percent and fraction differ — the
    Streamlit app's *Applied strain magnitude [%]* input.
  - **Tutorial** (`docs/tutorial.md`) — a worked walkthrough of the real
    workflow: measure *A*, *λ*, *w* off a micrograph → configure (CLI, app
    and Python side by side, with `--save-config`) → `wrinklefe analyze` →
    read the numbers → `wrinklefe critical` for the acceptance limit →
    `build_analysis_summary` / `export_summary` for the NCR attachment.
    Every command was run and every output block is real. This is the
    issue's "notebook **or equivalent worked walkthrough**": a Markdown page
    that `sphinx-build -W` builds and CI therefore keeps honest, rather than
    a `.ipynb` with no execution infrastructure behind it.
- Examples — **the four headline capabilities that had none** (issue #378 —
  Fixes #378): `11_penetration_gate.py` (~1 s) holds `theta_max` fixed while
  varying the laminate thickness so the angle-only model is flat by
  construction and every difference is the gate, then sweeps the position
  factor and runs one coupon gated vs ungated; `12_progressive_damage.py`
  (~13 s) prints the load-increment table with its peak and post-peak load
  drop and contrasts the ultimate-strength knockdown against the degenerate
  UD first-ply-failure retention; `13_crest_resin_pocket.py` (~20 s)
  compares the crest lens off / graded / binary and sweeps its height; and
  `14_stochastic_knockdown.py` (~14 s) propagates two measurement
  distributions to percentile knockdowns with a fixed seed, an LHS-vs-Monte
  Carlo check and a rank-correlation sensitivity screen. All four run on the
  minimal install the examples CI job uses (`pip install -e .`, no extras);
  `examples/README.md` gains their rows and every existing row was
  re-verified against the directory.
- App — **live progress for FE/CZM runs** (issue #377 — Fixes #377). The
  *Run analysis* status box now carries a real progress bar driven by the
  engine's existing `WrinkleAnalysis.run(progress_callback=...)` hook, so a
  long solve reports the phase it is in rather than sitting at a frozen
  10 %: *Building laminate* (0 %) → *Computing analytical predictions*
  (5 %) → *Assembling FE mesh* (10 %) → *Solving FE system* (25 %) →
  *Evaluating failure criteria* (75 %) → *Computing retention factors*
  (90 %) → *Analysis complete* (100 %), with the percentage in the bar text
  and the current phase echoed in the status label. Users on a 30–90 s
  Streamlit Cloud CZM run can now tell a slow solve from a hung one. The
  callback is defensive by construction — fractions are clamped and any
  widget error is swallowed and logged at debug — so a cosmetic update can
  never abort a solve that is already minutes in. **Coarse during CZM and
  progressive-damage solves:** `AnalysisResults`' CZM / progressive
  sub-paths do not emit progress of their own, so those runs hold at
  *Solving FE system* (25 %) for the whole Newton-Raphson / load-increment
  loop before completing; finer emits inside those loops would be an engine
  change and are left as a follow-up.
- Physics — **compaction-driven ply-thickness / local fibre-volume-fraction
  gradient** (issue #379, Part B — Fixes #379, completing the issue; builds
  on the Part A micromechanics). New `wrinklefe.core.compaction` derives a
  per-element local `Vf` from the deformed mesh with the kinematic rule
  `Vf_local = vf_nominal · h0/h` — fibre content conserved per element while
  the thickness change absorbs or expels resin — so a wrinkle constrained by
  rigid tooling thins and thickens its plies instead of being modelled at a
  constant ply thickness. Stretched trough elements turn resin-rich (low
  `Vf`, softer); compacted crest elements turn resin-starved (high `Vf`,
  stiffer), which is the treatment #371 left as a follow-up. Local materials
  are **ratio-anchored** on the measured preset,
  `P_local = P_preset · P_micro(Vf_local) / P_micro(vf_nominal)`, applied to
  the stiffnesses and CTEs, so the micromechanics model contributes only its
  `Vf` sensitivity and none of its 12–33 % absolute error; at
  `Vf_local == vf_nominal` every ratio is exactly 1.0 and the preset object
  itself is used. **Poisson ratios and all strengths stay at the preset
  values** — no mixing rule maps `Vf` to an allowable, and inventing one
  would be quietly wrong (local failure indices still move, because the
  local stiffness redistributes stress). `Vf` is quantized onto an
  `n_bins`-value grid anchored on the nominal value, so a mesh shares a few
  dozen material objects rather than one per element. Enabled with
  `AnalysisConfig(enable_vf_gradient=True)` (plus `vf_nominal`, `vf_fiber`,
  `vf_matrix`, `vf_max`) or `wrinklefe analyze --vf-gradient`; **opt-in and
  off by default**, FE-only, and restricted in v1 to
  `morphology="tool_flat"`, whose flat outer envelope is what keeps the
  per-column thickness — and hence the resin mass — conserved
  (`surface_pocket_side="both"` is the two-caul-plate case). Elements
  compacted past `vf_max` (default 0.75) saturate with a single counted
  warning: the rule carries no lateral resin flow along the ply. See
  `examples/10_vf_gradient_compaction.py`.
- Core — **constituent-based micromechanics: fibre + matrix + fibre volume
  fraction to a ply** (issue #379, Part A — the Vf-to-properties capability
  #379 names as its own blocker; the compaction kinematics that consume it
  are a separate change). New `wrinklefe.core.micromechanics` supplies
  `FiberProperties` (transversely isotropic) and `MatrixProperties`
  (isotropic, with `from_material` so an existing neat-resin card such as
  `EPOXY_S6C10` can be reused rather than duplicated), the mixing rules —
  Voigt rule of mixtures for `E1`/`nu12`/`nu23`, Halpin–Tsai for `E2`
  (ξ = 2) and `G12` (ξ = 1), transverse isotropy for `G23`, Schapery for
  the thermal expansion coefficients — and cited constituent presets:
  `FIBER_PRESETS` for the nine fibres behind the library systems (AS4,
  T300, T700S, IM7, IM10, IM6G, T800S, S-2 glass, Kevlar 49) and
  `MATRIX_PRESETS` for the two resins with published neat-resin data
  (3501-6, 8552). `OrthotropicMaterial.from_constituents(fiber, matrix,
  Vf)` builds the ply and validates it like any other card. **Strengths are
  not predicted** — no mixing rule here maps Vf to an allowable, so
  strengths, toughnesses and cohesive tractions are carried over unchanged
  from `strengths_from=` (or left at the defaults) and do not track Vf.
  Accuracy is stated from measurement, not aspiration: rebuilt at the Vf
  each preset documents in its own comment, the model lands within 12 % on
  `E1`, 26 % on `nu12`, 33 % on `E2`, 32 % on `G12` — except Kevlar-49 /
  epoxy `G12`, over-predicted by 84 % and recorded as an `xfail` rather
  than tuned away — and `nu23` is under-predicted throughout.
  `tests/test_micromechanics.py` pins every deviation. **No behaviour
  change:** nothing in the analysis pipeline consumes the module yet, no
  existing preset value moved, and the validation ledger is unaffected.
- App — **"Find acceptable limit" mode and the acceptance limit on the NCR
  attachment** (issue #280, app slice — *Fixes #280*, completing the issue).
  The Streamlit sidebar gains a **Find acceptable limit** button directly
  under *Run analysis* (it searches the same inputs) plus an
  *Acceptable-limit settings* expander: the searched parameter (amplitude
  or wavelength), the target as either a knockdown factor or an absolute
  MPa allowable, and — in Expert mode — the bracket, scan-point count and
  root tolerance. The cost is stated before the click: the search runs on
  the analytical path (~25 evaluations, well under a second) unless the
  config enables an FE-only feature, in which case the expander warns that
  it will spend ~25 FE solves. A converged search renders the limit, the
  objective achieved at it and the evaluation count, and says plainly that
  the reported value is the **conservative** one — backed off from the raw
  `brentq` root and verified by an extra forward run — with the raw root
  quoted only for transparency. A refusal renders the engine's own message
  **verbatim** (`no_crossing` as a warning, since it can simply mean
  nothing in range fails the criterion; `non_monotonic` / `flat` /
  `target_unreachable` as errors), and the scanned curve and evaluation
  ledger are shown either way because the curve is the diagnosis. **Apply
  this limit to the sidebar** seeds the geometry widget through the #375
  pending-seed mechanism, and the stored search is flagged stale by the
  same payload hash-compare the run results use (#374).
  `io.export.build_analysis_summary` gains an optional `critical_limit=`
  mapping (additive; existing callers unaffected) that adds a
  `critical_limit` block to the summary and an *Acceptance limit
  (goal-seek)* sub-block to the Markdown and PDF attachments, carrying the
  same non-binding MRB language as the rest of the summary. The app
  attaches it **only** when the stored search converged *and* its config
  payload equals the payload the displayed results were computed from — an
  acceptance limit derived for different inputs must never ride along on an
  NCR — and captions the omission when it does not.
- Analysis — **inverse goal-seek: maximum acceptable wrinkle amplitude**
  (issue #280). New `wrinklefe.goalseek.find_critical_value(base_config, *,
  parameter='amplitude', target_knockdown=…, …)` scans the resolved search
  range, brackets the single sign change of `objective(parameter) - target`,
  root-finds with `scipy.optimize.brentq`, and then **backs the answer off
  to the safe side** so the returned value satisfies the criterion under a
  real forward evaluation rather than merely to within the root tolerance.
  Returns a `CriticalValueResult` with the critical value, the achieved
  knockdown and strength, a ready-to-run `critical_config`, the full scan
  curve and the evaluation ledger (`GoalSeekEvaluation` rows, with
  `summary()` / `plot()`). The upper bound is derived from the laminate
  thickness for amplitude and clamped by the tool_flat validation bound and
  by `mesh_shear_diagnostics`' `amplitude_safe` on the FE path; search
  direction and monotonicity are *measured* from a log-spaced scan, never
  assumed; and the search refuses — with an actionable,
  measurement-quoting message, never a scipy traceback — when the target is
  never crossed (`no_crossing`), never met (`target_unreachable`), reached
  by more than one root (`non_monotonic`, e.g. graded morphology vs
  wavelength, whose measured curve is U-shaped with an interior minimum
  near λ ≈ 3 mm), or unresolvable because the parameter is inert for the
  configuration (`flat`). FE-only config features (CZM, resin pockets,
  progressive damage, non-uniform transverse mode) force or refuse the FE
  path rather than silently no-opping, and integer mesh/ply-count fields
  are refused by name. Targets may be a knockdown factor or an absolute
  strength in MPa. New `wrinklefe critical` subcommand
  (`--parameter`/`--target-knockdown`/`--target-strength`/`--objective`/
  `--bracket`/`--max-value`/`--scan-points`/`--rtol`/`--analytical-only`/
  `--config`/`--save-config`/`--save-plot`/`--output-json`/`--output-csv`
  plus the `converge`-style geometry flags), defaulting to the analytical
  path (~25 evaluations, well under a second); exit 2 for bad input, 1 when
  no acceptable limit exists in range. New `examples/09_acceptance_limit.py`.
  The forward model is untouched (validation-ledger zero drift).
- Solver — **iterative-solver controls reachable from `AnalysisConfig`**
  (issue #265). The CG/ILU knobs that were hardcoded in
  `StaticSolver._solve_iterative` are now config fields, plumbed through
  `StaticSolver` and reachable via `--config`: `iterative_rtol` (default
  `1e-10`), `iterative_maxiter` (`10000`), `ilu_drop_tol` (`1e-4`),
  `ilu_fill_factor` (`None` → SciPy default), and `preconditioner` ∈
  {`ilu`, `jacobi`, `none`} (`ilu`) so users can pick the low-memory
  diagonal preconditioner (or none) on huge meshes. Defaults reproduce
  the previous hardcoded values bit-for-bit, so an existing iterative
  solve is unchanged (validation-ledger zero drift). Round-trip through
  `to_dict`/`from_dict`.
- Solver — **Newton/CZM convergence-failure diagnostics + actionable hint**
  (issue #262). A failed nonlinear solve used to return a bare
  `converged: False`. `NewtonRaphsonSolver.solve()` now also returns
  `failure_diagnostics` (first failing increment: index + load fraction,
  iteration count, final `||R_phys||` / BC violation / `||du||`, the tail
  of the residual history, line-search status, tangent-singular flag, and a
  classified `failure_reason` ∈ {`tangent_singular`, `diverged`,
  `stagnated`, `iteration_cap`}) and `failure_hint`, a single string naming
  the knob to turn (`czm_n_load_increments`, `czm_newton_tol`, the applied
  strain, `max_newton_iter`, or the arc-length roadmap). The `_newton_step`
  tuple return is unchanged and the numerical path is untouched, so
  converged results are bit-identical. `AnalysisResults` gains
  `czm_failure_diagnostics` / `czm_failure_hint`; the CLI prints the hint to
  stderr and the Streamlit CZM section shows it as an error on a
  non-converged `--enable-czm` run.
- App — **config upload/download and through-width transverse controls**
  (issue #375, app slice — *Fixes #375*, completing the issue). The
  Streamlit sidebar gains a **Config file** section: **Download config
  (JSON)** serialises the current effective `AnalysisConfig` from the live
  sidebar state (works before any run) and round-trips with the CLI
  `--config` / `--save-config` flags, and **Load config (JSON/YAML)**
  reads a saved case back into the sidebar. Loading parses via
  `AnalysisConfig.from_dict` (bad files surface an `st.error`, never a
  crash) and *seeds* the widgets from the config — the seed is staged in
  `session_state` and applied at the top of the sidebar before the widgets
  instantiate (then a rerun), so Streamlit's set-after-instantiate error is
  avoided; custom materials route through the custom-editor keys and stale
  results are cleared. Expert mode also exposes the through-width
  **transverse** envelope (`transverse_mode` selectbox plus span/width
  inputs); a non-uniform mode forces the FE path and is threaded into the
  run config, while the incompatible transverse+CZM combo shows a sidebar
  warning and is not threaded into an invalid config.
- CLI — **config-first sweeps/compares, transverse + stochastic exposure,
  and ergonomics** (issue #375, CLI slice). `sweep` and `compare` gain
  `--config PATH`: the file supplies the base `AnalysisConfig` (laminate,
  material, mesh, penetration gate) so a **UD amplitude sweep through the
  penetration gate** — previously impossible from the CLI — is reachable
  with `wrinklefe sweep --config ud_gate.json --parameter amplitude ...`;
  explicitly-passed geometry flags override the file via the #259
  SUPPRESS-default precedence. `analyze` exposes the through-width
  transverse surface (`--transverse-mode {uniform,gaussian_decay,
  sinusoidal_y,elliptical}`, `--transverse-span`, `--transverse-width`; a
  non-uniform mode forces the FE path) and the previously config-only
  `--nz-per-ply`, `--ply-thickness` (sets the gate D/T) and `--output-csv`.
  A new **`wrinklefe stochastic`** subcommand wraps
  `stochastic.probabilistic_analysis`: a `--config` base plus repeatable
  `--distribution FIELD:DIST:P1:P2` specs (`normal`/`uniform`/`lognormal`),
  `--n-samples`/`--seed`/`--method`, printing percentile knockdowns and
  writing JSON/CSV. `wrinklefe --version` now reads the installed package
  metadata instead of a hardcoded string. (App config upload/download and
  the app transverse controls are a separate follow-up; this slice does not
  complete #375.)
- App / CLI — **`tool_flat` surface-pocket controls live in the morphology
  definition** (issue #371, Part B — *Fixes #371*, completing the issue).
  The Streamlit Morphology selectbox gains **`tool_flat`** (Expert mode)
  with its own schematic cartoon (flat pinned face, uniform-amplitude core,
  amber resin wedges at the troughs). Selecting it renders the pinned-side
  and **surface-transition-plies** controls *directly under* the morphology
  controls — no longer in the FE expert section — alongside a live
  inversion-bound caption (max safe amplitude = `0.8 · S · t / nz`) and a
  pre-run warning naming both remedies, so the config `ValueError` is never
  the user's first feedback. Surface pockets **auto-enable** for `tool_flat`
  (no checkbox). The Analyze-tab cross-section renders the *thick* pockets
  via the actual `tool_flat` decay (pinned plies flat, uniform core); a
  fidelity test binds the analytic preview gap area to the multi-layer FE
  tagged volume (~10 %). The CLI accepts `--morphology tool-flat` (alias for
  `tool_flat`) and `--surface-transition-plies N`, with help noting the
  pockets auto-enable and the amplitude bound.
- Analysis — **`tool_flat` morphology with significant surface resin
  pockets** (issue #371, Part A). A new through-thickness decay mode /
  morphology (`morphology="tool_flat"`) models a wrinkle cured against
  rigid tooling: a uniform-amplitude core, a short linear transition over
  `surface_transition_plies` plies (new config field, default 2), and an
  **exactly-flat pinned surface** on `surface_pocket_side`
  (`"top"`/`"bottom"`/`"both"`). This fixes the root cause of the reported
  bug that surface pockets were invisible and mechanically negligible: the
  linear-decay morphologies (`stack`/`convex`/`concave`, `graded` with
  `decay_floor=0`) spread the wave across the whole thickness, so at the
  24-ply defaults (t=0.183 mm, A=0.5 mm) the outermost undulating ply moved
  only ~0.045 mm — a trough gap of **~0.25 of one ply thickness**. Under
  `tool_flat` the mismatch collects at the flat surface, so the trough
  pocket is ≈ the full amplitude (**~2.7 ply thicknesses** at defaults).
  Surface pockets **auto-enable** for `tool_flat` (they are its defining
  physics; skipped only for `analytical_only`), and the analytical path
  equals `uniform` (M_f = 1.0; the pocket effect is FE-only). Toggling the
  pockets now moves `modulus_retention_global` by a clearly significant
  margin (measured **~3.1–3.4 %** at A=0.5–0.55 mm, `surface_transition_plies=4`,
  side `both`) versus a negligible **~0.3 %** for the legacy `stack`
  morphology — a ~10× larger effect, the original complaint resolved.
  `compute_surface_resin_blend` tags **all** stretched layers in the
  multi-ply transition zone (volume-conserving), and its height metric is
  now the tilt-invariant vertical stretch (top-face minus bottom-face
  mean-z) so a realistically-localized wrinkle's in-plane slope no longer
  corrupts the gap. A verified element-inversion bound (`amplitude ≤ 0.8 ·
  surface_transition_plies · ply_thickness / nz_per_ply`) is enforced at
  construction with a message naming both remedies; multi-wrinkle / CZM /
  transverse combinations raise `NotImplementedError`. App/UX exposure of
  the new morphology is a deliberate follow-up (Part B).
- Analysis — **through-width (transverse) wrinkle surfaces reachable from
  `AnalysisConfig`** (issue #300). The already-implemented, already-tested
  `WrinkleSurface3D` transverse modes are now selectable through three new
  config fields: `transverse_mode`
  (`"uniform"`/`"gaussian_decay"`/`"sinusoidal_y"`/`"elliptical"`, default
  `"uniform"`), `transverse_span` (→ `span_y`, `None` tracks
  `domain_width`), and `transverse_width` (→ `width_y`, `None` resolves to
  `span_y / 4` — a localized mid-width patch). With the default
  `"uniform"` the pipeline still builds the bare x-only
  `GaussianSinusoidal`, so results are bit-identical (regression-safe). A
  non-uniform mode wraps the profile in a `WrinkleSurface3D` on the FE
  single-wrinkle path so the crest amplitude varies across the specimen
  width; at the same crest amplitude a localized wrinkle predicts a milder
  knockdown than the uniform baseline (real manufacturing wrinkles are
  localized, and the uniform assumption overstates the defect volume).
  FE-only and single-wrinkle for now: analytical-only, multi-wrinkle
  (`wrinkles`), and `enable_czm` combinations are rejected at construction
  with actionable messages. New `examples/transverse_wrinkle_knockdown.py`
  demonstrates localized-vs-uniform knockdown; CLI/app exposure is a
  deliberate follow-up.
- CLI — **wrinkle-defect capabilities on `analyze`** (issue #346). The
  new defect models shipped for the scripting API are now reachable from
  the command line: `--wrinkle-z-position Z` (off-mid-plane wrinkle,
  validated to `[0, 1]`), `--gate {li2024-moulded,li2025-vacbag}` (the
  two-parameter (θ, D/T) penetration gate, selecting a calibrated
  `GateParameters` preset), `--resin-pocket` (crest resin lens),
  `--surface-resin-pockets` / `--surface-pocket-side {top,bottom,both}`
  (tool-flat surface pockets), and `--progressive` / `--increments N`
  (load-stepping ultimate strength). The flags inherit the #259
  config-file precedence, so any flag left off keeps the `--config`
  value and any flag passed overrides it (and is written by
  `--save-config`); the FE-only features force the FE path with the same
  precedence as `--enable-czm` rather than silently no-op'ing under
  `--analytical-only`. The result summary now prints the progressive
  knockdown when a progressive run happened. `sweep --parameter
  wrinkle_z_position` works end-to-end over a `--config` base. The
  surface-resin-pocket flags (issue #361, newer than #346) are included
  as part of the same CLI-reachability story. Buckling was **deferred**:
  the linearized microbuckling solver is standalone diagnostic
  infrastructure with no `AnalysisConfig` knob and is documented as not a
  usable knockdown predictor, so there is no clean field to surface.
- Streamlit app — **surface resin pockets in the through-thickness
  cross-section** (issue #361, Part 4 follow-up). The Analyze-tab
  cross-section now shades the neat-resin pockets that fill the wrinkle
  troughs under a tool-flat surface, and an expert-mode sidebar toggle
  (*Surface resin pockets*, with a top/bottom/both side selector) drives
  both the preview and the FE run so a solve models exactly what the
  picture shows. The zone is rendered *analytically* — the amber fill
  between the flat tool line and the deformed outermost undulating ply —
  rather than by deforming an FE mesh at render time, so the preview
  stays responsive; a fidelity test cross-checks that rendered gap area
  against the resin volume `compute_surface_resin_blend` tags on a coarse
  mesh (agree within ~4%, well inside #361's 10% conservation tolerance).
  Shown only for tool-flat morphologies (`stack`/`convex`/`concave`, or
  `graded` with `decay_floor=0`); an incompatible morphology (`uniform`,
  or `graded` with a non-zero floor) is withheld from the config and
  flagged with a sidebar note, so a run can never build an invalid
  config. Opt-in and off by default (feature-off preview unchanged).
- Tool-flat surfaces with surface resin pockets (issue #361).
  Parts cured against rigid tooling / a caul sheet keep perfectly flat
  outer surfaces while the fibres undulate internally; the wrinkle
  troughs fill with neat resin just under the flat surface. New
  `SurfacePocketSpec` / `compute_surface_resin_blend` in
  `wrinklefe.core.resin_pocket` tag, per column, the transition element
  that stretches to span the gap between the flat surface and the
  outermost undulating ply, weighting it by the excess-stretch fraction
  `max(0, (h - h0) / h)` — exactly volume-conserving (equal to the
  integrated kinematic gap `-w(x)·decay_last`). Enabled via
  `AnalysisConfig.enable_surface_resin_pockets`, `surface_pocket_side`
  (`top`/`bottom`/`both`) and `surface_pocket_min_gap`, reusing
  `resin_pocket_material` / `resin_pocket_graded`. FE-only effect
  (`modulus_retention_global` and first-ply failure in the isotropic
  resin zone); it composes with the crest lens (per-element maximum) and
  needs no solver changes. Requires a tool-flat morphology whose decay
  reaches 0 at the chosen surface (`stack`/`convex`/`concave`, or
  `graded` with `decay_floor=0`); `uniform` and `graded` with a non-zero
  floor are rejected with a message naming the fix. Disabled by default
  (bit-identical results when off).
- Save / load an `AnalysisConfig` (issue #259).
  `AnalysisConfig.to_dict()` / `from_dict()` provide a round-trippable,
  `config_version`-stamped serialisation; `save_json` / `load_json`
  (and extension-dispatching `save` / `load`, plus optional YAML when
  PyYAML is installed) read and write config files. Library materials
  serialise by preset name, custom materials inline, and penetration-gate
  presets by their registry name (new
  `wrinklefe.core.penetration_gate.GATE_PRESETS`). Loading rejects
  unknown keys and version mismatches loudly. The `wrinklefe analyze`
  CLI gains `--config PATH` (load a config, with explicitly-passed flags
  overriding the file) and `--save-config PATH` (write the effective
  config). A follow-up will surface the same config download/upload in
  the Streamlit app.
- Streamlit app — **Through-thickness cross-section** on the Configure
  tab: a new panel that draws the deformed ply stack in the (x, z) plane
  so users can see how the wrinkle manifests through the laminate
  thickness. It reuses the real `WrinkleConfiguration.apply_to_nodes`
  field the FE mesh uses, so the picture faithfully tracks the active
  morphology, its through-thickness amplitude decay, the dual-wrinkle
  phase offset, and any in-plane amplitude profile. Each ply is a band
  coloured by fibre angle (the same hue map as the layup visualizer),
  with the wrinkle-interface plies outlined for the dual morphologies.
- Monte-Carlo / Latin-hypercube uncertainty propagation (issue #301):
  `wrinklefe.stochastic.probabilistic_analysis(base_config,
  distributions, n_samples, seed, method="lhs"|"mc")` samples
  `AnalysisConfig` fields from user distributions (`("normal", m, s)`,
  `("uniform", lo, hi)`, `("lognormal", mu, sigma)`, or any frozen
  `scipy.stats` distribution), runs the analytical path per sample, and
  returns a `ProbabilisticResults` with percentile
  knockdowns/strengths, mean ± std, the input samples for sensitivity
  scatter, an optional histogram+scatter `plot()`, and a `summary()`
  that explicitly labels the output as model-input-propagation
  statistics — **not** CMH-17 A-/B-basis allowables. Fixed seeds are
  fully reproducible; degenerate (zero-variance) distributions
  reproduce the deterministic result exactly; invalid draws fail loudly
  instead of being clipped; `n_workers` reuses the #260 process pool
  for FE-path sampling. 1000 analytical samples run in ~0.7 s for
  UD/gate configs.
- Vectorized `_laminate_modulus_knockdown` (issue #301 enabler): the
  multidirectional analytical modulus knockdown ran a per-(ply,
  x-station) Python loop (12,000 6×6 rotations/condensations), making
  every multidirectional analytical run take ~1.2 s. The wrinkle-tilt
  rotation and plane-stress condensation are now batched over the whole
  (ply, x) grid (the in-plane rotation is computed once per ply, and
  `T_sigma(θ)^-1 = T_sigma(−θ)` replaces the batched solve) — 1.18 s →
  22 ms (~52×) per analytical run, results identical to the loop
  (regression-tested against it; ledger baselines zero-drift).
- Mesh-resolution warning (issue #306): `WrinkleMesh.generate` warns
  when the hex mesh samples the wrinkle wavelength with fewer than 4
  elements (element `dx` vs `lambda`), naming the offending spacing and
  the `nx` needed — under-sampled wrinkles previously produced silent
  aliasing noise.
- CZM-capable glass/aramid presets (issue #268): `S2_GLASS_EPOXY` and
  `KEVLAR49_EPOXY` now carry representative interlaminar toughness
  (`GIc`/`GIIc` in the published glass/aramid ranges), so
  `enable_czm=True` runs for every built-in material instead of raising
  for those two. Configurations that previously errored now produce
  cohesive-damage results.
- Coordinate-aware maxima (issue #297):
  `FieldResults.max_displacement_location()` and
  `max_stress_location()` return the physical `(x, y, z)` of the
  governing node / element centroid alongside the value.
- Failure-mode breakdown plot (issue #269, part 1):
  `wrinklefe.viz.plot_failure_mode_breakdown` with a stable
  per-failure-mode colour map (`MODE_COLORS`).
- Process-parallel parametric sweeps (issue #260):
  `WrinkleAnalysis.parametric_sweep(..., n_workers=N)`,
  `wrinklefe.sweep.run_sweep(..., n_workers=N)`, and
  `wrinklefe sweep --parallel N` fan the independent per-point solves
  out over a `ProcessPoolExecutor` (`N=0` uses all CPU cores; the
  default `N=1` keeps the exact sequential path). Results are identical
  to and ordered like the sequential run — measured 3.6× at 4 workers
  on an 8-value full-FE amplitude sweep. `run_sweep` progress becomes
  completion-based in parallel mode; `KeyboardInterrupt` (or a worker
  failure) cancels the queued futures instead of draining the pool.
  Peak memory scales with workers × per-solve footprint — size `N` by
  available RAM for fine meshes.
- Vectorized `FieldResults.max_principal_stress` (issue #295): the
  per-Gauss-point Python double loop (one `np.linalg.eigvalsh` call per
  point) is replaced by a single batched eigen-solve on the
  `(n_elem, n_gp, 3, 3)` tensor stack — ~6× on a 50k-element × 8-GP
  field (4.2 s → 0.7 s on first access; results identical to 1e-10,
  regression-tested against the old loop kept as the test oracle).
  Element centroids are now computed once per `FieldResults` (lazy
  `element_centers` property) instead of rebuilt on every
  `stress_through_thickness` call (372 ms → 5 ms per query on the same
  mesh).
- Vectorized `evaluate_field` for LaRC05, Puck, and Budiansky–Fleck
  (issue #299) — the three most expensive criteria were the last ones
  running the base class's per-Gauss-point Python loop. The
  fracture-plane / action-plane searches broadcast over a
  `(N, n_theta)` grid processed in cache-sized row blocks; measured on
  an 80,000-point field: LaRC05 12.0 s → 0.51 s (23×), Puck 13.8 s →
  0.67 s (21×), kink-band 93 ms → 7 ms (13×). Failure post-processing
  no longer dwarfs the linear solve when these criteria are enabled
  (the progressive-damage crack-band loop, which evaluates
  MaxStress+LaRC05 every equilibrium iteration, inherits the speedup).
  Outputs are **bit-identical** to per-point `evaluate()` — enforced by
  a per-criterion equivalence suite using exact array equality across
  randomized samples covering every branch regime. The base-class loop
  fallback now logs at DEBUG so future criteria authors notice.
- Penetration-gate validation harness (issue #161): the calibrated UD
  gate — the only strength path sensitive to wrinkle amplitude and
  through-thickness position independently of the peak angle — is now
  pinned in the reproducible ledger. `scripts/validate.py` scores a
  per-case gate column (with drift detection and `--update` re-pinning)
  for any dataset naming a `penetration_gate` preset; the Li 2025
  dataset carries `expected_gate_kd` baselines and a per-case `z_frac`
  through-thickness position (S-A-2 rides the gate's `P(z)` factor to
  its measured near-surface KD). The issue's acceptance criteria are
  permanent regression tests: the S-M-2/4/5 amplitude trio (identical
  20° angle, measured KD 0.629/0.943/1.000) lands at +2.2 %/−0.6 %/
  −0.3 % (±15 % band), orderings asserted monotonic, all six cases
  within the ±20 % parity band. README and VALIDATION.md updated to
  document the in-repo reproducible UD validation.
- Cohesive-zone delamination in multi-wrinkle FE (issue #283):
  `enable_czm=True` now runs with an `AnalysisConfig.wrinkles` list
  instead of raising `NotImplementedError`. Cohesive layers are inserted
  along the **full length** of every nominated interface, and
  `czm_interfaces="near_crest"` nominates the interface nearest *each*
  wrinkle (deduplicated) — wrinkles sharing an interface index get one
  continuous cohesive surface, so a delamination initiating at one crest
  can propagate toward its neighbour (crest-to-crest link-up, the Li 2025
  multi-wrinkle failure pattern; see
  `examples/08_multi_wrinkle_czm_linkup.py`). Regression anchors: a
  one-entry `wrinkles` list reproduces the scalar-config CZM solution
  bit-tight; far-separated wrinkles match independent single-wrinkle
  solves within a few percent with an intact interface between them;
  scalar (named-morphology) CZM interface resolution is unchanged.
- Analytical stiffness (axial-modulus) knockdown on the analytical path:
  `AnalysisResults.analytical_modulus_knockdown`, a closed-form CLT
  series-average of the off-axis lamina modulus over the wrinkle profile
  (`analysis._profile_modulus_knockdown`). Populated for unidirectional
  layups (loading-independent, zero FE cost) — the closed-form companion
  to the FE `modulus_retention`, which previously was the only stiffness
  knockdown (the analytical path reported none). Surfaced in
  `AnalysisResults.summary()`, the `analyze --output-json` payload,
  `results_to_dict`, and the Streamlit app; validated by
  `validation/validate_modulus.py` (analytical MAE 3.9 % / 1.2 % on the
  Li 2025 / Hsiao & Daniel UD datasets).
- Resin-pocket material zone (`wrinklefe.core.resin_pocket`:
  `ResinPocketSpec`, `compute_resin_mask`, `compute_resin_blend`):
  a graded neat-epoxy lens at the wrinkle crest, tagged into the FE mesh
  via `AnalysisConfig.enable_resin_pocket` /
  `resin_pocket_graded` / `resin_pocket_material` /
  `resin_pocket_height_scale` / `resin_pocket_length_scale`. The modulus
  and fibre-misalignment angle blend together so the wrinkle defect is
  counted once. Adds `OrthotropicMaterial.isotropic()` / `.blend()` and
  an isotropic neat-epoxy card `EPOXY_S6C10`.
- Progressive-damage FE solver
  (`wrinklefe.solver.progressive_damage`: `ProgressiveDamageSolver`,
  `ProgressiveDamageResult`), enabled via
  `AnalysisConfig.enable_progressive_damage` with
  `progressive_n_increments` / `progressive_residual_factor` /
  `progressive_max_strain`. Load-steps to ultimate load with optional
  crack-band (Bažant–Oh) regularization — the first FE route to a real
  UD compression knockdown.
- Two-parameter (θ, D/T, z) penetration gate
  (`wrinklefe.core.penetration_gate`: `GateParameters`,
  `penetration_gate_kd`, `angle_floor`, `position_factor`,
  `predict_from_geometry`, `calibrate_gate`, plus presets
  `GATE_LI2024_MOULDED` and `GATE_LI2025_VACBAG`), wired through
  `AnalysisConfig.penetration_gate`. A closed-form UD predictor
  `KD = 1 − (1 − KD_angle(θ))·S(D/T)·P(z)` at zero FE cost.
- Linear buckling / geometric stiffness
  (`wrinklefe.solver.buckling`: `LinearBucklingSolver`,
  `BucklingResult`, `microbuckling_knockdown`), backed by
  `Hex8Element.geometric_stiffness_matrix` and
  `assemble_geometric_stiffness`.
- Movable wrinkle through-thickness position
  (`AnalysisConfig.wrinkle_z_position`, 0.5 = mid-plane); the graded
  decay centres there.
- `AC318_S6C10_vacbag` material card — the Li 2025 vacuum-bag
  realization of the AC318 / S6C10-800 S-glass/epoxy prepreg (measured
  Xc = 335.5 MPa, E1 = 50.8 GPa).
- `IM6G_3501_6` material card — Hercules IM6G / 3501-6 carbon/epoxy
  (Vf 0.66) from Hsiao & Daniel (1996), the material behind validation
  **Dataset G** (UD carbon, measured stiffness *and* strength knockdown).
  This brings the built-in `MaterialLibrary` to 12 cards (11
  fibre-reinforced systems + the `EPOXY_S6C10` neat-epoxy card).
- Stiffness-only validation driver (`validation/validate_modulus.py`):
  compares WrinkleFE's axial Young's-modulus knockdown — the FE
  `modulus_retention` and a closed-form CLT series-average estimate of
  the off-axis lamina modulus over the wrinkle profile — against the
  measured modulus knockdown in the UD datasets E (Li 2024), F (Li 2025),
  and G (Hsiao & Daniel 1996). It is the first stiffness (as opposed to
  strength) validation in the repository; analytical MAE 3.9 % (F) /
  1.2 % (G), FE MAE 6.9 % (F) / 5.1 % (G).
- Stiffness validation chart (`validation/plot_modulus_validation.py` →
  `validation/fig_modulus_validation.png`): the modulus counterpart of
  the strength parity chart — predicted-vs-experimental modulus knockdown
  (analytical and FE) vs misalignment angle, plus a parity panel, across
  the UD modulus datasets E/F/G.
- Combined validation parity chart
  (`validation/plot_all_validation.py` →
  `validation/fig_all_validation_parity.png`): a predicted-vs-experimental
  parity plot of all single-wrinkle cases (Datasets A–F) inside a ±20%
  band, each predicted with the model that applies to it
  (Budiansky–Fleck / three-mechanism for multidirectional A–D, the
  penetration gate for UD E/F).
- `wrinklefe sweep` and `wrinklefe compare` gained `--output-json` and
  `--output-csv`: machine-readable batch results (a JSON array of
  per-run objects matching `analyze --output-json`, and a tidy
  one-row-per-run CSV with full float precision). The stdout tables are
  unchanged.
- `examples/` directory of runnable workflow scripts (basic knockdown,
  parametric sweep, morphology comparison, CZM delamination, export
  round-trip, custom material, mesh convergence), executed in CI.
- Sphinx documentation site (`docs/`) with an autogenerated API
  reference, published configuration for Read the Docs.
- **Theory: physics & mechanics** documentation page (`docs/theory.md`):
  a consolidated, code-accurate reference for the wrinkle kinematics,
  the CLT-weighted Budiansky–Fleck kink-band (incl. the confinement /
  block-penalty yield-strain model), the tension three-mechanism
  minimum, the unidirectional penetration gate, the resin-pocket and
  progressive-damage / crack-band routes, and the cohesive-zone law.
- `mesh_convergence_study()` helper and a `wrinklefe converge` CLI
  command for refinement studies.
- Structured `logging` across the analysis pipeline; `wrinklefe ... -v`
  attaches a DEBUG stderr handler.
- Multi-wrinkle finite-element solve (`AnalysisConfig.wrinkles`),
  including overlapping/interacting wrinkles.
- Contracted layup notation in the NCR validation summary
  (`to_contracted_layup`).
- Committed validation-ledger harness (`scripts/validate.py`,
  `tests/test_validation/ledger.json`).
- `provenance` block on JSON exports and the NCR summary recording the
  installed version, numerics stack, platform, and timestamp.
- GitHub issue forms, pull-request template, and this changelog.
- Docstring examples are now executed in CI (issue #296). A dedicated
  `doctests` job runs `pytest --doctest-modules src/wrinklefe`, and
  `doctest_optionflags = "NORMALIZE_WHITESPACE ELLIPSIS"` is set. The
  runnable `>>>` examples were made exact (real expected output,
  NumPy-version-stable reprs) so they act as regression guards; examples
  that need a generated mesh or a full FE solve are marked
  `# doctest: +SKIP`. Kept out of the default `addopts` so a doc example
  cannot block the core suite. Docstring-only change — no numeric
  results shift.
- Tooling — **pre-commit hooks, a coverage floor, and an armed benchmark
  gate** (issue #376 — Fixes #376). Six tooling gaps that let quality
  regress silently, closed together; no library code changes and the
  validation ledger is zero-drift.
  - **`.pre-commit-config.yaml`**, which the repo lacked entirely. Mirrors
    the `lint` CI job exactly: `ruff check` (full pyproject ruleset) on the
    commit stage, whole-tree `mypy src/wrinklefe app.py streamlit_viz.py`
    on the pre-push stage. Both hooks are `language: system`, so they run
    the tools already installed in your `[all,dev]` environment and cannot
    drift from CI's versions — a pinned upstream mirror would install a
    second ruff/mypy and, for mypy, one without streamlit/plotly, which is
    the failure mode issue #374 fixed in CI. Install with
    `pre-commit install` and `pre-commit install --hook-type pre-push`; see
    CONTRIBUTING.md.
  - **A coverage floor** — `[tool.coverage.report] fail_under` in
    `pyproject.toml`. Coverage was measured and uploaded but never
    enforced, so a new module could land at zero coverage with no signal.
    The floor is set from a measured full-suite total, deliberately a
    couple of points below it so ordinary variation does not flake the
    build. It gates the `test-full` job, whose run is what the floor was
    measured against; the OS/Python matrix job also passes `--cov` but
    deselects the slow integration solves, so its coverage is structurally
    lower and it explicitly opts out with `--cov-fail-under=0`. The floor
    is a ratchet: raise it as coverage rises, never lower it to go green.
  - **CLI tests for `converge` and `materials`**, which had zero
    arg-parsing or wiring coverage. Sixteen tests following the existing
    patched-engine pattern: flag→config mapping, `--refine` parsing,
    `--layup`/`--material` threading, `--save-plot`, clean one-line error
    handling, and bogus-argument exit-2 cases that assert the engine never
    runs.
  - **A requirements-consistency test.** `requirements.txt` had drifted to
    `numpy>=2.1` against pyproject's `numpy>=1.24` with only a
    "keep versions in sync" comment holding the line.
    `tests/test_requirements_consistency.py` now asserts the files cannot
    contradict pyproject (which stays the source of truth for install
    metadata) while still allowing the deploy file's tighter pins.
  - **The benchmark regression gate is armed.** The 2x median compare step
    had never executed, because `tests/test_benchmarks/baseline/` was never
    bootstrapped. A baseline is now committed. See **Fixed** for the
    artifact-upload bug that made the documented bootstrap impossible.
  - **The three `.claude/skills` docs are current again.**
    `expand-mypy-coverage` and `fix-ruff-violations` still described
    migrations that completed long ago (init-file-only mypy with ~102
    errors; a `--select E9,F63,F7,F82` ruff starter scope with ~840
    deferred violations), so an agent following them verbatim would have
    *narrowed* CI scope. Both are rewritten as regression guards over the
    real, finished scopes, and `pre-commit-hooks` now points at the actual
    config above.

### Fixed
- CLT — **ply stress was recovered from the total strain under a thermal
  load** (issue #273). `Laminate.ply_stresses_global` computed
  `sigma = Qbar @ eps` with the *total* strain, which is correct only
  when `delta_T == 0`; the free-expansion part has to be removed first
  (`sigma = Qbar (eps − alpha_global · delta_T)`). Without it an
  unrestrained single ply heated with no mechanical load reported
  non-zero stress, and a cross-ply cool-down came out with the matrix
  direction in *compression* where the physics puts it in tension: for
  IM7/8552 `[0/90]s` at `delta_T = -155` the 90° transverse stress went
  from −2.25 MPa to its correct **+34.3 MPa**. The sign of the
  matrix-cracking driver was inverted — the same class of error #133
  fixed at the resultant level. The correction is gated on
  `delta_T != 0`, so every purely mechanical call is unchanged and the
  validation ledger shows zero drift.
- Docs — the **stochastic module's Jensen argument had the sign
  backwards** (issue #394 — Fixes #394). `wrinklefe.stochastic` claimed
  the compressive knockdown `1 / (1 + theta_eff / gamma_Y)` is *concave*
  in the misalignment angle, so a deterministic run at the mean input
  would overestimate the mean knockdown. The law is **convex**
  (`d²KD/dθ² = 2 / (γ_Y²(1 + θ/γ_Y)³) > 0`), and the measured gap on
  `examples/14_stochastic_knockdown.py` runs the other way (`+0.0011`,
  sampling *above* the point value) — which is why the example
  deliberately reported the gap without explaining it. Measured over
  A ∈ [0.05, 1.0] mm at λ = 16 mm, the compressive path is convex at
  every sampled point, while the **tension path and the compressive path
  below the penetration gate's `D/T` clamp are locally concave**, so the
  sign is not global. The docstring now says so and anchors on the
  robust statement the example already makes: read the percentiles
  (P5/P10 sit ~0.038 below the deterministic answer), not the sign of
  the mean gap. `viz.plot_kinkband_concavity`'s caption, title and
  `plot_jensen_gap`'s docstring carried the same inverted claim and are
  corrected — the plot itself was always drawing the convex picture. New
  tests in `tests/test_stochastic.py` measure the curvature by second
  differences instead of restating it. No numerical behaviour changes.
- Typing — `viz.plots_2d` **no longer annotates against a package that
  does not exist**. `plot_strength_distribution` and `plot_jensen_gap`
  imported `MonteCarloResults` / `JensenGapResult` from
  `wrinklefe.statistics.*` under a blanket `import-not-found` ignore, so
  mypy checked nothing about either argument. They now use local
  structural `Protocol`s naming exactly the attributes each function
  reads. (`wrinklefe.stochastic.ProbabilisticResults` is *not* the right
  target: it exposes `knockdown` / `strength_MPa` and no per-sample
  morphology labels or Jensen-gap breakdown.)
- Packaging — the **source distribution no longer ships `figures/`**
  (issue #264). An unused `setuptools-scm` in `[build-system].requires`
  installed a git file-finder into the isolated build environment, which
  swept every tracked path into the sdist — including ~17.6 MB of
  rendered PNGs, two of them 8.4 MB each, plus `.github/`, `.claude/` and
  the regenerated validation plots and CSVs. `setuptools-scm` is dropped
  (nothing consumed it; the version is static) and `MANIFEST.in` now
  declares the sdist contents explicitly, with
  `scripts/check_sdist_contents.sh` failing the build in both CI and the
  release workflow if an excluded path reappears. The sdist drops from
  3,752,180 to 949,736 bytes compressed (23.2 MB to 3.4 MB uncompressed);
  the wheel is unchanged.
- Citation — **one canonical repository URL** (issue #284).
  `CITATION.cff`, the README's plain-text and BibTeX citations, the clone
  commands and `usage_tracking.REPO_URL` spelled the repository three
  different ways (`wrinklefe`, `WrinkleFE`, `wrinkleFE`), which is how a
  citation count fragments; all now use
  `https://github.com/elhajjar1/wrinkleFE`. The README also states plainly
  that the existing DOI badge is the *article's* and that a software DOI
  follows once Zenodo archiving is enabled.
- Tooling — **the `benchmark-timings` CI artifact was always empty**
  (issue #376 — Fixes #376). `.benchmarks/` is a dot directory and
  `actions/upload-artifact@v4` skips hidden files by default, so the
  `benchmarks` job logged "No files were found with the provided path" and
  uploaded nothing while still reporting success. That made CONTRIBUTING's
  documented baseline-bootstrap procedure — "download the
  `benchmark-timings` artifact from a green `main` run" — impossible to
  follow, which is why the 2x regression gate stayed dormant for its whole
  life. Fixed with `include-hidden-files: true`.
- Phase 0 correctness batch (issue #374, *Fixes #374*) — six small hazards
  from a full-codebase scan; none change numerical results for valid runs
  (ledger zero-drift):
  - **Silent modulus-retention fallback.** The local (σ₁₁ proxy) and global
    (reaction-based) FE modulus-retention blocks swallowed every exception
    and set the no-knockdown value `1.0` (the local block logged nothing),
    so a bug in the FE stiffness path read as a clean bill of health. Both
    now log a `WARNING` with `exc_info` and set a companion boolean flag
    (`AnalysisResults.modulus_retention_failed` /
    `modulus_retention_global_failed`) so a fallback `1.0` is
    distinguishable from a genuinely computed `1.0`. The value stays a
    `float` (its many consumers call `float()`/format it); the flag is
    surfaced in `summary()` and serialised only when set.
  - **Stale app results.** `reset_inputs()` now also drops the run-derived
    `results` / `cfg_payload`, and the Analyze tab renders an "inputs have
    changed since this run" banner when the live sidebar no longer matches
    the payload the shown results were computed from. The Reset button moved
    to an `on_click` callback, fixing a latent `StreamlitAPIException` from
    writing widget-keyed state after instantiation.
  - **CI mypy blind spot.** The lint job now installs the `streamlit` extra
    so mypy type-checks the app against real streamlit/plotly types (was
    `Any`), and the documented local gate passes on `main` (the `app.py`
    `reset_inputs` loop annotation is fixed). `check_untyped_defs` is now
    enabled repo-wide and clean; the `parametric_sweep.py` "str → float |
    None" suspect was a typing false-positive (a legitimately heterogeneous
    params dict inferred too narrowly), fixed by annotating `DEFAULTS`.
  - **`.gitignore` traps.** The wholesale `examples/` and `validation/*`
    ignores hid new scripts from `git status`; replaced by per-directory
    `.gitignore`s that ignore only generated outputs, so drivers stay
    visible. Added the previously phantom
    `examples/08_multi_wrinkle_czm_linkup.py` (crest-to-crest CZM link-up)
    that the README already listed.
  - **Dead scipy guard.** Removed the `except ImportError: pass` around
    `scipy.stats.gaussian_kde` in `viz/plots_2d.py` (scipy is a hard
    dependency; the guard would have silently skipped the KDE).
- Results-export schema drift (issue #345): the structured JSON export
  (`wrinklefe.io.results.results_to_dict` / `export_results_json`) and the
  NCR validation summary (`wrinklefe.io.export.build_analysis_summary`)
  silently dropped several `AnalysisResults` fields. `results_to_dict`
  now serialises `modulus_retention_global` and `analytical_onset_knockdown`,
  and — for progressive-damage runs — a gated `progressive` block
  (`strength_MPa`, `pristine_strength_MPa`, `knockdown`, `n_increments`,
  and the `(strain, stress)` load history). The NCR markdown/PDF renderers
  now surface the global coupon modulus retention and a progressive-damage
  section, wired end-to-end from the Streamlit Export tab. A new
  dataclass-walking drift-guard test asserts every `AnalysisResults` field
  is either exported or on an explicit allowlist, so a future field cannot
  silently go unexported. Analytical-only runs are unchanged (no
  `progressive` block, no empty NCR rows).
- Dual-wrinkle amplitude contract in the FE mesh (issue #305): the
  `stack`/`convex`/`concave` morphologies build the mesh by summing two
  through-thickness–decayed displacement fields, and each constituent
  previously carried the full amplitude `A`, so the in-phase `stack` mesh
  peaked at `2A` — double the intended geometry, meshed fibre angle and FE
  knockdown, and inconsistent with the analytical
  `theta_max = arctan(2*pi*A/lambda)`. Each constituent is now generated at
  half amplitude `A/2` (`morphology._profile_at_half_amplitude`), so the
  `stack` mesh composes to exactly `A` and its fibre angle matches the
  analytical profile. Explicit multi-wrinkle `WrinkleSpec` configurations
  are unaffected (each listed wrinkle keeps its specified amplitude). See
  Numerical results.
- Linear-buckling eigensolve correctness/robustness
  (`solver/buckling.py`): for a wrinkled (non-uniform) pre-stress the
  geometric "mass" matrix `M = -K_geo` is **indefinite**, which violated
  the SPD assumption of the previous `eigsh` shift-invert. It returned
  spurious, run-to-run-varying eigenvalues — and on macOS arm64 sometimes
  no surviving positive mode at all, so `critical_load_factor` came back
  `inf` and the buckling-knockdown test flaked in CI. The solve is now the
  symmetric-definite pencil `M φ = μ K φ` (the material stiffness `K` is
  SPD), with `λ = 1/μ` and a deterministic ARPACK start vector — finite,
  reproducible, and matching a dense reference. This is infrastructure
  only; `microbuckling_knockdown` is still **not** the production UD
  predictor (see Numerical results).
- Documentation accuracy (physics audit): the README "How It Works" and
  `ARCHITECTURE.md` confinement section now show the full three-parameter
  effective yield strain `gamma_Y_eff = max(0.032 + 0.050·f_conf −
  0.010·max(n_block−1,0), 0.016)` (the block-penalty and floor terms were
  missing), document `theta_eff = M_f·theta_max` and the optional
  Argon–Fleck quadratic term, and correct the graded through-thickness
  decay scale to `max(λ/2, A)`. The `MaterialLibrary` docstring (and its
  doctest) now lists all 12 registered cards (11 fibre-reinforced systems
  + the `EPOXY_S6C10` neat-epoxy card) instead of a stale list of nine.
- API-reference docstrings (physics audit, follow-up): the LaRC05 module
  docstring described "iterative φ_c computation" and a Ramberg-Osgood
  nonlinear-shear amplification that the code does not apply — corrected
  to the linear closed-form load-induced φ_c that is actually used (the
  `max_phi_c_iter` / `phi_c_tol` parameters are documented as reserved /
  unused). The linear-buckling module docstring now carries the item-D.4
  negative-finding note (the eigenvalue over-predicts the UD wrinkle
  knockdown; use the penetration gate instead), matching the README and
  the new theory page.
- `wrinklefe sweep` now validates its inputs: an unknown `--parameter`,
  `--min >= --max`, or `--steps < 2` print a one-line error and exit
  non-zero (code 2) before any solve, instead of a raw traceback or a
  silently degenerate sweep. The `--parameter` help no longer implies
  only amplitude/wavelength/width are accepted.
- Graded-morphology compression knockdown now honours `decay_floor`
  (previously inert in compression while honoured in tension).
- JSON export stamps the real installed version instead of a hardcoded
  `0.1.0` literal.
- Latent crashes surfaced by static typing: `np.trapz` removal under
  numpy 2.0 in the stress-resultants path; `WrinkleSurface3D` attribute
  access in `max_angle` / `fiber_angles_at_nodes`.

### Changed
- Viz — the **Plotly figure library is now part of the package** as
  `wrinklefe.viz.plotly_figs` (issue #286 — Fixes #286). It used to live in
  the repository-root `streamlit_viz.py`, which the src layout keeps out of
  the wheel: `pip install wrinklefe` shipped none of the interactive 3D
  layer, so it was reachable only from a checkout. Notebook users can now
  do
  ```python
  from wrinklefe.viz import mesh3d_figure, stress_contour_figure
  ```
  after `pip install 'wrinklefe[plotly]'` — a new lean extra, so the
  interactive figures no longer require the whole Streamlit server stack
  (the `streamlit` extra still carries plotly as before).
  - **Plotly stays optional.** The names are re-exported through a PEP 562
    module `__getattr__`, so `import wrinklefe.viz` never imports plotly;
    the import happens on first attribute access, and without plotly that
    access raises an `ImportError` naming the install command rather than a
    bare `No module named 'plotly'`. Both halves are verified against the
    built **wheel** in clean venvs by the `build` CI job.
  - **`streamlit_viz.py` remains at the root as a deprecated re-export
    shim**, so `import streamlit_viz` keeps working for the hosted
    Streamlit deployment and any external references. The functions
    themselves are unchanged — the shim re-exports the very same objects.
    Deprecation is stated in its docstring only: no `DeprecationWarning`,
    because the app re-imports it on every script rerun.
  - The Plotly figure API is now in the Sphinx reference (`api/viz`),
    which it could not be while the module lived outside the package.
- Docs — **`internal/ARCHITECTURE.md`'s data flow now matches the code**
  (issue #378 — Fixes #378). The old diagram described a 2024-era pipeline
  and omitted most of what has been added since. It now follows
  `WrinkleAnalysis.run` as written: the wrinkle-field construction
  (transverse surfaces, multi-wrinkle), the analytical branch through the
  penetration gate, the ordered per-element material overrides (crest resin
  lens → surface resin pockets → compaction `Vf` gradient), the three FE
  branches (CZM, progressive damage, linear), and the retention-factor
  baseline. Adds the three defect mechanisms the modelling section omitted
  (surface resin pockets, the compaction `Vf` gradient, through-width
  transverse surfaces) and module-table rows for the capabilities that had
  none (`core/compaction.py`, `core/micromechanics.py`, `core/layup.py`,
  `core/cohesive_mesh.py`, `elements/cohesive8.py`, the remaining failure
  criteria, `convergence.py`, `goalseek.py`, `stochastic.py`,
  `io/results.py`).
- Tooling — **the benchmark compare step now runs on every build**
  (issue #376 — Fixes #376), instead of skipping for want of a baseline.
  Because no runner-generated artifact existed to seed it (see **Fixed**),
  the committed baseline is container-generated: produced on Python 3.12 so
  its pytest-benchmark machine id matches the job's storage key, but on
  different hardware from a GitHub runner. Absolute timings are therefore
  not comparable, so the step carries `continue-on-error: true` and runs as
  a visible report rather than a hard gate. CONTRIBUTING.md documents the
  three steps to promote it to blocking once a runner artifact is
  available. Also relaxes `pytest==9.1.1` to `pytest>=9,<10` in
  `requirements-test.txt` so patch upgrades are not frozen.
- App — **the analysis run is no longer wrapped in `@st.cache_data`**
  (issue #377 — Fixes #377). That decorator was what blocked live progress
  widgets (issue #242: Streamlit refuses element calls made inside a
  cache-decorated function against a layout block created outside it). The
  solve now runs uncached and the *result dict* is cached by hand in
  `st.session_state`, keyed on exactly the same hashable `cfg_payload` the
  decorator used to hash — so re-running an identical configuration is
  still instant and returns the same object. Unlike `@st.cache_data`, the
  manual cache is **bounded to the 4 most recently used payloads**: FE
  results carry mesh, displacement, stress and failure-index arrays, and
  keeping every distinct run of a browsing session was needless memory.
  *Reset to defaults* clears the cache along with the displayed results.
  `run_analysis_cached(cfg_payload, progress_callback=None)` remains the
  entry point and returns an unchanged result-dict shape; the solve itself
  now lives in an undecorated `_run_analysis`.
- Solver — **the ILU→diagonal preconditioner fallback is now loud and
  narrow** (issue #265). When the iterative solver's ILU factorisation
  fails, `StaticSolver._solve_iterative` emits a `logging.WARNING`
  *unconditionally* (previously only a `print` under `verbose`), quoting
  the original error type and message and stating that the diagonal
  (Jacobi) fallback is in effect. Only the exceptions `spilu` raises for a
  genuine factorisation failure (`RuntimeError`, `MemoryError`,
  `ValueError`) are caught; any other exception type now propagates
  instead of masquerading as "ILU failed". The CG non-convergence
  `RuntimeError` now also names the active preconditioner, the iterations
  used, the cap, and the final relative residual. No change to converged
  numerics.
- App — **surface-pocket controls relocated into the morphology
  definition** (issue #371, Part B). The standalone *Surface resin pockets*
  expander in the Expert FE section is gone; its controls now live directly
  under the Morphology selector. For `tool_flat` the pockets are implicit
  (auto-enabled, shown as a caption) with the pinned-side and
  transition-ply controls inline; the legacy tool-flat morphologies
  (`stack`/`convex`/`concave`, `graded` with `decay_floor=0`) keep an
  *advanced* opt-in whose help explains their pockets are inherently small
  (~0.25 ply thickness — use `tool_flat` for significant pockets). The
  `sb_surface_transition_plies` widget key joins `DEFAULTS` so *Reset to
  defaults* round-trips it. Rationale: under the linear-decay morphologies
  the pockets were mechanically negligible by construction (measured
  ~0.25 ply-thickness trough gap), so a tooling-dominated wrinkle belongs
  in the morphology definition, not as an add-on toggle.
- Packaging — **PyVista/VTK moved to an optional `vtk` extra** (issue
  #302). Plain `pip install wrinklefe` no longer pulls in VTK (~150 MB
  lighter) and stays headless-safe; the 3D cohesive-zone plots
  (`plot_interface_damage_3d` / `plot_crack_front_3d`) now require
  `pip install "wrinklefe[vtk]"` (also included in `[all]`). PyVista was
  already imported lazily, so `import wrinklefe`, the CLI, the Streamlit
  app, and the docs build are unaffected when it is absent; the
  `_require_pyvista` error message now names the `wrinklefe[vtk]` extra.
- Streamlit app — acknowledgment gate and intro reworded to a
  professional-tool framing (issue #333): the gate leads with what
  WrinkleFE computes rather than "free academic software", and the
  supporting copy, email placeholder, and acknowledgment checkbox use
  neutral, work-agnostic wording. Gate mechanics, the
  `WRINKLEFE_DISABLE_GATE` off-switch, usage logging, and the
  MIT/attribution facts are unchanged.
- Streamlit app — the **Configure** and **Results** tabs are merged into
  a single default **Analyze** tab (issue #334), leaving three tabs
  `["Analyze", "Export", "Help"]`. The wrinkle/laminate preview lives in
  an expander that is open before the first run and auto-collapses once
  results exist, so results lead the view after a run. Export and Help
  are unchanged. (Supersedes the #358 tab-order entry below.)
- Results-export `SCHEMA_VERSION` bumped `1.0` → `1.1` (issue #345): the
  additive `modulus_retention_global`, `analytical_onset_knockdown`, and
  gated `progressive` fields in the structured JSON export. Additive only;
  existing consumers of 1.0 fields are unaffected.
- `AnalysisConfig` now validates ply angles at construction (issue #344)
  — **breaking only for previously-accepted invalid inputs**: a config
  with `|angle| > 90` (e.g. `angles=[900.0, 0.0, 452.0]`) now raises
  `ValueError` naming the offending index and value instead of silently
  flowing a non-canonical angle into CLT trig, where the tension-mechanism
  heuristic mis-classified it (a 900° ply read as a 90° ply). The check
  reuses the shared `validate_ply_angle` rule from `core.layup` (issue
  #343), so the parser and the config validator can never drift. Valid
  layups (decimals, `±90`, long stacks) are unaffected. Follow-up:
  `Laminate.from_angles` is intentionally left unvalidated for now.
- `parse_layup` input strictness (issue #308) — **breaking for inputs
  that previously parsed**: ply-angle tokens with `|angle| > 90` (e.g.
  the repeat-count-like `[02/902]s`, which silently parsed as 2° and
  902° plies) and leading-zero tokens now raise `ValueError` instead of
  building a wrong laminate; ASCII `+-45`/`-+45` are now accepted as
  the ± shorthand. Scripts feeding the newly-rejected forms must switch
  to explicit repeat syntax (e.g. `[0_2/90_2]s`).
- Mesh aspect-ratio warning re-baselined (issue #303): the warning now
  compares each element against the mesh's own median aspect ratio and
  flags only outliers, instead of a fixed 10:1 threshold that flagged
  71–100 % of elements on typical thin-ply meshes (default runs are now
  quiet; genuinely anomalous elements still warn).
- CI enforces the full Ruff ruleset and `mypy` over the whole tree.
- Streamlit app — the **Cohesive Zone Modeling** sidebar controls now
  render only in **Expert mode**. CZM requires the full nonlinear FE
  solve, which is itself expert-only (novice mode forces
  `analytical_only=True`), so the control now sits with the other expert
  FE settings instead of the simplified novice sidebar. The analytical
  path and CZM behaviour are unchanged.
- Streamlit app — tab layout: **Configure** is now the first (default)
  tab, so the app opens on the laminate/geometry view instead of the
  intro. The old **Overview** tab moved to the end and is renamed
  **Help**. Within the Configure tab the **Wrinkle mid-surface profile**
  plot moved to the bottom, so the through-thickness cross-section leads.

### Removed
- Docs — **the internal CZM execution plan is no longer published**
  (issue #378 — Fixes #378). `docs/czm_plan.md` was a one-line shim that
  pulled `internal/CZM_PLAN.md` — a *completed* internal execution plan,
  agent-orchestration tables and all — into the user-facing Sphinx site.
  The shim is deleted and `czm_plan` is dropped from the toctree; the note
  itself stays on disk for repository readers. Sphinx skips its
  "not included in any toctree" check for any file some page `include`s,
  which is why the remaining `internal/*.md` need no marker (their shims
  still exist); removing this one's shim would have tripped that check under
  `-W`, so `internal/CZM_PLAN.md` is listed in `exclude_patterns` in
  `docs/conf.py` with the reasoning recorded there. No user-facing content
  is lost: nothing in the published docs linked to it.
- The dead `export` optional-dependency extra (`meshio` was never
  imported; native `.inp`/VTK writers need no extra).

### Numerical results
- **Thermal residual stress (issue #273, Stage 1)**: opt-in and off by
  default (`delta_T = 0.0`), so every existing result — the validation
  ledger included — is byte-identical. When it is switched on the CLT
  numbers move, as intended and as measured: on an IM7/8552
  quasi-isotropic `[0/45/-45/90]s` wrinkled coupon (graded, A = 0.5 mm,
  λ = 16 mm, compression) a `delta_T = -155` cure cool-down adds
  **+34.3 MPa** of transverse tension to every ply — 55 % of
  `Yt = 62.3 MPa` — and takes the CLT first-ply-failure load factor from
  173.3 to 3.62 (**−97.9 %**), with the critical mode flipping from
  `matrix_compression` to `matrix_tension`. The reference mechanical
  resultant the pipeline uses is deliberately small, so the thermal term
  dominates; that is precisely why omitting it from a run that asks for
  it would be a wrong number rather than a small one. The closed-form
  `analytical_knockdown` is unchanged (Budiansky–Fleck has no
  temperature term).
- **Compaction Vf gradient (issue #379, Part B)**: opt-in and off by
  default, so every existing result — the ledger included — is unchanged.
  When it is switched on for a `tool_flat` run the FE numbers move, as
  intended and as measured: on a 24-ply UD two-caul case (`tool_flat`,
  `surface_pocket_side="both"`, A = 0.25 mm, IM7/8552, 40 × 4 mm domain)
  `modulus_retention_global` goes 0.937591 with the binary surface pockets
  to 0.944565 with the gradient (+0.006974, +0.744 %) — the compacted crest
  band stiffens more than the resin-rich trough softens, and the
  continuous field replaces a binary neat-resin tag. Against a wrinkle with
  no trough treatment at all (0.957389) the gradient is 0.012824 *softer*,
  so it lands between the two, as the physics implies. 64 of 1920 elements
  saturate at the `vf_max = 0.75` cap at that amplitude.
- **Penetration gate × multi-wrinkle (issue #342)**: with a
  `penetration_gate` preset and the geometry supplied via
  `AnalysisConfig.wrinkles`, the gate previously took its angle from the
  wrinkle specs but its penetration `D/T` from the leftover scalar
  `cfg.amplitude` (typically the unused 0.366 default) — a silently
  plausible wrong knockdown (0.98 instead of 0.64 on the issue's repro).
  The gate now evaluates per spec — `theta_i = arctan(2πA_i/λ_i)`,
  `D_i/T = A_i/T`, and `z_i = (ply_interface+1)/n_plies` through the
  position factor `P(z)` (`cfg.wrinkle_z_position` is a scalar-path
  parameter and is ignored when specs are present) — and returns the
  weakest-link (minimum) knockdown over the wrinkles. Scalar-config gate
  results are unchanged (pinned ledger baselines show zero drift); any
  gate × `wrinkles` configuration returns different (correct) values.
- **LaRC05 / Puck last-bit normalization (issue #299)**: the scalar
  `evaluate()` paths now square via an explicit product (`x * x`)
  instead of `x ** 2` — scalar `np.float64` pow routes through libm and
  could land 1 ULP away from the exact product the vectorized field
  path computes. Failure indices from these two criteria may therefore
  shift by at most one floating-point ULP (≈ 1e-16 relative) toward the
  exactly-rounded value; no physical or tolerance-visible change.
- **Linear-buckling microbuckling knockdown**: with the eigensolve
  corrected (indefinite `-K_geo` handled via the symmetric-definite
  pencil), the bifurcation load of the homogenised ply-mesh *rises* with
  the wrinkle (tilted fibres carry less destabilising axial pre-stress;
  e.g. the Li 20 mm coupon goes pristine λ ≈ 8.30 → amplitude-0.6 wrinkle
  λ ≈ 8.65), so `microbuckling_knockdown` returns ≈ 1.0 (no knockdown)
  rather than the spurious sub-1.0 values the old indefinite-`M` solve
  produced. This sharpens the documented negative finding (item D.4): the
  linear eigenvalue gets the wrinkle-knockdown *sign* wrong, which is why
  the UD wrinkle knockdown is taken from the penetration gate. No
  production prediction path consumes `microbuckling_knockdown`, so no
  user-facing knockdown changes.
- **Penetration gate**: when `AnalysisConfig.penetration_gate` is set to
  a `GateParameters` preset, `analytical_knockdown` (and
  `analytical_strength_MPa`) are computed from the two-parameter
  (θ, D/T, z) gate instead of Budiansky–Fleck, so UD configurations
  return different (calibrated) knockdowns. The default
  `penetration_gate=None` preserves previous results bit-for-bit.
- **Graded compression `decay_floor`**: configurations that set
  `decay_floor` under compression now produce different (correct)
  knockdowns. The default `decay_floor=0.0` preserves previous results
  bit-for-bit.
- **Multi-wrinkle fibre angles**: fibre-misalignment fields now derive
  from the slope of the composed displacement field ("compose then
  differentiate"). FE results shift for dual-wrinkle morphologies
  (`stack`/`convex`/`concave`) wherever through-thickness decay < 1 or
  wrinkles overlap; single-wrinkle results are unchanged. Analytical
  predictions are unaffected.
- **Dual-wrinkle mesh amplitude (issue #305)**: the `stack`/`convex`/
  `concave` FE meshes previously peaked at up to `2A` because each of the
  two summed constituents carried the full amplitude `A`. Each constituent
  is now half amplitude, so the in-phase `stack` mesh peaks at exactly `A`
  and its meshed fibre angle drops to the analytical
  `arctan(2*pi*A/lambda)` (roughly halved). FE-derived quantities (fibre
  angles, FE knockdown, stiffness retention) shift for these three
  morphologies; single-wrinkle (`uniform`/`graded`) and explicit
  `WrinkleSpec` multi-wrinkle configurations are unchanged, and analytical
  predictions (which already used the configured `A`) are unaffected.

## [1.0.0]

Initial public release: analytical Budiansky–Fleck knockdown plus a 3-D
finite-element pipeline with LaRC05/Hashin/Puck ply failure and
cohesive-zone delamination; five wrinkle morphologies; the material
library; JSON/CSV/Abaqus/VTK export; a command-line interface; and the
Streamlit web application.

[Unreleased]: https://github.com/elhajjar1/wrinkleFE/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/elhajjar1/wrinkleFE/releases/tag/v1.0.0
