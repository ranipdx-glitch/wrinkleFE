# WrinkleFE

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://wrinklefe.streamlit.app/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs41598--025--06693--4-blue.svg)](https://doi.org/10.1038/s41598-025-06693-4)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)
[![GitHub stars](https://img.shields.io/github/stars/elhajjar1/wrinkleFE?style=social)](https://github.com/elhajjar1/wrinkleFE)

An open-source Python finite element package for predicting strength and stiffness knockdown in composite laminates containing fiber waviness defects.

> ⭐ **Found WrinkleFE useful?** Please [star the repository](https://github.com/elhajjar1/wrinkleFE) and [cite it](#citation) — it's a free academic project, and stars and citations are what keep it supported.

## Try it in your browser

The fastest way to use WrinkleFE is the hosted Streamlit app — no install required:

### [Launch the WrinkleFE Streamlit app](https://wrinklefe.streamlit.app/)

Pick a material, set the wrinkle amplitude / wavelength / morphology, and the
app returns the analytical knockdown, plots, and (optionally) a full FE solve.
Public link, no account needed.

## Features

- **Compression model:** CLT-weighted Budiansky-Fleck kink-band with layup-dependent confinement
- **Tension model:** Three-mechanism (fiber cos^2 theta, Hashin matrix, curved-beam sigma_33 delamination) with thick-ply in-situ correction
- **3D finite element:** Structured hexahedral mesh with LaRC04/05 failure criteria
- **Five morphologies:** Stack, convex, concave, uniform, graded (with configurable decay floor)
- **Graded averaging:** Through-thickness ply-averaged knockdown for graded wrinkles
- **Movable wrinkle position:** Configurable through-thickness placement (`wrinkle_z_position`, 0.5 = mid-plane)
- **Multi-wrinkle configurations:** Arbitrary N-wrinkle layouts via `AnalysisConfig.wrinkles` — a list of `WrinkleSpec(amplitude, wavelength, width, ply_interface, phase_offset)` — through the analytical, FE, penetration-gate (per-spec, weakest-link) and CZM paths
- **Cohesive-zone delamination (CZM):** Bilinear traction–separation interface elements (`enable_czm=True`) with per-interface damage, energy and crack-length reporting — including continuous cohesive surfaces across adjacent wrinkles for crest-to-crest delamination link-up (see `examples/08_multi_wrinkle_czm_linkup.py`)
- **Resin-pocket material zone:** Graded neat-epoxy lens at the wrinkle crest (modulus + fibre-angle blend, counted once) via `wrinklefe.core.resin_pocket`
- **Tool-flat surfaces & surface resin pockets:** Parts cured against rigid tooling / a caul sheet keep perfectly flat outer surfaces while the fibres undulate internally; the wrinkle troughs fill with neat resin under the flat surface (`enable_surface_resin_pockets`, `surface_pocket_side` = `top`/`bottom`/`both`). FE-only, trough-following, volume-conserving; composes with the crest lens
- **Compaction Vf / ply-thickness gradient:** Under rigid tooling the resin is squeezed out of the compacted regions and pools where the geometry opens up, so ply thickness and local fibre volume fraction vary through the wrinkle (`enable_vf_gradient`, `tool_flat` only — `surface_pocket_side="both"` is the two-caul-plate case). `Vf_local = vf_nominal · h0/h` per element (fibre content conserved), and the local ply is the preset **scaled by the micromechanics Vf ratio** (`wrinklefe.core.compaction`); FE-only, opt-in, and the continuous generalization of the binary surface pockets
- **Progressive-damage FE:** Load-stepping `ProgressiveDamageSolver` to ultimate load with optional crack-band (Bažant–Oh) regularization — the first FE route to a UD compression knockdown
- **Penetration gate (θ, D/T, z):** Closed-form two-parameter UD predictor `KD = 1 − (1 − KD_angle(θ))·S(D/T)·P(z)` with calibrated presets; zero FE cost (`wrinklefe.core.penetration_gate`)
- **Linear buckling:** Geometric-stiffness eigenvalue solve (`LinearBucklingSolver`) with a microbuckling knockdown — verified structural-buckling infrastructure, but a homogenised-continuum eigenvalue does not capture the *fibre-scale* wrinkle knockdown (it gets the sign wrong: the bifurcation load *rises* with the wrinkle), so it is **not** the production UD predictor (the penetration gate is); see the [modelling findings](docs/wrinkle_modeling_findings.md)
- **11 built-in laminate materials** (AS4/3501-6, IM7/8552, T300/914, T700/2510, AC318/S6C10 S-glass/epoxy, T800S/M21, IM10/8552, IM6G/3501-6 carbon/epoxy — the Hsiao & Daniel 1996 wavy-UD study, S-2 glass/epoxy, Kevlar-49/epoxy, plus `AC318_S6C10_vacbag` — the Li 2025 vacuum-bag realization, measured Xc=335.5 MPa, E1=50.8 GPa) plus an isotropic neat-epoxy card (`EPOXY_S6C10`) for the resin-pocket zone
- **Inverse goal-seek (acceptance limits):** `wrinklefe.goalseek.find_critical_value` / `wrinklefe critical` answer the disposition question directly — the largest wrinkle amplitude (or any float config field) that still meets a target knockdown or allowable strength. The answer is backed off to the safe side and verified by a real forward run, and the search refuses with an actionable message when the curve admits no unique root
- **Process-parallel sweeps:** `n_workers=N` on both sweep APIs and `wrinklefe sweep --parallel N` fan the independent per-point solves across CPU cores (results identical to and ordered like the sequential run)
- **Comprehensive test suite** covering all modules (run `pytest` to see the current count)

## Developer / library install

If you want to script against the package or contribute to development:

```bash
pip install wrinklefe
```

Or install the latest source:

```bash
git clone https://github.com/elhajjar1/wrinkleFE.git
cd wrinklefe
pip install -e ".[all]"
```

The 3D cohesive-zone renders (`plot_interface_damage_3d` /
`plot_crack_front_3d`) use PyVista, which pulls in VTK (~150 MB). It is an
optional dependency, so plain `pip install wrinklefe` stays lean and
headless-safe; install `pip install "wrinklefe[vtk]"` when you need those
plots (it is already included in `[all]`).

Verify the install:

```bash
python -c "import wrinklefe; print('WrinkleFE installed successfully')"
```

Run the test suite:

```bash
pytest
```

## Quick Start

### Streamlit web app

The fastest path is the hosted Streamlit instance — no install, no Python
required: **<https://wrinklefe.streamlit.app/>**. Pick a material from the
sidebar (or **Custom…** to enter your own elastic constants and strength
allowables), enter a layup in contracted notation (e.g. `[0/45/-45/90]_3s`),
set the wrinkle geometry, and click **Run analysis**. The app ships with
per-morphology schematic cartoons (including the `tool_flat`
surface-resin-pocket morphology, whose pinned-side and transition-ply
controls appear inline with the Morphology selector), a live wrinkle
preview, and the same analytical + FE pipeline as the Python API. Expert
mode also exposes the through-width `transverse_mode` envelope, and the
sidebar **Config file** section can **download** the current settings as a
portable `AnalysisConfig` JSON — round-tripping with the CLI `--config` /
`--save-config` flags — or **load** a saved case (JSON/YAML) back into the
sidebar.

A run reports live progress while it solves: the status box tracks the
pipeline phase by phase (building the laminate → analytical predictions →
assembling the FE mesh → solving → evaluating failure → retention
factors), so a long FE or CZM solve is distinguishable from a hung one.
Re-running an identical configuration is served from a small in-session
result cache instead of re-solving; *Reset to defaults* clears it.

A second sidebar button, **Find acceptable limit**, runs the inverse
search on those same inputs — the browser form of
[`wrinklefe critical`](#finding-the-maximum-acceptable-defect). Set a
target knockdown (or an absolute MPa allowable) in the
*Acceptable-limit settings* expander and the app reports the largest
amplitude the laminate can carry and still meet it, together with the
scanned knockdown curve and the per-evaluation ledger. The reported
number is the conservative one — backed off from the root and verified
by a real forward run — and one click applies it back to the sidebar so
you can run the full analysis at the limit. When the search converged
for the *same* configuration a run was made on, the limit is carried
onto the NCR validation summary on the **Export** tab; a limit searched
against different inputs is deliberately left off. If the curve admits
no unique answer the app shows the engine's diagnosis and the measured
curve instead of a number. See
[`DEPLOYMENT_STREAMLIT.md`](docs/internal/DEPLOYMENT_STREAMLIT.md) for the full feature
tour and instructions for self-hosting.

To run the app locally:

```bash
pip install -r requirements.txt
pip install -e .
streamlit run app.py
```

The app's interactive 3D views are ordinary library functions, so you can
draw the same figures in a notebook — the mesh surface, the σ contour, the
deformed shape, the failure-index surface and the y-slice scatters:

```bash
pip install 'wrinklefe[plotly]'
```

```python
from wrinklefe.viz import mesh3d_figure, stress_contour_figure

fig = stress_contour_figure(nodes, elements, stress_per_elem)
fig.show()
```

Plotly is optional: `import wrinklefe.viz` works without it, and only
asking for one of these figures raises an error telling you what to
install. See the [`wrinklefe.viz` API reference](docs/api/viz.rst).

### Python API

```python
from wrinklefe.analysis import AnalysisConfig, WrinkleAnalysis

config = AnalysisConfig(
    amplitude=0.366, wavelength=16.0, width=12.0,
    morphology="stack", loading="compression",
)
result = WrinkleAnalysis(config).run()
print(result.summary())
```

Runnable scripts for the common workflows — parametric sweeps,
morphology comparison, CZM delamination, the penetration gate,
progressive damage, the crest resin lens, stochastic propagation, export
round-trips, custom materials — live in [`examples/`](examples/); each
states its expected runtime and output, and CI executes them all so they
stay current.

Three pages carry most of what a new user needs:

- [Tutorial](docs/tutorial.md) — a worked walkthrough from a micrograph
  measurement to the NCR validation attachment.
- [Interpreting results](docs/interpreting_results.md) — what each
  headline number means, what it does not mean, which one governs, and
  the severity bands (`io/export.py` remains authoritative for those).
- [Units & conventions](docs/units_conventions.md) — mm / MPa / N/mm,
  strain and `Vf` as fractions, angle and coordinate conventions, signs.

The full API reference and user guide are built from [`docs/`](docs/)
with Sphinx (`pip install -e ".[docs]" && sphinx-build -W docs
docs/_build`) and published at
<https://wrinklefe.readthedocs.io>.

`amplitude` (`A`) is the **half-amplitude** [mm]: the peak displacement
of the wrinkled mid-surface from the flat (unwrinkled) reference plane,
so `z(x) = A·cos(2πx/λ)` (modulated by the envelope) and the
peak-to-trough height is `2A`. For a measured wrinkle (e.g. from a
cross-section micrograph or CT slice), `A = (z_max − z_min) / 2`. The
peak fibre misalignment angle scales as `θ_max ≈ arctan(2πA/λ)`, which
drives the Budiansky-Fleck compressive knockdown.

### Wrinkle geometry parameters

All wrinkle length parameters use a single, consistent unit:
**millimetres (mm)** — the same unit as `ply_thickness` and
`domain_length` (the default `amplitude=0.366` mm is exactly two ply
thicknesses of 0.183 mm). Lengths are **not** normalized by thickness.
The longitudinal coordinate `x` runs along the laminate in the fibre
direction; out-of-plane displacement `z(x)` is measured from the flat
(undeformed) mid-surface. Angles are in **radians**.

This table is the canonical reference for every wrinkle-geometry
parameter exposed by `AnalysisConfig`, the CLI, and the Streamlit UI.
The `AnalysisConfig` docstrings, CLI `--help`, and Streamlit `help=`
tooltips all mirror these definitions; the
`tests/test_param_docs_match.py` regression test pins the defaults so
the docs and the dataclass cannot drift.

| Parameter | Units | Default | Definition | Constraint |
|-----------|-------|---------|------------|------------|
| `amplitude` (A) | mm | `0.366` | Half-amplitude: peak displacement of the wrinkled mid-surface from the flat reference, with `z(x) = A·cos(2π(x − x₀)/λ)` (modulated by the envelope) and peak-to-trough height `2A`. For a measured wrinkle, `A = (z_max − z_min)/2`. | ≥ 0 (`0` = flat / no wrinkle) |
| `wavelength` (λ) | mm | `16.0` | Spatial period of the `cos(2π(x − x₀)/λ)` carrier along the longitudinal x-direction (crest-to-crest distance). Wavenumber `k = 2π/λ`. | > 0 |
| `width` (w) | mm | `12.0` | Longitudinal envelope decay length about the centre `x₀`. Exact meaning is profile-dependent: Gaussian 1/e length scale in `exp(−(x−x₀)²/w²)`, tapered flat-top extent (`\|x−x₀\| < w/2`), or triangular half-base (`\|x−x₀\| < w`). Also used as the transverse (y-direction) extent of the wrinkle in 3-D dual-wrinkle / graded mesh deformation. | > 0 |
| `phase` (φ) | rad | `None` | Explicit dual-wrinkle phase offset between the two wrinkle centrelines. `None` derives φ from `morphology` via `MORPHOLOGY_PHASES` (stack φ=0, convex φ=+π/2, concave φ=−π/2). A float overrides the named-morphology phase so arbitrary offsets can be swept (e.g. `0` to `π`). Ignored for single-wrinkle morphologies (`uniform`, `graded`). | finite when set |
| `decay_floor` | dimensionless | `0.0` | Graded morphology only: minimum fraction of the wrinkle amplitude retained at the laminate outer surfaces. `0.0` = full decay to zero amplitude at the surfaces (pure graded); `1.0` = no decay (equivalent to `uniform`). | in `[0, 1]` |
| `amplitude_profile` | name | `"constant"` | Spatially varying in-plane modulation of the wrinkle amplitude `A`, applied on top of the wrinkle's own longitudinal envelope. `"constant"` (default) preserves the legacy uniform `A`; `"gaussian"` multiplies `A` by `exp(−(s/d)²)`; `"linear"` multiplies `A` by `max(0, 1 − \|s\|/d)` (clipped). `s` is the coordinate from the wrinkle centre along `amplitude_profile_axis` and `d` is `amplitude_profile_decay_length`. | one of `constant`, `gaussian`, `linear` |
| `amplitude_profile_decay_length` | mm | `None` | Decay length `d` (mm) for the Gaussian sigma or linear-decay extent. `None` falls back to the wrinkle profile's own `width`. Ignored when `amplitude_profile == "constant"`. | finite and > 0 when set |
| `amplitude_profile_axis` | axis | `"x"` | In-plane axis along which the amplitude modulation runs. Pick `"y"` for an independent transverse tapering of `A` that does not stack with the existing longitudinal envelope on `x`. | one of `x`, `y` |
| `transverse_mode` | name | `"uniform"` | Through-width (transverse, y-direction) wrinkle-surface envelope `f(y)`. `"uniform"` (default) builds the bare x-only wrinkle exactly as before (bit-identical). The non-uniform modes wrap the profile in a `WrinkleSurface3D` so the crest amplitude varies across the specimen width: `"gaussian_decay"` decays it toward the edges (localized mid-width defect), `"sinusoidal_y"` ripples it across the width, and `"elliptical"` confines it to a mid-width patch. **FE-only** — a non-uniform mode requires `analytical_only=False` and is not yet combinable with multi-wrinkle (`wrinkles`) or `enable_czm` (both rejected at construction). | one of `uniform`, `gaussian_decay`, `sinusoidal_y`, `elliptical` |
| `transverse_span` | mm | `None` | Specimen width `span_y` seen by the transverse envelope. `None` tracks `domain_width` so the envelope always spans the meshed y-extent. Ignored when `transverse_mode == "uniform"`. | finite and > 0 when set |
| `transverse_width` | mm | `None` | Transverse localization half-width `width_y`: the Gaussian 1/e length for `"gaussian_decay"` and the ellipse half-width for `"elliptical"` (ignored by `"uniform"`/`"sinusoidal_y"`). `None` resolves to `span_y / 4` — a localized mid-width patch whose amplitude has fallen to `exp(−4) ≈ 0.018` of the crest at the edges (`gaussian_decay`) or that occupies the central half of the width (`elliptical`). | finite and > 0 when set |
| `surface_transition_plies` | count | `2` | `tool_flat` morphology only: number of plies over which the amplitude ramps linearly from the full-amplitude core to zero at the pinned (tool-flat) surface. A short transition concentrates the full-amplitude trough mismatch into a thin surface band (significant resin pockets), but the crest-side transition elements compress by `amplitude / surface_transition_plies`, so the amplitude is bounded: `amplitude ≤ 0.8 · surface_transition_plies · ply_thickness / nz_per_ply` (enforced at construction; exceed it and the message names both remedies). | integer ≥ 1 |

Peak fibre misalignment: `θ_max ≈ arctan(2πA/λ)` (exact for a pure
cosine; dimensionless because A and λ share the mm length unit). See the
`WrinkleProfile` class docstring in `src/wrinklefe/core/wrinkle.py` for
the full per-profile geometric definitions.

#### Through-width (transverse) variation

By default a wrinkle is treated as uniform across the full specimen width
(`transverse_mode="uniform"`). Real manufacturing wrinkles are usually
*localized* — high amplitude mid-width, fading toward the edges — so a
uniform-width assumption overstates the defect volume and biases the
knockdown conservative. Set `transverse_mode` to `"gaussian_decay"`,
`"sinusoidal_y"`, or `"elliptical"` to run the FE path with a full
`z(x, y)` wrinkle surface (`transverse_span` and `transverse_width` tune
the width envelope). At the same crest amplitude a localized wrinkle
predicts a *milder* knockdown than the uniform baseline. This is FE-only
and single-wrinkle for now (analytical-only, multi-wrinkle, and CZM runs
are rejected at construction); the CLI/app knobs are a follow-up. See
`examples/transverse_wrinkle_knockdown.py` for a runnable comparison.

### Thermal / cure-residual loading (`delta_T`)

A laminate leaves the autoclave stress-free at cure temperature and is
used at room temperature. Cooling down locks in residual stress, because
each ply's transverse (matrix) contraction is restrained by its
neighbours' stiff fibres. That residual stress superposes directly onto
the matrix tension/compression states the failure criteria evaluate, so
it is not a second-order detail: on the IM7/8552 quasi-isotropic coupon
below a `delta_T = -155` cool-down adds **+34.3 MPa** of transverse
tension to every ply — 55 % of `Yt = 62.3 MPa`.

> **Sign convention.** `delta_T` is the temperature change **from the
> stress-free (cure) state**: `delta_T = T_service − T_stress_free`.
> **A cure cool-down is NEGATIVE.** A 177 °C cure taken to 22 °C service
> is `delta_T = -155`, *not* `+155` and *not* `22`. A positive value
> means the laminate is hotter than its stress-free state. Get this
> backwards and the residual matrix stress flips from tension to
> compression — a laminate prone to cure microcracking will look safe.

| Parameter | Units | Default | Definition | Constraint |
|-----------|-------|---------|------------|------------|
| `delta_T` | °C | `0.0` | Uniform temperature change from the stress-free (cure) state, `T_service − T_stress_free`. Negative for a cure cool-down. Adds the CLT thermal resultants `N^T`/`M^T` to the ABD solve on the analytical path and the element thermal initial-strain load vector `∫ Bᵀ C ε_th dV` on the FE path; stresses are recovered from the mechanical (not total) strain on both. `0.0` (default) leaves every result bit-identical. | finite, `\|delta_T\| ≤ 1000` |

```python
config = AnalysisConfig(
    amplitude=0.5, wavelength=16.0, width=12.0,
    morphology="graded", loading="compression",
    angles=[0, 45, -45, 90, 90, -45, 45, 0],
    ply_thickness=0.183,
    delta_T=-155.0,       # 177 °C cure -> 22 °C service
)
result = WrinkleAnalysis(config).run()
print(result.summary())          # states the ΔT and its sense
print(result.failure_report.critical_mode)   # -> matrix_tension
```

```bash
wrinklefe analyze --delta-T -155 --angles "[0/45/-45/90]s"
```

**Both paths honour it (issue #273).** The analytical path adds the CLT
thermal resultants to the ABD solve; the FE path assembles the element
thermal initial-strain load vector `∫ Bᵀ C ε_th dV` and subtracts the
thermal strain during stress recovery, so `σ = C̄ (B u − ε_th)`. Because
the wrinkle rotates the fibre frame, the CTE mismatch concentrates in
exactly the elements where the failure criteria are evaluated. A flat
laminate solved this way reproduces the closed-form CLT ply stresses to
better than 0.5 %.

Two deliberate asymmetries: the pristine **retention baseline** is
solved at the same `delta_T` (so the retention factor compares like with
like), while the **measured modulus** is solved at `delta_T = 0` — a
strain-independent reaction offset from a residual load is not a
stiffness change, and folding it in would report a spurious modulus
shift.

Moisture (`LoadState.delta_C`, and the `beta1/2/3` swelling coefficients
on every material preset) is deliberately **not** exposed on
`AnalysisConfig` yet: nothing in the CLT solve consumes `delta_C`, so a
config field would be a silent no-op — the exact failure mode this
change removes for temperature.

### Tension analysis

```python
config = AnalysisConfig(
    amplitude=0.366, wavelength=16.0, width=12.0,
    morphology="stack", loading="tension",
    angles=[0, 45, 90, -45, 0, 45, -45, 0, 0, -45, 45, 0, -45, 90, 45, 0],
    ply_thickness=0.152,
)
result = WrinkleAnalysis(config).run()
print(result.summary())
```

### Graded morphology (embedded wrinkle)

```python
config = AnalysisConfig(
    amplitude=0.5, wavelength=15.0, width=11.0,
    morphology="graded", decay_floor=0.0,
    loading="compression",
)
result = WrinkleAnalysis(config).run()
print(result.analytical_knockdown)
```

### Resin pocket + progressive-damage FE (UD compression)

For a unidirectional wrinkle, tag the soft neat-epoxy lens at the crest
and load-step to ultimate with the progressive-damage solver:

```python
from wrinklefe.analysis import AnalysisConfig, WrinkleAnalysis
from wrinklefe.core.material import MaterialLibrary

lib = MaterialLibrary()
config = AnalysisConfig(
    amplitude=0.366, wavelength=16.0, width=12.0,
    morphology="graded", loading="compression",
    material=lib.get("AC318_S6C10"),
    enable_resin_pocket=True,                       # graded epoxy lens at the crest
    resin_pocket_material=lib.get("EPOXY_S6C10"),   # default if left None
    enable_progressive_damage=True,                 # load-step to ultimate
    progressive_n_increments=15,
)
result = WrinkleAnalysis(config).run()
print(result.progressive_knockdown, result.progressive_strength_MPa)
```

### Tool-flat surfaces & surface resin pockets (FE)

Parts cured against rigid tooling (or under a caul sheet) keep
**perfectly flat outer surfaces**: the fibre undulation is confined to
the interior, and where the outermost undulating ply dips away from the
flat surface the gap fills with **neat resin** — surface-visible pockets
over the wrinkle troughs, thinning to nothing over the crests. WrinkleFE's
default through-thickness decay already leaves the outer surfaces exactly
flat (for `stack`/`convex`/`concave`, or `graded` with `decay_floor=0`);
`enable_surface_resin_pockets` supplies the missing *material* — it tags
the stretched transition elements as fibre-free isotropic resin (a
stiffness hole and a matrix-cracking site where fibre-misalignment
criteria are meaningless).

```python
config = AnalysisConfig(
    amplitude=0.354, wavelength=7.4, width=3.7,
    morphology="graded", loading="compression",   # graded, decay_floor=0 ⇒ flat surface
    material=lib.get("AC318_S6C10"), angles=[0.0] * 15, ply_thickness=0.42,
    enable_surface_resin_pockets=True,             # trough pockets under the flat surface
    surface_pocket_side="both",                    # "top" | "bottom" | "both"
    resin_pocket_material=lib.get("EPOXY_S6C10"),  # default if left None
)
result = WrinkleAnalysis(config).run()
print(result.modulus_retention_global)
```

The pocket geometry is volume-conserving (the tagged resin equals the
integrated kinematic gap between the flat surface and the outermost
undulating ply) and reuses the crest-lens material plumbing, so the two
zones **compose** (per-element maximum) when both are enabled. This is an
**FE-only** effect: the closed-form analytical path keeps using fibre
angles only. A `uniform` morphology (never flat) or `graded` with
`decay_floor > 0` (wavy surface) is rejected with a message naming the
fix.

#### The `tool_flat` morphology (significant pockets)

Under the smooth-decay morphologies above (`stack`/`convex`/`concave`,
or `graded` with `decay_floor=0`) the amplitude tapers **linearly** across
the whole thickness, so the outermost undulating ply barely moves — at the
24-ply defaults the trough gap is only ~0.25 of a *single* ply thickness,
invisible in the preview and mechanically negligible. Those morphologies
model a co-cured *gradual* wrinkle, not a tooling-dominated one.

`morphology="tool_flat"` models the tooling-dominated case directly: a
uniform-amplitude core (like `uniform`), a short linear transition over
`surface_transition_plies` (default 2), and an **exactly-flat pinned
surface** on `surface_pocket_side` (`"top"`/`"bottom"`/`"both"`). The
kinematic mismatch across the transition equals the *full* amplitude, so
the trough pocket is ≈ `amplitude` deep (≈2.7 ply thicknesses at
defaults) — significant and visible. Surface pockets **auto-enable** for
`tool_flat` (they are its defining physics), and the analytical path
equals `uniform` (M_f = 1.0; the pocket effect is FE-only).

```python
config = AnalysisConfig(
    morphology="tool_flat", amplitude=0.35, wavelength=16.0, width=12.0,
    surface_pocket_side="both",       # tool-flat face(s), also the pocket side
    surface_transition_plies=3,       # amplitude ramp width (>= 1)
    material=lib.get("IM7_8552"), angles=[0.0] * 24, ply_thickness=0.183,
)
result = WrinkleAnalysis(config).run()   # pockets auto-enabled
print(result.modulus_retention_global)
```

On the crest side the `surface_transition_plies` transition elements
*compress* by `amplitude / surface_transition_plies`, so the amplitude is
bounded — `amplitude ≤ 0.8 · surface_transition_plies · ply_thickness /
nz_per_ply`; beyond it the elements would invert and construction fails
with a message naming both remedies (more transition plies, or a smaller
amplitude).

Those compressed crest-side elements are exactly what the **compaction
Vf gradient** models (issue #379): set `enable_vf_gradient=True` and the
local fibre volume fraction follows the thickness change,
`Vf_local = vf_nominal · h0/h`, so the crest band compacts and stiffens
while the trough turns resin-rich and softens — one continuous field
instead of a binary neat-resin tag:

```python
config = AnalysisConfig(
    morphology="tool_flat", amplitude=0.25, surface_pocket_side="both",
    enable_vf_gradient=True,          # opt-in; tool_flat only in v1
    vf_nominal=None,                  # None -> the card's documented Vf
    vf_max=0.75,                      # compaction cap (square packing)
    material=lib.get("IM7_8552"), angles=[0.0] * 24, ply_thickness=0.183,
    analytical_only=False,
)
```

The local card is the **preset scaled by the micromechanics ratio**
`P_micro(Vf_local) / P_micro(vf_nominal)` for the stiffnesses and CTEs, so
the trend is modelled without importing the micromechanics model's
absolute error. **Poisson ratios and all strengths stay at the preset
values** — no mixing rule predicts a strength from `Vf` (local failure
indices still move, because the local stiffness redistributes stress).
Elements compacted past `vf_max` saturate, counted in a single warning;
the rule carries no lateral resin flow along the ply. When the gradient is
on it supersedes the binary surface-pocket tag (it is its continuous
generalization); the machined crest resin lens composes unchanged.

### Penetration gate (θ, D/T, z) — UD, zero FE cost

The closed-form two-parameter gate predicts a UD knockdown directly from
geometry. Call it on its own with a calibrated preset:

```python
from wrinklefe.core.penetration_gate import penetration_gate_kd, GATE_LI2024_MOULDED

kd = penetration_gate_kd(theta_deg=8.0, dt=0.10, params=GATE_LI2024_MOULDED)
print(kd)
```

Or drive it through `AnalysisConfig.penetration_gate` so
`analytical_knockdown` (and `analytical_strength_MPa`) come from the gate
instead of Budiansky–Fleck (use `GATE_LI2025_VACBAG` with
`AC318_S6C10_vacbag` for the vacuum-bag realization):

```python
from wrinklefe.analysis import AnalysisConfig, WrinkleAnalysis
from wrinklefe.core.material import MaterialLibrary
from wrinklefe.core.penetration_gate import GATE_LI2024_MOULDED

lib = MaterialLibrary()
config = AnalysisConfig(
    amplitude=0.366, wavelength=16.0, width=12.0,
    morphology="uniform", loading="compression",
    material=lib.get("AC318_S6C10"),
    penetration_gate=GATE_LI2024_MOULDED,
)
result = WrinkleAnalysis(config).run(analytical_only=True)
print(result.analytical_knockdown)
```

When `penetration_gate` is left unset (the default `None`), the
analytical knockdown is unchanged. With a multi-wrinkle configuration
(`AnalysisConfig.wrinkles`, below) the gate evaluates each wrinkle on
its own geometry — `theta_i = arctan(2πA_i/λ_i)`, penetration
`D_i/T = A_i/T`, and the through-thickness position factor `P(z_i)`
from the spec's ply interface — and returns the weakest-link (minimum)
knockdown.

### Multi-wrinkle configurations

Real laminates often carry several wrinkles. Passing a list of
`WrinkleSpec` entries overrides the single/dual-wrinkle dispatch and
places each wrinkle at its own ply interface with its own geometry and
longitudinal position (`phase_offset` shifts a crest by
`phase·λ/2π`):

```python
import numpy as np
from wrinklefe.analysis import AnalysisConfig, WrinkleAnalysis, WrinkleSpec

config = AnalysisConfig(
    morphology="graded", loading="compression",
    angles=[0.0] * 14, ply_thickness=0.44,
    wrinkles=[
        WrinkleSpec(amplitude=0.75, wavelength=12.9, width=6.45,
                    ply_interface=6, phase_offset=-2.0 * np.pi),
        WrinkleSpec(amplitude=0.75, wavelength=12.9, width=6.45,
                    ply_interface=6, phase_offset=+2.0 * np.pi),
    ],
)
result = WrinkleAnalysis(config).run()
```

The composed displacement and fibre-angle fields feed the FE solve
("compose then differentiate"), the penetration gate scores each
wrinkle and takes the weakest link, and `enable_czm=True` inserts
cohesive surfaces along the full length of every wrinkle-nominated
interface — wrinkles sharing an interface get one continuous surface,
so a delamination can propagate crest-to-crest between neighbours
(`examples/08_multi_wrinkle_czm_linkup.py` demonstrates the link-up
vs far-separated contrast).

#### When a CZM solve fails to converge

CZM solves can fail for physical reasons — the tangent goes indefinite
near peak load, the load steps are too coarse for the damage evolution,
or the residual stagnates. When that happens the run still returns
(partial) results with `czm_converged = False`, and the solver now
classifies *why* and hands back one actionable hint. `AnalysisResults`
carries `czm_failure_diagnostics` (increment index + load fraction,
iteration count, final residual, residual history, line-search status)
and `czm_failure_hint`; the CLI prints the hint to stderr and the
Streamlit CZM section shows it as an error. The hint names the knob to
turn:

- **still decreasing at the iteration cap** → raise `max_newton_iter` or
  loosen `czm_newton_tol`;
- **stagnated / diverged** → increase `czm_n_load_increments` or reduce
  the applied strain;
- **tangent singular (post-peak snap-back)** → more load increments or a
  smaller strain; displacement/arc-length control is the robust fix on
  the roadmap.

### Uncertainty propagation (probabilistic analysis)

Measured wrinkle geometry is uncertain — amplitude and wavelength come
from a micrograph or C-scan with error. `probabilistic_analysis`
samples `AnalysisConfig` fields from user-supplied distributions
(Latin-hypercube by default, plain Monte-Carlo optional) and runs the
analytical path per sample, turning "the model says 0.64" into "P5–P95
= 0.59–0.86 given my measurement uncertainty" — the form an NCR
disposition rationale actually needs:

```python
from wrinklefe.analysis import AnalysisConfig
from wrinklefe.core.penetration_gate import GATE_LI2025_VACBAG
from wrinklefe.stochastic import probabilistic_analysis

base = AnalysisConfig(
    amplitude=0.75, wavelength=12.9, width=6.45,
    angles=[0.0] * 14, ply_thickness=0.44, morphology="graded",
    penetration_gate=GATE_LI2025_VACBAG,
)
prob = probabilistic_analysis(
    base,
    {"amplitude": ("normal", 0.75, 0.08),
     "wavelength": ("normal", 12.9, 1.0)},
    n_samples=1000, seed=42,
)
print(prob.summary())                     # P5/P50/P95, mean ± std
print(prob.knockdown_percentile(5.0))     # 5th-percentile knockdown
prob.plot()                               # histogram + sensitivity scatter
```

Distributions accept `("normal", mean, std)`, `("uniform", lo, hi)`,
`("lognormal", mu, sigma)` or any frozen `scipy.stats` distribution; a
fixed `seed` makes the whole analysis reproducible, and `n_workers`
reuses the sweep process pool for FE-path sampling. 1000 analytical
samples run in under a second for UD/gate configs (~20 s for a 24-ply
multidirectional layup, or seconds with `n_workers`).

> **Not A-/B-basis values.** The reported percentiles are
> *model-input-propagation* statistics — the deterministic model driven
> by sampled geometry. They are **not** CMH-17 A-/B-basis allowables
> (one-sided tolerance bounds on physical test data with prescribed
> confidence) and must not be presented as basis values in
> certification paperwork.

### Batch parametric sweeps

For exploring how the knockdown varies across a parameter range, use
`WrinkleAnalysis.parametric_sweep` to sweep a single
`AnalysisConfig` field (any numeric field — `amplitude`, `wavelength`,
`width`, `phase`, `applied_strain`, ...):

```python
from wrinklefe.analysis import AnalysisConfig, WrinkleAnalysis

base = AnalysisConfig(
    amplitude=0.366, wavelength=16.0, width=12.0,
    morphology="stack", loading="compression",
)
results = WrinkleAnalysis.parametric_sweep(
    base, parameter="amplitude", values=[0.1, 0.2, 0.3, 0.4],
    analytical_only=True,
)
for r in results:
    print(f"A={r.config.amplitude:.3f}  KD={r.analytical_knockdown:.4f}")
```

Every sweep point is an independent analysis, so full-FE sweeps
parallelize across processes: pass `n_workers=N` (`0` = all CPU cores)
and the solves fan out over a process pool with results returned in
the same order as `values` — measured 3.6× at 4 workers on an 8-value
FE sweep. Peak memory scales with `n_workers` × the per-solve
footprint, so size the worker count by available RAM for fine meshes.

For multi-parameter cross-product sweeps with JSON output and plots,
use `wrinklefe.sweep.run_sweep`:

```python
import numpy as np
from wrinklefe.sweep import run_sweep, save_sweep_results, plot_sweep_results

sweep = run_sweep({
    "amplitude":  np.linspace(0.183, 0.549, 3),
    "wavelength": np.linspace(8.0, 24.0, 3),
})
save_sweep_results(sweep, "./sweep_output/")
plot_sweep_results(sweep, "./sweep_output/")
```

### Command line

```bash
wrinklefe --help

# Single-parameter sweep (analytical-only is the default, fast)
wrinklefe sweep --parameter amplitude --min 0.1 --max 0.5 --steps 5

# Full-FE sweep across 4 worker processes (--parallel 0 = all cores)
wrinklefe sweep --parameter amplitude --min 0.1 --max 0.5 --steps 8 \
    --no-analytical-only --parallel 4

# Sweep over a saved config's laminate/gate (the headline UD capability):
# --config supplies the base (material, layup, gate); --parameter/--min/
# --max/--steps drive the variation over it. `compare` takes --config too.
wrinklefe sweep --config ud_gate.json --parameter amplitude \
    --min 0.2 --max 0.9 --steps 6

# Inverse: the largest amplitude that still meets a 0.85 knockdown
wrinklefe critical --parameter amplitude --target-knockdown 0.85
```

Through-width (transverse) wrinkle surfaces and the mesh/thickness knobs
are exposed on `analyze` (a non-uniform `--transverse-mode` forces the FE
path; `--ply-thickness` sets the gate's D/T):

```bash
wrinklefe analyze --transverse-mode gaussian_decay --transverse-width 5 \
    --nz-per-ply 2 --ply-thickness 0.25
```

#### Wrinkle-defect capabilities

`analyze` also exposes the newer defect models. The FE-only features
(`--resin-pocket`, `--surface-resin-pockets`, `--vf-gradient`,
`--progressive`) force the
full FE solve, so they take precedence over `--analytical-only` / `--no-fe`:

```bash
# Two-parameter (theta, D/T) penetration gate — the best UD predictor —
# selecting a calibrated preset (UD-scoped; not for multidirectional laminates)
wrinklefe analyze --gate li2025-vacbag --morphology uniform \
    --angles 0,0,0,0,0,0,0,0 --interface-1 3 --interface-2 4 --amplitude 0.5

# Off-mid-plane wrinkle position (fraction of thickness, in [0, 1])
wrinklefe analyze --morphology graded --wrinkle-z-position 0.7

# Crest resin pocket + progressive-damage ultimate strength (FE)
wrinklefe analyze --resin-pocket --progressive --increments 15

# Surface resin pockets under a tool-flat surface (FE)
wrinklefe analyze --surface-resin-pockets --surface-pocket-side both

# tool_flat morphology — significant pockets (auto-enabled). The
# hyphenated --morphology tool-flat is an alias for tool_flat; amplitude
# is bounded by 0.8 * surface-transition-plies * ply_thickness / nz_per_ply.
wrinklefe analyze --morphology tool-flat --amplitude 0.35 \
    --surface-pocket-side both --surface-transition-plies 3

# Compaction Vf / ply-thickness gradient — two caul plates (FE, tool_flat)
wrinklefe analyze --morphology tool-flat --amplitude 0.25 \
    --surface-pocket-side both --vf-gradient
```

| Flag | Config field | Notes |
| --- | --- | --- |
| `--wrinkle-z-position Z` | `wrinkle_z_position` | Fraction of thickness in `[0, 1]` (0.5 = midplane); graded morphology |
| `--gate {li2024-moulded,li2025-vacbag}` | `penetration_gate` | Calibrated `GateParameters` preset; UD-scoped |
| `--resin-pocket` | `enable_resin_pocket` | Crest resin lens (FE-only) |
| `--surface-resin-pockets` / `--surface-pocket-side {top,bottom,both}` | `enable_surface_resin_pockets` / `surface_pocket_side` | Tool-flat surface pockets (FE-only). Auto-enabled for `--morphology tool-flat` |
| `--surface-transition-plies N` | `surface_transition_plies` | `tool_flat` only: amplitude-ramp width (≥ 1, default 2); bounds the amplitude |
| `--vf-gradient` | `enable_vf_gradient` | Compaction Vf / ply-thickness gradient (FE-only, `tool_flat` only); supersedes the binary surface pockets |
| `--vf-nominal V` / `--vf-fiber F` / `--vf-matrix M` | `vf_nominal` / `vf_fiber` / `vf_matrix` | Ratio anchor and constituents; omit to use the material card's documented values |
| `--progressive` / `--increments N` | `enable_progressive_damage` / `progressive_n_increments` | Load-stepping ultimate strength (FE-only) |

The long tail of finer knobs (custom `GateParameters`, resin-pocket
geometry scales, progressive load-ramp targets) stays reachable through
`--config` (see below). `sweep` sweeps any numeric `AnalysisConfig` field,
including `wrinkle_z_position`, over its `--config` base setup:

```bash
wrinklefe sweep --parameter wrinkle_z_position --min 0.2 --max 0.8 \
    --steps 4 --morphology graded
```

### Iterative solver controls

`solver="iterative"` runs conjugate gradient with an incomplete-LU (ILU)
preconditioner — memory-efficient for large meshes (>100 k DOFs). Its
knobs are `AnalysisConfig` fields, reachable through `--config` (their
defaults reproduce the historical hardcoded values, so an existing
iterative run is unchanged):

| Config field | Default | Meaning |
| --- | --- | --- |
| `iterative_rtol` | `1e-10` | CG relative-residual convergence tolerance (> 0) |
| `iterative_maxiter` | `10000` | CG iteration cap (≥ 1) |
| `ilu_drop_tol` | `1e-4` | ILU drop tolerance — the main quality/memory knob; larger = sparser, cheaper, weaker (≥ 0) |
| `ilu_fill_factor` | `None` | Upper bound on ILU fill; `None` keeps SciPy's default (≥ 1 when set) |
| `preconditioner` | `"ilu"` | `"ilu"`, `"jacobi"` (diagonal — much lower memory for huge meshes), or `"none"` |

If ILU construction fails (out of memory, or a structurally/numerically
singular factor), the solver now emits a **`WARNING`** and falls back to
the diagonal (Jacobi) preconditioner instead of switching silently — an
ill-conditioned matrix can make that fallback orders of magnitude slower,
so the warning names the original error and points at
`preconditioner="jacobi"` (to silence it) or the direct solver. Only the
narrow set of exceptions `spilu` uses to signal a real factorisation
failure is caught; any other error propagates. On non-convergence the
`RuntimeError` names the active preconditioner, the iterations used, the
cap, and the final relative residual.

### Saving and reusing a configuration

`analyze` can persist and reload a full `AnalysisConfig`. `--save-config`
writes the *effective* configuration (after any `--config` file and CLI
overrides are applied); `--config` reloads it. Any flag given on the same
command line as `--config` overrides the file value, while flags left off
keep the file's value:

```bash
# Save the effective config to a file, then reuse it verbatim
wrinklefe analyze --amplitude 0.4 --morphology concave --save-config case.json
wrinklefe analyze --config case.json

# Reuse the file but override one parameter (0.9 wins over the file's value)
wrinklefe analyze --config case.json --amplitude 0.9
```

The same round-trip is available programmatically via
`AnalysisConfig.to_dict()` / `from_dict()` and the
`save_json` / `load_json` (and extension-dispatching `save` / `load`)
helpers. The JSON pins a `config_version` field; loading a file with an
unknown key or a mismatched version fails loudly. YAML is supported when
PyYAML is installed (it is not a required dependency). Library materials
serialise by name, custom materials inline, and penetration-gate presets
serialise by their registry name.

### Probabilistic (stochastic) analysis

`wrinklefe stochastic` propagates input distributions through the analysis
to report percentile knockdowns — the "P5 knockdown under measurement
uncertainty" number — from a `--config` base plus one or more repeatable
`--distribution FIELD:DIST:P1:P2` specs (`DIST` is `normal`, `uniform`, or
`lognormal`). A fixed `--seed` makes the whole run reproducible; results
write to JSON (percentiles + per-sample arrays) and/or a per-sample CSV.

```bash
wrinklefe stochastic --config case.json \
    --distribution amplitude:normal:0.5:0.05 \
    --distribution wavelength:uniform:12:20 \
    --n-samples 1000 --seed 42 --method lhs \
    --output-json prob.json --output-csv prob.csv
```

These are model-**input-propagation** statistics, **not** CMH-17 A-/B-basis
allowables (the printed summary carries the disclaimer).

### Finding the maximum acceptable defect

`wrinklefe critical` inverts the analysis: instead of "given this wrinkle,
what is the knockdown?", it answers "given this allowable, what is the
largest wrinkle we can accept?" — the number that goes into an inspection
criterion or an NCR disposition. Pass `--target-knockdown` (a knockdown
factor) or `--target-strength` (an absolute MPa allowable), optionally
`--parameter` (any float `AnalysisConfig` field; default `amplitude`),
`--objective`, `--bracket LO HI`, `--max-value`, `--scan-points` and
`--rtol`. `--save-config` writes the config at the critical value — the
run to attach to the disposition — and `--output-json` / `--output-csv`
carry the scan curve and the full evaluation ledger.

```bash
# Compression: largest amplitude meeting a 0.85 knockdown
wrinklefe critical --parameter amplitude --target-knockdown 0.85 \
    --save-config limit.json

# Tension, expressed as an absolute allowable instead
wrinklefe critical --loading tension --target-strength 1020
```

**The returned value satisfies the criterion under a real forward
evaluation, not merely to within the root tolerance.** `brentq` converges
to a root, not to the side of it that satisfies the inequality, so the
search backs the answer off toward safety and verifies it with an extra
solve; the printed `Criterion satisfied:` line is that check.

The search runs on the analytical path by default (about 25 forward
evaluations, well under a second; `--bracket LO HI --scan-points 3` cuts
it to about 15). `--no-analytical-only` forces the full FE pipeline, where
each evaluation is 4 linear solves — minutes, not seconds. Root-finding is
inherently sequential, so there is no `--parallel`.

Neither monotonicity nor uniqueness is assumed: the direction is measured
from a log-spaced scan of the search range, and when the target is never
crossed (`no_crossing`), never met (`target_unreachable`), reached by more
than one root (`non_monotonic` — e.g. `morphology="graded"` against
`wavelength`, whose curve is U-shaped), or unresolvable because the
parameter is inert for the configuration (`flat`), the search exits 1 with
a message quoting the measurements behind the refusal — never a scipy
traceback. "No crossing in range" is not an error: it means no defect size
in range fails your criterion.

The same search is a button in the Streamlit app — **Find acceptable
limit**, directly under *Run analysis* — which additionally carries the
limit onto the NCR validation summary when the search and the run were
made on the same configuration.

### Exporting results to CSV / JSON

Numeric outputs (load factor, per-ply failure index, knockdown factors,
stress-field summary) can be written to a schema-versioned JSON or a
Pandas-friendly per-ply CSV for downstream comparison and plotting in
Excel, Pandas, or shared Jupyter notebooks:

```python
from wrinklefe.analysis import AnalysisConfig, WrinkleAnalysis
from wrinklefe.io.results import export_results_csv, export_results_json

result = WrinkleAnalysis(AnalysisConfig()).run()

export_results_json(result, "results.json")  # schema-versioned JSON
export_results_csv(result, "per_ply.csv")    # per-ply tabular CSV
```

The JSON output is deterministic (`sort_keys=True`), pins a top-level
`schema_version` field, and reduces large numpy arrays (e.g.
per-Gauss-point stress fields) to summary statistics so the file stays
compact. The CSV is one row per ply with columns `ply_index,
angle_deg, max_FI, min_RF, critical_mode, critical_criterion`, suitable
for `pandas.read_csv` or `csv.DictReader`.

Every JSON export carries a `provenance` block recording the installed
WrinkleFE version (never a hardcoded literal), the Python/numpy/scipy
versions, the platform, a UTC timestamp, and a solver snapshot — so a
result file can be audited and reproduced against the validation
ledger. The NCR validation summary (`build_analysis_summary`) embeds
the same block, and the top-level `wrinklefe_version` field reflects
the real installed version.

The Streamlit web app exposes the same exports as **Download results as
JSON** and **Download per-ply results as CSV** buttons on the Export
tab.

### Building a ply from its constituents (fibre volume fraction)

The 11 built-in systems above are fixed cards. When the question is
*"what if my fibre volume fraction is 0.55 rather than 0.60?"*,
`wrinklefe.core.micromechanics` builds the ply from its constituents —
a fibre, a neat resin, and a Vf — instead:

```python
from wrinklefe.core.material import MaterialLibrary, OrthotropicMaterial
from wrinklefe.core.micromechanics import FIBER_PRESETS, MATRIX_PRESETS

lib = MaterialLibrary()

ply = OrthotropicMaterial.from_constituents(
    FIBER_PRESETS["IM10"],            # Hexcel HexTow IM10 carbon
    MATRIX_PRESETS["EPOXY_8552"],     # HexPly 8552 neat resin
    0.55,                             # fibre volume fraction
    strengths_from=lib.get("IM10_8552"),   # allowables carried over as-is
)
ply.E1, ply.E2, ply.G12, ply.nu12
```

The rules are the standard ones, each cited in the module docstring:
Voigt rule of mixtures for `E1`, `nu12` and `nu23`, Halpin–Tsai for `E2`
(ξ = 2) and `G12` (ξ = 1), transverse isotropy for `G23`, and Schapery
for the thermal expansion coefficients. `FIBER_PRESETS` covers the nine
fibres behind the library systems (AS4, T300, T700S, IM7, IM10, IM6G,
T800S, S-2 glass, Kevlar 49) and `MATRIX_PRESETS` the two resins with
published neat-resin data (3501-6, 8552); any isotropic card already in
the library — `EPOXY_S6C10`, say — can be used as the matrix through
`MatrixProperties.from_material`. Every constituent constant is sourced
in a comment: `E1f`/`alpha1f` from the manufacturer data sheet, the
transverse set from Daniel & Ishai (2006) Table A.2/A.3.

**Two limitations, both deliberate.**

1. **Strengths are not predicted.** No mixing rule here maps Vf to an
   allowable, so `Xt`/`Xc`/`Yt`/`Yc`/`S12`/… are carried over unchanged
   from `strengths_from` (or left at the defaults) and do *not* track
   Vf. Longitudinal strength is set by fibre-strength statistics and
   misalignment, transverse and shear strengths by the matrix and the
   fibre–matrix interface; a Vf-scaled strength model would be quietly
   wrong.
2. **The elastic predictions are approximate.** Rebuilt from
   constituents at the Vf each preset documents in its own source
   comment, the model lands within 12 % on `E1`, 26 % on `nu12`, 33 %
   on `E2` and 32 % on `G12` — except Kevlar-49/epoxy `G12`, which
   Halpin–Tsai over-predicts by 84 % (published aramid `G12f` values
   scatter by an order of magnitude; that case is a recorded `xfail`
   rather than a tuned fibre constant). `nu23` is under-predicted
   throughout. `tests/test_micromechanics.py` pins every one of those
   deviations. Use the model for the *trend* — how properties move with
   Vf, anchored on a measured ply — not as a source of allowables.

Nothing in the analysis pipeline consumes this yet: it is the
prerequisite capability for modelling resin squeeze-out under a
constrained wrinkle (issue #379), where the local Vf varies from
element to element.

## Validation

### What the in-repo tests check

The integration tests under `tests/test_integration/` exercise the full
`WrinkleAnalysis` pipeline and assert physical-sanity properties of the
analytical knockdown rather than reproducing absolute experimental
strengths:

- `test_elhajjar_validation.py` (10 tests): zero-amplitude returns
  knockdown ≈ 1, knockdown decreases monotonically with amplitude,
  morphology ordering convex > stack > concave in compression, knockdown
  stays in `(0, 1]`, and strength equals `Xc * knockdown`.
- `test_tension_validation.py` (13 tests): tension pipeline completes,
  uses `Xt` (not `Xc`), three-mechanism (`kd_fiber`, `kd_matrix`,
  `kd_oop`) decomposition is populated with the controlling mode, and
  tension knockdown is no more severe than compression for the same
  defect.

These act as regression guards on the analytical model. For the
multidirectional datasets (Elhajjar 2025, Mukhopadhyay 2015, Li et al.
2026) the repository does not ship the experimental data points, so
case-level error statistics for them are not reproducible in-repo.
(Tracking issue: #22.)

For the **unidirectional** datasets the situation is different: the
committed validation ledger
([`tests/test_validation/ledger.json`](tests/test_validation/ledger.json))
carries the digitized measured knockdowns for **Li et al. (2025)**
(Dataset F, 6 single-wrinkle S-glass compression cases) and **Hsiao &
Daniel (1996)** (Dataset G, carbon), and one command regenerates the
full per-case predicted-vs-measured table with drift detection against
pinned baselines:

```bash
python scripts/validate.py
```

For Li (2025) the ledger scores three predictors per case: the plain
Budiansky–Fleck angle floor, the closed-form **modulus** knockdown, and
the calibrated **penetration gate** — the UD strength path (issue #161)
that is sensitive to amplitude and through-thickness position
independently of the peak angle. On the S-M-2/4/5 trio (identical 20°
peak angle, amplitude 1.5/1.0/0.5 mm, measured KD 0.63/0.94/1.00 — a
~60 % strength spread invisible to any angle-only model) the gate lands
within **2.2 % / 0.6 % / 0.3 %**; over all six cases (including the
near-surface S-A-2 via the position factor) the mean absolute KD error
is **0.035**, with every case inside the ±20 % parity band. These
acceptance criteria are pinned as permanent regression tests in
[`tests/test_validation/test_ledger.py`](tests/test_validation/test_ledger.py).

### Quantitative validation against experiment

Comparison of the analytical predictions against published experimental
data is documented in the accompanying paper:

- Elhajjar, R. (2025). *Fat-tailed failure strength distributions and
  manufacturing defects in advanced composites.* Scientific Reports,
  15:25977. https://doi.org/10.1038/s41598-025-06693-4

Additional datasets referenced by the model calibration (Mukhopadhyay
et al., 2015; Li et al., 2026) are cited in [References](#references)
below. Reproducing case-level pass/fail tables for those
multidirectional datasets from this repository alone is not currently
possible — their raw data are not included.

For a consolidated predicted-vs-experimental view, the script
[`validation/plot_all_validation.py`](validation/plot_all_validation.py)
regenerates `validation/fig_all_validation_parity.png`: a single parity
plot of every single-wrinkle case (Datasets A–F) inside a ±20% band,
with each dataset predicted by the model that physically applies to it
(Budiansky–Fleck / three-mechanism for the multidirectional cases A–D,
the penetration gate for the UD cases E/F).

### Stiffness (modulus) knockdown

Besides strength, WrinkleFE reports a **stiffness** knockdown of the axial
Young's modulus two ways: the FE `modulus_retention` (wrinkled vs pristine
from the linear static solve, any layup) and — for unidirectional layups —
a closed-form `analytical_modulus_knockdown` (a CLT series-average of the
off-axis lamina modulus over the wrinkle profile, no FE solve). The script
[`validation/validate_modulus.py`](validation/validate_modulus.py) scores
both against the UD datasets that report a *measured modulus* — **F**
(Li 2025, S-glass), **G** (Hsiao & Daniel 1996, carbon — the
`IM6G_3501_6` card), and the indicative **E** (Li 2024). The analytical estimate lands at
3.9 % MAE (F) / 1.2 % (G) and the FE at 6.9 % (F) / 5.1 % (G). The data
and both models agree that stiffness is far more wrinkle-tolerant than
strength: the modulus knockdown stays ≈0.81–0.98 for the S-glass cases
and only reaches ≈0.52–0.57 for a carbon uniform wrinkle at θ = 15°. The
script
[`validation/plot_modulus_validation.py`](validation/plot_modulus_validation.py)
renders the comparison as `validation/fig_modulus_validation.png` —
knockdown-vs-angle and a predicted-vs-experimental parity plot across all
three datasets.

## Supported morphologies

WrinkleFE ships five wrinkle morphologies (defined in
`src/wrinklefe/core/morphology.py`). They differ along *two* independent
axes: **how many wrinkles** are placed in the laminate and **how the
amplitude varies through the thickness**. The first three names below
are *dual-wrinkle* modes distinguished by the phase offset φ between
two adjacent wrinkle centrelines (the through-thickness amplitude
follows a linear taper from the wrinkle interface plies down to zero
at the laminate outer surfaces). The last three are *single-wrinkle*
modes that swap that taper for a different through-thickness profile.

| Morphology   | # wrinkles | Phase φ | Through-thickness amplitude            | M_f (compression) | When to use                                                                                                                                          |
|--------------|------------|---------|-----------------------------------------|-------------------|------------------------------------------------------------------------------------------------------------------------------------------------------|
| `stack`      | 2          | 0       | Linear decay, 1 at interface → 0 at surfaces | 1.0 (baseline)    | Two aligned wrinkles, peaks-over-peaks. The dual-wrinkle reference case used to scale `convex` / `concave`.                                          |
| `convex`     | 2          | +π/2    | Linear decay, 1 at interface → 0 at surfaces | < 1               | Two phase-shifted wrinkles whose interface bulges outward. *Least* damaging dual-wrinkle case in compression.                                       |
| `concave`    | 2          | −π/2    | Linear decay, 1 at interface → 0 at surfaces | > 1               | Two phase-shifted wrinkles whose interface pinches inward. *Most* damaging dual-wrinkle case in compression — design-driving.                       |
| `uniform`    | 1          | n/a     | Full amplitude on **every** ply (no decay)   | 1.0 (no pairing)  | A single through-thickness-wide wrinkle — every ply wavy with the same A. Conservative bound and sanity-check baseline.                              |
| `graded`     | 1          | n/a     | Linear decay from mid-ply to surfaces, with floor `decay_floor` ∈ [0, 1] | 1.0 (no pairing) | An embedded wrinkle that fades toward the surface plies. `decay_floor=0` is pure graded; `decay_floor=1` collapses to `uniform`.                    |
| `tool_flat`  | 1          | n/a     | Full-amplitude core, linear ramp over `surface_transition_plies`, **flat** at the pinned surface | 1.0 (no pairing) | A wrinkle cured against rigid tooling / a caul sheet: the pinned-flat surface fills the wrinkle troughs with **significant** resin pockets (auto-enabled). Shares `uniform`'s analytical knockdown; amplitude is bounded (see below). |

### `stack` vs `uniform` — what's the difference?

These two get conflated because both have `M_f = 1.0`, but they model
very different defects:

- **`stack`** places **two** wrinkles at adjacent interfaces with
  φ = 0 (aligned crests). Through the thickness the wrinkle decays
  linearly from the interface plies to zero at the outer surfaces —
  surface plies are flat.
- **`uniform`** places a **single** wrinkle and disables the
  through-thickness decay — every ply, including the outer surfaces,
  is displaced by the full profile.

For the same `amplitude` / `wavelength`, `apply_to_nodes` therefore
produces *different* deformed meshes: `stack` has a wrinkle
concentrated near the interface plies (and flat top/bottom plies),
while `uniform` has a wrinkle of the same amplitude at every single
ply. The `M_f = 1.0` coincidence is purely the analytical knockdown
parameter — the FE geometry, the per-ply fibre-angle field, and the
predicted ply-by-ply failure are not the same.

## Supported failure criteria

The criteria below live in `src/wrinklefe/failure/` and can be selected
through `FailureEvaluator` or used independently:

- **LaRC04/05** (`larc05.py`) — Pinho/Camanho 3-D criterion with
  fibre-kinking under compression, in-situ matrix strengths, and a
  fracture-plane search. Default for the FE solve.
- **Tsai-Wu** (`tsai_wu.py`) — 3-D tensor-polynomial criterion with a
  configurable interaction coefficient.
- **Tsai-Hill** (`tsai_hill.py`) — 3-D extension of the classical
  quadratic Tsai-Hill index.
- **Hashin** (`hashin.py`) — 3-D Hashin criterion with separate
  fibre-tension/-compression and matrix-tension/-compression modes.
- **Puck** (`puck.py`) — action-plane (Mode A/B/C) inter-fibre-failure
  criterion with simplified fibre failure.
- **Maximum Stress** (`max_stress.py`) and **Maximum Strain**
  (`max_strain.py`) — non-interactive checks against the principal
  material-frame allowables.
- **Budiansky-Fleck kink-band** (`kinkband.py`) — analytical compression
  knockdown with an optional interlaminar damage coupling
  (`InterlaminarDamage`); this is the model exposed in the
  `analytical_knockdown` field of `AnalysisResults`.
- **Progressive damage** (`progressive.py`) — `PlyDiscount` and
  `ContinuumDamage` post-failure stiffness reduction models that wrap
  any of the criteria above.

## How It Works

The full mechanism-by-mechanism derivation lives on the
[Theory: physics & mechanics](https://wrinklefe.readthedocs.io/en/latest/theory.html)
page (`docs/theory.md`). The essentials:

### Wrinkle kinematics

A wrinkle is reduced to its peak fibre-misalignment angle

```
theta_max = arctan(2*pi*A / lambda)
```

(half-amplitude `A`, wavelength `λ`) and — for unidirectional laminates —
its through-thickness penetration `D/T = A/T` (`T` = laminate thickness).
`theta_eff = M_f * theta_max` folds in the morphology factor `M_f`
(`stack` = 1, `convex` < 1, `concave` > 1).

### Compression — CLT-weighted Budiansky–Fleck kink-band

```
KD_lam      = f_0 * KD_BF + (1 - f_0)
KD_BF       = 1 / (1 + r + c_AF * r^2),   r = theta_eff / gamma_Y_eff
gamma_Y_eff = max(0.032 + 0.050 * f_confined
                        - 0.010 * max(n_block_max - 1, 0),  0.016)
```

`f_0` is the axial-stiffness fraction carried by the 0° plies (the plies
that kink); the `(1 - f_0)` term is the off-axis plies riding through at
full strength. The matrix shear-yield strain `gamma_Y_eff` **rises** with
the confinement `f_confined` (off-axis neighbours bracing the 0° plies
against kink-band rotation) and **falls** with the longest run of
consecutive 0° plies `n_block_max` (blocked 0° plies kink more easily),
floored at half the UD value. So a dispersed `[0/45/90/-45]s` resists
wrinkle knockdown far better than a blocked `[0_4/90_4]s`. The optional
Argon–Fleck quadratic term `c_AF` (`kink_band_quadratic_coeff`) defaults
to `0` — the pure linear Budiansky–Fleck floor.

### Tension — three-mechanism minimum, CLT-weighted

```
KD_lam = f_0 * min(cos^2(theta), KD_matrix, KD_oop) + (1 - f_0)
```

The 0° ply knockdown is the *most severe* of three competing mechanisms:
fibre load-rotation (`cos^2(theta)`), in-situ matrix cracking (a
Hashin/LaRC `σ22`–`τ12` interaction with a thick-ply in-situ strength
correction, `KD_matrix`), and a curved-beam out-of-plane delamination
check (`KD_oop`: the wrinkle curvature drives an interlaminar `σ33` at the
crest and `τ13` at the flanks). A Benzeggagh–Kenane mixed-mode
delamination-*onset* knockdown is reported alongside, and the tension
knockdown is floored by the compression value for the same defect
("tension is never worse than compression").

### Graded morphology

For the `graded` morphology the knockdown is averaged over the wrinkle
profile in both the longitudinal (`x`) and through-thickness (`z`)
directions. The **compression** path weights each ply by a Gaussian
through-thickness envelope centred at `wrinkle_z_position` (decay scale
`max(λ/2, A)`); the **tension** path uses an analogous linear taper. In
both, `decay_floor` sets the surface-ply amplitude: `0` is a fully
embedded wrinkle that fades to flat at the surfaces, `1` collapses to
`uniform`.

## References

- Elhajjar, R. (2025). Scientific Reports, 15:25977.
- Li, Y. et al. (2026). Composites Part A, 205:109719.
- Li, X. et al. (2024). Composites Science and Technology, 256:110762.
- Li, Y. et al. (2025). Polymer Composites, 46:15176-15187.
- Hsiao, H.M. & Daniel, I.M. (1996). Composites Science and Technology, 56(5), 581-593.
- Budiansky, B. & Fleck, N.A. (1993). J. Mech. Phys. Solids, 41(1), 183-211.
- Pinho, S.T. et al. (2005). NASA-TM-2005-213530.
- Camanho, P.P. et al. (2006). Composites Part A, 37(2), 165-176.
- Jin, L. et al. (2026). Thin-Walled Structures, 219:114237.

## License

MIT - see [LICENSE](LICENSE)

## Changelog

Notable changes between versions — including any that shift predictions,
flagged under a **Numerical results** heading — are recorded in
[CHANGELOG.md](CHANGELOG.md).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md)

## Citation

If you use WrinkleFE in your research, please cite it. The quickest way is the
**"Cite this repository"** button on the
[GitHub page](https://github.com/elhajjar1/wrinkleFE) — it's generated from
[`CITATION.cff`](CITATION.cff) and exports APA or BibTeX. The full software
citation:

> Elhajjar, R. (2025). WrinkleFE: An open-source finite element package for strength prediction of wrinkled composite laminates (Version 1.0.0) [Computer software]. University of Wisconsin-Milwaukee. https://github.com/elhajjar1/wrinkleFE

```bibtex
@software{elhajjar2025wrinklefe,
  author = {Elhajjar, Rani},
  title = {{WrinkleFE}: An Open-Source Finite Element Package for Strength
           Prediction of Wrinkled Composite Laminates},
  year = {2025},
  version = {1.0.0},
  publisher = {GitHub},
  url = {https://github.com/elhajjar1/wrinkleFE},
  note = {University of Wisconsin-Milwaukee}
}
```

> **Software DOI:** WrinkleFE does not have one yet. The DOI badge above is
> the *article's*, not the software's. A citable software DOI will be added
> here — and to `CITATION.cff` — once Zenodo archiving is enabled for the
> repository, at which point every tagged release is archived automatically
> and a concept DOI resolves to the latest one
> ([issue #284](https://github.com/elhajjar1/wrinkleFE/issues/284)). Until
> then the URL above is the citation target; please cite the version you
> actually used.

Please also cite the underlying experimental validation data:

> Elhajjar, R. (2025). Fat-tailed failure strength distributions and manufacturing defects in advanced composites. *Scientific Reports*, 15, 25977. https://doi.org/10.1038/s41598-025-06693-4
