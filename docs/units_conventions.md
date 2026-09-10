# Units & conventions

WrinkleFE carries no unit system of its own: every input is a plain
float in the units below, and every output is in the units implied by
them. The rule of thumb is **mm / MPa / N**, which makes moduli and
stresses MPa, force resultants N/mm and fracture toughnesses N/mm.

Every row below is sourced from the code or docstring that defines it;
that definition, not this page, is authoritative.

## Quantities

| Quantity | Unit | Convention | Defined in |
|----------|------|-----------|------------|
| Lengths (`amplitude`, `wavelength`, `width`, `domain_length`, `domain_width`, `ply_thickness`) | mm | "All lengths are millimetres" | `AnalysisConfig`, [getting started](getting_started.md) |
| Wrinkle amplitude *A* | mm | **Half**-amplitude: peak displacement of the wrinkled mid-surface from the flat reference, so `z(x) = A·cos(2π(x−x₀)/λ)` and the peak-to-trough height is `2A`. For a measured wrinkle, `A = (z_max − z_min) / 2`. | `AnalysisConfig.amplitude` |
| Wavelength *λ* | mm | Spatial period of the cosine carrier — the crest-to-crest distance along *x*. Wavenumber `k = 2π/λ` (1/mm). Must be > 0. | `AnalysisConfig.wavelength` |
| Envelope width *w* | mm | Longitudinal envelope decay length about the wrinkle centre `x₀`: the Gaussian 1/e length in `exp(−(x−x₀)²/w²)`. Also the transverse (*y*) extent in the 3-D dual-wrinkle / graded mesh. Must be > 0. | `AnalysisConfig.width` |
| Young's and shear moduli (`E1`, `E2`, `E3`, `G12`, `G13`, `G23`) | MPa | Fibre (1), in-plane transverse (2), through-thickness (3) directions. | `OrthotropicMaterial` |
| Strength allowables (`Xt`, `Xc`, `Yt`, `Yc`, `Zt`, `Zc`, `S12`, `S13`, `S23`) | MPa | Longitudinal / transverse / through-thickness tensile and compressive strengths, and shear strengths. Stored as **positive magnitudes**. | `OrthotropicMaterial` |
| Poisson's ratios (`nu12`, `nu13`, `nu23`) | — | Dimensionless. The reciprocal ratios (`nu21`, `nu31`, `nu32`) are derived from compliance-matrix symmetry, not supplied. | `OrthotropicMaterial` |
| Fracture toughnesses (`GIc`, `GIIc`, `czm_GIc`, `czm_GIIc`) | N/mm | Mode-I / Mode-II interlaminar fracture toughness. Numerically identical to kJ/m² (1 N/mm = 1 kJ/m²), so a card quoted in kJ/m² needs no conversion. | `OrthotropicMaterial`, `AnalysisConfig` |
| Cohesive peak tractions (`sigma_max`, `tau_max`, `czm_sigma_max`, `czm_tau_max`) | MPa | Mode-I peak normal traction and Mode-II peak shear traction on the ply interface. | `OrthotropicMaterial`, `AnalysisConfig` |
| Cohesive penalty stiffness (`czm_penalty`) | N/mm³ | Initial interface stiffness of the bilinear traction–separation law. | `AnalysisConfig.czm_penalty` |
| Cohesive energy (`czm_energy_dissipated`) | N·mm | Total dissipated cohesive energy across all interfaces. | `AnalysisResults` |
| Strain (`applied_strain`, `progressive_max_strain`, `gamma_Y`) | — | A **fraction**, not a percentage: the default `applied_strain = -0.01` is 1 % compression. | `AnalysisConfig.applied_strain` |
| Ply and misalignment angles | degrees | See [Angles](#angles) below. | `core/layup.py`, `AnalysisConfig.angles` |
| Phase offset (`phase`, `WrinkleSpec.phase_offset`) | radians | Dual-wrinkle phase φ between the two wrinkle centrelines (`stack` = 0, `convex` = +π/2, `concave` = −π/2). | `AnalysisConfig.phase`, `MORPHOLOGY_PHASES` |
| Fibre volume fraction (`vf_nominal`, `vf_max`, `Vf_local`) | — | A **fraction** in `[0, 1]`, not a percentage. `vf_max` defaults to 0.75 (just under square packing). | `core/micromechanics.py`, `core/compaction.py` |
| Through-thickness position (`wrinkle_z_position`) | — | A **fraction of the laminate thickness *T***: `0.0` = bottom surface, `0.5` = midplane (default), `1.0` = top surface. Must be in `[0, 1]`. | `AnalysisConfig.wrinkle_z_position` |
| Penetration `D/T` | — | Dimensionless through-thickness penetration `D/T = A/T` used by the penetration gate. | `core/penetration_gate.py`, [getting started](getting_started.md) |
| Decay floor (`decay_floor`) | — | Dimensionless fraction in `[0, 1]`: the minimum fraction of the wrinkle amplitude retained at the laminate outer surfaces. `0.0` = full decay to zero, `1.0` = no decay. | `AnalysisConfig.decay_floor` |
| Knockdowns / retentions / morphology factor | — | Dimensionless ratios; `1.0` = no loss relative to the pristine laminate. | `AnalysisResults` |
| Damage index *D* | — | Dimensionless, `0` = pristine, `1` = full loss of load-carrying capacity (capped at 0.999). Reporting only. | `AnalysisResults.damage_index` |
| Force resultants (`Nx`, `Ny`, `Nxy`, `Qx`, `Qy`) | N/mm | Per unit width. | `LoadState` |
| Moment resultants (`Mx`, `My`, `Mxy`) | N·mm/mm | Per unit width. | `LoadState` |
| Coefficients of thermal expansion (`alpha1`, `alpha2`, `alpha3`) | 1/K | Per kelvin; a temperature *change* is numerically the same in K and °C. | `OrthotropicMaterial` |
| Temperature change (`delta_T`) | °C | Uniform temperature change **from the stress-free (cure) state**: `delta_T = T_service − T_stress_free`. **A cure cool-down is negative** — a 177 °C cure taken to 22 °C service is `delta_T = -155`, not `+155` and not `22`. A positive value means the laminate is hotter than its stress-free state. Numerically identical in K and °C because it is a *change*. On `AnalysisConfig` the value must satisfy `\|delta_T\| ≤ 1000`. Both solution paths honour it (issue #273): CLT thermal resultants on the analytical path, the element thermal initial-strain load vector on the FE path. | `AnalysisConfig.delta_T`, `LoadState.delta_T` |
| Coefficients of moisture expansion (`beta1`, `beta2`, `beta3`) | 1/%M | Hygroscopic swelling per percent moisture. | `OrthotropicMaterial` |
| Moisture change (`delta_C`) | % | Uniform moisture concentration change. | `LoadState.delta_C` |
| Nonlinear shear coefficient (`beta_shear`) | 1/MPa³ | Ramberg–Osgood coefficient in the LaRC05 shear response. | `OrthotropicMaterial.beta_shear` |

## Angles

- **Ply angles** (`AnalysisConfig.angles`, `--angles` / `--layup`) are in
  **degrees**. The canonical range is `|angle| <= 90`; `0`, `±90`,
  decimals and negatives are all valid, and anything outside the range
  is rejected at construction by `validate_ply_angle` — the single
  source of truth shared by the layup parser and `AnalysisConfig`.
  Contracted notation (`[0/45/-45/90]_3s`, `[0/±45/90]s`) is accepted
  wherever an explicit list is.
- **Fibre misalignment angles** are stored on `AnalysisResults` in
  **radians** — `max_angle_rad` (peak misalignment θ_max),
  `effective_angle_rad` (effective angle θ_eff), `mesh_max_angle_rad`
  (peak angle read off the FE mesh, so it accounts for the
  through-thickness decay). `AnalysisResults.summary()` and the
  Streamlit app both *display* them in degrees; the attributes
  themselves stay in radians.
- **Penetration-gate inputs** are in degrees: `penetration_gate_kd`
  takes `theta_deg`.
- **`alpha_0`** — the LaRC04/05 fracture-plane angle under pure
  transverse compression — is in degrees (default 53.0).
- **Dual-wrinkle phase** is the one angular quantity in radians on the
  input side (see the table above).

## Coordinate system

The mesh is a structured hexahedral grid with

- **x** — the longitudinal direction: the wrinkle carrier runs along it
  (`z(x) = A·cos(2πx/λ)`), the domain spans `domain_length`, and it is
  the **load direction**. `modulus_retention_global` sums the axial
  reaction on the loaded (`x_max`) face; `applied_strain` is the nominal
  strain along *x*.
- **y** — the transverse (width) direction, spanning `domain_width`.
  The through-width envelope modes (`transverse_mode`) vary the wrinkle
  amplitude along it.
- **z** — the through-thickness direction. Mesh z-coordinates align with
  ply boundaries; `wrinkle_z_position` places the wrinkle along it.

Node ordering runs *i* (x) fastest, then *j* (y), then *k* (z), so the
node index for grid point `(i, j, k)` is
`k·(ny+1)·(nx+1) + j·(nx+1) + i`. Hex8 element node numbering follows
the VTK/Abaqus convention (bottom face counter-clockwise, then top face
counter-clockwise).

The **material frame** is `1` = fibre direction, `2` = in-plane
transverse, `3` = through-thickness. Voigt ordering for exported stress
components is `[σ₁₁, σ₂₂, σ₃₃, τ₂₃, τ₁₃, τ₁₂]`, with engineering shear
strains in the compliance matrix.

## Signs

- **`applied_strain` is signed**: negative is compression, positive is
  tension. The default is `-0.01`.
- **CLT force resultants are signed the same way**: `LoadState(Nx=-1000.0)`
  is uniaxial compression.
- **Strength allowables are unsigned magnitudes.** `Xc`, `Yc`, `Zc` are
  stored positive even though they are compressive strengths; the sign
  convention lives on the load, not on the allowable.
- **`loading`** (`'compression'` / `'tension'`) selects the physics
  (kink-band vs the three-mechanism tension model). It is a separate
  switch from the sign of `applied_strain` — set both consistently.
- **`delta_T` is a change from the stress-free (cure) state, not an
  absolute temperature.** `delta_T = T_service − T_stress_free`, so a
  cure cool-down is **negative**: a 177 °C cure taken to 22 °C service
  is `delta_T = -155` — not `+155`, and not `22`. A positive value means
  the laminate is *hotter* than its stress-free state. This is the
  single easiest sign to invert, and inverting it flips the residual
  matrix stress of a cross-ply from tension to compression, which is the
  difference between predicting cure microcracking and missing it. On
  `AnalysisConfig` the field drives **both** solution paths (issue
  #273): the CLT thermal resultants on the analytical path, and the
  element thermal initial-strain load vector `∫ Bᵀ C ε_th dV` on the FE
  path. The pristine retention baseline is solved at the same
  `delta_T`; the measured modulus is deliberately solved at
  `delta_T = 0`, because a strain-independent reaction offset from a
  residual load is not a stiffness change.

## Percent vs fraction: the one place they differ

Strain is a fraction everywhere in the library and the CLI
(`applied_strain=-0.01`, `--strain -0.01`). The **Streamlit app is the
exception**: its sidebar input is *Applied strain magnitude [%]* and it
converts on the way in (`applied_strain = applied_strain_pct / 100.0`).
The app takes the magnitude only —

> Magnitude only — the sign is taken from the loading mode
> (compression → negative, tension → positive). Editing this value is
> preserved when you toggle the loading mode.

So `1.0` in the app with loading = compression is `applied_strain =
-0.01` in the library. Fibre volume fraction is a fraction in both.

## See also

- [Getting started](getting_started.md) — installation and a first run.
- [Interpreting results](interpreting_results.md) — what each output
  number means.
- The {class}`~wrinklefe.analysis.AnalysisConfig` API page for the
  per-field reference, including defaults and validation bounds.
