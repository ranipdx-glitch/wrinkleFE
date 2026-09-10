"""FE solve for multi-wrinkle configurations (issue #252).

Stage 1: N wrinkles with disjoint longitudinal supports (the Li 2025
specimen layout). Stage 2: overlapping/interacting wrinkles — both the
displacement and the fiber-angle field derive from the composed wrinkle
field ("compose then differentiate"), so two coincident half-amplitude
wrinkles reproduce the single full-amplitude result exactly. The CZM
pathway is covered separately in ``test_multi_wrinkle_czm.py``
(issue #283).
"""

from __future__ import annotations

import numpy as np
import pytest

from wrinklefe.analysis import AnalysisConfig, WrinkleAnalysis, WrinkleSpec

pytestmark = pytest.mark.integration

_ANGLES_8 = [0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0]


def _two_wrinkle_config(**overrides):
    """Two identical wrinkles centred at x = L/2 ± 8 mm (disjoint supports).

    Each spec: lambda = 8, width = 2 -> support = center ± 6 mm.
    phase_offset = ±2*pi shifts the centre by ±lambda = ±8 mm.
    """
    spec = dict(amplitude=0.15, wavelength=8.0, width=2.0, ply_interface=3)
    defaults = dict(
        amplitude=0.15, wavelength=8.0, width=2.0,
        morphology="graded", loading="compression",
        angles=list(_ANGLES_8),
        wrinkles=[
            WrinkleSpec(**spec, phase_offset=-2.0 * np.pi),
            WrinkleSpec(**spec, phase_offset=+2.0 * np.pi),
        ],
        nx=32, ny=2, nz_per_ply=1,
    )
    defaults.update(overrides)
    return AnalysisConfig(**defaults)


def _element_x_centroids(mesh) -> np.ndarray:
    return mesh.nodes[mesh.elements].mean(axis=1)[:, 0]


class TestMultiWrinkleFE:

    def test_two_wrinkles_fe_runs_with_local_fi_peaks(self):
        """FE solve runs for 2 non-overlapping wrinkles and produces a
        local max-FI peak near each wrinkle centre."""
        cfg = _two_wrinkle_config()
        result = WrinkleAnalysis(cfg).run()

        assert result.field_results is not None
        assert result.failure_indices, "failure evaluation must run"

        # Combined per-element max FI across criteria and Gauss points.
        fi = np.max(
            np.stack([
                np.asarray(arr).max(axis=1)
                for arr in result.failure_indices.values()
            ]),
            axis=0,
        )
        x_c = _element_x_centroids(result.mesh)
        L = cfg.domain_length
        centers = (L / 2.0 - 8.0, L / 2.0 + 8.0)

        # The global FI maximum must sit inside one of the two wrinkle
        # supports (the wrinkles, not the far field, govern).
        x_peak = x_c[int(np.argmax(fi))]
        assert any(abs(x_peak - c) < 6.0 for c in centers), (
            f"global FI peak at x={x_peak:.1f} mm is outside both "
            f"wrinkle supports {centers}"
        )

        # Each wrinkle produces a local elevation above the far field.
        # (LaRC05 carries a uniform baseline FI under compression, so
        # the peaks are elevations on that baseline, not zero-to-one.)
        #
        # MARGIN, RE-PINNED: this threshold was 1.1, calibrated while
        # ``recover_element_results`` composed the local-frame transform in
        # the wrong order (ply and wrinkle factors reversed, and the stress
        # matrix reused for engineering strain).  That defect inflated the
        # in-wrinkle stresses of off-axis plies -- and this case is the
        # worst one for it: a [0/45/-45/90]s layup with the wrinkle at the
        # 90/90 interface, where the error peaked at ~8 %.  Measured on this
        # exact config, the fix moves the peak/far ratios
        #     1.1882 -> 1.0652   (x = 6 mm)
        #     1.1740 -> 1.0722   (x = 22 mm)
        # while the far-field baseline is unchanged to 5e-8 (those elements
        # are pristine, so the transform order cannot touch them).  The
        # elevation the test exists to detect is still unambiguous; only the
        # inflated margin is gone.  1.05 keeps the assertion meaningful with
        # room for mesh/solver noise.
        far = np.ones_like(x_c, dtype=bool)
        for center in centers:
            far &= np.abs(x_c - center) > 6.0
        assert far.any()
        peaks = []
        for center in centers:
            near = np.abs(x_c - center) < 4.0      # within lambda/2
            assert near.any()
            assert fi[near].max() > 1.05 * fi[far].max(), (
                f"expected a local FI peak near x={center:.1f} mm"
            )
            peaks.append(fi[near].max())
        # Identical wrinkles -> the two local peaks agree closely.
        assert peaks[0] == pytest.approx(peaks[1], rel=0.05)

    def test_domain_sized_from_union_of_extents(self):
        cfg = _two_wrinkle_config()
        # half-span = |dx| + 3*width = 8 + 6 = 14 -> L = 28 (> 3*lambda).
        assert cfg.domain_length == pytest.approx(28.0)

    def test_single_spec_matches_scalar_graded_path(self):
        """A one-entry wrinkles list is tolerance-equal to the scalar
        graded config it denotes (same placement, same mesh density,
        same domain)."""
        common = dict(
            amplitude=0.15, wavelength=16.0, width=8.0,
            morphology="graded", loading="compression",
            angles=list(_ANGLES_8),
            nx=16, ny=2, nz_per_ply=1,
            domain_length=48.0,
        )
        scalar_cfg = AnalysisConfig(**common, interface_1=3)
        multi_cfg = AnalysisConfig(
            **common,
            wrinkles=[
                WrinkleSpec(
                    amplitude=0.15, wavelength=16.0, width=8.0,
                    ply_interface=3, phase_offset=0.0,
                )
            ],
        )
        r_scalar = WrinkleAnalysis(scalar_cfg).run()
        r_multi = WrinkleAnalysis(multi_cfg).run()

        # modulus_retention is an FP-reduced scalar (mean sigma_11 over the
        # mesh); the scalar-graded and one-entry-wrinkles paths reduce in a
        # different order, agreeing only to ~1e-8, so rel=1e-9 flakes across
        # platforms (e.g. macOS arm64). 1e-6 is well within any physical
        # meaning; the FI fields below stay strict.
        assert r_multi.modulus_retention == pytest.approx(
            r_scalar.modulus_retention, rel=1e-6
        )
        for crit, arr in r_scalar.failure_indices.items():
            np.testing.assert_allclose(
                r_multi.failure_indices[crit], arr, rtol=1e-9,
                err_msg=f"{crit} FI field diverged between the scalar "
                "graded path and the one-entry wrinkles path",
            )

    def test_coincident_half_wrinkles_equal_full_amplitude(self):
        """Issue #252 Stage 2 regression: two coincident half-amplitude
        wrinkles reproduce the single full-amplitude result — exact
        under compose-then-differentiate, where displacement and fiber
        angle both derive from the composed field."""
        common = dict(
            wavelength=16.0, width=8.0,
            morphology="graded", loading="compression",
            angles=list(_ANGLES_8),
            nx=16, ny=2, nz_per_ply=1,
            domain_length=48.0,
        )
        spec = dict(wavelength=16.0, width=8.0, ply_interface=3,
                    phase_offset=0.0)
        full = AnalysisConfig(
            **common, amplitude=0.15,
            wrinkles=[WrinkleSpec(amplitude=0.15, **spec)],
        )
        halves = AnalysisConfig(
            **common, amplitude=0.15,
            wrinkles=[
                WrinkleSpec(amplitude=0.075, **spec),
                WrinkleSpec(amplitude=0.075, **spec),
            ],
        )
        r_full = WrinkleAnalysis(full).run()
        r_halves = WrinkleAnalysis(halves).run()

        np.testing.assert_allclose(
            r_halves.mesh.nodes, r_full.mesh.nodes, rtol=0, atol=1e-12,
            err_msg="coincident halves must deform the mesh identically",
        )
        np.testing.assert_allclose(
            r_halves.mesh.fiber_angles, r_full.mesh.fiber_angles,
            rtol=0, atol=1e-12,
            err_msg="coincident halves must produce the same angle field",
        )
        # FP-reduced scalar (see test above): the full vs coincident-halves
        # paths reduce sigma_11 in a different order and agree only to ~1e-8,
        # even though the mesh nodes and angle field are bit-identical
        # (asserted at 1e-12 above). rel=1e-9 flakes on macOS arm64; 1e-6 is
        # far inside any physical tolerance.
        assert r_halves.modulus_retention == pytest.approx(
            r_full.modulus_retention, rel=1e-6
        )
        # The FI field, like the scalar modulus_retention above, is a
        # post-solve quantity that the two-spec composition reduces in a
        # different summation order than the single full-amplitude spec.
        # The mesh nodes and angle field are bit-identical (1e-12 above),
        # so the inputs match; only this FP-reduction order differs, by
        # ~7e-6 on macOS arm64 — rel=1e-9 flakes there. atol=5e-5 is far
        # inside any physical FI tolerance yet survives the platform noise.
        for crit, arr in r_full.failure_indices.items():
            np.testing.assert_allclose(
                r_halves.failure_indices[crit], arr, rtol=1e-5, atol=5e-5,
                err_msg=f"{crit} FI field diverged between coincident "
                "halves and the full-amplitude wrinkle",
            )

    def test_overlapping_wrinkles_run_through_fe(self):
        """Stage 2: longitudinally overlapping wrinkles solve end-to-end."""
        cfg = _two_wrinkle_config()
        cfg.wrinkles[1].phase_offset = 0.5 * np.pi  # dx = +2 mm: overlap
        r = WrinkleAnalysis(cfg).run()
        assert r.field_results is not None
        assert r.failure_indices
        assert 0.0 < r.modulus_retention <= 1.0

    def test_overlapping_wrinkles_still_run_analytically(self):
        cfg = _two_wrinkle_config()
        cfg.wrinkles[1].phase_offset = 0.5 * np.pi
        r = WrinkleAnalysis(cfg).run(analytical_only=True)
        assert 0.0 < r.analytical_knockdown <= 1.0

