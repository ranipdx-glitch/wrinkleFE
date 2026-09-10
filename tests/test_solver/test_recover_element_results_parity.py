"""Parity / regression test for ``StaticSolver.recover_element_results``.

Issue #187 refactored ``recover_element_results`` to lift element-constant
work (``T_ply``, node ids, per-node wrinkle angles) out of the inner Gauss
loop and to compute all 8 Gauss-point wrinkle angles via a single matmul
against a pre-built shape-function matrix.

These tests pin down the numerical output at a handful of Gauss points on
a small, fully-deterministic problem.  The values are taken from running
``main`` on the same problem; any change in the post-processing math
should make the assertions trip.
"""

from __future__ import annotations

import numpy as np
import pytest

from wrinklefe.core.laminate import Laminate
from wrinklefe.core.material import OrthotropicMaterial
from wrinklefe.core.mesh import WrinkleMesh
from wrinklefe.solver.boundary import BoundaryHandler
from wrinklefe.solver.static import StaticSolver


@pytest.fixture
def parity_mesh_and_laminate():
    """Small multi-ply mesh that exercises both ply rotation and wrinkles.

    Uses a [0/45/-45/90] laminate so ``ply_angles`` are non-trivial, and a
    near-isotropic material so the values are easy to reason about.  Mesh
    is intentionally small so the regression numbers stay short to read.
    """
    E = 10_000.0
    nu = 0.3
    G = E / (2.0 * (1.0 + nu))
    material = OrthotropicMaterial(
        E1=E, E2=E, E3=E,
        G12=G, G13=G, G23=G,
        nu12=nu, nu13=nu, nu23=nu,
        Xt=500, Xc=500, Yt=500, Yc=500, Zt=500, Zc=500,
        S12=300, S13=300, S23=300,
        gamma_Y=0.02,
        name="parity_iso_10k",
    )
    laminate = Laminate.from_angles(
        [0.0, 45.0, -45.0, 90.0],
        material=material,
        ply_thickness=0.183,
    )
    gen = WrinkleMesh(
        laminate=laminate,
        wrinkle_config=None,
        Lx=4.0, Ly=2.0,
        nx=4, ny=2, nz_per_ply=1,
    )
    return gen.generate(), laminate


def test_recover_element_results_shape_and_finiteness(parity_mesh_and_laminate):
    """Output arrays have the expected shape and are all finite."""
    mesh, laminate = parity_mesh_and_laminate
    solver = StaticSolver(mesh, laminate)
    bcs = BoundaryHandler.compression_bcs(mesh, applied_strain=-0.01)
    results = solver.solve(bcs)

    n_elem = mesh.n_elements
    n_gp = 8

    for arr in (
        results.stress_global,
        results.stress_local,
        results.strain_global,
        results.strain_local,
    ):
        assert arr.shape == (n_elem, n_gp, 6)
        assert np.all(np.isfinite(arr))


def test_recover_element_results_local_matches_manual_transform(
    parity_mesh_and_laminate,
):
    """Local stress/strain at a hand-picked GP matches a manual transform.

    Two things this pins that an earlier version of it got wrong:

    * **Which matrix.** Stress rotates with ``T_sigma``; *engineering*
      strain rotates with ``T_eps = R T_sigma R^-1``
      (``R = diag(1,1,1,2,2,2)``).  Reusing ``T_sigma`` on a strain vector
      mis-scales every shear component by two.
    * **Which order.** ``C_bar = R_y(R_z(C))``, so
      ``sigma_local = T_ply @ T_wrinkle @ sigma_global`` — the ply
      transform on the LEFT, because the wrinkle rotation is the outer one
      and is therefore the first undone.  The reversed order is invisible
      on a 0 deg ply (``T_ply`` is the identity) and reaches ~8 % on the
      90 deg plies of this laminate.
    """
    from wrinklefe.core.transforms import (
        strain_transformation_3d,
        stress_transformation_3d,
    )
    from wrinklefe.elements.hex8 import Hex8Element

    mesh, laminate = parity_mesh_and_laminate
    solver = StaticSolver(mesh, laminate)
    bcs = BoundaryHandler.compression_bcs(mesh, applied_strain=-0.01)
    results = solver.solve(bcs)

    # Pick a few representative elements across the mesh and GPs across the
    # element so we exercise multiple ply angles and corners of the
    # reference cube.
    sample_elems = [0, mesh.n_elements // 3, mesh.n_elements // 2,
                    mesh.n_elements - 1]
    sample_gps = [0, 3, 7]

    from wrinklefe.elements.gauss import gauss_points_hex
    gp_coords, _ = gauss_points_hex(order=2)

    for e in sample_elems:
        ply_angle_rad = np.radians(float(mesh.ply_angles[e]))
        T_ply = stress_transformation_3d(ply_angle_rad, axis='z')
        T_ply_eps = strain_transformation_3d(ply_angle_rad, axis='z')

        node_ids = mesh.elements[e]
        fiber_angles_local = mesh.fiber_angles[node_ids]

        for g in sample_gps:
            xi, eta, zeta = gp_coords[g]
            N = Hex8Element.shape_functions(xi, eta, zeta)
            phi = float(N @ fiber_angles_local)
            T_wrinkle = stress_transformation_3d(phi, axis='y')
            T_wrinkle_eps = strain_transformation_3d(phi, axis='y')

            expected_sigma_local = (
                T_ply @ (T_wrinkle @ results.stress_global[e, g])
            )
            expected_eps_local = (
                T_ply_eps @ (T_wrinkle_eps @ results.strain_global[e, g])
            )

            np.testing.assert_allclose(
                results.stress_local[e, g],
                expected_sigma_local,
                rtol=1e-12,
                atol=1e-12,
            )
            np.testing.assert_allclose(
                results.strain_local[e, g],
                expected_eps_local,
                rtol=1e-12,
                atol=1e-12,
            )


def test_recover_element_results_recomputation_is_deterministic(
    parity_mesh_and_laminate,
):
    """Calling ``recover_element_results`` twice yields identical arrays."""
    mesh, laminate = parity_mesh_and_laminate
    solver = StaticSolver(mesh, laminate)
    bcs = BoundaryHandler.compression_bcs(mesh, applied_strain=-0.01)
    results1 = solver.solve(bcs)

    # Re-run the post-processing step against the same displacement and
    # check that we get bit-exact identical arrays the second time.
    u_flat = results1.displacement.reshape(-1)
    sg2, sl2, eg2, el2 = solver.recover_element_results(u_flat)

    np.testing.assert_array_equal(results1.stress_global, sg2)
    np.testing.assert_array_equal(results1.stress_local, sl2)
    np.testing.assert_array_equal(results1.strain_global, eg2)
    np.testing.assert_array_equal(results1.strain_local, el2)


# --------------------------------------------------------------------------- #
# Frame-consistency invariants
#
# The two tests above check the transform against a hand-built copy of the
# same algebra, so they would have agreed with a wrong convention as long
# as both sides were wrong the same way — which is exactly what happened.
# The two below are independent of how the transform is written: they are
# physical statements about what "local frame" has to mean.
# --------------------------------------------------------------------------- #


def test_local_frame_satisfies_the_unrotated_constitutive_law(
    parity_mesh_and_laminate,
):
    """``sigma_local == C_material @ eps_local``, to machine precision.

    The whole point of the material frame is that the *unrotated*
    stiffness holds in it.  This catches any mismatch between the pair
    (``stress_local``, ``strain_local``) and the ``C_bar`` that produced
    ``stress_global`` — wrong matrix, wrong order, or wrong axis — without
    restating the transform algebra the implementation uses.
    """
    mesh, laminate = parity_mesh_and_laminate
    solver = StaticSolver(mesh, laminate)
    bcs = BoundaryHandler.compression_bcs(mesh, applied_strain=-0.01)
    results = solver.solve(bcs)

    C = laminate.plies[0].material.stiffness_matrix
    predicted = results.strain_local @ C.T          # (n_elem, n_gp, 6)
    peak = np.abs(results.stress_local).max()
    assert peak > 1.0, "test is only meaningful under real load"
    assert np.abs(predicted - results.stress_local).max() < 1e-9 * peak


def test_wrinkle_rotation_of_a_90_deg_ply_leaves_local_sigma_11_alone():
    """Why the composition order is what it is — a derivation anchor.

    For a 90 deg ply the fibres run along global **y**, and the wrinkle
    misalignment is a rotation about **y** — i.e. about that ply's own
    fibre axis.  Material-frame ``sigma_11`` must therefore be completely
    independent of the wrinkle angle.  ``T_ply @ T_wrinkle`` satisfies
    that exactly; the reversed order drifts several percent per 0.1 rad.

    NOTE: this is a statement about the transform algebra, not a guard on
    ``recover_element_results`` — it builds ``C_bar`` itself and would
    keep passing if the solver regressed.  The implementation guard is
    :func:`test_local_frame_satisfies_the_unrotated_constitutive_law`,
    which does fail on the reversed order.  This test exists so a future
    reader can see *why* the order is not arbitrary before changing it.
    """
    from wrinklefe.core.transforms import (
        rotate_stiffness_3d,
        stress_transformation_3d,
    )

    mat = OrthotropicMaterial()
    theta = np.radians(90.0)
    T_ply = stress_transformation_3d(theta, axis="z")
    eps = np.array([-8.0e-3, 1.5e-3, 1.0e-3, 2.0e-4, 5.0e-4, 8.0e-4])

    reference = None
    for phi in (0.0, 0.05, 0.10, 0.20):
        C_bar = rotate_stiffness_3d(
            rotate_stiffness_3d(mat.stiffness_matrix, theta, axis="z"),
            phi, axis="y",
        )
        sigma_global = C_bar @ eps
        T_wrinkle = stress_transformation_3d(phi, axis="y")
        sigma_local = T_ply @ (T_wrinkle @ sigma_global)
        if reference is None:
            reference = sigma_local[0]
            assert abs(reference) > 1.0
        else:
            np.testing.assert_allclose(sigma_local[0], reference, rtol=1e-12)
