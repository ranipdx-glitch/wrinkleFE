"""Tests for the CLT LoadState -> 3-D boundary-condition conversion.

``BoundaryHandler.load_state_to_bcs`` is the single definition of this
mapping; ``StaticSolver._load_state_to_bcs`` delegates to it. (It used to
carry a second, divergent implementation — issue #96 — which clamped the
whole ``x_min`` face and silently ignored ``Ny``/``Mx``/``My``.)

**Membrane states are applied as self-equilibrated tractions on all four
in-plane faces**, with rigid-body motion removed at three points. That
matters and is not cosmetic: the previous per-resultant scheme loaded only
``x_max`` and clamped or symmetry-fixed the opposing faces, which turns a
uniform state into a cantilever reaction. Its ``Nxy`` shear came out 37 %
high, and ``Nx + Nxy`` came out with the *wrong sign* on the shear strain,
because the ``symmetry_y`` plane the ``Nx`` branch added restrains exactly
the deformation the shear traction is trying to produce.

The tests that let that through asserted only the *shape* of the BC list —
which face, which DOF, what total force, and that force vectors superpose.
None of them solved. So the classes below are in two layers:

1. **Structure**: per-face totals, the traction pairing, the anchor set.
2. **Physics** (``TestMembraneAgainstCLT``): solve a flat laminate and
   compare the recovered strains against ``Laminate.midplane_strains``.
   That layer is what actually pins the mapping; it is independent of how
   the BCs are spelled.
"""

import numpy as np
import pytest

from wrinklefe.core.laminate import Laminate, LoadState
from wrinklefe.core.material import OrthotropicMaterial
from wrinklefe.core.mesh import WrinkleMesh
from wrinklefe.solver.boundary import BoundaryHandler
from wrinklefe.solver.static import StaticSolver

# ======================================================================
# Fixtures (mirror tests/test_solver/test_static.py conventions)
# ======================================================================

@pytest.fixture
def x850_material():
    """Default CYCOM X850/T800 material."""
    return OrthotropicMaterial()


@pytest.fixture
def single_ply_laminate(x850_material):
    """Single-ply [0] laminate with 0.183 mm thickness."""
    return Laminate.from_angles([0.0], material=x850_material,
                                ply_thickness=0.183)


@pytest.fixture
def two_ply_laminate(x850_material):
    """Two-ply [0/0] laminate so the mesh has a true midplane node row."""
    return Laminate.from_angles([0.0, 0.0], material=x850_material,
                                ply_thickness=0.183)


@pytest.fixture
def small_mesh(single_ply_laminate):
    """Small 3x2x1 mesh, domain (Lx=3, Ly=2, Lz=0.183)."""
    gen = WrinkleMesh(
        laminate=single_ply_laminate,
        wrinkle_config=None,
        Lx=3.0, Ly=2.0,
        nx=3, ny=2, nz_per_ply=1,
    )
    return gen.generate()


@pytest.fixture
def bending_mesh(two_ply_laminate):
    """3x2x2 mesh (nz=2) so x_max has an exact midplane node (z=0)."""
    gen = WrinkleMesh(
        laminate=two_ply_laminate,
        wrinkle_config=None,
        Lx=3.0, Ly=2.0,
        nx=3, ny=2, nz_per_ply=1,
    )
    return gen.generate()


# ======================================================================
# Helpers
# ======================================================================

def _force_components(mesh, bcs):
    """Return (sum Fx, sum Fy, sum Fz) of the assembled global force."""
    F = BoundaryHandler(mesh).get_force_dofs(bcs)
    return F[0::3].sum(), F[1::3].sum(), F[2::3].sum()


def _bc_kinds(bcs):
    return [b.bc_type for b in bcs]


def _node_fix_bcs(bcs):
    """BCs that pin explicit node_ids (the rigid-body anchors)."""
    return [b for b in bcs if b.node_ids is not None and b.bc_type == "fixed"]


def _face_force(mesh, bcs, face, dof):
    """Total applied force on one face in one direction.

    The membrane load set is self-equilibrated, so the *net* force over the
    whole body is zero by construction; the per-face total is the quantity
    that has to equal the CLT resultant times the loaded edge length.
    """
    sel = [b for b in bcs
           if b.bc_type == "pressure" and b.face == face and b.dofs == [dof]]
    return BoundaryHandler(mesh).get_force_dofs(sel)[dof::3].sum()


def _anchor_dofs(bcs):
    """The rigid-body anchor pattern as a sorted list of dof-lists."""
    return sorted(sorted(b.dofs) for b in _node_fix_bcs(bcs))


# ======================================================================
# Zero load / thermal-only -> no mechanical BCs
# ======================================================================

class TestNoLoad:
    """LoadState with no mechanical resultant maps to an empty BC list."""

    def test_zero_load_returns_empty(self, small_mesh):
        bcs = BoundaryHandler.load_state_to_bcs(LoadState(), small_mesh)
        assert bcs == []

    def test_thermal_only_returns_empty(self, small_mesh):
        """delta_T-only LoadState is mechanical-only here -> []. (#133:
        thermal expansion is handled through CLT, not these BCs.)"""
        bcs = BoundaryHandler.load_state_to_bcs(
            LoadState(delta_T=-100.0), small_mesh
        )
        assert bcs == []

    def test_zero_load_no_applied_force(self, small_mesh):
        bcs = BoundaryHandler.load_state_to_bcs(LoadState(), small_mesh)
        fx, fy, fz = _force_components(small_mesh, bcs)
        assert (fx, fy, fz) == (0.0, 0.0, 0.0)


# ======================================================================
# Structure: rigid-body anchors
# ======================================================================

class TestRigidBodyAnchors:
    """Six point constraints, and nothing that restrains a strain."""

    def test_anchor_pattern(self, small_mesh):
        bcs = BoundaryHandler.load_state_to_bcs(LoadState(Nx=100.0), small_mesh)
        assert _anchor_dofs(bcs) == [[0, 1, 2], [1, 2], [2]]

    def test_no_face_fixes_at_all(self, small_mesh):
        """A face fix would restrain the deformation being applied.

        This is the assertion that would have prevented the shear defect:
        the old scheme's ``symmetry_y`` on ``y_min`` is precisely what made
        ``Nx + Nxy`` come out with the wrong sign.
        """
        bcs = BoundaryHandler.load_state_to_bcs(
            LoadState(Nx=100.0, Ny=50.0, Nxy=30.0), small_mesh
        )
        for bc in bcs:
            assert bc.bc_type not in ("symmetry_x", "symmetry_y", "symmetry_z")
            assert not (bc.bc_type == "fixed" and bc.face is not None), (
                f"membrane BCs must not fix a whole face, got {bc.face}"
            )

    def test_anchors_are_three_distinct_nodes(self, small_mesh):
        bcs = BoundaryHandler.load_state_to_bcs(LoadState(Nxy=30.0), small_mesh)
        ids = {int(b.node_ids[0]) for b in _node_fix_bcs(bcs)}
        assert len(ids) == 3

    def test_degenerate_footprint_is_rejected(self, single_ply_laminate):
        """A mesh with no y extent cannot host the anchor set."""
        mesh = WrinkleMesh(
            laminate=single_ply_laminate, wrinkle_config=None,
            Lx=3.0, Ly=2.0, nx=3, ny=2, nz_per_ply=1,
        ).generate()
        mesh.nodes[:, 1] = 0.0          # collapse the y extent
        with pytest.raises(ValueError, match="rigid-body"):
            BoundaryHandler.load_state_to_bcs(LoadState(Nx=1.0), mesh)


# ======================================================================
# Structure: traction pairing and per-face totals
# ======================================================================

class TestMembraneTractions:
    """Every resultant loads BOTH opposing faces, equal and opposite."""

    def test_nx_loads_both_x_faces(self, small_mesh):
        Lx, Ly, _ = small_mesh.domain_size
        bcs = BoundaryHandler.load_state_to_bcs(LoadState(Nx=100.0), small_mesh)
        assert _face_force(small_mesh, bcs, "x_max", 0) == pytest.approx(100.0 * Ly)
        assert _face_force(small_mesh, bcs, "x_min", 0) == pytest.approx(-100.0 * Ly)

    def test_ny_loads_both_y_faces(self, small_mesh):
        Lx, Ly, _ = small_mesh.domain_size
        bcs = BoundaryHandler.load_state_to_bcs(LoadState(Ny=50.0), small_mesh)
        assert _face_force(small_mesh, bcs, "y_max", 1) == pytest.approx(50.0 * Lx)
        assert _face_force(small_mesh, bcs, "y_min", 1) == pytest.approx(-50.0 * Lx)

    def test_nxy_is_a_complementary_shear_pair(self, small_mesh):
        """Shear acts on the y-faces too — that is what makes it uniform.

        Loading only the x-faces (the old behaviour) is a tip-loaded
        cantilever, and gave a shear strain 37 % above the CLT value.
        """
        Lx, Ly, _ = small_mesh.domain_size
        bcs = BoundaryHandler.load_state_to_bcs(LoadState(Nxy=30.0), small_mesh)
        assert _face_force(small_mesh, bcs, "x_max", 1) == pytest.approx(30.0 * Ly)
        assert _face_force(small_mesh, bcs, "x_min", 1) == pytest.approx(-30.0 * Ly)
        assert _face_force(small_mesh, bcs, "y_max", 0) == pytest.approx(30.0 * Lx)
        assert _face_force(small_mesh, bcs, "y_min", 0) == pytest.approx(-30.0 * Lx)

    def test_load_set_is_self_equilibrated(self, small_mesh):
        """No net force on the body — the anchors carry no reaction."""
        bcs = BoundaryHandler.load_state_to_bcs(
            LoadState(Nx=100.0, Ny=50.0, Nxy=30.0), small_mesh
        )
        fx, fy, fz = _force_components(small_mesh, bcs)
        assert fx == pytest.approx(0.0, abs=1e-9)
        assert fy == pytest.approx(0.0, abs=1e-9)
        assert fz == pytest.approx(0.0, abs=1e-9)

    def test_force_sign_follows_resultant_sign(self, small_mesh):
        Lx, Ly, _ = small_mesh.domain_size
        bcs = BoundaryHandler.load_state_to_bcs(LoadState(Nx=-250.0), small_mesh)
        assert _face_force(small_mesh, bcs, "x_max", 0) == pytest.approx(-250.0 * Ly)

    def test_combined_state_superposes(self, small_mesh):
        """The combined force vector is the sum of the individual ones."""
        handler = BoundaryHandler(small_mesh)
        ls = LoadState(Nx=100.0, Ny=50.0, Nxy=30.0)
        combined = handler.get_force_dofs(
            BoundaryHandler.load_state_to_bcs(ls, small_mesh))
        parts = sum(
            handler.get_force_dofs(
                BoundaryHandler.load_state_to_bcs(one, small_mesh))
            for one in (LoadState(Nx=100.0), LoadState(Ny=50.0),
                        LoadState(Nxy=30.0))
        )
        np.testing.assert_allclose(combined, parts, rtol=1e-12, atol=1e-12)


# ======================================================================
# Physics: the layer that actually pins the mapping
# ======================================================================

@pytest.fixture(scope="module")
def flat_qi_mesh_and_laminate():
    """Flat quasi-isotropic laminate, fine enough to compare against CLT."""
    lam = Laminate.from_angles(
        [0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0],
        material=OrthotropicMaterial(), ply_thickness=0.125,
    )
    mesh = WrinkleMesh(
        laminate=lam, wrinkle_config=None,
        Lx=20.0, Ly=10.0, nx=10, ny=4, nz_per_ply=1,
    ).generate()
    return mesh, lam


class TestMembraneAgainstCLT:
    """Solve a flat laminate and compare with ``midplane_strains``.

    A flat laminate under a uniform membrane state has a closed-form CLT
    answer, and a correct BC mapping must reproduce it. Strains are
    averaged over interior elements, away from the traction faces where the
    discrete load introduces a boundary layer.

    Measured agreement with the corrected mapping (and, for contrast, what
    the previous per-resultant scheme gave):

    ==================  ==============  ==============
    state               corrected       previous
    ==================  ==============  ==============
    uniaxial ``Nx``     0.4 %           0.4 %
    biaxial             0.3 %           0.3 %
    pure shear          **0.4 %**       **37 %**
    compression+shear   **0.2 %**       **wrong sign**
    ==================  ==============  ==============
    """

    TOL = 0.06   # 6 % — a physics tolerance, not machine precision

    @staticmethod
    def _solve(mesh, lam, load):
        bcs = BoundaryHandler.load_state_to_bcs(load, mesh)
        res = StaticSolver(mesh, lam).solve(bcs, solver="direct")
        centres = res.element_centers
        Lx, Ly, _ = mesh.domain_size
        interior = (
            (np.abs(centres[:, 0] - Lx / 2.0) < 0.30 * Lx)
            & (np.abs(centres[:, 1] - Ly / 2.0) < 0.35 * Ly)
        )
        assert interior.any()
        eps = res.strain_global[interior].mean(axis=(0, 1))
        # Voigt [11, 22, 33, 23, 13, 12] -> CLT [eps_x, eps_y, gamma_xy]
        return np.array([eps[0], eps[1], eps[5]])

    def _check(self, mesh, lam, load):
        fe = self._solve(mesh, lam, load)
        clt = lam.midplane_strains(load)[:3]
        scale = np.abs(clt).max()
        assert scale > 1e-6, "test load must produce a real strain"
        for i, name in enumerate(("eps_x", "eps_y", "gamma_xy")):
            assert abs(fe[i] - clt[i]) < self.TOL * scale, (
                f"{name}: FE {fe[i]:.6e} vs CLT {clt[i]:.6e} "
                f"(tolerance {self.TOL:.0%} of {scale:.3e})"
            )

    def test_uniaxial(self, flat_qi_mesh_and_laminate):
        mesh, lam = flat_qi_mesh_and_laminate
        self._check(mesh, lam, LoadState(Nx=-500.0))

    def test_biaxial(self, flat_qi_mesh_and_laminate):
        mesh, lam = flat_qi_mesh_and_laminate
        self._check(mesh, lam, LoadState(Nx=-500.0, Ny=-250.0))

    def test_pure_shear(self, flat_qi_mesh_and_laminate):
        """The case the previous mapping got 37 % wrong."""
        mesh, lam = flat_qi_mesh_and_laminate
        self._check(mesh, lam, LoadState(Nxy=200.0))

    def test_compression_plus_shear(self, flat_qi_mesh_and_laminate):
        """The case the previous mapping got the *sign* wrong on."""
        mesh, lam = flat_qi_mesh_and_laminate
        self._check(mesh, lam, LoadState(Nx=-500.0, Nxy=200.0))

    def test_full_membrane_state(self, flat_qi_mesh_and_laminate):
        mesh, lam = flat_qi_mesh_and_laminate
        self._check(mesh, lam, LoadState(Nx=-500.0, Ny=-250.0, Nxy=200.0))

    def test_shear_strain_sign_follows_nxy(self, flat_qi_mesh_and_laminate):
        """Reversing Nxy reverses the shear strain — and nothing else."""
        mesh, lam = flat_qi_mesh_and_laminate
        pos = self._solve(mesh, lam, LoadState(Nx=-500.0, Nxy=+200.0))
        neg = self._solve(mesh, lam, LoadState(Nx=-500.0, Nxy=-200.0))
        assert pos[2] > 0.0 and neg[2] < 0.0
        assert pos[2] == pytest.approx(-neg[2], rel=0.02)
        # The axial response is essentially unchanged by the shear sign:
        # this layup is balanced, so A16 = A26 = 0 and CLT decouples shear
        # from axial exactly. The FE interior average still differs by
        # ~0.5 % between the two signs — a discretization effect from the
        # traction boundary layer, which the finite averaging window does
        # not sample symmetrically under sign reversal, not real coupling.
        assert pos[0] == pytest.approx(neg[0], rel=0.01)


# ======================================================================
# Pure Mx bending
# ======================================================================

class TestPureMx:
    """Mx -> linear through-thickness ux prescribed on every x_max node."""

    def test_one_displacement_bc_per_xmax_node(self, bending_mesh):
        bcs = BoundaryHandler.load_state_to_bcs(LoadState(Mx=2.0), bending_mesh)
        disp = [b for b in bcs if b.bc_type == "displacement"]
        xmax = bending_mesh.nodes_on_face("x_max")
        assert len(disp) == len(xmax)
        assert all(b.dofs == [0] for b in disp)

    def test_linear_through_thickness_profile(self, bending_mesh,
                                              two_ply_laminate):
        """ux(z) = kappa_x * (z - z_mid) * Lx with the physically-correct
        CLT contract kappa_x = Mx / D11 (see #149).  Midplane node ux == 0;
        top & bottom fibers equal and opposite."""
        Lx, Ly, Lz = bending_mesh.domain_size
        Mx = 2.0
        D11 = float(two_ply_laminate.D[0, 0])
        bcs = BoundaryHandler.load_state_to_bcs(LoadState(Mx=Mx), bending_mesh)
        disp = {int(b.node_ids[0]): b.value
                for b in bcs if b.bc_type == "displacement"}
        xmax = bending_mesh.nodes_on_face("x_max")
        z = bending_mesh.nodes[xmax, 2]
        z_mid = 0.5 * (z.min() + z.max())
        for nid in xmax:
            zz = float(bending_mesh.nodes[nid, 2])
            assert disp[nid] == pytest.approx((Mx / D11) * (zz - z_mid) * Lx)
        mid = [nid for nid in xmax
               if abs(float(bending_mesh.nodes[nid, 2]) - z_mid) < 1e-9]
        assert mid, "expected an exact midplane node row on the bending mesh"
        for nid in mid:
            assert disp[nid] == pytest.approx(0.0)
        top = max(xmax, key=lambda n: bending_mesh.nodes[n, 2])
        bot = min(xmax, key=lambda n: bending_mesh.nodes[n, 2])
        assert disp[top] == pytest.approx(-disp[bot])

    def test_curvature_scales_with_D11(self, bending_mesh, two_ply_laminate):
        """Physically kappa_x must be Mx / D11 (fix for #149)."""
        Lx, Ly, Lz = bending_mesh.domain_size
        Mx = 2.0
        D11 = float(two_ply_laminate.D[0, 0])
        bcs = BoundaryHandler.load_state_to_bcs(LoadState(Mx=Mx), bending_mesh)
        disp = {int(b.node_ids[0]): b.value
                for b in bcs if b.bc_type == "displacement"}
        xmax = bending_mesh.nodes_on_face("x_max")
        z = bending_mesh.nodes[xmax, 2]
        z_mid = 0.5 * (z.min() + z.max())
        top = max(xmax, key=lambda n: bending_mesh.nodes[n, 2])
        z_top = float(bending_mesh.nodes[top, 2])
        assert disp[top] == pytest.approx((Mx / D11) * (z_top - z_mid) * Lx)

    def test_bending_uses_the_same_anchor_set(self, bending_mesh):
        bcs = BoundaryHandler.load_state_to_bcs(LoadState(Mx=2.0), bending_mesh)
        assert _anchor_dofs(bcs) == [[0, 1, 2], [1, 2], [2]]
        for bc in bcs:
            assert bc.bc_type not in ("symmetry_x", "symmetry_y", "symmetry_z")


# ======================================================================
# Membrane + curvature is refused, not silently superposed
# ======================================================================

class TestMembranePlusBendingRejected:
    """Prescribed displacement and traction on the same face don't mix.

    The curvature terms prescribe ``ux`` on every ``x_max`` node while a
    membrane state applies traction to that same face. Superposing them
    describes neither load, so the mapping refuses rather than returning a
    BC list whose solve means nothing.
    """

    @pytest.mark.parametrize("load", [
        LoadState(Nx=100.0, Mx=2.0),
        LoadState(Ny=50.0, My=2.0),
        LoadState(Nxy=30.0, Mx=2.0),
    ])
    def test_combination_raises(self, bending_mesh, load):
        with pytest.raises(ValueError, match="membrane.*curvature"):
            BoundaryHandler.load_state_to_bcs(load, bending_mesh)

    def test_message_names_both_groups(self, bending_mesh):
        with pytest.raises(ValueError) as exc:
            BoundaryHandler.load_state_to_bcs(
                LoadState(Nx=100.0, Mx=2.0), bending_mesh)
        msg = str(exc.value)
        assert "Nx/Ny/Nxy" in msg and "Mx/My" in msg


# ======================================================================
# StaticSolver._load_state_to_bcs now delegates (was a second converter)
# ======================================================================

class TestStaticSolverConverterDelegates:
    """Issue #96's "two parallel converters" are now one.

    The private converter used to clamp the whole ``x_min`` face and read
    only ``Nx``/``Nxy`` — so a state carrying ``Ny``, ``Mx`` or ``My`` was
    solved as though those components were zero, silently. It now delegates,
    so there is nothing left to diverge.
    """

    @pytest.mark.parametrize("load", [
        LoadState(Nx=100.0),
        LoadState(Nxy=30.0),
        LoadState(Nx=100.0, Ny=50.0, Nxy=30.0),
        LoadState(),
    ])
    def test_identical_bc_lists(self, small_mesh, single_ply_laminate, load):
        solver = StaticSolver(small_mesh, single_ply_laminate)
        handler = BoundaryHandler(small_mesh)
        ss = solver._load_state_to_bcs(load)
        bh = BoundaryHandler.load_state_to_bcs(load, small_mesh)
        assert handler.get_constrained_dofs(ss) == handler.get_constrained_dofs(bh)
        np.testing.assert_allclose(
            handler.get_force_dofs(ss), handler.get_force_dofs(bh),
            rtol=1e-12, atol=1e-12,
        )

    def test_ny_is_no_longer_silently_dropped(self, small_mesh,
                                              single_ply_laminate):
        """The regression that motivated removing the duplicate."""
        Lx, Ly, _ = small_mesh.domain_size
        solver = StaticSolver(small_mesh, single_ply_laminate)
        bcs = solver._load_state_to_bcs(LoadState(Ny=50.0))
        assert _face_force(small_mesh, bcs, "y_max", 1) == pytest.approx(
            50.0 * Lx
        )

    def test_zero_load_gives_no_bcs(self, small_mesh, single_ply_laminate):
        solver = StaticSolver(small_mesh, single_ply_laminate)
        assert solver._load_state_to_bcs(LoadState()) == []
