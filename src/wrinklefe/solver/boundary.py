"""Boundary condition handling and CLT-to-3D load mapping.

Provides:

- :class:`BoundaryCondition` — a dataclass describing a single BC
  (fixed, displacement, force, pressure, or symmetry).
- :class:`BoundaryHandler` — resolves BCs against a mesh and applies
  them to the global system ``K u = F`` via the penalty method or the
  elimination (partitioning) method.
- Convenience functions for common loading scenarios (uniaxial
  compression, pure bending) and a mapping from CLT
  :class:`~wrinklefe.core.laminate.LoadState` to 3-D BCs.

References
----------
Bathe, K.-J. (2006). Finite Element Procedures, Chapter 4.
Cook, R.D. et al. (2002). Concepts and Applications of Finite Element
    Analysis, 4th ed., Chapter 2.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import sparse

from wrinklefe.core.laminate import LoadState
from wrinklefe.core.mesh import MeshData

# ======================================================================
# Penalty method for displacement BCs
# ======================================================================

# Textbook penalty scaling (Bathe, Finite Element Procedures, sec. 4.2.2):
# alpha = _PENALTY_SCALE * max(|diag(K)|).  For typical FE stiffness diagonals
# (~1e3 - 1e6 N/mm), this puts alpha in the 1e11 - 1e14 range, which is
# large enough to enforce u_i ~= u_prescribed to ~8 significant digits but
# small enough to keep cond(K) well below the float64 limit (~1e16).  A
# fixed mega-penalty like 1e20 silently destroys accuracy in the
# unconstrained DOFs through ill-conditioning.
_PENALTY_SCALE = 1.0e8


def apply_penalty_bcs(
    K: sparse.csc_matrix,
    F: np.ndarray,
    constrained_dofs: dict[int, float],
    *,
    in_place: bool = False,
    penalty: float | None = None,
) -> tuple[sparse.csc_matrix, np.ndarray]:
    """Apply displacement boundary conditions via the penalty method.

    For each constrained DOF *i* with prescribed value *u_i*:

    .. math::

        K_{ii} \\mathrel{+}= \\alpha, \\qquad F_i \\mathrel{+}= \\alpha \\, u_i

    where the penalty factor scales with the problem,
    :math:`\\alpha = \\text{scale} \\cdot \\max(|\\mathrm{diag}(K)|, 1)`,
    with ``scale = _PENALTY_SCALE = 1e8`` (Bathe sec. 4.2.2).  This
    preserves matrix symmetry and sparsity, and keeps cond(K) bounded.

    Parameters
    ----------
    K : scipy.sparse.csc_matrix
        Global stiffness matrix.
    F : np.ndarray
        Global force vector (length ``n_dof``).
    constrained_dofs : dict[int, float]
        Mapping ``{dof_index: prescribed_value}``.
    in_place : bool, keyword-only, optional
        If ``True``, mutate ``K`` and ``F`` in place (faster, no copy).
        Default ``False``: both inputs are copied, so the caller's data
        is preserved.  Pass ``in_place=True`` only when the caller owns
        ``K`` and ``F`` exclusively (e.g. inside ``StaticSolver.solve``).
    penalty : float or None, optional
        Override the computed penalty value.  Normally ``None`` (the
        default), in which case ``alpha`` is derived from
        ``max(|diag(K)|)`` per the convention above.  Provided as an
        escape hatch for legacy callers and tests; prefer the default.

    Returns
    -------
    K_modified : scipy.sparse.csc_matrix
        Stiffness matrix with diagonal augmented at constrained DOFs
        (same object as input if ``in_place=True``, else a copy).
    F_modified : np.ndarray
        Force vector with penalty contribution added at constrained DOFs
        (same object as input if ``in_place=True``, else a copy).

    Notes
    -----
    Uses ``scipy.sparse.diags`` plus a sparse add to inject the diagonal
    contribution — no LIL round-trip is needed.  This is both faster and
    avoids changing the sparsity pattern outside the diagonal.
    """
    if not constrained_dofs:
        if in_place:
            return K, F
        return K.copy(), F.copy()

    if penalty is None:
        diag_max = float(np.abs(K.diagonal()).max())
        alpha = _PENALTY_SCALE * max(diag_max, 1.0)
    else:
        alpha = float(penalty)

    n_dof = K.shape[0]
    dofs = np.fromiter(constrained_dofs.keys(), dtype=np.intp,
                        count=len(constrained_dofs))
    vals = np.fromiter(constrained_dofs.values(), dtype=np.float64,
                        count=len(constrained_dofs))

    # Build a sparse diagonal contribution: alpha at each constrained DOF.
    diag_data = np.zeros(n_dof, dtype=np.float64)
    diag_data[dofs] = alpha
    K_pen = sparse.diags(diag_data, 0, format="csc")

    if in_place:
        # K is csc; (csc + csc) -> csc, but it returns a new object.  We
        # cannot truly modify K's data array in place without touching
        # its sparsity pattern, so re-bind the local name and let the
        # caller pick up the returned reference.
        K_out = (K + K_pen).tocsc()
        F[dofs] += alpha * vals
        return K_out, F

    K_out = (K + K_pen).tocsc()
    F_out = F.copy()
    F_out[dofs] += alpha * vals
    return K_out, F_out


# ======================================================================
# Internal helpers
# ======================================================================

def _quad_areas(nodes: np.ndarray, quads: np.ndarray) -> np.ndarray:
    """Return the area of each Q4 quadrilateral, by triangle split.

    Splits each quad ``(n0, n1, n2, n3)`` along the ``n0-n2`` diagonal
    into two triangles and sums ``|d1 x d2| / 2`` for each.  This is
    exact for planar quads and a sensible bilinear approximation for
    mildly non-planar (e.g. wrinkled) faces.

    Parameters
    ----------
    nodes : np.ndarray
        Shape ``(n_nodes, 3)`` coordinate array.
    quads : np.ndarray
        Shape ``(n_quads, 4)`` corner-node-index array.

    Returns
    -------
    np.ndarray
        Shape ``(n_quads,)`` array of areas.
    """
    p0 = nodes[quads[:, 0]]
    p1 = nodes[quads[:, 1]]
    p2 = nodes[quads[:, 2]]
    p3 = nodes[quads[:, 3]]
    # Triangle 0-1-2
    cross1 = np.cross(p1 - p0, p2 - p0)
    # Triangle 0-2-3
    cross2 = np.cross(p2 - p0, p3 - p0)
    a1 = 0.5 * np.linalg.norm(cross1, axis=1)
    a2 = 0.5 * np.linalg.norm(cross2, axis=1)
    return np.asarray(a1 + a2)


# ======================================================================
# BoundaryCondition dataclass
# ======================================================================

@dataclass
class BoundaryCondition:
    """A single boundary condition applied to mesh nodes.

    Exactly one of *face* or *node_ids* must be specified to identify
    the affected nodes.  For face-based BCs the node IDs are resolved
    at application time via :meth:`MeshData.nodes_on_face`.

    Parameters
    ----------
    bc_type : str
        Type of boundary condition:

        - ``"fixed"`` — constrain specified DOFs to zero.
        - ``"displacement"`` — prescribe a non-zero displacement.
        - ``"force"`` — apply a point force to each listed node.
        - ``"pressure"`` — apply a distributed force over a face,
          divided equally among face nodes.
        - ``"symmetry_x"`` — symmetry about the yz-plane (ux = 0).
        - ``"symmetry_y"`` — symmetry about the xz-plane (uy = 0).
        - ``"symmetry_z"`` — symmetry about the xy-plane (uz = 0).

    face : str or None
        Mesh face identifier: ``'x_min'``, ``'x_max'``, ``'y_min'``,
        ``'y_max'``, ``'z_min'``, ``'z_max'``.  Mutually exclusive
        with *node_ids*.
    node_ids : np.ndarray or None
        Explicit array of node indices (0-based).  Mutually exclusive
        with *face*.
    dofs : list[int]
        Which translational DOFs to affect (0 = ux, 1 = uy, 2 = uz).
        Default is ``[0, 1, 2]`` (all three).
    value : float
        Prescribed displacement (mm) or force magnitude (N).
        For ``"pressure"`` BCs, this is the total force (N) applied to
        the face, distributed equally among face nodes.

    Examples
    --------
    >>> bc = BoundaryCondition(bc_type="fixed", face="x_min", dofs=[0])
    >>> bc = BoundaryCondition(bc_type="displacement", face="x_max",
    ...                        dofs=[0], value=-0.5)
    """

    bc_type: str
    face: str | None = None
    node_ids: np.ndarray | None = None
    dofs: list[int] = field(default_factory=lambda: [0, 1, 2])
    value: float = 0.0

    _VALID_TYPES = frozenset({
        "fixed", "displacement", "force", "pressure",
        "symmetry_x", "symmetry_y", "symmetry_z",
    })

    def __post_init__(self) -> None:
        if self.bc_type not in self._VALID_TYPES:
            raise ValueError(
                f"Unknown bc_type '{self.bc_type}'. "
                f"Must be one of {sorted(self._VALID_TYPES)}."
            )
        if self.face is None and self.node_ids is None:
            raise ValueError(
                "Either 'face' or 'node_ids' must be specified."
            )

    def resolve_nodes(self, mesh: MeshData) -> np.ndarray:
        """Return the node indices affected by this BC.

        If *face* is set, queries ``mesh.nodes_on_face(face)``.
        Otherwise returns *node_ids* directly.

        Parameters
        ----------
        mesh : MeshData
            The mesh to resolve face names against.

        Returns
        -------
        np.ndarray
            1-D array of node indices (0-based).
        """
        if self.node_ids is not None:
            return np.asarray(self.node_ids, dtype=np.intp)
        if self.face is None:
            raise ValueError(
                "BoundaryCondition needs either node_ids or face to "
                "resolve its node set"
            )
        return mesh.nodes_on_face(self.face)

    def effective_dofs(self) -> list[int]:
        """Return the DOF indices affected, accounting for symmetry types.

        For ``"symmetry_x"`` returns ``[0]``, for ``"symmetry_y"`` returns
        ``[1]``, for ``"symmetry_z"`` returns ``[2]``.  Otherwise returns
        ``self.dofs``.

        Returns
        -------
        list[int]
        """
        symmetry_map = {
            "symmetry_x": [0],
            "symmetry_y": [1],
            "symmetry_z": [2],
        }
        return symmetry_map.get(self.bc_type, self.dofs)


# ======================================================================
# BoundaryHandler
# ======================================================================

class BoundaryHandler:
    """Applies boundary conditions to the global system K u = F.

    Supports two methods for imposing prescribed displacements:

    - **Penalty method** (:meth:`apply_penalty`): adds a large diagonal
      value to the stiffness matrix.  Simple, preserves matrix size,
      works well with direct solvers.
    - **Elimination method** (:meth:`apply_elimination`): partitions the
      system into free and constrained DOFs, producing a reduced system.
      More accurate, but changes the matrix size.

    Parameters
    ----------
    mesh : MeshData
        The finite element mesh.

    Examples
    --------
    Illustrative usage (needs a generated ``mesh``, so the snippet is skipped
    under ``--doctest-modules``).

    >>> handler = BoundaryHandler(mesh)  # doctest: +SKIP
    >>> bcs = BoundaryHandler.compression_bcs(mesh, applied_strain=-0.01)  # doctest: +SKIP
    >>> constrained = handler.get_constrained_dofs(bcs)  # doctest: +SKIP
    >>> F = handler.get_force_dofs(bcs)  # doctest: +SKIP
    >>> K_mod, F_mod = handler.apply_penalty(K, F, constrained)  # doctest: +SKIP
    """

    def __init__(self, mesh: MeshData) -> None:
        self.mesh = mesh

    # ------------------------------------------------------------------
    # Resolve BCs to DOF-level prescriptions
    # ------------------------------------------------------------------

    def get_constrained_dofs(
        self, bcs: list[BoundaryCondition]
    ) -> dict[int, float]:
        """Convert displacement-type BCs to a DOF-value mapping.

        Processes BCs of type ``"fixed"``, ``"displacement"``,
        ``"symmetry_x"``, ``"symmetry_y"``, and ``"symmetry_z"``.

        Parameters
        ----------
        bcs : list[BoundaryCondition]
            All boundary conditions.  Force/pressure BCs are ignored.

        Returns
        -------
        dict[int, float]
            Mapping ``{global_dof_index: prescribed_value}``.
            For fixed and symmetry BCs the value is 0.0.
        """
        disp_types = {"fixed", "displacement",
                       "symmetry_x", "symmetry_y", "symmetry_z"}
        constrained: dict[int, float] = {}

        for bc in bcs:
            if bc.bc_type not in disp_types:
                continue

            nodes = bc.resolve_nodes(self.mesh)
            dof_list = bc.effective_dofs()
            val = 0.0 if bc.bc_type in ("fixed", "symmetry_x",
                                         "symmetry_y", "symmetry_z") else bc.value

            for nid in nodes:
                for d in dof_list:
                    global_dof = 3 * int(nid) + d
                    constrained[global_dof] = val

        return constrained

    def get_force_dofs(
        self, bcs: list[BoundaryCondition]
    ) -> np.ndarray:
        """Assemble the global force vector from force/pressure BCs.

        For ``"force"`` BCs, the value is applied directly to each
        specified node in the specified DOFs.

        For ``"pressure"`` BCs the value is the *total* force over the
        face and is distributed via **consistent face integration**
        (issue #50).  For each face quad with area :math:`A_e` we
        contribute :math:`t \\cdot A_e / 4` to each of its four corner
        nodes, where the traction :math:`t = \\text{value} / A_\\text{face}`
        is uniform over the face:

        .. math::

            F_i = \\int_\\Gamma N_i(\\xi, \\eta)\\, t\\, dA
                = t \\sum_{e \\ni i} \\frac{A_e}{4}.

        On a uniform grid this reproduces the corner-1/4, edge-1/2,
        interior-1 tributary-area weights; on a wrinkled mesh the
        per-quad areas vary so the weights are area-correct.  When a
        pressure BC is supplied via explicit ``node_ids`` (no face
        topology available) we fall back to equal per-node distribution.

        Parameters
        ----------
        bcs : list[BoundaryCondition]
            All boundary conditions.  Only ``"force"`` and ``"pressure"``
            types contribute; others are ignored.

        Returns
        -------
        np.ndarray
            Shape ``(n_dof,)`` global force vector.
        """
        n_dof = self.mesh.n_dof
        F = np.zeros(n_dof, dtype=np.float64)

        for bc in bcs:
            if bc.bc_type == "force":
                nodes = bc.resolve_nodes(self.mesh)
                for nid in nodes:
                    for d in bc.dofs:
                        F[3 * int(nid) + d] += bc.value

            elif bc.bc_type == "pressure":
                self._apply_pressure_bc(bc, F)

        return F

    def _apply_pressure_bc(
        self, bc: BoundaryCondition, F: np.ndarray,
    ) -> None:
        """Accumulate consistent nodal forces from a pressure BC.

        See :meth:`get_force_dofs` for the underlying math.  This method
        mutates ``F`` in place.
        """
        if bc.face is not None:
            face_elems = self.mesh.face_elements(bc.face)
            if face_elems.shape[0] == 0:
                return
            quad_areas = _quad_areas(self.mesh.nodes, face_elems)
            total_area = float(quad_areas.sum())
            if total_area <= 0.0:
                # Degenerate face — fall back to equal split to stay
                # consistent with the legacy behaviour.
                nodes = bc.resolve_nodes(self.mesh)
                n_nodes = len(nodes)
                if n_nodes == 0:
                    return
                share = bc.value / n_nodes
                for nid in nodes:
                    for d in bc.dofs:
                        F[3 * int(nid) + d] += share
                return

            traction = bc.value / total_area
            # Each face quad contributes t * A_e / 4 to each of its
            # 4 corner nodes (Q4 bilinear shape functions, uniform t).
            contrib = traction * quad_areas * 0.25  # shape (n_face_e,)
            for d in bc.dofs:
                # Scatter-add into F[3*nid + d] for each corner.
                np.add.at(
                    F, 3 * face_elems.ravel() + d,
                    np.repeat(contrib, 4),
                )
            return

        # No face topology — explicit node list.  Fall back to equal
        # division (legacy behaviour for ``node_ids``-based pressure).
        nodes = bc.resolve_nodes(self.mesh)
        n_nodes = len(nodes)
        if n_nodes == 0:
            return
        share = bc.value / n_nodes
        for nid in nodes:
            for d in bc.dofs:
                F[3 * int(nid) + d] += share

    # ------------------------------------------------------------------
    # Penalty method
    # ------------------------------------------------------------------

    def apply_penalty(
        self,
        K: sparse.csc_matrix,
        F: np.ndarray,
        constrained_dofs: dict[int, float],
        penalty: float | None = None,
    ) -> tuple[sparse.csc_matrix, np.ndarray]:
        """Apply the penalty method for prescribed displacements.

        Thin wrapper around :func:`apply_penalty_bcs` (module level) that
        preserves the historical signature.  Always returns *copies* of
        ``K`` and ``F`` (the input arrays are not mutated).

        For each constrained DOF *i* with prescribed value *v*:

        .. math::
            K_{ii} \\mathrel{+}= \\alpha, \\qquad F_i \\mathrel{+}= \\alpha \\, v

        Parameters
        ----------
        K : scipy.sparse.csc_matrix
            Global stiffness matrix (``n_dof x n_dof``).
        F : np.ndarray
            Global force vector (``n_dof,``).  Not mutated.
        constrained_dofs : dict[int, float]
            Mapping ``{dof_index: prescribed_value}`` from
            :meth:`get_constrained_dofs`.
        penalty : float or None, optional
            Override the penalty value.  ``None`` (the default) computes
            ``alpha = _PENALTY_SCALE * max(|diag(K)|, 1.0)`` per Bathe
            sec. 4.2.2 — the recommended setting; a fixed mega-penalty
            (e.g. ``1e20``) pushes cond(K) past the float64 limit and
            silently corrupts the unconstrained DOFs.

        Returns
        -------
        K_modified : scipy.sparse.csc_matrix
            Modified stiffness matrix in CSC format (copy).
        F_modified : np.ndarray
            Modified force vector (copy).
        """
        return apply_penalty_bcs(
            K, F, constrained_dofs, in_place=False, penalty=penalty,
        )

    # ------------------------------------------------------------------
    # Elimination method
    # ------------------------------------------------------------------

    def apply_elimination(
        self,
        K: sparse.csc_matrix,
        F: np.ndarray,
        constrained_dofs: dict[int, float],
    ) -> tuple[sparse.csc_matrix, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Apply the elimination (partitioning) method for prescribed displacements.

        Partitions the system into free (f) and constrained (c) DOFs:

        .. math::
            K_{ff} \\, u_f = F_f - K_{fc} \\, u_c

        This produces a reduced system of size ``n_free x n_free``.

        Parameters
        ----------
        K : scipy.sparse.csc_matrix
            Global stiffness matrix (``n_dof x n_dof``).
        F : np.ndarray
            Global force vector (``n_dof,``).
        constrained_dofs : dict[int, float]
            Mapping ``{dof_index: prescribed_value}``.

        Returns
        -------
        K_ff : scipy.sparse.csc_matrix
            Reduced stiffness matrix (``n_free x n_free``).
        F_reduced : np.ndarray
            Reduced force vector (``n_free,``), accounting for the
            effect of prescribed displacements.
        free_dofs : np.ndarray
            Sorted array of free (unconstrained) DOF indices.
        constrained_dof_indices : np.ndarray
            Sorted array of constrained DOF indices.
        constrained_values : np.ndarray
            Prescribed values corresponding to ``constrained_dof_indices``.
        """
        n_dof = K.shape[0]
        all_dofs = np.arange(n_dof, dtype=np.intp)

        # Sort constrained DOFs for consistent ordering
        c_dofs_sorted = np.array(sorted(constrained_dofs.keys()), dtype=np.intp)
        c_vals = np.array([constrained_dofs[d] for d in c_dofs_sorted],
                          dtype=np.float64)

        # Free DOFs = all DOFs minus constrained DOFs
        constrained_set = set(c_dofs_sorted)
        free_dofs = np.array(
            [d for d in all_dofs if d not in constrained_set], dtype=np.intp
        )

        # Extract sub-matrices using sparse indexing
        # K_ff: free x free, K_fc: free x constrained
        K_csr = K.tocsr()
        # scipy-stubs does not model ndarray fancy indexing on sparse
        # matrices; the runtime supports it.
        K_ff = K_csr[np.ix_(free_dofs, free_dofs)].tocsc()  # type: ignore[index]
        K_fc = K_csr[np.ix_(free_dofs, c_dofs_sorted)]  # type: ignore[index]

        # Build reduced RHS: F_f - K_fc * u_c
        F_free = F[free_dofs].copy()
        if c_vals.size > 0 and np.any(c_vals != 0.0):
            F_free -= K_fc @ c_vals

        return K_ff, F_free, free_dofs, c_dofs_sorted, c_vals

    # ------------------------------------------------------------------
    # CLT LoadState to 3-D boundary conditions
    # ------------------------------------------------------------------

    @staticmethod
    def _corner_node(mesh: MeshData, x: float, y: float, z: float) -> np.ndarray:
        """Index of the mesh node nearest the point ``(x, y, z)``."""
        d = np.sum((mesh.nodes - np.array([x, y, z])) ** 2, axis=1)
        return np.array([int(np.argmin(d))], dtype=np.intp)

    @staticmethod
    def _rigid_body_bcs(mesh: MeshData) -> list[BoundaryCondition]:
        """Six point constraints that remove rigid-body motion — and nothing else.

        A membrane state must be applied with **tractions only**: fixing a
        whole face would restrain the very deformation being applied (a
        clamped ``x_min`` suppresses the Poisson contraction, and a
        ``symmetry_y`` plane suppresses shear outright).  So rigid-body
        motion is removed at three points instead, statically determinately:

        ==========================  ====================================
        node                        constrained
        ==========================  ====================================
        ``(x_min, y_min, z_min)``   ``ux, uy, uz``  (3 translations)
        ``(x_max, y_min, z_min)``   ``uy, uz``      (rot. about z and x)
        ``(x_min, y_max, z_min)``   ``uz``          (rot. about y)
        ==========================  ====================================

        Six constraints in total.  Fixing ``uy`` at the second node picks
        the *simple-shear* gauge (``u_x = gamma y``, ``u_y = 0``) out of the
        family of fields that differ only by a rigid rotation — the strain
        state, which is what the solver reports, is identical either way.

        Raises
        ------
        ValueError
            If the three anchor points are not distinct (a degenerate mesh
            with no extent in x or y), since the constraint set would then
            leave a rigid mode and the stiffness matrix singular.
        """
        nodes = mesh.nodes
        x0, y0, z0 = nodes.min(axis=0)
        x1, y1 = nodes[:, 0].max(), nodes[:, 1].max()

        p0 = BoundaryHandler._corner_node(mesh, x0, y0, z0)
        p1 = BoundaryHandler._corner_node(mesh, x1, y0, z0)
        p2 = BoundaryHandler._corner_node(mesh, x0, y1, z0)
        if len({int(p0[0]), int(p1[0]), int(p2[0])}) != 3:
            raise ValueError(
                "Cannot remove rigid-body motion: the mesh has no extent in "
                f"x or y (anchor nodes {int(p0[0])}, {int(p1[0])}, "
                f"{int(p2[0])} are not distinct). A membrane load state "
                "needs a two-dimensional footprint."
            )
        return [
            BoundaryCondition(bc_type="fixed", node_ids=p0, dofs=[0, 1, 2]),
            BoundaryCondition(bc_type="fixed", node_ids=p1, dofs=[1, 2]),
            BoundaryCondition(bc_type="fixed", node_ids=p2, dofs=[2]),
        ]

    @staticmethod
    def _membrane_bcs(
        load: LoadState, mesh: MeshData
    ) -> list[BoundaryCondition]:
        """Self-equilibrated tractions for a uniform membrane state.

        Every non-zero in-plane resultant is applied on **both** opposing
        faces, with opposite sign — that is what makes the state uniform
        rather than a cantilever reaction:

        - ``Nx``  : x-traction on ``x_max`` (+) and ``x_min`` (-)
        - ``Ny``  : y-traction on ``y_max`` (+) and ``y_min`` (-)
        - ``Nxy`` : the **complementary** shear pair — y-traction on the
          x-faces *and* x-traction on the y-faces.  Applying it on the
          x-faces alone produces a tip-loaded cantilever, not pure shear.

        Resultants are per unit width, so a face force is the resultant
        times the length of the edge it acts on (``Ly`` for the x-faces,
        ``Lx`` for the y-faces); ``get_force_dofs`` then distributes each
        face total by consistent Q4 face integration.

        The load set is self-equilibrated, so the only kinematic
        constraints needed are the six rigid-body anchors from
        :meth:`_rigid_body_bcs`.
        """
        Lx, Ly, _Lz = mesh.domain_size
        bcs = BoundaryHandler._rigid_body_bcs(mesh)

        def traction(face: str, dof: int, value: float) -> None:
            if value != 0.0:
                bcs.append(BoundaryCondition(
                    bc_type="pressure", face=face, dofs=[dof], value=value,
                ))

        for face, sign in (("x_max", 1.0), ("x_min", -1.0)):
            traction(face, 0, sign * load.Nx * Ly)
            traction(face, 1, sign * load.Nxy * Ly)
        for face, sign in (("y_max", 1.0), ("y_min", -1.0)):
            traction(face, 1, sign * load.Ny * Lx)
            traction(face, 0, sign * load.Nxy * Lx)

        return bcs

    @staticmethod
    def _bending_bcs(
        load: LoadState, mesh: MeshData
    ) -> list[BoundaryCondition]:
        """Curvature (prescribed-displacement) BCs for ``Mx`` / ``My``.

        ``kappa_x = Mx / D11`` and ``kappa_y = My / D22`` from the
        decoupled CLT moment-curvature relation, imposed as a linear
        through-thickness displacement on the far face.
        """
        bcs: list[BoundaryCondition] = []
        Lx, Ly, _Lz = mesh.domain_size
        bcs.extend(BoundaryHandler._rigid_body_bcs(mesh))

        if abs(load.Mx) > 0:
            xmax_nodes = mesh.nodes_on_face("x_max")
            z_coords = mesh.nodes[xmax_nodes, 2]
            z_mid = 0.5 * (z_coords.min() + z_coords.max())
            if mesh.laminate is None:
                raise ValueError(
                    "Cannot map Mx to a curvature boundary condition: the "
                    "mesh has no attached laminate.  Use WrinkleMesh.generate() "
                    "or attach a Laminate to MeshData.laminate."
                )
            D11 = float(mesh.laminate.D[0, 0])
            if D11 == 0.0:
                raise ValueError(
                    "Laminate bending stiffness D11 is zero; cannot map Mx "
                    "to a curvature boundary condition."
                )
            kappa_x = load.Mx / D11
            for nid in xmax_nodes:
                z = float(mesh.nodes[nid, 2])
                bcs.append(BoundaryCondition(
                    bc_type="displacement",
                    node_ids=np.array([nid], dtype=np.intp),
                    dofs=[0],
                    value=kappa_x * (z - z_mid) * Lx,
                ))

        if abs(load.My) > 0:
            ymax_nodes = mesh.nodes_on_face("y_max")
            z_coords = mesh.nodes[ymax_nodes, 2]
            z_mid = 0.5 * (z_coords.min() + z_coords.max())
            if mesh.laminate is None:
                raise ValueError(
                    "Cannot map My to a curvature boundary condition: the "
                    "mesh has no attached laminate.  Use WrinkleMesh.generate() "
                    "or attach a Laminate to MeshData.laminate."
                )
            D22 = float(mesh.laminate.D[1, 1])
            if D22 == 0.0:
                raise ValueError(
                    "Laminate bending stiffness D22 is zero; cannot map My "
                    "to a curvature boundary condition."
                )
            kappa_y = load.My / D22
            for nid in ymax_nodes:
                z = float(mesh.nodes[nid, 2])
                bcs.append(BoundaryCondition(
                    bc_type="displacement",
                    node_ids=np.array([nid], dtype=np.intp),
                    dofs=[1],
                    value=kappa_y * (z - z_mid) * Ly,
                ))

        return bcs

    @staticmethod
    def load_state_to_bcs(
        load: LoadState, mesh: MeshData
    ) -> list[BoundaryCondition]:
        """Convert a CLT load state to 3-D boundary conditions.

        **Membrane states** (``Nx``, ``Ny``, ``Nxy``, in any combination)
        are applied as self-equilibrated tractions on all four in-plane
        faces, with rigid-body motion removed at three points — see
        :meth:`_membrane_bcs`.  Verified against the closed-form CLT
        solution for a flat laminate: uniaxial, biaxial, pure shear and
        combined compression-shear all reproduce ``midplane_strains`` to
        better than 1 % on this package's default mesh density.

        **Curvature states** (``Mx``, ``My``) are applied as prescribed
        linear through-thickness displacements on the far face, with
        ``kappa = M / D`` from the decoupled CLT relation.

        Mixing the two is **rejected**: a curvature BC prescribes
        displacement on ``x_max``/``y_max`` while a membrane state applies
        traction to those same faces, and superposing the two would report
        a stress state corresponding to neither load.

        Parameters
        ----------
        load : LoadState
            CLT load state with force/moment resultants.
        mesh : MeshData
            The mesh to resolve face node IDs.

        Returns
        -------
        list[BoundaryCondition]
            Boundary conditions suitable for :meth:`get_constrained_dofs`
            and :meth:`get_force_dofs`.  Empty when the state carries no
            mechanical load.

        Raises
        ------
        ValueError
            If membrane and curvature resultants are combined, or if the
            mesh cannot support the required rigid-body anchors.
        """
        has_membrane = (
            abs(load.Nx) > 0 or abs(load.Ny) > 0 or abs(load.Nxy) > 0
        )
        has_bending = abs(load.Mx) > 0 or abs(load.My) > 0

        if has_membrane and has_bending:
            raise ValueError(
                "LoadState combines membrane (Nx/Ny/Nxy) and curvature "
                "(Mx/My) resultants, which this BC mapping cannot apply "
                "together: the curvature terms prescribe displacement on "
                "the same faces the membrane terms load with traction, so "
                "superposing them would describe neither load. Apply them "
                "in separate runs, or use a membrane-only / curvature-only "
                "state."
            )

        if has_membrane:
            return BoundaryHandler._membrane_bcs(load, mesh)
        if has_bending:
            return BoundaryHandler._bending_bcs(load, mesh)
        return []

    # ------------------------------------------------------------------
    # Convenience BC generators
    # ------------------------------------------------------------------

    @staticmethod
    def compression_bcs(
        mesh: MeshData,
        applied_strain: float = -0.01,
    ) -> list[BoundaryCondition]:
        """Standard uniaxial compression boundary conditions.

        Sets up a displacement-controlled compression test:

        - ``x_min``: ux = 0 (fixed in loading direction).
        - ``x_max``: ux = applied_strain * Lx (prescribed displacement).
        - ``y_min``: uy = 0 (symmetry about xz-plane).
        - One corner node fully fixed (rigid body suppression).

        Parameters
        ----------
        mesh : MeshData
            The finite element mesh.
        applied_strain : float, optional
            Applied nominal strain (negative for compression).
            Default is ``-0.01`` (1% compressive strain).

        Returns
        -------
        list[BoundaryCondition]
            List of BCs ready for the handler.

        Examples
        --------
        Illustrative usage (needs a generated ``mesh``, so the snippet is
        skipped under ``--doctest-modules``).

        >>> bcs = BoundaryHandler.compression_bcs(mesh, applied_strain=-0.005)  # doctest: +SKIP
        >>> handler = BoundaryHandler(mesh)  # doctest: +SKIP
        >>> constrained = handler.get_constrained_dofs(bcs)  # doctest: +SKIP
        """
        Lx = mesh.domain_size[0]
        prescribed_disp = applied_strain * Lx

        # Identify one corner node for full rigid body suppression
        xmin_nodes = mesh.nodes_on_face("x_min")
        zmin_nodes = mesh.nodes_on_face("z_min")
        # Corner node at (x_min, y_min, z_min)
        corner_candidates = np.intersect1d(xmin_nodes, zmin_nodes)
        if corner_candidates.size > 0:
            ymin_nodes = mesh.nodes_on_face("y_min")
            full_corner = np.intersect1d(corner_candidates, ymin_nodes)
            corner_node = np.array(
                [full_corner[0] if full_corner.size > 0 else corner_candidates[0]],
                dtype=np.intp,
            )
        else:
            corner_node = np.array([xmin_nodes[0]], dtype=np.intp)

        bcs = [
            # Fix ux on x_min (loading face support)
            BoundaryCondition(bc_type="fixed", face="x_min", dofs=[0]),
            # Prescribe ux on x_max (loading face)
            BoundaryCondition(
                bc_type="displacement", face="x_max",
                dofs=[0], value=prescribed_disp,
            ),
            # Symmetry: uy = 0 on y_min
            BoundaryCondition(bc_type="symmetry_y", face="y_min"),
            # Fully fix corner node for rigid body suppression (uy, uz)
            BoundaryCondition(
                bc_type="fixed", node_ids=corner_node, dofs=[1, 2],
            ),
        ]

        return bcs

    @staticmethod
    def bending_bcs(
        mesh: MeshData,
        curvature: float = 0.001,
    ) -> list[BoundaryCondition]:
        """Pure bending boundary conditions.

        Applies a linear through-thickness displacement on ``x_max``
        to produce a bending deformation:

        .. math::
            u_x(z) = \\kappa \\, (z - z_{\\text{mid}}) \\, L_x

        where ``z_mid`` is the laminate midplane z-coordinate.

        Support conditions:

        - ``x_min``: ux = 0, uy = 0.
        - ``y_min``: uy = 0 (symmetry).
        - One corner node fully fixed.

        Parameters
        ----------
        mesh : MeshData
            The finite element mesh.
        curvature : float, optional
            Applied curvature (1/mm).  Positive curvature produces
            tension on the z_max surface.  Default is ``0.001``.

        Returns
        -------
        list[BoundaryCondition]
            List of BCs for bending analysis.
        """
        Lx = mesh.domain_size[0]

        # Identify corner node
        xmin_nodes = mesh.nodes_on_face("x_min")
        zmin_nodes = mesh.nodes_on_face("z_min")
        ymin_nodes = mesh.nodes_on_face("y_min")
        corner_candidates = np.intersect1d(xmin_nodes, zmin_nodes)
        full_corner = np.intersect1d(corner_candidates, ymin_nodes)
        corner_node = np.array(
            [full_corner[0] if full_corner.size > 0 else xmin_nodes[0]],
            dtype=np.intp,
        )

        bcs = [
            # Fix ux and uy on x_min
            BoundaryCondition(bc_type="fixed", face="x_min", dofs=[0, 1]),
            # Symmetry: uy = 0 on y_min
            BoundaryCondition(bc_type="symmetry_y", face="y_min"),
            # Fully fix corner node (uz for rigid body)
            BoundaryCondition(
                bc_type="fixed", node_ids=corner_node, dofs=[2],
            ),
        ]

        # Linear through-thickness displacement on x_max
        xmax_nodes = mesh.nodes_on_face("x_max")
        z_coords_xmax = mesh.nodes[xmax_nodes, 2]
        z_mid = 0.5 * (z_coords_xmax.min() + z_coords_xmax.max())

        for nid in xmax_nodes:
            z = float(mesh.nodes[nid, 2])
            ux = curvature * (z - z_mid) * Lx
            bcs.append(BoundaryCondition(
                bc_type="displacement",
                node_ids=np.array([nid], dtype=np.intp),
                dofs=[0],
                value=ux,
            ))

        return bcs
