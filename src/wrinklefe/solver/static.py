"""Linear static finite element solver for composite laminates.

Solves the equilibrium equation :math:`K u = F` for the displacement field,
then recovers stresses and strains at Gauss points in both global and local
(material) coordinate systems.

The solver supports:

- **Direct** solution via ``scipy.sparse.linalg.spsolve`` (robust for any size).
- **Iterative** solution via conjugate gradient with ILU preconditioner
  (memory-efficient for large problems with >100 k DOFs).
- Automatic conversion from a CLT :class:`~wrinklefe.core.laminate.LoadState`
  to 3-D boundary conditions.

Workflow
--------
1. Create ``StaticSolver(mesh, laminate)``
2. Define boundary conditions (list of ``BoundaryCondition``)
3. Call ``solver.solve(bcs)`` to get :class:`~wrinklefe.solver.results.FieldResults`

References
----------
Bathe, K.-J. (2006). Finite Element Procedures.
Zienkiewicz, O.C. & Taylor, R.L. (2000). The Finite Element Method, Vol. 1.
"""

from __future__ import annotations

import logging
import time
from contextlib import contextmanager
from typing import TYPE_CHECKING

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as spla

from wrinklefe.core.laminate import Laminate, LoadState
from wrinklefe.core.mesh import MeshData
from wrinklefe.core.transforms import stress_transformation_3d
from wrinklefe.solver.assembler import GlobalAssembler
from wrinklefe.solver.results import FieldResults

logger = logging.getLogger(__name__)


@contextmanager
def _suppress_assembly_warnings():
    """Suppress divide-by-zero and invalid-value warnings during assembly."""
    with np.errstate(divide='ignore', invalid='ignore'):
        yield

if TYPE_CHECKING:
    from wrinklefe.solver.boundary import BoundaryCondition


class StaticSolver:
    """Linear static finite element solver.

    Solves :math:`K \\cdot u = F` for the displacement field, then
    recovers stresses and strains at Gauss points.

    Parameters
    ----------
    mesh : MeshData
        Finite element mesh.
    laminate : Laminate
        Laminate definition with ply materials and orientations.
    element_type : str, optional
        Element formulation: ``'hex8'`` (standard 2x2x2 integration).
        Default is ``'hex8'``.
    iterative_rtol : float, optional
        Relative-residual convergence tolerance for the iterative (CG)
        solver.  Default ``1e-10``.
    iterative_maxiter : int, optional
        Maximum CG iterations for the iterative solver.  Default ``10000``.
    ilu_drop_tol : float, optional
        Drop tolerance for the ILU preconditioner (``spilu``).  Default
        ``1e-4``.
    ilu_fill_factor : float or None, optional
        Upper bound on the ILU fill (``spilu``'s ``fill_factor``).
        ``None`` (default) leaves SciPy's own default in place.
    preconditioner : str, optional
        Preconditioner for the iterative solver: ``'ilu'`` (default),
        ``'jacobi'`` (diagonal), or ``'none'`` (unpreconditioned).

    Notes
    -----
    The iterative-solver controls default to the values previously
    hardcoded in :meth:`_solve_iterative`, so an iterative solve with the
    defaults is bit-for-bit unchanged.  They are normally supplied by the
    :class:`~wrinklefe.analysis.WrinkleAnalysis` pipeline from the
    matching :class:`~wrinklefe.analysis.AnalysisConfig` fields.
    """

    def __init__(
        self,
        mesh: MeshData,
        laminate: Laminate,
        element_type: str = "hex8",
        *,
        iterative_rtol: float = 1e-10,
        iterative_maxiter: int = 10000,
        ilu_drop_tol: float = 1e-4,
        ilu_fill_factor: float | None = None,
        preconditioner: str = "ilu",
        delta_T: float = 0.0,
    ) -> None:
        self.mesh = mesh
        self.laminate = laminate
        self.element_type = element_type
        # Uniform temperature change from the stress-free (cure) state
        # (issue #273 Stage 2).  ``T_service - T_stress_free``, so a cure
        # cool-down is negative.  Threaded into the assembler (thermal
        # load vector) and used by :meth:`recover_element_results` to
        # subtract the thermal initial strain before computing stress.
        self.delta_T = float(delta_T)
        self.assembler = GlobalAssembler(
            mesh, laminate, element_type, delta_T=delta_T
        )

        # Iterative-solver controls (issue #265). Defaults reproduce the
        # previously hardcoded values in ``_solve_iterative`` bit-for-bit.
        self.iterative_rtol = iterative_rtol
        self.iterative_maxiter = iterative_maxiter
        self.ilu_drop_tol = ilu_drop_tol
        self.ilu_fill_factor = ilu_fill_factor
        self.preconditioner = preconditioner

        # Populated after solve
        self._displacement: np.ndarray | None = None
        self._K: sparse.csc_matrix | None = None
        self._constrained_dofs: dict[int, float] = {}

    # ------------------------------------------------------------------
    # Main solve interface
    # ------------------------------------------------------------------

    def solve(
        self,
        boundary_conditions: list[BoundaryCondition],
        solver: str = "direct",
        verbose: bool = False,
        keep_stiffness: bool = False,
    ) -> FieldResults:
        """Solve the static problem.

        Steps
        -----
        1. Assemble global stiffness matrix K.
        2. Assemble force vector F from boundary conditions, plus the
           thermal initial-strain load when ``delta_T != 0``.
        3. Apply displacement BCs via the penalty method.
        4. Solve K u = F.
        5. Post-process: recover stresses and strains.

        Parameters
        ----------
        boundary_conditions : list[BoundaryCondition]
            List of boundary conditions (displacement and force BCs).
        solver : str, optional
            ``'direct'`` uses ``spsolve``; ``'iterative'`` uses CG with
            ILU preconditioner. Default is ``'direct'``.
        verbose : bool, optional
            Deprecated and ignored. Progress is reported through the
            ``wrinklefe.solver.static`` logger (milestones at INFO,
            per-element detail at DEBUG).
        keep_stiffness : bool, optional
            If True, retain a copy of the unmodified (pre-penalty) global
            stiffness matrix on ``self._K`` after solve. Default is False,
            which avoids the memory cost of holding a full sparse matrix
            copy on every solver instance (relevant for parametric sweeps
            and analyses that keep multiple solvers alive). Enable only
            when callers need post-solve access to K (e.g., to compute
            reaction forces without re-assembly).

        Returns
        -------
        FieldResults
            Complete solution with displacement, stress, and strain fields.

        Raises
        ------
        RuntimeError
            If the iterative solver fails to converge.
        """
        from wrinklefe.solver.boundary import BoundaryHandler

        t0 = time.perf_counter()

        # 1. Assemble global stiffness
        logger.debug("Assembling global stiffness matrix...")
        with _suppress_assembly_warnings():
            K = self.assembler.assemble_stiffness(verbose=verbose)
        # Opt-in retention of the unmodified K (see ``keep_stiffness``).
        # Default: do not store a copy to avoid doubling peak FE memory.
        self._K = K.copy() if keep_stiffness else None

        t1 = time.perf_counter()
        logger.info("Assembly time: %.2f s", t1 - t0)

        # 2. Assemble force vector
        bc_handler = BoundaryHandler(self.mesh)
        F = bc_handler.get_force_dofs(boundary_conditions)

        # 2b. Thermal initial-strain load (issue #273 Stage 2).  Zero when
        # ``delta_T == 0``, so a purely mechanical run is unchanged.
        if self.delta_T != 0.0:
            F = F + self.assembler.assemble_thermal_force()

        # 3. Apply displacement BCs via penalty method
        self._constrained_dofs = bc_handler.get_constrained_dofs(
            boundary_conditions
        )
        K, F = self._apply_penalty_bcs(K, F, self._constrained_dofs, verbose)

        # 4. Solve
        logger.info(
            "Solving system (%d DOFs, solver=%s)...", self.mesh.n_dof, solver
        )

        if solver == "direct":
            u = self._solve_direct(K, F, verbose=verbose)
        elif solver == "iterative":
            u = self._solve_iterative(K, F, verbose=verbose)
        else:
            raise ValueError(
                f"Unknown solver '{solver}'. Use 'direct' or 'iterative'."
            )

        t2 = time.perf_counter()
        logger.info("Solve time: %.2f s", t2 - t1)
        t1 = t2

        # 5. Post-process
        logger.debug("Recovering element stresses and strains...")

        stress_g, stress_l, strain_g, strain_l = self.recover_element_results(
            u, verbose=verbose
        )

        # Reshape displacement to (n_nodes, 3)
        displacement = u.reshape(-1, 3)
        self._displacement = displacement

        t3 = time.perf_counter()
        logger.info("Post-processing time: %.2f s", t3 - t1)
        logger.info("Total solve time: %.2f s", t3 - t0)

        return FieldResults(
            displacement=displacement,
            stress_global=stress_g,
            stress_local=stress_l,
            strain_global=strain_g,
            strain_local=strain_l,
            mesh=self.mesh,
            laminate=self.laminate,
        )

    def solve_load_state(
        self,
        load: LoadState,
        solver: str = "direct",
        verbose: bool = False,
    ) -> FieldResults:
        """Convenience method: solve from a CLT LoadState.

        Converts the CLT force and moment resultants into 3-D boundary
        conditions on the mesh faces:

        - **x_min** face: fully clamped (ux = uy = uz = 0).
        - **x_max** face: uniform traction derived from Nx, Ny, Nxy
          distributed over the face area.
        - **y_min** / **y_max** faces: free (natural BC) unless Ny or Nxy
          are non-zero (handled through x_max traction).
        - **z_min** / **z_max** faces: free.

        Parameters
        ----------
        load : LoadState
            CLT load state (Nx, Ny, Nxy, Mx, My, Mxy).
        solver : str, optional
            ``'direct'`` or ``'iterative'``. Default is ``'direct'``.
        verbose : bool, optional
            Print progress. Default is ``False``.

        Returns
        -------
        FieldResults
            Complete solution data.
        """

        bcs = self._load_state_to_bcs(load)
        return self.solve(bcs, solver=solver, verbose=verbose)

    # ------------------------------------------------------------------
    # Linear algebra solvers
    # ------------------------------------------------------------------

    def _solve_direct(
        self,
        K: sparse.csc_matrix,
        F: np.ndarray,
        verbose: bool = False,
    ) -> np.ndarray:
        """Direct sparse solver using ``scipy.sparse.linalg.spsolve``.

        Parameters
        ----------
        K : scipy.sparse.csc_matrix
            Global stiffness matrix with BCs applied.
        F : np.ndarray
            Shape ``(n_dof,)`` global force vector.
        verbose : bool, optional
            Print solver info.

        Returns
        -------
        np.ndarray
            Shape ``(n_dof,)`` displacement vector.
        """
        u = spla.spsolve(K, F)
        if logger.isEnabledFor(logging.DEBUG):
            residual = np.linalg.norm(K @ u - F)
            logger.debug("Direct solver residual: %.4e", residual)
        return u

    def _diagonal_preconditioner(
        self, K: sparse.csc_matrix
    ) -> spla.LinearOperator:
        """Build a diagonal (Jacobi) ``LinearOperator`` from ``K``.

        Zero diagonal entries are replaced by 1.0 so the inverse is
        well-defined.
        """
        n = K.shape[0]
        diag = K.diagonal()
        diag[diag == 0] = 1.0
        return spla.LinearOperator(
            shape=(n, n),
            matvec=lambda x: x / diag,
            dtype=K.dtype,
        )

    def _build_preconditioner(
        self, K: sparse.csc_matrix
    ) -> tuple[spla.LinearOperator | None, str]:
        """Build the CG preconditioner and report which one is active.

        Honours ``self.preconditioner`` (``'ilu'`` | ``'jacobi'`` |
        ``'none'``).  For ``'ilu'``, a *narrow* fallback to the diagonal
        (Jacobi) preconditioner is engaged when ``spilu`` raises one of
        the errors it uses to signal a memory/singular/structural failure
        (``RuntimeError``, ``MemoryError``, ``ValueError``); that fallback
        is logged unconditionally at WARNING level.  Any other exception
        type propagates unchanged so genuinely unexpected bugs are not
        masked as an ILU failure.

        Returns
        -------
        (M_op, active) : tuple
            ``M_op`` is the ``LinearOperator`` (or ``None`` for the
            unpreconditioned case); ``active`` is the name of the
            preconditioner actually in effect (``'ilu'``, ``'jacobi'`` or
            ``'none'``) — used in the non-convergence diagnostics.
        """
        n = K.shape[0]
        precond = self.preconditioner.lower().strip()

        if precond == "none":
            logger.debug("Using unpreconditioned CG (preconditioner='none').")
            return None, "none"

        if precond == "jacobi":
            logger.debug("Building diagonal (Jacobi) preconditioner...")
            return self._diagonal_preconditioner(K), "jacobi"

        # Default: incomplete-LU.
        logger.debug("Building ILU preconditioner...")
        spilu_kwargs: dict = {"drop_tol": self.ilu_drop_tol}
        if self.ilu_fill_factor is not None:
            spilu_kwargs["fill_factor"] = self.ilu_fill_factor
        try:
            ilu = spla.spilu(K, **spilu_kwargs)
        except (RuntimeError, MemoryError, ValueError) as err:
            # Narrow fallback: only the failure modes ``spilu`` uses to
            # signal that the factorisation itself could not be built
            # (out of memory, structurally/numerically singular input).
            # Any other exception type is a bug and must not be masked.
            logger.warning(
                "ILU preconditioner failed (%s: %s); falling back to a "
                "diagonal (Jacobi) preconditioner — the solve may be much "
                "slower on an ill-conditioned matrix. Set "
                "preconditioner='jacobi' to silence this, or use the "
                "direct solver.",
                type(err).__name__,
                err,
            )
            return self._diagonal_preconditioner(K), "jacobi"

        M_op = spla.LinearOperator(
            shape=(n, n),
            matvec=ilu.solve,
            dtype=K.dtype,
        )
        return M_op, "ilu"

    def _solve_iterative(
        self,
        K: sparse.csc_matrix,
        F: np.ndarray,
        tol: float | None = None,
        maxiter: int | None = None,
        verbose: bool = False,
    ) -> np.ndarray:
        """Iterative CG solver with a configurable preconditioner.

        Uses ``scipy.sparse.linalg.cg`` with the preconditioner selected
        by ``self.preconditioner`` (an incomplete-LU factorisation by
        default), wrapped in a ``LinearOperator``.  Suitable for large
        problems (>100 k DOFs).

        The convergence tolerance, iteration cap, and preconditioner
        controls come from the instance attributes set on
        :class:`StaticSolver` (``iterative_rtol``, ``iterative_maxiter``,
        ``ilu_drop_tol``, ``ilu_fill_factor``, ``preconditioner``), which
        the :class:`~wrinklefe.analysis.WrinkleAnalysis` pipeline plumbs
        from :class:`~wrinklefe.analysis.AnalysisConfig`.  Their defaults
        reproduce the previously hardcoded values bit-for-bit.

        Parameters
        ----------
        K : scipy.sparse.csc_matrix
            Global stiffness matrix with BCs applied.
        F : np.ndarray
            Shape ``(n_dof,)`` global force vector.
        tol : float or None, optional
            Override for the CG relative-residual tolerance.  ``None``
            (default) uses ``self.iterative_rtol``.
        maxiter : int or None, optional
            Override for the CG iteration cap.  ``None`` (default) uses
            ``self.iterative_maxiter``.
        verbose : bool, optional
            Deprecated and ignored; progress is reported through the
            module logger.

        Returns
        -------
        np.ndarray
            Shape ``(n_dof,)`` displacement vector.

        Raises
        ------
        RuntimeError
            If the CG solver fails to converge.  The message names the
            active preconditioner, the iterations used, the iteration cap,
            and the final relative residual.
        """
        if tol is None:
            tol = self.iterative_rtol
        if maxiter is None:
            maxiter = self.iterative_maxiter

        M_op, active_preconditioner = self._build_preconditioner(K)

        # Iteration counter for logging and error reporting
        iter_count = [0]

        def _callback(xk: np.ndarray) -> None:
            iter_count[0] += 1

        callback = _callback

        # SciPy >=1.12 deprecated ``tol=`` in favour of ``rtol=``.
        u, info = spla.cg(K, F, rtol=tol, maxiter=maxiter, M=M_op,
                          callback=callback)

        if info != 0:
            b_norm = float(np.linalg.norm(F))
            rel_residual = (
                float(np.linalg.norm(K @ u - F)) / b_norm
                if b_norm > 0.0
                else float(np.linalg.norm(K @ u - F))
            )
            raise RuntimeError(
                f"CG solver did not converge: info={info} "
                f"(preconditioner={active_preconditioner!r}, "
                f"iterations={iter_count[0]}, maxiter={maxiter}, "
                f"rtol={tol}, final relative residual="
                f"{rel_residual:.4e})"
            )

        if logger.isEnabledFor(logging.INFO):
            residual = np.linalg.norm(K @ u - F)
            logger.info(
                "CG converged in %d iterations (preconditioner=%s), "
                "residual: %.4e",
                iter_count[0], active_preconditioner, residual,
            )

        return u

    # ------------------------------------------------------------------
    # Boundary condition helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _apply_penalty_bcs(
        K: sparse.csc_matrix,
        F: np.ndarray,
        constrained_dofs: dict[int, float],
        verbose: bool = False,
    ) -> tuple[sparse.csc_matrix, np.ndarray]:
        """Apply displacement boundary conditions via the penalty method.

        Thin wrapper around
        :func:`wrinklefe.solver.boundary.apply_penalty_bcs` that uses
        ``in_place=True`` (this solver owns ``K`` and ``F`` exclusively
        within :meth:`solve`).  See the helper docstring for the math
        and the choice of penalty scaling.
        """
        from wrinklefe.solver.boundary import _PENALTY_SCALE, apply_penalty_bcs

        if not constrained_dofs:
            return K, F

        K_out, F_out = apply_penalty_bcs(
            K, F, constrained_dofs, in_place=True,
        )

        if logger.isEnabledFor(logging.DEBUG):
            diag_max = float(np.abs(K.diagonal()).max())
            alpha = _PENALTY_SCALE * max(diag_max, 1.0)
            logger.debug(
                "Applied %d displacement BCs (penalty alpha=%.2e)",
                len(constrained_dofs), alpha,
            )

        return K_out, F_out

    def _load_state_to_bcs(self, load: LoadState) -> list:
        """Convert a CLT LoadState to 3-D boundary conditions.

        Parameters
        ----------
        load : LoadState
            CLT-level load state.

        Returns
        -------
        list[BoundaryCondition]
            List of BCs suitable for ``self.solve()``.
        """
        from wrinklefe.solver.boundary import BoundaryCondition

        bcs: list[BoundaryCondition] = []

        # Clamp x_min face: ux = uy = uz = 0
        x_min_nodes = self.mesh.nodes_on_face("x_min")
        bcs.append(
            BoundaryCondition(
                bc_type="fixed",
                node_ids=x_min_nodes,
                dofs=[0, 1, 2],
                value=0.0,
            )
        )

        # Apply traction on x_max face from Nx via a pressure BC so the
        # total face force is distributed by consistent face integration
        # (see boundary.BoundaryHandler.get_force_dofs and issue #50).
        # CLT Nx has units of force per unit width (N/mm), so the total
        # nodal force across the face is Nx * Ly.
        _, Ly, Lz = self.mesh.domain_size
        x_max_has_elements = self.mesh.nx > 0 and self.mesh.ny > 0 and self.mesh.nz > 0

        if x_max_has_elements and not np.isclose(load.Nx, 0.0):
            total_force_x = load.Nx * Ly
            bcs.append(
                BoundaryCondition(
                    bc_type="pressure",
                    face="x_max",
                    dofs=[0],
                    value=total_force_x,
                )
            )

        if x_max_has_elements and not np.isclose(load.Nxy, 0.0):
            total_force_y = load.Nxy * Ly
            bcs.append(
                BoundaryCondition(
                    bc_type="pressure",
                    face="x_max",
                    dofs=[1],
                    value=total_force_y,
                )
            )

        return bcs

    # ------------------------------------------------------------------
    # Post-processing
    # ------------------------------------------------------------------

    # Cached (gp_coords, N_gp, dN_dxi_gp) for the standard hex8 / 2x2x2
    # quadrature.  Populated lazily by :meth:`_hex8_gauss_shape_functions`
    # and reused across every solve because shape functions in natural
    # coordinates are identical for every hex8 element (issue #187).
    _hex8_gp_shape_cache: (
        tuple[np.ndarray, np.ndarray, np.ndarray] | None
    ) = None

    @classmethod
    def _hex8_gauss_shape_functions(
        cls,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return cached 2x2x2 Gauss-point coords, shape functions, and derivatives.

        Returns
        -------
        gp_coords : np.ndarray
            Shape ``(n_gp, 3)`` natural coordinates of the Gauss points
            (lexicographic order, matching ``Hex8Element._gauss_points``).
        N_gp : np.ndarray
            Shape ``(n_gp, 8)`` — row ``i`` holds the 8 trilinear shape
            functions evaluated at Gauss point ``i``.  For a per-element
            nodal field ``f`` (shape ``(8,)``), ``N_gp @ f`` returns the
            interpolated value at each Gauss point in one matmul.
        dN_dxi_gp : np.ndarray
            Shape ``(n_gp, 3, 8)`` — natural-coordinate derivatives of the
            8 shape functions at each Gauss point.  Constant across all
            hex8 elements; used to build per-element Jacobians via a
            single batched matmul.
        """
        if cls._hex8_gp_shape_cache is not None:
            return cls._hex8_gp_shape_cache

        from wrinklefe.elements.gauss import gauss_points_hex
        from wrinklefe.elements.hex8 import Hex8Element

        gp_coords, _ = gauss_points_hex(order=2)
        n_gp = gp_coords.shape[0]
        N_gp = np.empty((n_gp, 8))
        dN_dxi_gp = np.empty((n_gp, 3, 8))
        for i, (xi, eta, zeta) in enumerate(gp_coords):
            N_gp[i] = Hex8Element.shape_functions(xi, eta, zeta)
            dN_dxi_gp[i] = Hex8Element.shape_derivatives(xi, eta, zeta)
        cls._hex8_gp_shape_cache = (gp_coords, N_gp, dN_dxi_gp)
        return cls._hex8_gp_shape_cache

    def recover_element_results(
        self,
        displacement: np.ndarray,
        verbose: bool = False,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Post-process element-level stresses and strains.

        For each element:

        1. Extract element displacements from the global vector.
        2. Compute global-frame stress and strain at each Gauss point.
           When ``delta_T != 0`` the thermal initial strain is subtracted
           from the total strain before the constitutive law is applied
           (issue #273 Stage 2), so the reported stress is the mechanical
           (residual) stress; the reported strain stays the *total*
           kinematic strain ``B u``.
        3. Transform stress and strain to local material coordinates
           using the ply angle and wrinkle misalignment.  Stress uses the
           stress transformation ``T_sigma``; **engineering strain uses
           the strain transformation** ``T_eps = R T_sigma R^-1``
           (``R = diag(1,1,1,2,2,2)``), which is a different matrix
           wherever the rotation couples a normal component to a shear.

        Quantities that are constant per element (``T_ply``, the element's
        node ids and per-node wrinkle angles) are computed once outside
        the Gauss loop, and all wrinkle angles at the 8 Gauss points are
        obtained via a single ``N_gp @ fiber_angles_local`` matmul.  See
        issue #187 for the performance motivation.

        Parameters
        ----------
        displacement : np.ndarray
            Shape ``(n_dof,)`` global displacement vector.
        verbose : bool, optional
            Print progress. Default is ``False``.

        Returns
        -------
        stress_global : np.ndarray
            Shape ``(n_elements, n_gauss, 6)`` stress in global coordinates.
        stress_local : np.ndarray
            Shape ``(n_elements, n_gauss, 6)`` stress in local material coordinates.
        strain_global : np.ndarray
            Shape ``(n_elements, n_gauss, 6)`` strain in global coordinates.
        strain_local : np.ndarray
            Shape ``(n_elements, n_gauss, 6)`` strain in local material coordinates.
        """
        from wrinklefe.core.transforms import (
            rotate_stiffness_3d,
            strain_transformation_3d,
        )

        thermal = self.delta_T != 0.0
        n_elem = self.mesh.n_elements
        # Shape functions and natural-coord derivatives at the 8 Gauss
        # points are constant across all hex8 elements; build once and
        # reuse (cached on the class).
        _gp_coords, N_gp, dN_dxi_gp = self._hex8_gauss_shape_functions()
        n_gp = N_gp.shape[0]

        stress_global = np.empty((n_elem, n_gp, 6))
        stress_local = np.empty((n_elem, n_gp, 6))
        strain_global = np.empty((n_elem, n_gp, 6))
        strain_local = np.empty((n_elem, n_gp, 6))

        # Pre-fetch mesh arrays so attribute lookups don't happen in the
        # hot loop.
        mesh_elements = self.mesh.elements
        mesh_nodes = self.mesh.nodes
        mesh_fiber_angles = self.mesh.fiber_angles
        mesh_ply_angles = self.mesh.ply_angles
        mesh_ply_ids = self.mesh.ply_ids

        # Per-element material is resolved through the mesh's single
        # decision point (progressive-damage override > graded resin blend
        # > binary resin pocket > host ply), and the wrinkle-misalignment
        # retention through ``resin_angle_scale`` (0 at a fibre-free resin
        # centre, 1 in the bulk).  Stiffness is cached by material identity
        # (shared ply objects hit the cache; per-element blended/degraded
        # materials are distinct).
        plies = self.laminate.plies
        mat_C_cache: dict[int, np.ndarray] = {}
        # Material-axes CTE vectors, cached by material identity alongside
        # the stiffness (issue #273 Stage 2).  Only populated for a
        # thermally loaded run.
        mat_alpha_cache: dict[int, np.ndarray] = {}

        # Reusable scratch B-matrix template (zeros are stable between GPs
        # at the slots we never touch; we overwrite the populated slots
        # for every GP via vector assignment).
        B_scratch = np.zeros((n_gp, 6, 24))

        for e in range(n_elem):
            if e % 1000 == 0:
                logger.debug(
                    "Post-processing element %d/%d (%.1f%%)",
                    e, n_elem, 100.0 * e / n_elem,
                )

            node_ids = mesh_elements[e]
            node_coords = mesh_nodes[node_ids]  # (8, 3)

            # Batched Jacobians across all 8 GPs in one matmul.
            J_all = dN_dxi_gp @ node_coords  # (n_gp, 3, 3)
            J_inv_all = np.linalg.inv(J_all)  # (n_gp, 3, 3)
            dN_dx_all = J_inv_all @ dN_dxi_gp  # (n_gp, 3, 8) — physical-coord derivs

            # Build the 6x24 strain-displacement matrix at every GP via
            # vectorised slot assignment (matches Hex8Element.B_matrix).
            B = B_scratch
            # eps_11 = du/dx -> rows 0, cols 0,3,6,...,21
            B[:, 0, 0::3] = dN_dx_all[:, 0, :]
            # eps_22 = dv/dy -> rows 1, cols 1,4,...,22
            B[:, 1, 1::3] = dN_dx_all[:, 1, :]
            # eps_33 = dw/dz -> rows 2, cols 2,5,...,23
            B[:, 2, 2::3] = dN_dx_all[:, 2, :]
            # gamma_23 = dv/dz + dw/dy -> rows 3, cols (1,2,4,5,...)
            B[:, 3, 1::3] = dN_dx_all[:, 2, :]
            B[:, 3, 2::3] = dN_dx_all[:, 1, :]
            # gamma_13 = du/dz + dw/dx -> rows 4
            B[:, 4, 0::3] = dN_dx_all[:, 2, :]
            B[:, 4, 2::3] = dN_dx_all[:, 0, :]
            # gamma_12 = du/dy + dv/dx -> rows 5
            B[:, 5, 0::3] = dN_dx_all[:, 1, :]
            B[:, 5, 1::3] = dN_dx_all[:, 0, :]

            # Element nodal displacements (24,) via flat node-id expansion.
            u_elem = displacement[
                (3 * node_ids[:, None] + np.arange(3)).ravel()
            ]

            # Strain at each GP — one batched matmul.
            eps_g = B @ u_elem  # (n_gp, 6)

            # Rotated stiffness per GP: ply rotation (about z) once per
            # element, then wrinkle rotation (about y) per GP.
            ply_material = plies[int(mesh_ply_ids[e])].material
            mat_e = self.mesh.element_material(e, ply_material)
            mid = id(mat_e)
            C_material = mat_C_cache.get(mid)
            if C_material is None:
                C_material = mat_e.stiffness_matrix
                mat_C_cache[mid] = C_material

            ply_angle_rad = np.radians(float(mesh_ply_angles[e]))
            if abs(ply_angle_rad) > 1.0e-15:
                C_ply = rotate_stiffness_3d(C_material, ply_angle_rad, axis='z')
            else:
                C_ply = C_material

            # Ply-frame CTE vector, rotated to the laminate frame with the
            # INVERSE strain transformation (issue #273 Stage 2) so that
            # ``sigma = C_bar (eps - alpha dT)`` holds in global axes.  The
            # wrinkle rotation is applied per Gauss point below, mirroring
            # the stiffness.
            if thermal:
                alpha_material = mat_alpha_cache.get(mid)
                if alpha_material is None:
                    alpha_material = np.array([
                        float(getattr(mat_e, "alpha1", 0.0)),
                        float(getattr(mat_e, "alpha2", 0.0)),
                        float(getattr(mat_e, "alpha3", 0.0)),
                        0.0, 0.0, 0.0,
                    ])
                    mat_alpha_cache[mid] = alpha_material
                if abs(ply_angle_rad) > 1.0e-15:
                    alpha_ply = (
                        strain_transformation_3d(-ply_angle_rad, axis='z')
                        @ alpha_material
                    )
                else:
                    alpha_ply = alpha_material

            angle_scale = self.mesh.resin_angle_scale(e)
            if angle_scale == 0.0:
                wrinkle_angles_gp = np.zeros(n_gp)
            else:
                fiber_angles_local = mesh_fiber_angles[node_ids]  # (8,)
                # One matmul gives the interpolated wrinkle angle per GP,
                # scaled by the resin retention factor.
                wrinkle_angles_gp = angle_scale * (N_gp @ fiber_angles_local)

            T_ply = stress_transformation_3d(ply_angle_rad, axis='z')
            # Engineering STRAIN needs its own transformation: T_eps =
            # R T_sig R^-1 with R = diag(1,1,1,2,2,2).  Reusing T_sig on a
            # strain vector silently mis-scales every component that the
            # rotation couples to a shear.
            T_ply_eps = strain_transformation_3d(ply_angle_rad, axis='z')
            sig_g = np.empty((n_gp, 6))
            for g in range(n_gp):
                phi = float(wrinkle_angles_gp[g])
                if abs(phi) > 1.0e-15:
                    C_bar = rotate_stiffness_3d(C_ply, phi, axis='y')
                else:
                    C_bar = C_ply
                if thermal:
                    if abs(phi) > 1.0e-15:
                        alpha_g = (
                            strain_transformation_3d(-phi, axis='y')
                            @ alpha_ply
                        )
                    else:
                        alpha_g = alpha_ply
                    # Only the mechanical strain carries stress.
                    sig_g[g] = C_bar @ (eps_g[g] - alpha_g * self.delta_T)
                else:
                    sig_g[g] = C_bar @ eps_g[g]

                # Composition order is fixed by how ``C_bar`` was built:
                #     C_bar = R_y(R_z(C)) = Ts(phi)^-1 Ts(th)^-1 C Te(th) Te(phi)
                # so   sigma_g = C_bar eps_g   implies
                #     Ts(th) Ts(phi) sigma_g = C [Te(th) Te(phi) eps_g].
                # The PLY transform therefore sits on the LEFT: the wrinkle
                # rotation is the outer one taking the intermediate frame to
                # global, so it is the first one undone.  (Reversing this is
                # invisible for a 0 deg ply or a pristine element, where one
                # factor is the identity, but reaches ~10 % on a 90 deg ply
                # at phi = 0.1 rad.)
                T_wrinkle = stress_transformation_3d(phi, axis='y')
                stress_local[e, g] = T_ply @ (T_wrinkle @ sig_g[g])

                T_wrinkle_eps = strain_transformation_3d(phi, axis='y')
                strain_local[e, g] = T_ply_eps @ (T_wrinkle_eps @ eps_g[g])

            stress_global[e] = sig_g
            strain_global[e] = eps_g

        logger.debug(
            "Post-processing element %d/%d (100.0%%) -- done.", n_elem, n_elem
        )

        return stress_global, stress_local, strain_global, strain_local

    # Cached (extrapolation_matrix, node_to_gp) for hex8 / 2x2x2.
    # Populated lazily by :meth:`_extrapolate_to_nodes`.
    _hex8_extrap_cache: tuple[np.ndarray, np.ndarray] | None = None

    @classmethod
    def _build_hex8_extrapolation(cls) -> tuple[np.ndarray, np.ndarray]:
        """Build (and cache) the hex8 / 2x2x2 Gauss-to-node extrapolation matrix.

        Two ordering conventions are at play and **must not be conflated**:

        - ``gauss_points_hex(order=2)`` orders the 8 Gauss points
          **lexicographically** in ``(xi, eta, zeta)`` via
          ``np.meshgrid(..., indexing="ij")`` — zeta varies fastest.
        - :data:`wrinklefe.elements.hex8._NODE_COORDS` orders the 8 nodes by
          the **VTK / Abaqus** convention (bottom face CCW, then top face CCW).

        The i-th Gauss point is **not** at the natural coordinate of the
        i-th node.  Pairing them by index alone (a tempting but wrong
        simplification) scrambles the extrapolated nodal field — this is
        the trap reported in issue #51.

        We build the relationship explicitly: for each VTK node *j*, find
        the lexicographic Gauss-point index ``node_to_gp[j]`` whose natural
        coordinates have the same sign pattern as node *j*.  That gives
        the 1-to-1 nearest-neighbour mapping the caller usually wants for
        diagnostics, and is used by the regression test in
        ``tests/test_solver/test_extrapolate.py``.

        The extrapolation itself uses the full inverse shape-function
        matrix (not just the nearest-neighbour pairing), which is exact
        for a tri-linear field:

        .. math::

            \\mathbf{N}_{gp}[i, j] = N_j(\\xi_{gp,i})
            \\quad\\Rightarrow\\quad
            \\mathbf{f}_{nodes} = \\mathbf{N}_{gp}^{-1} \\, \\mathbf{f}_{gp}

        where row *i* uses **lex** Gauss-point ordering and column *j* uses
        **VTK** node ordering — so the output is naturally indexed by VTK
        node order, matching ``mesh.elements`` connectivity.

        Returns
        -------
        N_inv : np.ndarray
            Shape ``(8, 8)`` — extrapolation matrix; ``f_nodes = N_inv @ f_gp``
            where ``f_gp`` is in **lex** order and ``f_nodes`` is in **VTK**
            order.
        node_to_gp : np.ndarray
            Shape ``(8,)`` integer array; ``node_to_gp[j]`` is the lex GP
            index nearest to VTK node *j*.  For the default conventions
            this is ``[0, 4, 6, 2, 1, 5, 7, 3]``.
        """
        if cls._hex8_extrap_cache is not None:
            return cls._hex8_extrap_cache

        from wrinklefe.elements.gauss import gauss_points_hex
        from wrinklefe.elements.hex8 import _NODE_COORDS, Hex8Element

        gp_coords, _ = gauss_points_hex(order=2)  # lex order, shape (8, 3)

        # Build node->GP nearest mapping by matching sign patterns.  Each
        # GP sits at the centroid of one of the 8 sub-cubes; for a hex8
        # the unambiguous pairing is by sign of (xi, eta, zeta).
        node_to_gp = np.empty(8, dtype=int)
        for j in range(8):
            target = _NODE_COORDS[j]  # (+/-1, +/-1, +/-1)
            # Use sign-match: argmin over distance is equivalent here and
            # is robust to any future change in the GP magnitudes.
            sign_match = np.all(
                np.sign(gp_coords) == np.sign(target)[None, :], axis=1
            )
            matches = np.flatnonzero(sign_match)
            if matches.size != 1:
                raise RuntimeError(
                    f"Failed to pair hex8 VTK node {j} with a unique 2x2x2 "
                    f"Gauss point (matches={matches.tolist()}). Did the GP or node ordering "
                    "convention change?"
                )
            node_to_gp[j] = int(matches[0])

        # Sanity check: mapping must be a permutation of 0..7.
        if sorted(node_to_gp.tolist()) != list(range(8)):
            raise RuntimeError(
                "hex8 node->GP mapping is not a permutation: "
                f"{node_to_gp.tolist()}"
            )

        # Build shape-function matrix: row i uses lex GP i, col j is VTK node j.
        N_gp = np.empty((8, 8))
        for i in range(8):
            xi, eta, zeta = gp_coords[i]
            N_gp[i] = Hex8Element.shape_functions(xi, eta, zeta)
        N_inv = np.linalg.inv(N_gp)

        cls._hex8_extrap_cache = (N_inv, node_to_gp)
        return cls._hex8_extrap_cache

    def _extrapolate_to_nodes(self, gauss_values: np.ndarray) -> np.ndarray:
        """Extrapolate values from 2x2x2 Gauss points to 8 hex nodes.

        Uses the inverse of the shape function matrix evaluated at the
        Gauss points.  For a hex8 element with 2x2x2 Gauss quadrature
        this is exact for any tri-linear field.

        **Ordering contract** (see :meth:`_build_hex8_extrapolation` for
        the full discussion of issue #51):

        - ``gauss_values`` is indexed by **lexicographic** Gauss-point
          order — exactly what ``Hex8Element._gauss_points`` and
          ``stress_at_gauss_points`` produce.
        - The returned nodal array is indexed by **VTK / Abaqus** node
          order — exactly what ``mesh.elements`` connectivity expects.

        These two orders are *not* the same; a naive "pair GP i with node
        i" simplification would scramble the result spatially.

        Parameters
        ----------
        gauss_values : np.ndarray
            Shape ``(8,)`` or ``(8, n_components)`` values at the 8 Gauss
            points, in lexicographic order.

        Returns
        -------
        np.ndarray
            Same shape as ``gauss_values`` — extrapolated nodal values in
            VTK node order.

        Notes
        -----
        The Gauss points for the 2-point rule are at
        ``xi = +/- 1/sqrt(3) ~ +/- 0.57735``.  The extrapolation matrix
        is constant across all hex8 elements and is cached on the class.
        """
        gauss_values = np.asarray(gauss_values, dtype=float)
        squeeze = False
        if gauss_values.ndim == 1:
            gauss_values = gauss_values[:, np.newaxis]
            squeeze = True
        if gauss_values.shape[0] != 8:
            raise ValueError(
                "gauss_values must have 8 rows (one per 2x2x2 Gauss point), "
                f"got shape {gauss_values.shape}."
            )

        N_inv, _node_to_gp = self._build_hex8_extrapolation()
        nodal = np.asarray(N_inv @ gauss_values)
        return nodal[:, 0] if squeeze else nodal
