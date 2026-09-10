"""High-level analysis pipeline for wrinkled composite laminates.

Provides :class:`WrinkleAnalysis`, a one-stop orchestrator that chains:

1. Laminate definition (material + stacking sequence)
2. Wrinkle geometry and morphology configuration
3. Mesh generation
4. Static FE solve
5. Failure evaluation
6. Optional buckling analysis
7. Optional Monte Carlo / Jensen gap statistics

This module is the primary user-facing entry point for typical workflows.
All lower-level modules (core, elements, solver, failure, statistics) are
accessed through this pipeline.

Examples
--------
Minimal compression analysis::

    >>> from wrinklefe.analysis import WrinkleAnalysis, AnalysisConfig
    >>> config = AnalysisConfig(
    ...     amplitude=0.366, wavelength=16.0, width=12.0,
    ...     morphology="concave", loading="compression",
    ... )
    >>> analysis = WrinkleAnalysis(config)
    >>> result = analysis.run()  # doctest: +SKIP
    >>> print(result.summary())  # doctest: +SKIP

References
----------
- Elhajjar, R. (2025). Scientific Reports, 15:25977 (fat-tail statistics).
- Jin, L. et al. (2026). Thin-Walled Structures, 219:114237 (wrinkle geometry).
- Budiansky, B. & Fleck, N.A. (1993). J. Mech. Phys. Solids, 41(1), 183-211.
"""

from __future__ import annotations

import itertools
import json
import logging
import math
import os
from collections.abc import Callable, Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass, fields, replace
from pathlib import Path
from typing import Any, Literal, cast

import numpy as np

from wrinklefe.core.cohesive_mesh import insert_cohesive_interface
from wrinklefe.core.laminate import Laminate, LoadState
from wrinklefe.core.layup import validate_ply_angle
from wrinklefe.core.material import MaterialLibrary, OrthotropicMaterial
from wrinklefe.core.mesh import MeshData, WrinkleMesh
from wrinklefe.core.morphology import (
    MORPHOLOGY_PHASES,
    SINGLE_WRINKLE_MODES,
    WrinkleConfiguration,
    WrinklePlacement,
)
from wrinklefe.core.penetration_gate import (
    GATE_PRESETS,
    GateParameters,
    penetration_gate_kd,
)
from wrinklefe.core.transforms import rotate_stiffness_3d
from wrinklefe.core.wrinkle import (
    GaussianSinusoidal,
    WrinkleProfile,
    WrinkleSurface3D,
)
from wrinklefe.elements.cohesive8 import (
    Cohesive8Element,
    CohesiveProperties,
    make_initial_state,
)
from wrinklefe.failure.delamination import build_delamination_report
from wrinklefe.failure.evaluator import FailureEvaluator, LaminateFailureReport
from wrinklefe.solver.assembler import GlobalAssembler
from wrinklefe.solver.boundary import BoundaryHandler
from wrinklefe.solver.nonlinear import NewtonRaphsonSolver
from wrinklefe.solver.results import FieldResults
from wrinklefe.solver.static import StaticSolver

logger = logging.getLogger(__name__)

# Analytical damage model constants (Section 6 of CLAUDE.md)
_D0 = 0.15       # Base damage coefficient
_BETA_ANGLE = 3.0  # Angle sensitivity
_THETA_CRIT = 0.1  # Critical angle (rad)
_A_REF = 0.183    # Reference amplitude (1 ply thickness, mm)

# Sanity bound on the cure-residual temperature change (deg C, issue
# #273).  ``delta_T`` is a *change from the stress-free state*, not an
# absolute temperature, so realistic magnitudes are O(100-200 deg C);
# anything past this is almost always an absolute temperature typed into
# the wrong field.
_DELTA_T_MAX = 1000.0

# Number of x-integration points for profile-proportional knockdown
_N_PROFILE_PTS = 500

# Confinement model constants
# Calibrated with CLT-weighted BF against Elhajjar (2025), T700/2510, and
# Mukhopadhyay (2015) blocked-layup compression cases.
_GAMMA_Y_UD = 0.032   # UD matrix yield strain (no confinement)
_ALPHA_CONF = 0.050   # confinement boost coefficient (per off-axis-neighbour score)
# Block-size penalty: each additional 0-deg ply in a consecutive run beyond
# the first reduces gamma_Y_eff by this amount. Captures the empirical
# observation that blocked layups such as Mukhopadhyay's [0_2] (effectively
# [0_4] across the symmetry plane) kink more easily than the neighbour-
# counting confinement score alone predicts: inner 0-deg faces of a block
# are bracketed by another 0-deg ply that does not constrain lateral
# expansion of the kink band. Only applied when at least one off-axis ply
# exists so pure UD ([0]_n) remains at the UD calibration point.
_BETA_BLOCK = 0.010   # per-extra-ply block penalty on gamma_Y_eff
# Lower bound on gamma_Y_eff so a long 0-deg block cannot drive it negative
# or arbitrarily close to zero (which would otherwise produce a degenerate
# Budiansky-Fleck knockdown).
_GAMMA_Y_FLOOR = _GAMMA_Y_UD / 2.0  # = 0.016 (UD half-strain floor)


def _confined_fraction(angles: list[float], tol: float = 5.0) -> float:
    """Weighted confinement fraction for 0-degree plies.

    Each 0-degree ply is scored by how many of its neighbors are off-axis:
        - Both neighbors off-axis (or free surface): score = 1.0
        - One neighbor off-axis, one neighbor 0-deg:  score = 0.5
        - Both neighbors 0-deg (block interior):      score = 0.0

    This partial-confinement model correctly handles both dispersed
    layups (e.g., [0/45/90/-45]) and blocked layups (e.g., [0_2/90_2]).

    Used to compute the effective matrix yield strain:
        gamma_Y_eff = gamma_Y_UD + alpha * f_confined

    Calibrated (with CLT-weighted compression) against:
        - UD [0]_n:       f=0.000, gamma_Y=0.032
        - Mukhopadhyay (2015): f=0.417, gamma_Y=0.053  (blocked [0_2])
        - Elhajjar (2025): f=0.833, gamma_Y=0.074  (dispersed)
    """
    n_0 = sum(1 for a in angles if abs(a) < tol)
    if n_0 == 0:
        return 0.0
    score = 0.0
    for i, a in enumerate(angles):
        if abs(a) < tol:
            left_ok = (abs(angles[i - 1]) > tol) if i > 0 else True
            right_ok = (abs(angles[i + 1]) > tol) if i < len(angles) - 1 else True
            if left_ok and right_ok:
                score += 1.0
            elif left_ok or right_ok:
                score += 0.5
            # else: both neighbors are 0-deg → score += 0.0
    return score / n_0


def _effective_gamma_Y(angles: list[float]) -> float:
    """Compute layup-dependent effective matrix yield strain.

    Three-parameter model::

        gamma_Y_eff = max(
            gamma_Y_UD + alpha_conf * f_confined
                       - beta_block * max(n_block_max - 1, 0),
            gamma_Y_floor,
        )

    where:

    * ``f_confined`` is the weighted confinement fraction of 0-degree
      plies (0 = unconfined, 1 = fully interspersed; see
      :func:`_confined_fraction`).  The linear ``alpha_conf`` term
      captures the constraint that off-axis plies impose on kink-band
      lateral expansion in multidirectional laminates.
    * ``n_block_max`` is the longest run of consecutive 0-degree plies
      (see :func:`_max_consecutive_zero_plies`).  The ``beta_block``
      term penalises long 0-deg blocks: each additional ply inside a
      block beyond the first contributes another increment of
      lateral-expansion freedom that the neighbour-counting confinement
      score does not capture.  Inner 0-deg faces of a block are
      bracketed by another 0-deg ply that does not constrain kink-band
      lateral expansion, so the matrix yields at a lower applied shear
      strain in blocked layups than in dispersed layups with the same
      ``f_confined``.

    The block penalty is only applied when at least one off-axis ply
    exists in the layup.  Pure UD ``[0]_n`` would otherwise be driven
    below the calibration point by the penalty term; with the guard, UD
    retains ``gamma_Y_eff = gamma_Y_UD = 0.032`` regardless of ``n``.

    The result is floored at ``_GAMMA_Y_FLOOR`` (= gamma_Y_UD / 2) so
    very thick 0-blocks cannot drive ``gamma_Y_eff`` arbitrarily close
    to zero, which would otherwise produce a degenerate Budiansky-Fleck
    knockdown.

    With CLT-weighted compression (``KD_lam = f0 * KD_BF + (1 - f0)``),
    the confinement effect is separated from load redistribution.

    Calibration anchors (three-parameter model, beta_block = 0.010):

    ======================================  =======  ===============  ==========
    Layup                                   ``f``    ``n_block_max``  ``gamma_Y``
    ======================================  =======  ===============  ==========
    UD ``[0]_n``                            ~0.13    n (guard skips)  0.032
    Mukhopadhyay ``[..../0_2]_3s``          ~0.42    4 (block of 4    ~0.023
                                                     at symmetry
                                                     plane)
    Elhajjar ``[0/45/90/-45/0/45/-45/0]_s`` ~0.83    2 (only at the   ~0.064
                                                     symmetry plane)
    ======================================  =======  ===============  ==========
    """
    fc = _confined_fraction(angles)
    # Guard: pure UD has no off-axis plies; the block penalty would
    # otherwise drive its gamma_Y below the calibration point.
    n_off_axis = sum(1 for a in angles if abs(a) >= 5.0)
    if n_off_axis == 0:
        return _GAMMA_Y_UD
    n_block_max = _max_consecutive_zero_plies(angles)
    block_penalty = _BETA_BLOCK * max(n_block_max - 1, 0)
    gamma_Y = _GAMMA_Y_UD + _ALPHA_CONF * fc - block_penalty
    return max(gamma_Y, _GAMMA_Y_FLOOR)


def _profile_proportional_kd(
    amplitude: float,
    wavelength: float,
    width: float,
    domain_length: float,
    ply_thickness: float,
    n_plies: int,
    gamma_Y: float,
    theta_max: float,
    *,
    morphology_factor: float = 1.0,
    through_thickness_decay: bool = True,
    z_position_fraction: float = 0.5,
    decay_scale: float | None = None,
    decay_floor: float = 0.0,
    kink_band_quadratic_coeff: float = 0.0,
) -> float:
    """Budiansky-Fleck knockdown averaged over the wrinkle profile.

    Instead of applying a single peak-angle knockdown to all plies,
    this function computes the local fibre angle at every (x, z) point
    in the laminate and averages the Budiansky-Fleck response:

        KD_lam = (1/N) * sum_p [ (1/L_s) * int KD(x, z_p) dx ]

    where the local angle at position (x, z_p) is:

        theta(x, z_p) = |dz_w/dx| * M_f * Phi(z_p)

    with:
        |dz_w/dx|  = slope of the GaussianSinusoidal wrinkle profile
        M_f        = morphology factor (accounts for dual-wrinkle interaction)
        Phi(z_p)   = exp(-(z_p - z_c)^2 / (2 * sigma^2))  (through-thickness
                     decay, standard Gaussian convention with sigma =
                     ``decay_scale``).  This differs from the legacy form
                     ``exp(-(z_p - T/2)^2 / A^2)`` (no factor of 2 — implicit
                     ``sqrt(2)*A`` scale): the new form makes the decay
                     scale an explicit standard deviation.
        z_c        = z_position_fraction * T (laminate thickness)

    The decay scale defaults to ``max(wavelength / 2, amplitude)`` when
    not provided.  The wrinkle's longitudinal extent (set by the
    wavelength) is the physical scale over which a buried wrinkle
    perturbs the through-thickness fibre orientation field: a short-
    wavelength wrinkle decays through only a few plies; a long-wavelength
    wrinkle reaches further.  The legacy ``A``-based scale almost always
    falls inside this default for the calibrated datasets, but for thick
    UD laminates with long wavelengths (e.g. Li 2024) the legacy form
    confined the wrinkle effect to just the midplane plies, leaving
    the laminate KD near 1 even at high amplitudes.

    When *through_thickness_decay* is False, Phi(z_p) = 1 for all plies
    (all plies see the same longitudinal profile).  This is appropriate
    for dual-wrinkle morphologies (stack/convex/concave) where the wrinkle
    extends through the full thickness.

    The per-point knockdown uses the Argon-Fleck quadratic extension of
    the Budiansky-Fleck closed form::

        r = theta_eff / gamma_Y
        KD = 1.0 / (1.0 + r + c_AF * r^2)

    With ``c_AF = 0`` (default) the legacy linear form is recovered.
    Non-zero ``c_AF`` improves the high-angle response (theta > ~20 deg)
    where the linear form systematically over-predicts strength.

    Parameters
    ----------
    amplitude : float
        Wrinkle amplitude A [mm].
    wavelength : float
        Full sinusoidal wavelength lambda [mm].
    width : float
        Gaussian envelope half-width w [mm].
    domain_length : float
        Specimen / domain length L_s [mm].
    ply_thickness : float
        Ply thickness [mm].
    n_plies : int
        Total number of plies.
    gamma_Y : float
        Effective matrix yield shear strain.
    theta_max : float
        Maximum unattenuated fibre angle [rad].
    morphology_factor : float
        Morphology factor M_f that scales the effective angle to account
        for dual-wrinkle interaction (convex < 1.0, concave > 1.0,
        stack = 1.0, graded = 1.0).  Default 1.0.
    through_thickness_decay : bool
        If True (default), apply Gaussian through-thickness decay centred
        at ``z_position_fraction * T`` with scale ``decay_scale`` (see
        below).  If False, all plies see the full wrinkle angle profile
        (Phi = 1).
    z_position_fraction : float
        Fraction of the laminate thickness at which the wrinkle through-
        thickness decay is centred.  ``0.5`` (default) centres the decay
        at the midplane, reproducing the legacy behaviour; ``0.0`` and
        ``1.0`` place the decay centre at the bottom and top surfaces,
        respectively.  Only consulted when ``through_thickness_decay`` is
        True.
    decay_scale : float or None
        Through-thickness Gaussian standard deviation [mm].  When None
        (default) the auto formula ``max(wavelength / 2, amplitude)`` is
        used.  Must be strictly positive when provided.
    decay_floor : float
        Minimum fraction of the wrinkle angle retained at any ply
        (issue #254): the through-thickness term becomes
        ``decay_floor + (1 - decay_floor) * raw`` with ``raw`` the
        Gaussian above — the same floor semantics the tension graded
        path applies, so a sign-flipped load sees the same envelope.
        ``0.0`` (default) reproduces the legacy pure-Gaussian decay
        bit-for-bit; ``1.0`` disables the decay (every ply sees the
        full angle).  Caller is responsible for the [0, 1] range
        (``AnalysisConfig`` validates it).
    kink_band_quadratic_coeff : float
        Argon-Fleck quadratic coefficient ``c_AF`` (dimensionless).
        Default 0.0 (legacy linear BF).  Must be >= 0.

    Returns
    -------
    float
        Profile-averaged BF knockdown factor (0, 1].
    """
    T = n_plies * ply_thickness
    z_center = z_position_fraction * T
    L_s = domain_length

    # Resolve the through-thickness decay scale.  The auto default uses
    # the wrinkle's longitudinal extent (lambda / 2) so long-wavelength
    # wrinkles reach further through the thickness; falls back to A so
    # short-wavelength wrinkles do not collapse to an unphysically small
    # decay scale.
    if decay_scale is None:
        sigma = max(wavelength / 2.0, amplitude)
    else:
        sigma = float(decay_scale)
    sigma_sq2 = 2.0 * sigma * sigma
    c_AF = float(kink_band_quadratic_coeff)

    # Longitudinal profile: compute |dz/dx| at each x-point
    x = np.linspace(-L_s / 2.0, L_s / 2.0, _N_PROFILE_PTS)

    # z(x) = A * exp(-x^2/w^2) * cos(2*pi*x/lambda)
    gauss_env = np.exp(-(x ** 2) / (width ** 2))
    cos_term = np.cos(2.0 * np.pi * x / wavelength)

    # Slope via analytical derivative (more accurate than np.gradient)
    sin_term = np.sin(2.0 * np.pi * x / wavelength)
    dz_dx = amplitude * gauss_env * (
        (-2.0 * x / (width ** 2)) * cos_term
        - (2.0 * np.pi / wavelength) * sin_term
    )
    theta_x = np.abs(np.arctan(dz_dx)) * morphology_factor  # M_f-scaled angle

    # Average over plies and x-positions
    kd_sum = 0.0
    for p in range(n_plies):
        z_p = (p + 0.5) * ply_thickness
        if through_thickness_decay:
            raw_p = np.exp(-((z_p - z_center) ** 2) / sigma_sq2)
            phi_p = decay_floor + (1.0 - decay_floor) * raw_p
        else:
            phi_p = 1.0
        theta_xz = theta_x * phi_p  # local angle at (x, z_p)
        r = theta_xz / gamma_Y
        kd_xz = 1.0 / (1.0 + r + c_AF * r * r)
        kd_sum += np.mean(kd_xz)

    return kd_sum / n_plies


def _is_unidirectional(angles: Sequence[float], tol: float = 5.0) -> bool:
    """True when every ply is a 0-degree (axial) ply within ``tol`` degrees.

    0 and 180 degrees are equivalent fibre directions. Used to dispatch
    the analytical modulus knockdown to its scalar unidirectional fast path
    (:func:`_profile_modulus_knockdown`); multidirectional layups take the
    laminate generalization (:func:`_laminate_modulus_knockdown`) instead.
    """
    if not angles:
        return False
    for a in angles:
        off = abs(float(a)) % 180.0
        if min(off, 180.0 - off) > tol:
            return False
    return True


def _profile_modulus_knockdown(
    amplitude: float,
    wavelength: float,
    width: float,
    domain_length: float,
    ply_thickness: float,
    n_plies: int,
    E1: float,
    E2: float,
    G12: float,
    nu12: float,
    *,
    morphology_factor: float = 1.0,
    through_thickness_decay: bool = True,
    z_position_fraction: float = 0.5,
    decay_scale: float | None = None,
    decay_floor: float = 0.0,
) -> float:
    r"""Axial Young's-modulus knockdown ``E_x / E_x0`` for a wavy UD laminate.

    A Classical-Lamination-Theory series-average of the off-axis lamina
    modulus over the wrinkle profile — the same off-axis-compliance
    integration as Hsiao & Daniel (1996, Compos. Sci. Technol. 56:581).
    At every ``(x, z_p)`` the local fibre tilt is the same field used by
    :func:`_profile_proportional_kd`,
    ``theta(x, z_p) = |dz_w/dx| * M_f * Phi(z_p)``, and the off-axis axial
    modulus of a 0-degree ply tilted by ``theta`` is::

        1/E_x(theta) = cos^4/E1 + (1/G12 - 2 nu12/E1) cos^2 sin^2 + sin^4/E2

    The plies share the membrane strain, so the section modulus at a
    station is the through-thickness mean
    ``E_sec(x) = <E_x(theta(x, z_p))>_p``; along the load direction the
    compliances add, so the effective modulus is the series (harmonic)
    average ``E_eff = 1 / <1/E_sec(x)>_x``. The knockdown is
    ``E_eff / E1`` (the pristine UD axial modulus is ``E1``).

    Linear-elastic, hence loading-independent. This closed form assumes a
    0-degree base ply, so it covers single-wrinkle **unidirectional** axial
    layups exactly (see :func:`_is_unidirectional`). Multidirectional and
    multi-wrinkle configurations are handled by the laminate generalization
    :func:`_laminate_modulus_knockdown`, which reduces to this result for
    ``[0]_n``; this scalar form is retained as the UD fast path so the
    pinned UD baselines stay numerically identical.

    Returns
    -------
    float
        Axial-modulus knockdown in ``(0, 1]``.
    """
    T = n_plies * ply_thickness
    z_center = z_position_fraction * T
    L_s = domain_length
    if decay_scale is None:
        sigma = max(wavelength / 2.0, amplitude)
    else:
        sigma = float(decay_scale)
    sigma_sq2 = 2.0 * sigma * sigma

    # Same longitudinal angle field as the strength profile-average.
    x = np.linspace(-L_s / 2.0, L_s / 2.0, _N_PROFILE_PTS)
    gauss_env = np.exp(-(x ** 2) / (width ** 2))
    cos_term = np.cos(2.0 * np.pi * x / wavelength)
    sin_term = np.sin(2.0 * np.pi * x / wavelength)
    dz_dx = amplitude * gauss_env * (
        (-2.0 * x / (width ** 2)) * cos_term
        - (2.0 * np.pi / wavelength) * sin_term
    )
    theta_x = np.abs(np.arctan(dz_dx)) * morphology_factor

    shear_coupling = 1.0 / G12 - 2.0 * nu12 / E1
    Ex_sum = np.zeros_like(x)
    for p in range(n_plies):
        z_p = (p + 0.5) * ply_thickness
        if through_thickness_decay:
            raw_p = np.exp(-((z_p - z_center) ** 2) / sigma_sq2)
            phi_p = decay_floor + (1.0 - decay_floor) * raw_p
        else:
            phi_p = 1.0
        theta = theta_x * phi_p
        c2 = np.cos(theta) ** 2
        s2 = np.sin(theta) ** 2
        inv_Ex = c2 * c2 / E1 + shear_coupling * c2 * s2 + s2 * s2 / E2
        Ex_sum += 1.0 / inv_Ex
    E_sec = Ex_sum / n_plies              # through-thickness mean at each x
    E_eff = 1.0 / np.mean(1.0 / E_sec)    # series-average along x
    return float(E_eff / E1)


def _plane_stress_qbar_tilted(
    stiffness_3d: np.ndarray, phi_rad: float, theta_rad: float
) -> np.ndarray:
    """Plane-stress reduced stiffness for a ply rotated in-plane *and* tilted.

    The ply's full 3D stiffness ``[C]`` (material axes) is rotated by the
    in-plane orientation ``phi`` about the through-thickness ``z`` axis and
    by the out-of-plane wrinkle tilt ``theta`` about the transverse ``y``
    axis (the combination of the ply angle with the local wrinkle slope is
    a composition of these two principal-axis rotations,
    :func:`wrinklefe.core.transforms.rotate_stiffness_3d`). The rotated 3D
    stiffness is then statically condensed to plane stress
    (``sigma_33 = tau_23 = tau_13 = 0``) onto the in-plane Voigt indices
    ``(11, 22, 12)`` so it can enter the laminate membrane (A) matrix like a
    standard ``Q-bar``.

    For ``theta = 0`` this returns the ordinary CLT ``Q-bar(phi)`` and for
    ``phi = 0`` the axial term ``1/inv(Q-bar)[0, 0]`` equals the off-axis
    formula used by :func:`_profile_modulus_knockdown`, so the laminate
    knockdown reduces exactly to the UD scalar result for ``[0]_n``.
    """
    c_rot = rotate_stiffness_3d(stiffness_3d, phi_rad, axis="z")
    c_rot = rotate_stiffness_3d(c_rot, theta_rad, axis="y")
    # Voigt split: in-plane (11, 22, 12) vs out-of-plane (33, 23, 13).
    ip = [0, 1, 5]
    op = [2, 3, 4]
    c_ii = c_rot[np.ix_(ip, ip)]
    c_io = c_rot[np.ix_(ip, op)]
    c_oi = c_rot[np.ix_(op, ip)]
    c_oo = c_rot[np.ix_(op, op)]
    return np.asarray(c_ii - c_io @ np.linalg.solve(c_oo, c_oi))


def _laminate_modulus_knockdown(
    slope_field: np.ndarray,
    ply_decays: np.ndarray,
    angles: Sequence[float],
    stiffness_3d: np.ndarray,
    ply_thickness: float,
    E_x0: float,
) -> float:
    r"""Axial-modulus knockdown ``E_x / E_x0`` for an arbitrary wavy laminate.

    Generalizes :func:`_profile_modulus_knockdown` from a single 0-degree
    base ply to an arbitrary stacking sequence and to an already-composed
    multi-wrinkle slope field, via a CLT membrane (A-matrix) series-average
    — the laminate form of the Hsiao & Daniel (1996) off-axis-compliance
    integration.

    At every longitudinal station ``x`` each ply ``p`` carries a local fibre
    tilt ``theta(x, z_p) = arctan|sum_w (dz_w/dx) Phi_w(z_p)|`` (the composed
    slope field, decayed through the thickness exactly as the FE composes it
    in :meth:`WrinkleConfiguration.fiber_angles_at_nodes` — "compose then
    differentiate"). The ply's plane-stress stiffness is rotated to its
    in-plane angle ``phi_p`` plus that tilt
    (:func:`_plane_stress_qbar_tilted`) and summed into the membrane
    stiffness ``A(x) = sum_p Q-bar_p(phi_p, theta) * t``. The plies share the
    membrane strain, so the section axial modulus is
    ``E_x_section(x) = 1 / (a11(x) * T)`` with ``a11 = inv(A)[0, 0]`` and
    ``T`` the total thickness; along ``x`` the compliances add, giving the
    series (harmonic) average ``E_eff = 1 / <1/E_x_section(x)>_x``.

    The knockdown is ``E_eff / E_x0`` where ``E_x0`` is the **flat** laminate
    axial modulus (``Laminate.Ex``). Because the flat-laminate A-matrix is
    recovered exactly when every tilt is zero, the knockdown is exactly
    ``1.0`` for a degenerate (zero-amplitude) wrinkle and lies in ``(0, 1]``
    otherwise. Off-axis plies, already carrying little axial load, are
    insensitive to the axial misalignment, so a multidirectional layup is
    knocked down less than the same wrinkle in pure UD.

    Linear-elastic, hence loading-independent.

    The per-ply tilt is built as ``sum_w slope_field[w] * ply_decays[p, w]``
    — the signed per-wrinkle slopes, each scaled by that wrinkle's
    through-thickness decay, summed *before* the angle is taken. A single
    wrinkle is just the one-entry (``N_w == 1``) case of this composition.

    Parameters
    ----------
    slope_field : np.ndarray
        Shape ``(N_w, N_x)`` per-wrinkle longitudinal slope ``dz_w/dx``
        *before* the through-thickness decay, evaluated along the wrinkle.
    ply_decays : np.ndarray
        Shape ``(n_plies, N_w, N_x)`` through-thickness decay ``Phi_w(z_p)``
        applied to each wrinkle's slope for each ply.
    angles : Sequence[float]
        Ply in-plane orientations ``phi_p`` in **degrees**.
    stiffness_3d : np.ndarray
        The common ply ``6x6`` 3D stiffness ``[C]`` (material axes).
    ply_thickness : float
        Uniform ply thickness [mm].
    E_x0 : float
        Pristine flat-laminate axial modulus ``Laminate.Ex`` [MPa].

    Returns
    -------
    float
        Axial-modulus knockdown in ``(0, 1]``.
    """
    phis = [math.radians(float(a)) for a in angles]
    n_plies = len(phis)
    total_thickness = n_plies * ply_thickness

    # Compose the per-wrinkle slopes (each decayed for this ply) then
    # take the local tilt angle: theta(p, x) = arctan|sum_w dz_w/dx * Phi_w|.
    slope_px = np.einsum("wx,pwx->px", slope_field, ply_decays)
    theta_px = np.arctan(np.abs(slope_px))  # (n_plies, n_x)

    # The in-plane (z-axis) rotation depends only on the ply angle, so
    # rotate the base stiffness once per ply; the wrinkle tilt (y-axis)
    # rotation and the plane-stress condensation are batched over the
    # whole (ply, x) grid (issue #301: the previous per-(ply, x) Python
    # loop made every multidirectional analytical run take seconds).
    c_phi = np.stack(
        [rotate_stiffness_3d(stiffness_3d, phi, axis="z") for phi in phis]
    )  # (n_plies, 6, 6)

    # Batched y-axis stress-transformation matrices T_sigma(theta) and
    # the engineering-strain counterparts T_eps = R T_sigma R^-1
    # (row/column Reuter scaling), mirroring
    # wrinklefe.core.transforms.stress/strain_transformation_3d.
    def _tsigma_y(cos_t: np.ndarray, sin_t: np.ndarray) -> np.ndarray:
        c2, s2, sc = cos_t * cos_t, sin_t * sin_t, sin_t * cos_t
        t = np.zeros(cos_t.shape + (6, 6), dtype=float)
        t[..., 0, 0] = c2
        t[..., 0, 2] = s2
        t[..., 0, 4] = -2.0 * sc
        t[..., 1, 1] = 1.0
        t[..., 2, 0] = s2
        t[..., 2, 2] = c2
        t[..., 2, 4] = 2.0 * sc
        t[..., 3, 3] = cos_t
        t[..., 3, 5] = sin_t
        t[..., 4, 0] = sc
        t[..., 4, 2] = -sc
        t[..., 4, 4] = c2 - s2
        t[..., 5, 3] = -sin_t
        t[..., 5, 5] = cos_t
        return t

    cos_t = np.cos(theta_px)
    sin_t = np.sin(theta_px)
    t_sig = _tsigma_y(cos_t, sin_t)
    # Rotation transformations invert by angle negation:
    # T_sigma(theta)^-1 == T_sigma(-theta) — a matmul instead of a
    # batched 6x6 solve.
    t_sig_inv = _tsigma_y(cos_t, -sin_t)
    reuter = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0])
    t_eps = t_sig * (reuter[:, None] / reuter[None, :])

    # C_rot = T_sigma^-1 @ C_phi @ T_eps, batched over (ply, x).
    c_rot = t_sig_inv @ (c_phi[:, None, :, :] @ t_eps)

    # Plane-stress condensation onto the in-plane Voigt indices
    # (11, 22, 12): Q-bar = C_ii - C_io @ C_oo^-1 @ C_oi.
    ip = np.array([0, 1, 5])
    op = np.array([2, 3, 4])
    c_ii = c_rot[..., ip[:, None], ip[None, :]]
    c_io = c_rot[..., ip[:, None], op[None, :]]
    c_oi = c_rot[..., op[:, None], ip[None, :]]
    c_oo = c_rot[..., op[:, None], op[None, :]]
    q_bar = c_ii - c_io @ np.linalg.solve(c_oo, c_oi)

    # Membrane A(x) = sum_p Q-bar_p * t; shared membrane strain gives
    # E_section(x) = 1 / (inv(A)[0,0] * T).
    a_mat = q_bar.sum(axis=0) * ply_thickness  # (n_x, 3, 3)
    inv_E_section = np.linalg.inv(a_mat)[:, 0, 0] * total_thickness

    # E_section(x) = 1 / inv_E_section; series-average the compliance.
    e_eff = 1.0 / float(np.mean(inv_E_section))
    return float(e_eff / E_x0)


def _max_consecutive_zero_plies(angles: list[float], tol: float = 5.0) -> int:
    """Find maximum number of consecutive 0-degree plies in a layup.

    Used by the tension OOP model to determine the effective curved-beam
    thickness h_eff = n_adj × t_ply for interlaminar stress prediction.

    Returns the true maximum consecutive count, which is ``0`` when the
    layup contains no 0-degree plies.  Callers must guard against the
    zero case (e.g. for the curved-beam OOP path, ``n_adj == 0`` means
    there is no continuous 0-degree block to develop interlaminar stress).
    """
    max_count = 0
    count = 0
    for a in angles:
        if abs(a) < tol:
            count += 1
            max_count = max(max_count, count)
        else:
            count = 0
    return max_count


# ======================================================================
# Public helper: amplitude → wavelength estimation
# ======================================================================

def estimate_wavelength_from_amplitude(
    amplitude: float,
    *,
    K_lambda: float = 19.9,
    lambda_min: float = 8.2,
    lambda_ref: float = 0.366,
    scaling: str = "sqrt",
) -> float:
    """Estimate wrinkle wavelength from amplitude when lambda is not measured.

    Some validation datasets report only the wrinkle amplitude ``A`` and
    require an external rule to recover the wavelength ``lambda`` needed
    by the Budiansky-Fleck peak-fibre-angle model
    ``theta_max = arctan(2*pi*A/lambda)``.  Two scaling rules are
    supported:

    * ``"linear"`` (legacy):
        ``lambda = K_lambda * A``.  Reproduces the original convention
        documented in §1.4 of ``VALIDATION_DATA.md`` and used by older
        validation harnesses.  Because ``theta_max`` reduces to
        ``arctan(2*pi/K_lambda)``, the predicted peak fibre angle is
        *constant* in ``A`` under this rule.  That is unphysical for
        severe wrinkles, where larger amplitudes produce steeper local
        fibre rotations and stronger compressive knockdowns.

    * ``"sqrt"`` (default, recommended):
        ``lambda = K_lambda * sqrt(A * lambda_ref)``.  At the reference
        amplitude ``A = lambda_ref`` the sqrt rule matches the legacy
        linear rule exactly (``lambda = K_lambda * lambda_ref``), so
        the calibration of mild wrinkles is preserved.  For
        ``A > lambda_ref`` the wavelength grows sub-linearly with
        amplitude, so ``theta_max = arctan(2*pi*A/lambda)`` increases
        monotonically with ``A`` and the model captures the experimental
        knockdown collapse seen at high D/T in the Elhajjar (2025)
        compression dataset.

    Both rules are clamped from below at ``lambda_min`` so that tiny
    amplitudes do not produce vanishingly small wavelengths.

    Parameters
    ----------
    amplitude : float
        Wrinkle amplitude *A* [mm].  Must be non-negative.
    K_lambda : float, optional
        Slope coefficient. Defaults to ``19.9`` (Elhajjar 2025
        T700/2510 calibration).
    lambda_min : float, optional
        Lower clamp on the returned wavelength [mm].  Defaults to
        ``8.2`` mm (Elhajjar 2025).
    lambda_ref : float, optional
        Reference amplitude [mm] at which the sqrt scaling is anchored
        to the legacy linear rule.  Defaults to ``0.366`` mm (two ply
        thicknesses in the Elhajjar T700/2510 layup).  Ignored when
        ``scaling == "linear"``.
    scaling : {"sqrt", "linear"}, optional
        Scaling rule selector.  Defaults to ``"sqrt"``.

    Returns
    -------
    float
        Wavelength ``lambda`` in mm, lower-bounded by ``lambda_min``.

    Raises
    ------
    ValueError
        If ``scaling`` is not one of ``"sqrt"`` or ``"linear"``.

    Examples
    --------
    Legacy linear rule (constant peak angle):

    >>> estimate_wavelength_from_amplitude(0.5, scaling="linear")
    9.95

    Sub-linear sqrt rule (default):  at the reference amplitude the
    raw rule matches the legacy linear rule exactly
    (``K_lambda * lambda_ref = 19.9 * 0.366 = 7.2834``), which then
    clamps to ``lambda_min = 8.2``:

    >>> estimate_wavelength_from_amplitude(0.366)
    8.2
    """
    if scaling not in ("sqrt", "linear"):
        raise ValueError(
            f"estimate_wavelength_from_amplitude: scaling must be "
            f"'sqrt' or 'linear', got {scaling!r}"
        )
    if scaling == "linear":
        lam = K_lambda * amplitude
    else:  # "sqrt"
        # lambda = K_lambda * sqrt(A * lambda_ref) -- matches linear at
        # A = lambda_ref, grows sub-linearly for A > lambda_ref.  Guard
        # against negative amplitude reaching the sqrt.
        lam = K_lambda * math.sqrt(max(amplitude, 0.0) * lambda_ref)
    return max(lam, lambda_min)


# ======================================================================
# Configuration
# ======================================================================

@dataclass
class WrinkleSpec:
    """Single-wrinkle specification used to assemble multi-wrinkle configs.

    A list of :class:`WrinkleSpec` instances passed to
    :class:`AnalysisConfig` via the ``wrinkles`` field overrides the
    single/dual-wrinkle dispatch in :meth:`WrinkleAnalysis.run`, allowing
    arbitrary N-wrinkle layouts at arbitrary ply interfaces with arbitrary
    phase offsets to be analysed (see Dataset F / Li et al. 2025).

    Parameters
    ----------
    amplitude : float
        Wrinkle half-amplitude *A* [mm], strictly positive.
    wavelength : float
        Wavelength lambda [mm], strictly positive.
    width : float
        Gaussian envelope half-width *w* [mm], strictly positive.
    ply_interface : int
        Ply interface index passed to :class:`WrinklePlacement`. For a
        laminate with *N* plies valid indices are 0 through N-2 inclusive.
    phase_offset : float, optional
        Phase offset phi [rad] relative to the reference wrinkle.
        Default 0.0.
    """

    amplitude: float
    wavelength: float
    width: float
    ply_interface: int
    phase_offset: float = 0.0


#: Schema version stamped into :meth:`AnalysisConfig.to_dict` output and
#: verified by :meth:`AnalysisConfig.from_dict`.  Bump this whenever the
#: serialised config layout changes in a non-round-trippable way so old
#: files fail loudly rather than load into a mismatched schema.
CONFIG_VERSION = 1


def _material_to_jsonable(mat: OrthotropicMaterial | None) -> dict | None:
    """Serialise a material as a library-preset reference or inline custom.

    A material equal (field-for-field) to the like-named library preset is
    written as ``{"preset": name}`` so the file stays compact and tracks
    library updates; any other material — including one that reuses a
    library name but tweaks a property — is written as
    ``{"custom": material.to_dict()}`` so no information is lost.
    """
    if mat is None:
        return None
    try:
        preset = MaterialLibrary().get(mat.name)
    except KeyError:
        preset = None
    if preset is not None and preset.to_dict() == mat.to_dict():
        return {"preset": mat.name}
    return {"custom": mat.to_dict()}


def _material_from_jsonable(
    value: dict | None, *, field: str
) -> OrthotropicMaterial | None:
    """Inverse of :func:`_material_to_jsonable`."""
    if value is None:
        return None
    if (
        not isinstance(value, dict)
        or len(value) != 1
        or next(iter(value)) not in ("preset", "custom")
    ):
        raise ValueError(
            f"AnalysisConfig.from_dict: {field} must be null, "
            f"{{'preset': name}}, or {{'custom': {{...}}}}, got {value!r}"
        )
    if "preset" in value:
        try:
            return MaterialLibrary().get(value["preset"])
        except KeyError as exc:
            raise ValueError(
                f"AnalysisConfig.from_dict: {field} references unknown "
                f"material preset {value['preset']!r}"
            ) from exc
    return OrthotropicMaterial.from_dict(value["custom"])


def _gate_to_jsonable(gate: GateParameters | None) -> dict | None:
    """Serialise a penetration gate as a registered-preset reference.

    Only the calibrated presets in
    :data:`wrinklefe.core.penetration_gate.GATE_PRESETS` are serialisable
    (the parameters are material-realization specific and are not meant to
    be hand-authored inline).  An unregistered custom gate raises so the
    lossy write surfaces loudly instead of silently dropping data.
    """
    if gate is None:
        return None
    preset = GATE_PRESETS.get(gate.name)
    if preset is not None and preset == gate:
        return {"preset": gate.name}
    raise ValueError(
        f"AnalysisConfig.to_dict: penetration_gate {gate.name!r} is not a "
        f"registered preset and inline custom gates are not serialisable; "
        f"use one of {sorted(GATE_PRESETS)} or None"
    )


def _gate_from_jsonable(value: dict | None) -> GateParameters | None:
    """Inverse of :func:`_gate_to_jsonable`."""
    if value is None:
        return None
    if not isinstance(value, dict) or set(value) != {"preset"}:
        raise ValueError(
            "AnalysisConfig.from_dict: penetration_gate must be null or "
            f"{{'preset': name}}, got {value!r}"
        )
    name = value["preset"]
    if name not in GATE_PRESETS:
        raise ValueError(
            f"AnalysisConfig.from_dict: unknown penetration_gate preset "
            f"{name!r}; available {sorted(GATE_PRESETS)}"
        )
    return GATE_PRESETS[name]


@dataclass
class AnalysisConfig:
    """Configuration for a wrinkle analysis run.

    Collects all user-specified parameters in a single object that
    can be serialised, compared, and passed to :class:`WrinkleAnalysis`.

    Parameters
    ----------
    amplitude : float
        Wrinkle half-amplitude *A* [mm]: the peak displacement of the
        wrinkled mid-surface from the flat (unwrinkled) reference plane,
        so ``z(x) = A * cos(2*pi*(x - x0) / lambda)`` (modulated by the
        envelope) and the peak-to-trough height is ``2A``. For a
        measured wrinkle, ``A = (z_max - z_min) / 2``. Units: mm.
        Default 0.366 (two ply thicknesses).

        Effect on knockdown: A enters the maximum fibre misalignment
        angle through the closed-form ``theta_max = arctan(2*pi*A /
        lambda)`` used in ``_compute_analytical``, so for small A/lambda
        the peak fibre angle scales nearly linearly with A and amplifies
        the Budiansky-Fleck compressive knockdown.
    wavelength : float
        Spatial period of the cosine carrier *lambda* [mm]: the
        crest-to-crest distance of the underlying
        ``cos(2*pi*(x - x0) / lambda)`` carrier along the longitudinal
        x-direction.  The wavenumber is ``k = 2*pi/lambda`` (1/mm).
        Default 16.0.  Must be > 0.
    width : float
        Longitudinal envelope decay length *w* [mm] about the wrinkle
        centre ``x0``.  Exact meaning is profile-dependent: Gaussian
        1/e length scale in ``exp(-(x - x0)**2 / w**2)``, tapered
        flat-top extent (``|x - x0| < w/2``), or triangular half-base
        (``|x - x0| < w``).  Also used as the transverse (y-direction)
        extent of the wrinkle in the 3-D dual-wrinkle / graded mesh
        deformation.  Default 12.0.  Must be > 0.
    morphology : str
        Morphology name. Five values are accepted; the first three are
        *dual-wrinkle* modes distinguished by phase, the last two are
        *single-wrinkle* modes distinguished by their through-thickness
        amplitude profile:

        - ``'stack'`` (default) — two aligned wrinkles, φ = 0. Linear
          through-thickness decay from the interface plies to zero at
          the outer surfaces. M_f = 1.0 (dual-wrinkle baseline).
        - ``'convex'`` — two wrinkles, φ = +π/2 (interface bulges
          outward). Same through-thickness decay. M_f < 1, least
          damaging in compression.
        - ``'concave'`` — two wrinkles, φ = −π/2 (interface pinches
          inward). Same through-thickness decay. M_f > 1, most
          damaging in compression.
        - ``'uniform'`` — *one* wrinkle, full amplitude on every ply
          (no through-thickness decay). M_f = 1.0 because there is no
          pairwise interaction, but the deformed mesh and per-ply
          fibre-angle field differ from ``'stack'``: every ply
          (including the outer surfaces) carries the full profile.
        - ``'graded'`` — one wrinkle, linear decay from mid-ply to the
          surfaces with the ``decay_floor`` knob (0 = full decay to
          zero, 1 = same as ``'uniform'``).
    phase : float or None
        Explicit dual-wrinkle phase offset phi [radians] between the two
        wrinkle centrelines.  When ``None`` (default), the phase is
        derived from ``morphology`` via :data:`MORPHOLOGY_PHASES`
        (stack=0, convex=+pi/2, concave=-pi/2).  When set to a float, it
        overrides the named-morphology phase, allowing arbitrary
        dual-wrinkle phase offsets to be analysed or swept (e.g.
        between 0 and pi).  Ignored for single-wrinkle morphologies
        (``'uniform'``, ``'graded'``).  Must be finite when set.
    decay_floor : float
        Graded morphology only (dimensionless, in ``[0, 1]``): minimum
        fraction of the wrinkle amplitude retained at the laminate outer
        surfaces.  ``0.0`` (default) means full decay to zero amplitude
        at the surfaces (pure graded); ``1.0`` means no decay
        (equivalent to ``uniform``).  Values outside ``[0, 1]`` are
        rejected by ``__post_init__``.
    wrinkle_z_position : float
        Through-thickness position of the (single-wrinkle) decay centre,
        expressed as a fraction of the laminate thickness *T*.  ``0.0``
        places the wrinkle at the bottom surface, ``0.5`` (default) at
        the midplane (legacy behaviour), and ``1.0`` at the top surface.
        Used by the graded morphology path to shift the Gaussian through-
        thickness decay centre off the midplane and to bias the linear
        per-ply tension grading; mirrors the Above / Middle / Below
        wrinkle locations of Li et al. (2025) Dataset F.  Ignored for the
        ``stack``, ``convex``, ``concave`` and ``uniform`` morphologies,
        whose through-thickness behaviour is set by the morphology itself
        or by the ``interface_1`` / ``interface_2`` interface indices.
        Must be a finite float in ``[0.0, 1.0]``.
    amplitude_profile : {"constant", "gaussian", "linear"}
        Spatially varying in-plane amplitude modulation applied on top
        of the wrinkle's own longitudinal envelope (see
        :class:`~wrinklefe.core.morphology.WrinkleConfiguration`).
        ``"constant"`` (default) preserves the legacy behaviour --
        the wrinkle amplitude *A* is used uniformly across the in-plane
        domain.  ``"gaussian"`` multiplies *A* by ``exp(-(s/d)**2)`` and
        ``"linear"`` by ``max(0, 1 - |s|/d)`` (clipped), where *s* is
        the in-plane coordinate (relative to the wrinkle centre) along
        ``amplitude_profile_axis`` and *d* is
        ``amplitude_profile_decay_length``.
    amplitude_profile_decay_length : float or None
        Length scale *d* (mm) controlling the Gaussian sigma or the
        linear-decay extent.  ``None`` (default) falls back to the
        wrinkle profile's own ``width``, so the amplitude tapers on the
        same length scale as the envelope.  Must be positive when
        provided.  Ignored when ``amplitude_profile == "constant"``.
    amplitude_profile_axis : {"x", "y"}
        In-plane axis along which the amplitude modulation runs.
        Default ``"x"``.  Pick ``"y"`` (transverse axis) for an
        independent in-plane tapering of *A* that does not stack with
        the existing longitudinal envelope on *x*.
    transverse_mode : {"uniform", "gaussian_decay", "sinusoidal_y", "elliptical"}
        Through-width (transverse *y*) wrinkle-surface envelope (#300).
        Default ``"uniform"`` builds the bare x-only wrinkle exactly as
        before (bit-identical). The non-uniform modes wrap the profile in
        a :class:`~wrinklefe.core.wrinkle.WrinkleSurface3D` so the
        amplitude varies across the specimen width: ``"gaussian_decay"``
        decays it toward the edges, ``"sinusoidal_y"`` ripples it across
        the width, and ``"elliptical"`` confines it to a mid-width patch.
        FE-only — a non-uniform mode requires ``analytical_only=False`` and
        is not combinable with ``wrinkles`` (multi-wrinkle) or
        ``enable_czm`` (both rejected at construction).
    transverse_span : float or None
        Specimen width ``span_y`` (mm) seen by the transverse envelope.
        ``None`` (default) tracks ``domain_width`` so the envelope always
        spans the meshed y-extent. Must be > 0 when set. Ignored when
        ``transverse_mode == "uniform"``.
    transverse_width : float or None
        Transverse localization half-width ``width_y`` (mm): the Gaussian
        1/e length for ``"gaussian_decay"`` and the ellipse half-width for
        ``"elliptical"`` (ignored by ``"uniform"``/``"sinusoidal_y"``).
        ``None`` (default) resolves to ``span_y / 4`` — a localized
        mid-width patch. Must be > 0 when set.
    loading : str
        Loading mode: ``'compression'`` or ``'tension'``.
        Default is ``'compression'``.
    material : OrthotropicMaterial or None
        Composite material.  ``None`` uses the default IM7/8552.
    angles : list[float] or None
        Ply angles in degrees.  ``None`` uses a quasi-isotropic
        ``[0/45/-45/90]_3s`` layup (24 plies).
    interface_1 : int or None
        Ply interface for the first wrinkle.  ``None`` (default) auto-
        derives an interior interface from the resolved layup so small
        laminates (< 13 plies) work out of the box: ``interface_1 =
        max(0, n_plies // 2 - 1)``.  For the default 24-ply layup this
        resolves to ``11`` (backwards-compatible).
    interface_2 : int or None
        Ply interface for the second wrinkle.  ``None`` (default) auto-
        derives ``interface_2 = min(n_plies - 1, n_plies // 2)``.  For
        the default 24-ply layup this resolves to ``12``
        (backwards-compatible).
    nx : int
        Mesh divisions in x.  Default 12.
    ny : int
        Mesh divisions in y.  Default 6.
    nz_per_ply : int
        Mesh divisions per ply in z.  Default 1.
    domain_length : float
        Domain length in x [mm].  Default ``3 * wavelength``.
    domain_width : float
        Domain width in y [mm].  Default 20.0.
    applied_strain : float
        Applied nominal strain for displacement-controlled loading.
        Default ``-0.01`` (1 % compression).
    delta_T : float
        Uniform temperature change **from the stress-free (cure)
        state**, in deg C.  Default ``0.0`` (no thermal load).

        **Sign convention** — ``delta_T`` is *T_service - T_stress_free*,
        so a cure cool-down is **negative**: a 177 deg C cure taken to a
        22 deg C service temperature is ``delta_T = -155``.  A positive
        value means the laminate is *hotter* than its stress-free state.
        Getting this sign backwards flips the residual matrix stress from
        tension to compression, so state it explicitly whenever the value
        is reported.

        Both solution paths honour the temperature (issue #273).  The
        CLT / analytical path adds the thermal resultants to the ABD
        solve and recovers ply stresses from the mechanical strain
        (Stage 1); the FE path assembles the element thermal
        initial-strain load vector ``int B^T C eps_th dV`` and subtracts
        the thermal strain during stress recovery (Stage 2), so the
        reported FE stresses, failure indices and retention factors are
        residual-stress inclusive.  The retention baseline is solved at
        the same ``delta_T``, so the retention factor compares like with
        like; the measured-modulus solve deliberately runs at
        ``delta_T = 0`` because a thermal offset in the reaction force
        would corrupt a stiffness measurement.  Must be finite with
        ``|delta_T| <= 1000``.
    load_state : LoadState or None
        General CLT load state (issue #275).  ``None`` (default) keeps
        the legacy ``loading`` / ``applied_strain`` pair as the entire
        load vocabulary and leaves every result bit-identical.

        When set, its resultants define the FE mechanical load directly,
        so biaxial (``Ny``) and in-plane-shear (``Nxy``) states become
        reachable — the states whose transverse and shear components
        drive the matrix failure modes, and which a uniaxial-only surface
        cannot express.  ``applied_strain`` then no longer sets the FE
        boundary conditions; it still drives the separate uniaxial probe
        that measures modulus retention, which is a stiffness property
        and deliberately load-state-independent.

        **Membrane components only.**  ``Mx``/``My``/``Mxy`` are
        rejected: the BC mapping applies curvature as a prescribed
        displacement on the same faces a membrane state loads with
        traction, so the two cannot be superposed.  ``Qx``/``Qy`` have no
        mapping at all.  ``delta_T`` must be zero on the load state — set
        the temperature on :attr:`delta_T` instead, so one quantity has
        one owner.

        Strength under a combined state is reported as a **proportional
        load factor** (:attr:`AnalysisResults.load_state_factor`): the scalar
        the whole state is multiplied by to reach first failure.  That
        reduces to the usual definition for a uniaxial state.
    solver : str
        Linear solver: ``'direct'`` or ``'iterative'``.  Default ``'direct'``.
    iterative_rtol : float
        Relative-residual convergence tolerance for the iterative (CG)
        solver.  Default ``1e-10``.  Must be > 0.  Ignored on the direct
        path.
    iterative_maxiter : int
        Maximum CG iterations for the iterative solver.  Default
        ``10000``.  Must be >= 1.  Ignored on the direct path.
    ilu_drop_tol : float
        Drop tolerance for the ILU preconditioner (``scipy`` ``spilu``);
        the main preconditioner quality/memory knob.  Default ``1e-4``.
        Must be >= 0.  Ignored unless ``preconditioner='ilu'``.
    ilu_fill_factor : float or None
        Upper bound on the ILU fill (``spilu``'s ``fill_factor``).
        ``None`` (default) leaves SciPy's own default in place and is
        only forwarded when set.  Must be >= 1 when given.  Ignored
        unless ``preconditioner='ilu'``.
    preconditioner : str
        Preconditioner for the iterative solver: ``'ilu'`` (default,
        incomplete-LU), ``'jacobi'`` (diagonal — much lower memory for
        very large meshes), or ``'none'`` (unpreconditioned CG).  Ignored
        on the direct path.
    verbose : bool
        Print progress information.  Default ``False``.
    through_thickness_decay_scale : float or None
        Optional override (mm) for the through-thickness Gaussian
        standard deviation used by the graded morphology's profile-
        proportional KD path and the tension graded-averaging block.
        ``None`` (default) triggers the auto formula
        ``max(wavelength / 2, amplitude)``.  Must be > 0 when set.
    kink_band_quadratic_coeff : float
        Argon-Fleck quadratic coefficient ``c_AF`` (dimensionless) in
        the extended Budiansky-Fleck closed form
        ``KD = 1 / (1 + r + c_AF * r**2)`` with
        ``r = theta_eff / gamma_Y_eff``.  Default ``0.0`` recovers the
        legacy linear BF response.  Must be >= 0.
    """

    # Wrinkle geometry
    amplitude: float = 0.366
    wavelength: float = 16.0
    width: float = 12.0

    # Morphology
    morphology: str = "stack"
    # Explicit dual-wrinkle phase offset phi [rad]. None → derive from
    # `morphology` (stack=0, convex=+pi/2, concave=-pi/2). A float
    # overrides the named-morphology phase so arbitrary phases can be
    # analysed/swept. Ignored for single-wrinkle modes (uniform/graded).
    phase: float | None = None
    decay_floor: float = 0.0  # graded mode: min amplitude fraction at surfaces (0–1)
    # Single-wrinkle through-thickness position as a fraction of the laminate
    # thickness (0 = bottom surface, 0.5 = midplane, 1 = top surface). Only
    # consulted by the graded morphology path; ignored for stack/convex/
    # concave/uniform.
    wrinkle_z_position: float = 0.5

    # Spatially varying in-plane amplitude profile (#178 follow-up). Defaults
    # mirror WrinkleConfiguration so the legacy "constant" behaviour is
    # preserved when callers leave them unset.
    amplitude_profile: str = "constant"
    amplitude_profile_decay_length: float | None = None
    amplitude_profile_axis: str = "x"

    # Through-width (transverse / y-direction) wrinkle surface (#300).
    # Selects a WrinkleSurface3D transverse envelope f(y) so the wrinkle
    # amplitude can vary across the specimen width instead of being uniform
    # across it (the default). ``"uniform"`` builds the bare, x-only
    # GaussianSinusoidal exactly as before (bit-identical), so the default
    # path is untouched. The non-uniform modes decay the amplitude toward
    # the edges (``"gaussian_decay"``), ripple it across the width
    # (``"sinusoidal_y"``), or confine it to a mid-width elliptical patch
    # (``"elliptical"``) — matching real, localized manufacturing wrinkles.
    # FE-only: the transverse variation only manifests in the mesh, so a
    # non-uniform mode requires the FE path (analytical_only=False) and is
    # not yet combinable with multi-wrinkle (``wrinkles``) or ``enable_czm``
    # (rejected at validation; see __post_init__). CLI/app exposure is a
    # deliberate follow-up.
    transverse_mode: str = "uniform"
    # Total specimen width span_y (mm) seen by the transverse envelope.
    # ``None`` (default) tracks the mesh width ``domain_width`` so the
    # envelope always spans the meshed y-extent. Must be > 0 when set.
    transverse_span: float | None = None
    # Transverse localization half-width width_y (mm): the Gaussian 1/e
    # length scale for ``"gaussian_decay"`` and the ellipse half-width for
    # ``"elliptical"`` (ignored by ``"uniform"`` / ``"sinusoidal_y"``).
    # ``None`` (default) resolves to ``span_y / 4`` — a localized mid-width
    # patch whose amplitude has fallen to exp(-4) ≈ 0.018 of the crest at the
    # edges (gaussian_decay) or that occupies the central half of the width
    # (elliptical). Must be > 0 when set.
    transverse_width: float | None = None

    # Loading
    loading: str = "compression"

    # Material & laminate
    material: OrthotropicMaterial | None = None
    angles: list[float] | None = None

    # Ply thickness
    ply_thickness: float = 0.183  # mm (1 ply thickness for CYCOM X850/T800)

    # Wrinkle placement. ``None`` triggers auto-derivation in
    # ``__post_init__`` from ``len(angles)`` so small laminates work out
    # of the box (issues #154/#156). For the default 24-ply layup the
    # auto-derived pair is (11, 12), preserving backwards compatibility.
    interface_1: int | None = None
    interface_2: int | None = None

    # Multi-wrinkle override. When non-empty, this overrides the named
    # single/dual-wrinkle dispatch in WrinkleAnalysis.run, allowing
    # arbitrary N-wrinkle configurations (Li et al. 2025 Dataset F).
    # Each spec contributes one WrinklePlacement; FE solve is currently
    # out of scope for this path (set analytical_only=True).
    wrinkles: list[WrinkleSpec] | None = None

    # Mesh
    nx: int = 12
    ny: int = 6
    nz_per_ply: int = 1
    domain_length: float = 0.0  # 0 → auto = 3 * wavelength
    domain_width: float = 20.0

    # Loading parameters
    applied_strain: float = -0.01

    # Thermal / cure-residual loading (issue #273).
    #
    # ``delta_T`` is the uniform temperature change FROM the stress-free
    # (cure) state, in deg C: ``T_service - T_stress_free``.  Cool-down
    # from cure is therefore NEGATIVE — a 177 C cure taken to 22 C is
    # ``delta_T = -155``.  Default 0.0 leaves every result bit-identical.
    #
    # Consumed by both paths since #273 Stage 2: the CLT path through
    # ``Laminate.thermal_resultants`` -> ``midplane_strains`` -> ply
    # stress recovery, and the FE path through the element thermal
    # initial-strain load vector ``int B^T C eps_th dV``.
    delta_T: float = 0.0

    # General load state (issue #275).  ``None`` (default) keeps the
    # legacy ``loading`` / ``applied_strain`` pair as the whole load
    # vocabulary and leaves every result bit-identical.  When set, the
    # CLT resultants it carries define the FE mechanical load directly
    # (via ``BoundaryHandler.load_state_to_bcs``), which is the only way
    # to reach biaxial and in-plane-shear states — exactly the ones whose
    # transverse/shear components drive the matrix failure modes.
    #
    # Membrane components only for now: ``Mx``/``My`` are rejected
    # because the BC mapping cannot superpose a prescribed-displacement
    # curvature with membrane tractions, and ``Qx``/``Qy`` have no
    # mapping at all.  ``delta_T`` stays on this config rather than on
    # the load state, so temperature has exactly one owner.
    load_state: LoadState | None = None

    # Solver
    solver: str = "direct"

    # Iterative-solver controls (issue #265). Only consulted on the
    # ``solver="iterative"`` path; inert (but validated) for the direct
    # solver. Defaults reproduce the previously hardcoded values in
    # ``StaticSolver._solve_iterative`` bit-for-bit, so the default
    # iterative solve is unchanged.
    #
    # ``iterative_rtol``   — CG relative-residual convergence tolerance.
    # ``iterative_maxiter``— CG iteration cap.
    # ``ilu_drop_tol``     — ILU drop tolerance (the main preconditioner
    #                        quality/memory knob; larger = sparser/cheaper
    #                        factor, weaker preconditioner).
    # ``ilu_fill_factor``  — ILU fill upper bound; ``None`` leaves SciPy's
    #                        ``spilu`` default (only forwarded when set).
    # ``preconditioner``   — ``"ilu"`` (default), ``"jacobi"`` (diagonal,
    #                        low memory for huge meshes), or ``"none"``.
    iterative_rtol: float = 1e-10
    iterative_maxiter: int = 10000
    ilu_drop_tol: float = 1e-4
    ilu_fill_factor: float | None = None
    preconditioner: str = "ilu"

    # Analytical-only mode (skip FE assembly)
    analytical_only: bool = False

    # Verbosity
    verbose: bool = False

    # Through-thickness Gaussian decay scale [mm] used by the graded
    # morphology's profile-proportional KD path.  ``None`` (default)
    # triggers the auto formula ``max(wavelength / 2, amplitude)`` so the
    # decay scale tracks the wrinkle's longitudinal extent rather than
    # its (much smaller) amplitude.  Provide an explicit positive float
    # to pin the decay scale (e.g. to reproduce the pre-PR amplitude-
    # based behaviour for regression testing).  Must be > 0 when set.
    through_thickness_decay_scale: float | None = None

    # Argon-Fleck quadratic coefficient ``c_AF`` (dimensionless) for the
    # extended Budiansky-Fleck closed form ``KD = 1/(1 + r + c_AF*r^2)``
    # where ``r = theta_eff / gamma_Y_eff``.  Default 0.0 recovers the
    # legacy linear BF response; positive values are required to match
    # the high-angle response of thick UD wrinkled coupons (Li 2024 /
    # Li 2025 high-amplitude cases) where the linear form systematically
    # over-predicts strength.  Must be >= 0.
    kink_band_quadratic_coeff: float = 0.0

    # ------------------------------------------------------------------
    # Two-parameter (theta, D/T) penetration-gate knockdown (item D.3).
    # ------------------------------------------------------------------
    # When set, the analytical knockdown is computed from the
    # penetration-gate model instead of the Budiansky-Fleck kink-band
    # path: KD = 1 - (1 - KD_angle(theta)) * min(1, (D/T / dt0)**p), which
    # reproduces both the angle and the through-thickness penetration
    # dependence the Li UD grids expose (E MAE 2.8 %, F MAE 6.0 % vs the
    # angle-only/FE models' ~20-30 %).  Material-realization specific —
    # use ``wrinklefe.core.penetration_gate.GATE_LI2024_MOULDED`` /
    # ``GATE_LI2025_VACBAG`` or calibrate your own.  UD-scoped: do NOT set
    # for multidirectional/blocked laminates.
    penetration_gate: GateParameters | None = None

    # ------------------------------------------------------------------
    # Resin-pocket material zone (Li et al. 2024/2025 UD glass datasets).
    # ------------------------------------------------------------------
    # When ``enable_resin_pocket`` is True, the FE path tags the hex
    # elements inside the cosine resin lens at the wrinkle crest and
    # assigns them an isotropic epoxy material (``resin_pocket_material``,
    # default the built-in ``EPOXY_S6C10`` card) instead of the host ply's
    # fibre-direction material, with the fibre-misalignment angle zeroed.
    # This captures the soft, fibre-free inclusion the machined cosine
    # insert leaves at the crest — a real compressive-knockdown mechanism
    # the homogenised-ply mesh otherwise misses.  No effect on the
    # analytical path or when ``analytical_only=True``.
    enable_resin_pocket: bool = False
    resin_pocket_material: OrthotropicMaterial | None = None
    # Graded transition (default): the pocket modulus blends smoothly from
    # neat resin at the lens centre to the host fibre material at the
    # boundary, and the fibre-misalignment angle is scaled by (1 - weight),
    # so the wrinkle defect is counted once.  A binary fibre/resin jump
    # (``False``) over-weakens via a spurious stress concentration that
    # double-counts the misaligned-fibre crest knockdown.
    resin_pocket_graded: bool = True
    # Crest half-height of the lens as a multiple of the wrinkle
    # half-amplitude A (mm at center, tapering to 0 at the longitudinal
    # edges).  Must be > 0 when the pocket is enabled.
    resin_pocket_height_scale: float = 1.0
    # Longitudinal half-extent of the lens as a multiple of wavelength/2
    # (the cosine insert support).  Must be > 0 when enabled.
    resin_pocket_length_scale: float = 1.0

    # ------------------------------------------------------------------
    # Surface resin pockets (tool-flat outer surface; issue #361).
    # ------------------------------------------------------------------
    # Parts cured against rigid tooling / a caul sheet have perfectly flat
    # outer surfaces: the fibre undulation is confined to the interior and
    # the wrinkle troughs fill with neat resin just under the flat surface.
    # When ``enable_surface_resin_pockets`` is True the FE path tags, in the
    # chosen outer band(s), the transition element that stretches to span
    # the gap between the flat surface and the outermost undulating ply, and
    # blends in the isotropic resin card by the excess-stretch fraction
    # (fibre-free, misalignment angle suppressed).  This is the *complement*
    # of the crest lens (surface-bonded, trough-following) and reuses
    # ``resin_pocket_material`` / ``resin_pocket_graded`` for the material
    # and blend behaviour.  FE-only effect (``modulus_retention_global`` and
    # FPF in the resin zone); no analytical-path change.  Requires a
    # tool-flat morphology whose decay reaches 0 at the chosen surface
    # (``stack``/``convex``/``concave`` or ``graded`` with
    # ``decay_floor == 0``); ``uniform`` and ``graded`` with a floor leave
    # wavy surfaces and are rejected.
    enable_surface_resin_pockets: bool = False
    # Which tool-flat surface(s) to tag: "top" (+z), "bottom" (-z), "both".
    # For the ``tool_flat`` morphology this doubles as ``tool_side`` — the
    # surface(s) pinned exactly flat by the tooling.
    surface_pocket_side: str = "top"
    # Absolute gap (mm) below which an element is treated as flat, so
    # numerically-flat regions do not produce resin "dust".  ``None`` uses
    # the ply-thickness-scaled default (0.01 * ply_thickness).
    surface_pocket_min_gap: float | None = None
    # ``tool_flat`` morphology only: number of plies over which the wrinkle
    # amplitude ramps linearly from the full-amplitude core to zero at the
    # pinned surface.  A short transition concentrates the full-amplitude
    # trough mismatch into a thin surface band (significant pockets), but the
    # crest-side transition elements compress by ``amplitude /
    # surface_transition_plies`` — so the amplitude is bounded (see
    # ``_validate``): ``A <= 0.8 * surface_transition_plies * ply_thickness /
    # nz_per_ply``.  Must be >= 1.  Default 2.
    surface_transition_plies: int = 2

    # ------------------------------------------------------------------
    # Compaction-driven fibre-volume-fraction gradient (issue #379).
    # ------------------------------------------------------------------
    # A wrinkle constrained by rigid tooling does not keep a constant ply
    # thickness: resin is squeezed out of the compacted regions and pools
    # where the geometry opens up.  When ``enable_vf_gradient`` is True the
    # FE path derives a per-element local fibre volume fraction from the
    # deformed element height (``Vf_local = vf_nominal * h0 / h``, fibre
    # content conserved) and installs per-element materials obtained by
    # scaling the preset card with the micromechanics ``Vf`` ratio
    # (:mod:`~wrinklefe.core.compaction`): stiffnesses and CTEs move,
    # Poisson ratios and ALL strengths stay at the preset values (the
    # mixing rules do not predict strengths — documented limitation).
    #
    # Opt-in and FE-only.  v1 requires ``morphology="tool_flat"``, whose
    # flat outer envelope makes the per-column thickness (and therefore the
    # resin mass) conserved by construction; ``surface_pocket_side="both"``
    # is the two-caul-plate case.  When enabled, the gradient SUPERSEDES the
    # binary surface-pocket tagging (it is the continuous generalization of
    # it — the resin-rich trough is the low-Vf end of the same field), while
    # the machined crest resin lens (``enable_resin_pocket``) composes
    # unchanged.
    enable_vf_gradient: bool = False
    # Fibre volume fraction the preset material card represents (the ratio
    # anchor).  ``None`` resolves from
    # :data:`~wrinklefe.core.compaction.CONSTITUENT_DEFAULTS` for the
    # documented library systems; supply it for any other card.
    vf_nominal: float | None = None
    # Constituent presets used for the ``Vf`` ratio.  ``None`` resolves from
    # ``CONSTITUENT_DEFAULTS``; explicit values always win.  ``vf_fiber`` is
    # a :data:`~wrinklefe.core.micromechanics.FIBER_PRESETS` key,
    # ``vf_matrix`` a ``MATRIX_PRESETS`` key or an isotropic material-library
    # card name.
    vf_fiber: str | None = None
    vf_matrix: str | None = None
    # Upper clamp on the local fibre volume fraction (default 0.75, just
    # under square packing).  Elements compacted past it saturate, which is
    # counted and warned once: the rule carries no lateral resin flow.
    vf_max: float = 0.75

    # ------------------------------------------------------------------
    # Progressive-damage FE path (load-stepping ply-discount to ultimate
    # load).  When ``enable_progressive_damage`` is True the FE solve runs
    # the :class:`~wrinklefe.solver.progressive_damage.ProgressiveDamageSolver`
    # on the wrinkled mesh and a pristine baseline, populating
    # ``progressive_strength_MPa`` and ``progressive_knockdown`` on the
    # result.  This is the only path that carries UD compression past
    # first-ply failure (the linear LaRC05 index never activates for
    # pristine UD).  No effect on ``analytical_only`` runs; not combinable
    # with ``enable_czm``.
    enable_progressive_damage: bool = False
    progressive_n_increments: int = 15
    progressive_residual_factor: float = 0.1
    # Target nominal strain magnitude for the load-stepping ramp.  ``None``
    # auto-sizes it to ~1.8x the fibre failure strain (Xc / E1) so the
    # ramp brackets the peak load.  Must be > 0 when set.
    progressive_max_strain: float | None = None

    # ------------------------------------------------------------------
    # Cohesive zone modelling (delamination prediction).
    # ------------------------------------------------------------------
    # v1: bilinear intrinsic CZM with Benzeggagh-Kenane mode-mixity.
    # Off by default; when ``enable_czm=True`` the FE solve switches
    # from :class:`StaticSolver` (linear) to
    # :class:`NewtonRaphsonSolver` and inserts zero-thickness cohesive
    # elements at the requested ply interfaces.  None of the other
    # czm_* fields have any effect when ``enable_czm=False``.
    enable_czm: bool = False
    # Which ply interfaces to insert cohesive elements at:
    #   * ``"all"`` — every interior ply interface,
    #   * ``"near_crest"`` — for a scalar (named-morphology) config, the
    #     single interface whose z-coordinate is closest to the wrinkle
    #     peak; for a multi-wrinkle config (``wrinkles`` set), the
    #     interface nearest *each* wrinkle, deduplicated — wrinkles
    #     sharing an interface index get one continuous cohesive
    #     surface, enabling crest-to-crest delamination link-up
    #     (issue #283),
    #   * ``list[int]`` — explicit list of interface indices in
    #     ``[0, n_plies-1)`` (0 = bottom-most interior interface).
    czm_interfaces: list[int] | str = "near_crest"
    czm_law: str = "bilinear"
    czm_GIc: float | None = None      # N/mm; None -> material default
    czm_GIIc: float | None = None
    czm_sigma_max: float | None = None  # MPa; None -> material default
    czm_tau_max: float | None = None
    czm_penalty: float = 1.0e6           # N/mm^3 initial interface stiffness
    czm_BK_eta: float = 1.45
    # Default bumped 20 -> 100 after the Phase 7 NASA TM DCB validation
    # showed that the original 20-increment default was too coarse for
    # accurate post-peak crack-propagation tracking; coarse increments
    # also amplified the cohesive law's d=1 corner artefact (post-peak
    # see-saw oscillations) by ~33 %.  100 fixed equal increments lands
    # the predicted peak load within experimental scatter and integrated
    # energy within 7 % for the IM7/8552 DCB benchmark.  Users who need
    # faster turnaround can lower this; users running publication-grade
    # validations should bump to 200.
    czm_n_load_increments: int = 100
    czm_newton_tol: float = 1.0e-4

    def __post_init__(self) -> None:
        if self.domain_length <= 0:
            if self.wrinkles:
                # Multi-wrinkle: size the domain from the union of the
                # wrinkle extents (per-spec center offset from its phase
                # plus the 3*width Gaussian support), not from the
                # scalar wavelength field. The 3*wavelength floor keeps
                # a single centred spec consistent with the scalar path.
                half_span = max(
                    abs(s.phase_offset) * s.wavelength / (2.0 * math.pi)
                    + 3.0 * s.width
                    for s in self.wrinkles
                )
                self.domain_length = max(
                    2.0 * half_span, 3.0 * self.wavelength
                )
            else:
                self.domain_length = 3.0 * self.wavelength
        if self.material is None:
            self.material = MaterialLibrary().get("IM7_8552")
        if self.angles is None:
            # Quasi-isotropic [0/45/-45/90]_3s → 24 plies
            base: list[float] = [0, 45, -45, 90]
            self.angles = (base * 3) + list(reversed(base * 3))

        # Auto-derive interior interface indices when the user did not
        # specify them. The two interfaces sit symmetrically about the
        # mid-thickness, so for the default 24-ply layup this resolves
        # to (11, 12) — i.e. backwards-compatible with the previous
        # hard-coded dataclass defaults. For small layups (< 13 plies)
        # it picks valid in-range indices instead of crashing in
        # ``_validate`` (issues #154 / #156).
        n_plies = len(self.angles) if self.angles is not None else 0
        if n_plies > 0:
            mid = n_plies // 2
            if self.interface_1 is None:
                self.interface_1 = max(0, mid - 1)
            if self.interface_2 is None:
                self.interface_2 = min(n_plies - 1, mid)

        # Surface resin pockets ARE the physics of the ``tool_flat``
        # morphology (the flat pinned surface would otherwise silently
        # stretch the fibre elements over the troughs — the old #371 bug).
        # Auto-enable them on the FE path.  Skipped for analytical_only,
        # where the pocket has no effect and the run equals ``uniform``.
        morph_name = (
            self.morphology.lower().strip()
            if isinstance(self.morphology, str)
            else self.morphology
        )
        if (
            morph_name == "tool_flat"
            and not self.analytical_only
            and not self.enable_surface_resin_pockets
        ):
            logger.warning(
                "Surface resin pockets auto-enabled for the 'tool_flat' "
                "morphology (they are its defining physics: the flat pinned "
                "surface fills the wrinkle troughs with neat resin). Set "
                "enable_surface_resin_pockets=True to silence this message."
            )
            self.enable_surface_resin_pockets = True

        self._validate()

    def _validate(self) -> None:
        """Fail fast on physically invalid configuration.

        Called from :meth:`__post_init__` after defaults (``domain_length``,
        ``material``, ``angles``) have been resolved.  Each check raises a
        :class:`ValueError` naming the offending field and value so that
        misconfiguration surfaces at construction time rather than as an
        obscure traceback deep in the solver/mesh path.
        """
        # ``angles`` is filled in __post_init__ before _validate runs.
        assert self.angles is not None
        # --- Ply angles ------------------------------------------------
        # Enforce the canonical fibre-angle range (|angle| <= 90) using
        # the shared parser rule so a mis-typed layup (e.g. 900°) fails
        # loudly at construction instead of flowing into CLT trig and
        # being silently mis-classified by the tension-mechanism heuristic.
        for i, angle in enumerate(self.angles):
            validate_ply_angle(
                float(angle), context=f"AnalysisConfig.angles[{i}] = "
            )
        # --- Wrinkle geometry -----------------------------------------
        # amplitude == 0 is a legitimate "no wrinkle" (flat) case: the
        # mid-surface profile z(x) = A * envelope reduces to 0, so the
        # closed-form misalignment angle arctan(2*pi*A/lambda) is 0.
        # Only a negative amplitude is physically meaningless.
        if self.amplitude < 0:
            raise ValueError(
                f"AnalysisConfig.amplitude must be >= 0 "
                f"(0 = flat / no wrinkle), got {self.amplitude}"
            )
        if not (self.wavelength > 0):
            # Strictly positive: lambda divides into the closed-form
            # slope (2*pi*A/lambda) and the auto-derived domain_length.
            raise ValueError(
                f"AnalysisConfig.wavelength must be > 0, "
                f"got {self.wavelength}"
            )
        if not (self.width > 0):
            raise ValueError(
                f"AnalysisConfig.width must be > 0, got {self.width}"
            )
        if not (self.domain_length > 0):
            raise ValueError(
                f"AnalysisConfig.domain_length must be > 0, "
                f"got {self.domain_length}"
            )
        if not (self.domain_width > 0):
            raise ValueError(
                f"AnalysisConfig.domain_width must be > 0, "
                f"got {self.domain_width}"
            )
        if not (self.ply_thickness > 0):
            raise ValueError(
                f"AnalysisConfig.ply_thickness must be > 0, "
                f"got {self.ply_thickness}"
            )

        # --- Morphology / phase ---------------------------------------
        # run() resolves morphology by name (possibly lower/strip'd) via
        # WrinkleConfiguration.from_morphology_name; an explicit numeric
        # ``phase`` overrides the named-morphology phase for dual-wrinkle
        # modes but the name is still consumed (single-wrinkle modes and
        # the from_morphology_name fallback both require a known name).
        valid_morphologies = sorted(
            set(MORPHOLOGY_PHASES) | set(SINGLE_WRINKLE_MODES)
        )
        if (
            not isinstance(self.morphology, str)
            or self.morphology.lower().strip() not in valid_morphologies
        ):
            raise ValueError(
                f"AnalysisConfig.morphology must be one of "
                f"{valid_morphologies}, got {self.morphology!r}"
            )
        if self.phase is not None and not math.isfinite(float(self.phase)):
            raise ValueError(
                f"AnalysisConfig.phase must be finite when set, "
                f"got {self.phase}"
            )
        if not (0.0 <= self.decay_floor <= 1.0):
            raise ValueError(
                f"AnalysisConfig.decay_floor must be in [0, 1], "
                f"got {self.decay_floor}"
            )

        # Single-wrinkle through-thickness position as a fraction of T.
        # Must be a finite float in [0, 1]; NaN and infinities are rejected
        # so they cannot silently shift the decay centre off the laminate.
        try:
            wz_value = float(self.wrinkle_z_position)
        except (TypeError, ValueError):
            raise ValueError(
                f"AnalysisConfig.wrinkle_z_position must be a finite float "
                f"in [0, 1], got {self.wrinkle_z_position!r}"
            )
        if not math.isfinite(wz_value) or not (0.0 <= wz_value <= 1.0):
            raise ValueError(
                f"AnalysisConfig.wrinkle_z_position must be a finite float "
                f"in [0, 1], got {self.wrinkle_z_position!r}"
            )

        # --- Amplitude profile (spatially varying A modulation) -------
        valid_amplitude_profiles = ("constant", "gaussian", "linear")
        if (
            not isinstance(self.amplitude_profile, str)
            or self.amplitude_profile.lower().strip() not in valid_amplitude_profiles
        ):
            raise ValueError(
                f"AnalysisConfig.amplitude_profile must be one of "
                f"{list(valid_amplitude_profiles)}, got {self.amplitude_profile!r}"
            )
        valid_amplitude_profile_axes = ("x", "y")
        if (
            not isinstance(self.amplitude_profile_axis, str)
            or self.amplitude_profile_axis.lower().strip()
            not in valid_amplitude_profile_axes
        ):
            raise ValueError(
                f"AnalysisConfig.amplitude_profile_axis must be one of "
                f"{list(valid_amplitude_profile_axes)}, "
                f"got {self.amplitude_profile_axis!r}"
            )
        if self.amplitude_profile_decay_length is not None and not (
            self.amplitude_profile_decay_length > 0.0
            and math.isfinite(self.amplitude_profile_decay_length)
        ):
            raise ValueError(
                f"AnalysisConfig.amplitude_profile_decay_length must be "
                f"a finite positive float when set, "
                f"got {self.amplitude_profile_decay_length}"
            )

        # --- Through-width (transverse) wrinkle surface (#300) --------
        # Reuse WrinkleSurface3D's own mode set so the two never drift.
        if (
            not isinstance(self.transverse_mode, str)
            or self.transverse_mode not in WrinkleSurface3D._VALID_MODES
        ):
            raise ValueError(
                f"AnalysisConfig.transverse_mode must be one of "
                f"{sorted(WrinkleSurface3D._VALID_MODES)}, "
                f"got {self.transverse_mode!r}"
            )
        # span_y / width_y overrides must be positive and finite when set;
        # ``None`` resolves to domain_width and span_y/4 at build time and is
        # therefore always valid.
        if self.transverse_span is not None and not (
            self.transverse_span > 0.0 and math.isfinite(self.transverse_span)
        ):
            raise ValueError(
                f"AnalysisConfig.transverse_span must be a finite positive "
                f"float when set, got {self.transverse_span}"
            )
        if self.transverse_width is not None and not (
            self.transverse_width > 0.0 and math.isfinite(self.transverse_width)
        ):
            raise ValueError(
                f"AnalysisConfig.transverse_width must be a finite positive "
                f"float when set, got {self.transverse_width}"
            )
        if self.transverse_mode != "uniform":
            # The transverse envelope only manifests in the FE mesh, so it is
            # meaningless on the x-only analytical path — fail fast instead of
            # silently ignoring the requested surface.
            if self.analytical_only:
                raise ValueError(
                    "AnalysisConfig.transverse_mode="
                    f"{self.transverse_mode!r} requires the FE path but "
                    "analytical_only=True. The transverse surface only "
                    "manifests in the mesh; set analytical_only=False or "
                    "transverse_mode='uniform'."
                )
            # Multi-wrinkle FE and CZM composition with a 3-D surface are out
            # of scope for this first cut — reject the combination up front
            # (issue #300) rather than surprising the user mid-solve.
            if self.wrinkles is not None:
                raise NotImplementedError(
                    "AnalysisConfig.transverse_mode="
                    f"{self.transverse_mode!r} is not yet supported together "
                    "with a multi-wrinkle 'wrinkles' list. Use a single-"
                    "wrinkle config (wrinkles=None) or transverse_mode="
                    "'uniform'."
                )
            if self.enable_czm:
                raise NotImplementedError(
                    "AnalysisConfig.transverse_mode="
                    f"{self.transverse_mode!r} is not yet supported together "
                    "with enable_czm=True. Disable CZM or use "
                    "transverse_mode='uniform'."
                )

        # --- Loading --------------------------------------------------
        valid_loadings = ("compression", "tension")
        if (
            not isinstance(self.loading, str)
            or self.loading.lower().strip() not in valid_loadings
        ):
            raise ValueError(
                f"AnalysisConfig.loading must be one of "
                f"{list(valid_loadings)}, got {self.loading!r}"
            )
        if not math.isfinite(self.applied_strain):
            raise ValueError(
                f"AnalysisConfig.applied_strain must be finite, "
                f"got {self.applied_strain}"
            )

        # --- Thermal / cure-residual load (issue #273) ----------------
        if not math.isfinite(self.delta_T):
            raise ValueError(
                f"AnalysisConfig.delta_T must be finite, got {self.delta_T}"
            )
        if abs(self.delta_T) > _DELTA_T_MAX:
            raise ValueError(
                f"AnalysisConfig.delta_T = {self.delta_T} deg C is outside "
                f"the supported range |delta_T| <= {_DELTA_T_MAX:.0f}. "
                f"delta_T is the temperature change FROM the stress-free "
                f"(cure) state, not an absolute temperature: a 177 C cure "
                f"taken to 22 C service is delta_T = -155, not -273 or 22."
            )
        # --- General load state (issue #275) --------------------------
        if self.load_state is not None:
            ls = self.load_state
            if not isinstance(ls, LoadState):
                raise ValueError(
                    "AnalysisConfig.load_state must be a LoadState or None, "
                    f"got {type(ls).__name__}."
                )
            for name in ("Nx", "Ny", "Nxy", "Mx", "My", "Mxy",
                         "Qx", "Qy", "delta_T", "delta_C"):
                value = float(getattr(ls, name))
                if not math.isfinite(value):
                    raise ValueError(
                        f"AnalysisConfig.load_state.{name} must be finite, "
                        f"got {value}."
                    )
            # Membrane-only: the BC mapping refuses to superpose a
            # prescribed-displacement curvature with membrane tractions
            # (they act on the same faces), and Qx/Qy have no mapping.
            for name in ("Mx", "My", "Mxy"):
                if float(getattr(ls, name)) != 0.0:
                    raise ValueError(
                        f"AnalysisConfig.load_state.{name}="
                        f"{getattr(ls, name)} is not supported yet: "
                        "curvature resultants are applied as prescribed "
                        "displacements on the same faces a membrane state "
                        "loads with traction, so the two cannot be "
                        "superposed (issue #275). Use a membrane-only "
                        "state (Nx/Ny/Nxy)."
                    )
            for name in ("Qx", "Qy"):
                if float(getattr(ls, name)) != 0.0:
                    raise ValueError(
                        f"AnalysisConfig.load_state.{name}="
                        f"{getattr(ls, name)} is not supported: transverse "
                        "shear resultants have no 3-D boundary-condition "
                        "mapping. Use a membrane-only state (Nx/Ny/Nxy)."
                    )
            if not any(abs(float(getattr(ls, n))) > 0.0
                       for n in ("Nx", "Ny", "Nxy")):
                raise ValueError(
                    "AnalysisConfig.load_state carries no membrane "
                    "resultant (Nx, Ny and Nxy are all zero), which maps "
                    "to an empty boundary-condition set and a singular "
                    "system. Set at least one, or leave load_state=None "
                    "to use the applied_strain path."
                )
            # One owner per quantity: temperature lives on the config.
            for name in ("delta_T", "delta_C"):
                if float(getattr(ls, name)) != 0.0:
                    raise ValueError(
                        f"AnalysisConfig.load_state.{name}="
                        f"{getattr(ls, name)} must be zero. Environmental "
                        "loading is owned by the config, not the load "
                        "state: set AnalysisConfig.delta_T instead (a "
                        "second place to set the same quantity is how "
                        "sign errors get in). Moisture is not exposed at "
                        "all - nothing in the solve consumes delta_C."
                    )
            # FE-only: the closed-form analytical knockdown is uniaxial.
            if self.analytical_only:
                raise ValueError(
                    "AnalysisConfig.load_state requires the FE path but "
                    "analytical_only=True. The closed-form knockdown is "
                    "defined for a uniaxial state only; a general load "
                    "state is applied through 3-D boundary conditions. "
                    "Set analytical_only=False or load_state=None."
                )
            if self.enable_czm:
                raise ValueError(
                    "AnalysisConfig.load_state is not yet combinable with "
                    "enable_czm: the CZM path builds its own uniaxial "
                    "compression boundary conditions, so the load state "
                    "would be silently ignored (issue #275)."
                )
            if self.enable_progressive_damage:
                raise ValueError(
                    "AnalysisConfig.load_state is not yet combinable with "
                    "enable_progressive_damage: that solver ramps a "
                    "prescribed displacement to ultimate load, and a "
                    "force-controlled load state has no applied_strain to "
                    "ramp (issue #275)."
                )

        # --- Wrinkle placement (interface indices) --------------------
        n_plies = len(self.angles)
        for name in ("interface_1", "interface_2"):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool):
                raise ValueError(
                    f"AnalysisConfig.{name} must be an int, got {value!r}"
                )
            if not (0 <= value < n_plies):
                raise ValueError(
                    f"AnalysisConfig.{name} must be in [0, {n_plies}) "
                    f"(0 <= interface < number of plies), got {value}"
                )

        # --- Mesh resolution (structural integers) --------------------
        for name in ("nx", "ny", "nz_per_ply"):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool):
                raise ValueError(
                    f"AnalysisConfig.{name} must be an int, got {value!r}"
                )
            if value < 1:
                raise ValueError(
                    f"AnalysisConfig.{name} must be >= 1, got {value}"
                )

        # --- Solver ---------------------------------------------------
        valid_solvers = ("direct", "iterative")
        if (
            not isinstance(self.solver, str)
            or self.solver.lower().strip() not in valid_solvers
        ):
            raise ValueError(
                f"AnalysisConfig.solver must be one of "
                f"{list(valid_solvers)}, got {self.solver!r}"
            )

        # --- Iterative-solver controls (issue #265) -------------------
        # Positive finite CG tolerance.
        if not (
            isinstance(self.iterative_rtol, (int, float))
            and not isinstance(self.iterative_rtol, bool)
            and math.isfinite(self.iterative_rtol)
            and self.iterative_rtol > 0.0
        ):
            raise ValueError(
                f"AnalysisConfig.iterative_rtol must be a finite positive "
                f"float, got {self.iterative_rtol!r}"
            )
        # Positive integer iteration cap.
        if (
            not isinstance(self.iterative_maxiter, int)
            or isinstance(self.iterative_maxiter, bool)
            or self.iterative_maxiter < 1
        ):
            raise ValueError(
                f"AnalysisConfig.iterative_maxiter must be an int >= 1, "
                f"got {self.iterative_maxiter!r}"
            )
        # Non-negative finite ILU drop tolerance.
        if not (
            isinstance(self.ilu_drop_tol, (int, float))
            and not isinstance(self.ilu_drop_tol, bool)
            and math.isfinite(self.ilu_drop_tol)
            and self.ilu_drop_tol >= 0.0
        ):
            raise ValueError(
                f"AnalysisConfig.ilu_drop_tol must be a finite float >= 0, "
                f"got {self.ilu_drop_tol!r}"
            )
        # Optional ILU fill factor: None or a finite float >= 1.
        if self.ilu_fill_factor is not None and not (
            isinstance(self.ilu_fill_factor, (int, float))
            and not isinstance(self.ilu_fill_factor, bool)
            and math.isfinite(self.ilu_fill_factor)
            and self.ilu_fill_factor >= 1.0
        ):
            raise ValueError(
                f"AnalysisConfig.ilu_fill_factor must be a finite float "
                f">= 1 or None, got {self.ilu_fill_factor!r}"
            )
        valid_preconditioners = ("ilu", "jacobi", "none")
        if (
            not isinstance(self.preconditioner, str)
            or self.preconditioner.lower().strip() not in valid_preconditioners
        ):
            raise ValueError(
                f"AnalysisConfig.preconditioner must be one of "
                f"{list(valid_preconditioners)}, got {self.preconditioner!r}"
            )

        # --- Through-thickness decay scale ----------------------------
        # Optional positive float (mm) overriding the auto formula
        # ``max(wavelength / 2, amplitude)`` used by the graded
        # morphology profile-proportional path.  Must be strictly
        # positive and finite when provided.
        if self.through_thickness_decay_scale is not None and not (
            self.through_thickness_decay_scale > 0.0
            and math.isfinite(self.through_thickness_decay_scale)
        ):
            raise ValueError(
                f"AnalysisConfig.through_thickness_decay_scale must be a "
                f"finite positive float when set, "
                f"got {self.through_thickness_decay_scale}"
            )

        # --- Argon-Fleck quadratic coefficient ------------------------
        # Dimensionless non-negative float (0.0 = legacy linear BF).
        if (
            not isinstance(self.kink_band_quadratic_coeff, (int, float))
            or isinstance(self.kink_band_quadratic_coeff, bool)
            or not math.isfinite(self.kink_band_quadratic_coeff)
            or self.kink_band_quadratic_coeff < 0.0
        ):
            raise ValueError(
                f"AnalysisConfig.kink_band_quadratic_coeff must be a "
                f"finite float >= 0, "
                f"got {self.kink_band_quadratic_coeff!r}"
            )

        # --- Penetration gate -----------------------------------------
        if self.penetration_gate is not None and not isinstance(
            self.penetration_gate, GateParameters
        ):
            raise ValueError(
                "AnalysisConfig.penetration_gate must be a GateParameters "
                f"or None, got {type(self.penetration_gate).__name__}"
            )

        # --- Resin-pocket zone ----------------------------------------
        # Fields are no-ops when ``enable_resin_pocket`` is False; still
        # type-checked so misconfigurations surface at construction time.
        if not isinstance(self.enable_resin_pocket, bool):
            raise ValueError(
                f"AnalysisConfig.enable_resin_pocket must be a bool, "
                f"got {self.enable_resin_pocket!r}"
            )
        for name in ("resin_pocket_height_scale", "resin_pocket_length_scale"):
            val = getattr(self, name)
            if not (
                isinstance(val, (int, float))
                and not isinstance(val, bool)
                and math.isfinite(val)
                and val > 0.0
            ):
                raise ValueError(
                    f"AnalysisConfig.{name} must be a finite positive "
                    f"float, got {val!r}"
                )

        # --- Surface resin pockets (tool-flat surface; issue #361) ----
        if not isinstance(self.enable_surface_resin_pockets, bool):
            raise ValueError(
                f"AnalysisConfig.enable_surface_resin_pockets must be a "
                f"bool, got {self.enable_surface_resin_pockets!r}"
            )
        if self.surface_pocket_side not in ("top", "bottom", "both"):
            raise ValueError(
                f"AnalysisConfig.surface_pocket_side must be 'top', "
                f"'bottom' or 'both', got {self.surface_pocket_side!r}"
            )
        if self.surface_pocket_min_gap is not None and not (
            isinstance(self.surface_pocket_min_gap, (int, float))
            and not isinstance(self.surface_pocket_min_gap, bool)
            and math.isfinite(self.surface_pocket_min_gap)
            and self.surface_pocket_min_gap >= 0.0
        ):
            raise ValueError(
                f"AnalysisConfig.surface_pocket_min_gap must be a "
                f"non-negative finite float or None, got "
                f"{self.surface_pocket_min_gap!r}"
            )
        if self.enable_surface_resin_pockets:
            # Surface pockets are an FE-only material effect; they have no
            # analytical-path knockdown, so an analytical-only run would
            # silently drop them (mirrors the FE-only nature of the crest
            # pocket, but rejected explicitly so it cannot no-op unnoticed).
            if self.analytical_only:
                raise ValueError(
                    "AnalysisConfig: enable_surface_resin_pockets requires "
                    "the FE path (an isotropic resin zone in the mesh); it "
                    "has no effect under analytical_only=True. Set "
                    "analytical_only=False or disable surface resin pockets."
                )
            # The chosen surface must be tool-flat, i.e. the through-thickness
            # decay must reach 0 there. ``uniform`` never decays (wavy
            # surfaces); ``graded`` with a non-zero floor leaves residual
            # surface waviness.  Name the conflict and the fix.
            morph = (
                self.morphology.lower().strip()
                if isinstance(self.morphology, str)
                else self.morphology
            )
            if morph == "uniform":
                raise ValueError(
                    "AnalysisConfig: uniform morphology has wavy outer "
                    "surfaces (displacement never decays through the "
                    "thickness); surface resin pockets require a tool-flat "
                    "morphology. Use 'stack'/'convex'/'concave', or 'graded' "
                    "with decay_floor=0."
                )
            if morph == "graded" and self.decay_floor > 0.0:
                raise ValueError(
                    "AnalysisConfig: graded morphology with decay_floor="
                    f"{self.decay_floor} leaves the outer surfaces wavy "
                    "(the floor is a residual surface amplitude); surface "
                    "resin pockets require a tool-flat surface. Set "
                    "decay_floor=0 (or use 'stack'/'convex'/'concave')."
                )

        # --- tool_flat morphology (issue #371) ------------------------
        # ``surface_transition_plies`` is validated for every config (it is a
        # dataclass field), then the tool_flat-specific geometry bound and
        # combination rules are enforced only when the morphology is selected.
        if (
            not isinstance(self.surface_transition_plies, int)
            or isinstance(self.surface_transition_plies, bool)
            or self.surface_transition_plies < 1
        ):
            raise ValueError(
                f"AnalysisConfig.surface_transition_plies must be an int "
                f">= 1, got {self.surface_transition_plies!r}"
            )
        morph = (
            self.morphology.lower().strip()
            if isinstance(self.morphology, str)
            else self.morphology
        )
        if morph == "tool_flat":
            # Element-inversion bound (verified in the #371 spike): on the
            # crest side each of the ``surface_transition_plies`` transition
            # elements compresses by ``amplitude / surface_transition_plies``
            # over an element height ``ply_thickness / nz_per_ply``, so it
            # inverts once ``amplitude >= surface_transition_plies *
            # ply_thickness / nz_per_ply``.  Reject above a safe 0.8 fraction,
            # naming BOTH remedies.
            a_safe = (
                0.8
                * self.surface_transition_plies
                * self.ply_thickness
                / self.nz_per_ply
            )
            if self.amplitude > a_safe:
                raise ValueError(
                    "AnalysisConfig: tool_flat amplitude "
                    f"{self.amplitude:.4g} mm exceeds the safe bound "
                    f"{a_safe:.4g} mm — the crest-side transition elements "
                    "would invert (negative Jacobian). Either increase "
                    "surface_transition_plies to >= "
                    f"{math.ceil(self.amplitude * self.nz_per_ply / (0.8 * self.ply_thickness))} "
                    f"or reduce amplitude below {a_safe:.4g} mm "
                    "(the bound is 0.8 * surface_transition_plies * "
                    "ply_thickness / nz_per_ply)."
                )
            # tool_flat is a single-wrinkle FE morphology whose pockets and
            # decay are not verified in combination with the multi-wrinkle,
            # CZM, or through-width paths — reject rather than silently drop
            # the tool-flat decay (which would revert to the old linear taper).
            if self.wrinkles is not None:
                raise NotImplementedError(
                    "AnalysisConfig: tool_flat morphology is single-wrinkle; "
                    "it is not yet combinable with the multi-wrinkle "
                    "'wrinkles' override (which would ignore the tool-flat "
                    "decay). Use a single wrinkle or a different morphology."
                )
            if self.enable_czm:
                raise NotImplementedError(
                    "AnalysisConfig: tool_flat morphology is not yet "
                    "combinable with enable_czm (the surface-pocket resin "
                    "zone and cohesive elements are unverified together)."
                )
            if self.transverse_mode != "uniform":
                raise NotImplementedError(
                    "AnalysisConfig: tool_flat morphology is not yet "
                    "combinable with a non-uniform transverse_mode "
                    f"({self.transverse_mode!r}); the tool-flat decay is "
                    "verified for the x-only wrinkle only."
                )

        # --- Compaction Vf gradient (issue #379) ----------------------
        if not isinstance(self.enable_vf_gradient, bool):
            raise ValueError(
                f"AnalysisConfig.enable_vf_gradient must be a bool, "
                f"got {self.enable_vf_gradient!r}"
            )
        for name in ("vf_nominal",):
            val = getattr(self, name)
            if val is not None and not (
                isinstance(val, (int, float))
                and not isinstance(val, bool)
                and math.isfinite(val)
                and 0.0 < val < 1.0
            ):
                raise ValueError(
                    f"AnalysisConfig.{name} must be a float in (0, 1) or "
                    f"None, got {val!r}"
                )
        if not (
            isinstance(self.vf_max, (int, float))
            and not isinstance(self.vf_max, bool)
            and math.isfinite(self.vf_max)
            and 0.0 < self.vf_max <= 0.9
        ):
            raise ValueError(
                f"AnalysisConfig.vf_max must be a float in (0, 0.9], got "
                f"{self.vf_max!r}"
            )
        for name in ("vf_fiber", "vf_matrix"):
            val = getattr(self, name)
            if val is not None and not isinstance(val, str):
                raise ValueError(
                    f"AnalysisConfig.{name} must be a preset name (str) or "
                    f"None, got {val!r}"
                )
        if self.enable_vf_gradient:
            # FE-only: the gradient is a per-element material field, so an
            # analytical-only run would silently drop it (same rule as the
            # resin-pocket / surface-pocket features).
            if self.analytical_only:
                raise ValueError(
                    "AnalysisConfig: enable_vf_gradient requires the FE "
                    "path (it installs per-element materials in the mesh); "
                    "it has no effect under analytical_only=True. Set "
                    "analytical_only=False or disable the Vf gradient."
                )
            # v1 restriction: only the tool_flat morphology has a flat outer
            # envelope, which is what makes the per-column thickness — and
            # therefore the resin mass — conserved under the kinematic rule.
            # Every other morphology leaves an undulating outer surface, so
            # the compaction would create or destroy material.
            if morph != "tool_flat":
                raise ValueError(
                    "AnalysisConfig: enable_vf_gradient is restricted to "
                    f"morphology='tool_flat' in v1, got {self.morphology!r}. "
                    "Only a tool-flat outer envelope conserves the "
                    "per-column thickness (and hence the resin mass) under "
                    "the compaction rule Vf = vf_nominal * h0/h; an "
                    "undulating outer surface would create or destroy "
                    "material. Use morphology='tool_flat' (with "
                    "surface_pocket_side='both' for the two-caul-plate "
                    "case) or disable the Vf gradient."
                )
            if self.enable_czm:
                raise NotImplementedError(
                    "AnalysisConfig: enable_vf_gradient is not yet "
                    "combinable with enable_czm (the per-element Vf "
                    "materials and cohesive elements are unverified "
                    "together)."
                )
            # Constituents must resolve now, not deep in the mesh path.
            from wrinklefe.core.compaction import VfGradientSpec

            material_name = (
                self.material.name if self.material is not None else ""
            )
            try:
                VfGradientSpec.for_material(
                    material_name,
                    fiber=self.vf_fiber,
                    matrix=self.vf_matrix,
                    vf_nominal=self.vf_nominal,
                    vf_max=self.vf_max,
                )
            except ValueError as exc:
                raise ValueError(
                    f"AnalysisConfig: enable_vf_gradient cannot resolve the "
                    f"micromechanics anchor: {exc}"
                ) from exc

        # --- Progressive-damage path ----------------------------------
        if not isinstance(self.enable_progressive_damage, bool):
            raise ValueError(
                f"AnalysisConfig.enable_progressive_damage must be a bool, "
                f"got {self.enable_progressive_damage!r}"
            )
        if self.enable_progressive_damage and self.enable_czm:
            raise ValueError(
                "AnalysisConfig: enable_progressive_damage and enable_czm "
                "cannot both be True (separate nonlinear FE paths)."
            )
        if (
            not isinstance(self.progressive_n_increments, int)
            or isinstance(self.progressive_n_increments, bool)
            or self.progressive_n_increments < 1
        ):
            raise ValueError(
                f"AnalysisConfig.progressive_n_increments must be an int "
                f">= 1, got {self.progressive_n_increments!r}"
            )
        if not (
            isinstance(self.progressive_residual_factor, (int, float))
            and not isinstance(self.progressive_residual_factor, bool)
            and 0.0 < self.progressive_residual_factor < 1.0
        ):
            raise ValueError(
                f"AnalysisConfig.progressive_residual_factor must be in "
                f"(0, 1), got {self.progressive_residual_factor!r}"
            )
        if self.progressive_max_strain is not None and not (
            isinstance(self.progressive_max_strain, (int, float))
            and not isinstance(self.progressive_max_strain, bool)
            and math.isfinite(self.progressive_max_strain)
            and self.progressive_max_strain > 0.0
        ):
            raise ValueError(
                f"AnalysisConfig.progressive_max_strain must be a finite "
                f"positive float when set, got {self.progressive_max_strain!r}"
            )

        # --- Cohesive zone modelling ----------------------------------
        # All CZM fields are no-ops when ``enable_czm`` is False; we
        # still type-check them so misconfigurations surface at
        # construction time rather than mid-solve.
        if not isinstance(self.enable_czm, bool):
            raise ValueError(
                f"AnalysisConfig.enable_czm must be a bool, "
                f"got {self.enable_czm!r}"
            )
        if self.czm_law != "bilinear":
            raise ValueError(
                f"AnalysisConfig.czm_law: only 'bilinear' is supported "
                f"in v1, got {self.czm_law!r}"
            )
        if isinstance(self.czm_interfaces, str):
            if self.czm_interfaces not in ("all", "near_crest"):
                raise ValueError(
                    "AnalysisConfig.czm_interfaces must be 'all', "
                    "'near_crest', or a list of int, got "
                    f"{self.czm_interfaces!r}"
                )
        elif isinstance(self.czm_interfaces, list):
            n_plies = len(self.angles)
            for i, idx in enumerate(self.czm_interfaces):
                if not isinstance(idx, int) or isinstance(idx, bool):
                    raise ValueError(
                        f"AnalysisConfig.czm_interfaces[{i}] must be an "
                        f"int, got {idx!r}"
                    )
                # Valid interior ply interfaces are 0 .. n_plies-2.
                if not (0 <= idx <= n_plies - 2):
                    raise ValueError(
                        f"AnalysisConfig.czm_interfaces[{i}] must be in "
                        f"[0, {n_plies - 2}], got {idx}"
                    )
        else:
            raise ValueError(
                "AnalysisConfig.czm_interfaces must be a string ('all' "
                "or 'near_crest') or a list[int], got "
                f"{type(self.czm_interfaces).__name__}"
            )
        for name in ("czm_GIc", "czm_GIIc", "czm_sigma_max", "czm_tau_max"):
            val = getattr(self, name)
            if val is None:
                continue
            if not (
                isinstance(val, (int, float))
                and not isinstance(val, bool)
                and math.isfinite(val)
                and val > 0.0
            ):
                raise ValueError(
                    f"AnalysisConfig.{name} must be a finite positive "
                    f"float when set, got {val!r}"
                )
        if not (
            isinstance(self.czm_penalty, (int, float))
            and not isinstance(self.czm_penalty, bool)
            and math.isfinite(self.czm_penalty)
            and self.czm_penalty > 0.0
        ):
            raise ValueError(
                f"AnalysisConfig.czm_penalty must be a finite positive "
                f"float, got {self.czm_penalty!r}"
            )
        if not (
            isinstance(self.czm_BK_eta, (int, float))
            and not isinstance(self.czm_BK_eta, bool)
            and math.isfinite(self.czm_BK_eta)
            and self.czm_BK_eta > 0.0
        ):
            raise ValueError(
                f"AnalysisConfig.czm_BK_eta must be a finite positive "
                f"float, got {self.czm_BK_eta!r}"
            )
        if (
            not isinstance(self.czm_n_load_increments, int)
            or isinstance(self.czm_n_load_increments, bool)
            or self.czm_n_load_increments < 1
        ):
            raise ValueError(
                f"AnalysisConfig.czm_n_load_increments must be an int "
                f">= 1, got {self.czm_n_load_increments!r}"
            )
        if not (
            isinstance(self.czm_newton_tol, (int, float))
            and not isinstance(self.czm_newton_tol, bool)
            and math.isfinite(self.czm_newton_tol)
            and self.czm_newton_tol > 0.0
        ):
            raise ValueError(
                f"AnalysisConfig.czm_newton_tol must be a finite positive "
                f"float, got {self.czm_newton_tol!r}"
            )

        # --- Multi-wrinkle override -----------------------------------
        # When ``wrinkles`` is provided, it overrides the named
        # single/dual-wrinkle dispatch.  An empty list is rejected
        # because the intent is ambiguous (use None for the default).
        if self.wrinkles is not None:
            if not isinstance(self.wrinkles, list):
                raise ValueError(
                    "AnalysisConfig.wrinkles must be a list of WrinkleSpec "
                    f"or None, got {type(self.wrinkles).__name__}"
                )
            if len(self.wrinkles) == 0:
                raise ValueError(
                    "AnalysisConfig.wrinkles must contain at least one "
                    "WrinkleSpec; pass None to use the default dispatch."
                )
            for i, spec in enumerate(self.wrinkles):
                if not isinstance(spec, WrinkleSpec):
                    raise ValueError(
                        f"AnalysisConfig.wrinkles[{i}] must be a "
                        f"WrinkleSpec, got {type(spec).__name__}"
                    )
                if not (spec.amplitude > 0 and math.isfinite(spec.amplitude)):
                    raise ValueError(
                        f"AnalysisConfig.wrinkles[{i}].amplitude must be "
                        f"a positive finite float, got {spec.amplitude}"
                    )
                if not (spec.wavelength > 0 and math.isfinite(spec.wavelength)):
                    raise ValueError(
                        f"AnalysisConfig.wrinkles[{i}].wavelength must be "
                        f"a positive finite float, got {spec.wavelength}"
                    )
                if not (spec.width > 0 and math.isfinite(spec.width)):
                    raise ValueError(
                        f"AnalysisConfig.wrinkles[{i}].width must be "
                        f"a positive finite float, got {spec.width}"
                    )
                if (
                    not isinstance(spec.ply_interface, int)
                    or isinstance(spec.ply_interface, bool)
                ):
                    raise ValueError(
                        f"AnalysisConfig.wrinkles[{i}].ply_interface must "
                        f"be an int, got {spec.ply_interface!r}"
                    )
                # WrinklePlacement valid range is [0, n_plies - 2]
                if not (0 <= spec.ply_interface <= n_plies - 2):
                    raise ValueError(
                        f"AnalysisConfig.wrinkles[{i}].ply_interface must "
                        f"be in [0, {n_plies - 2}] (n_plies={n_plies}), "
                        f"got {spec.ply_interface}"
                    )
                if not math.isfinite(spec.phase_offset):
                    raise ValueError(
                        f"AnalysisConfig.wrinkles[{i}].phase_offset must "
                        f"be finite, got {spec.phase_offset}"
                    )

    # ------------------------------------------------------------------
    # Serialisation (round-trippable save / load) — issue #259
    # ------------------------------------------------------------------
    #: Object-valued fields serialised through dedicated encoders rather
    #: than emitted as-is.  Every remaining field is a JSON-native scalar
    #: or list, so this list plus the plain fields covers the dataclass in
    #: full — the completeness guard in the test-suite pins that invariant.
    _NON_PRIMITIVE_FIELDS = frozenset(
        {"material", "resin_pocket_material", "angles", "wrinkles",
         "penetration_gate"}
    )

    def to_dict(self) -> dict:
        """Serialise this config to a plain, JSON-round-trippable dict.

        The dict carries a ``config_version`` key and one entry per
        dataclass field.  Non-primitive fields are handled explicitly:
        materials become ``{"preset": name}`` or ``{"custom": {...}}``,
        ``angles`` the resolved ply list, ``wrinkles`` a list of plain
        dicts, and ``penetration_gate`` a ``{"preset": name}`` reference.
        Because the values are the *resolved* ones (post ``__post_init__``),
        the file is self-contained and ``load`` -> ``to_dict`` is
        idempotent.

        See Also
        --------
        from_dict : Inverse constructor.
        save_json, load_json : File helpers.
        """
        out: dict = {"config_version": CONFIG_VERSION}
        for f in fields(self):
            value = getattr(self, f.name)
            if f.name in ("material", "resin_pocket_material"):
                out[f.name] = _material_to_jsonable(value)
            elif f.name == "wrinkles":
                out[f.name] = (
                    None if value is None else [asdict(s) for s in value]
                )
            elif f.name == "load_state":
                out[f.name] = None if value is None else asdict(value)
            elif f.name == "penetration_gate":
                out[f.name] = _gate_to_jsonable(value)
            elif f.name == "angles":
                out[f.name] = (
                    None if value is None else [float(a) for a in value]
                )
            elif isinstance(value, np.ndarray):
                out[f.name] = value.tolist()
            elif isinstance(value, (list, tuple)):
                out[f.name] = [
                    float(v) if isinstance(v, (np.integer, np.floating)) else v
                    for v in value
                ]
            elif isinstance(value, np.integer):
                out[f.name] = int(value)
            elif isinstance(value, np.floating):
                out[f.name] = float(value)
            else:
                out[f.name] = value
        return out

    @classmethod
    def from_dict(cls, data: dict) -> AnalysisConfig:
        """Reconstruct an :class:`AnalysisConfig` from :meth:`to_dict` output.

        The ``config_version`` must match :data:`CONFIG_VERSION` and every
        remaining key must name a dataclass field — an unknown key or a
        version mismatch raises :class:`ValueError` naming the offender
        rather than being silently ignored.  Construction runs the normal
        ``__post_init__`` / ``_validate`` path, so a malformed file fails
        with the same messages a bad in-code config would.

        Raises
        ------
        ValueError
            On version mismatch, unknown keys, or invalid field values.
        """
        if not isinstance(data, dict):
            raise ValueError(
                f"AnalysisConfig.from_dict expects a dict, got "
                f"{type(data).__name__}"
            )
        payload = dict(data)
        version = payload.pop("config_version", None)
        if version != CONFIG_VERSION:
            raise ValueError(
                f"AnalysisConfig.from_dict: unsupported config_version "
                f"{version!r} (expected {CONFIG_VERSION})"
            )
        valid = {f.name for f in fields(cls)}
        unknown = set(payload) - valid
        if unknown:
            raise ValueError(
                f"AnalysisConfig.from_dict: unknown key(s) "
                f"{sorted(unknown)}; valid fields are {sorted(valid)}"
            )
        kwargs = dict(payload)
        for mkey in ("material", "resin_pocket_material"):
            if mkey in kwargs:
                kwargs[mkey] = _material_from_jsonable(kwargs[mkey], field=mkey)
        if kwargs.get("load_state") is not None:
            raw = kwargs["load_state"]
            if not isinstance(raw, dict):
                raise ValueError(
                    "AnalysisConfig.from_dict: load_state must be null or "
                    f"an object, got {type(raw).__name__}"
                )
            valid_ls = {f.name for f in fields(LoadState)}
            unknown_ls = set(raw) - valid_ls
            if unknown_ls:
                raise ValueError(
                    f"AnalysisConfig.from_dict: unknown load_state key(s) "
                    f"{sorted(unknown_ls)}; valid fields are "
                    f"{sorted(valid_ls)}"
                )
            kwargs["load_state"] = LoadState(**raw)
        if kwargs.get("wrinkles") is not None:
            specs = kwargs["wrinkles"]
            if not isinstance(specs, list):
                raise ValueError(
                    "AnalysisConfig.from_dict: wrinkles must be null or a "
                    f"list, got {type(specs).__name__}"
                )
            rebuilt: list[WrinkleSpec] = []
            for i, entry in enumerate(specs):
                if not isinstance(entry, dict):
                    raise ValueError(
                        f"AnalysisConfig.from_dict: wrinkles[{i}] must be a "
                        f"dict, got {type(entry).__name__}"
                    )
                spec_fields = {f.name for f in fields(WrinkleSpec)}
                extra = set(entry) - spec_fields
                if extra:
                    raise ValueError(
                        f"AnalysisConfig.from_dict: wrinkles[{i}] has "
                        f"unknown key(s) {sorted(extra)}"
                    )
                rebuilt.append(WrinkleSpec(**entry))
            kwargs["wrinkles"] = rebuilt
        if "penetration_gate" in kwargs:
            kwargs["penetration_gate"] = _gate_from_jsonable(
                kwargs["penetration_gate"]
            )
        return cls(**kwargs)

    def to_json(self) -> str:
        """Return the config as a deterministic JSON string."""
        return json.dumps(self.to_dict(), indent=2, sort_keys=True)

    def save_json(self, path: str | Path) -> None:
        """Write the config to ``path`` as JSON (parent dirs created)."""
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(self.to_json(), encoding="utf-8")

    @classmethod
    def load_json(cls, path: str | Path) -> AnalysisConfig:
        """Load a config from a JSON file written by :meth:`save_json`."""
        text = Path(path).read_text(encoding="utf-8")
        return cls.from_dict(json.loads(text))

    def save_yaml(self, path: str | Path) -> None:
        """Write the config as YAML.

        Optional: requires PyYAML (not a hard dependency).  Raises
        :class:`ImportError` with an actionable message when it is absent.
        """
        try:
            import yaml  # type: ignore[import-untyped]
        except ImportError as exc:  # pragma: no cover - env-dependent
            raise ImportError(
                "PyYAML is required for YAML config I/O; install it "
                "(`pip install pyyaml`) or use the JSON helpers."
            ) from exc
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(
            yaml.safe_dump(self.to_dict(), sort_keys=True), encoding="utf-8"
        )

    @classmethod
    def load_yaml(cls, path: str | Path) -> AnalysisConfig:
        """Load a config from a YAML file (requires PyYAML)."""
        try:
            import yaml  # type: ignore[import-untyped]
        except ImportError as exc:  # pragma: no cover - env-dependent
            raise ImportError(
                "PyYAML is required for YAML config I/O; install it "
                "(`pip install pyyaml`) or use the JSON helpers."
            ) from exc
        text = Path(path).read_text(encoding="utf-8")
        return cls.from_dict(yaml.safe_load(text))

    def save(self, path: str | Path) -> None:
        """Save the config, dispatching to YAML for ``.yaml`` / ``.yml``."""
        if Path(path).suffix.lower() in (".yaml", ".yml"):
            self.save_yaml(path)
        else:
            self.save_json(path)

    @classmethod
    def load(cls, path: str | Path) -> AnalysisConfig:
        """Load a config, dispatching to YAML for ``.yaml`` / ``.yml``."""
        if Path(path).suffix.lower() in (".yaml", ".yml"):
            return cls.load_yaml(path)
        return cls.load_json(path)


# ======================================================================
# Results
# ======================================================================

@dataclass
class AnalysisResults:
    """Aggregated results from a complete wrinkle analysis.

    Attributes
    ----------
    config : AnalysisConfig
        The configuration used for this analysis.
    mesh : MeshData
        The generated finite element mesh.
    wrinkle_config : WrinkleConfiguration
        The wrinkle configuration object.
    laminate : Laminate
        The laminate definition.
    morphology_factor : float
        Aggregate morphology factor M_f.
    max_angle_rad : float
        Maximum fibre misalignment angle (radians).
    effective_angle_rad : float
        Effective fibre angle theta_eff (radians).
    damage_index : float
        Interlaminar damage index D.
    analytical_knockdown : float
        Analytical combined knockdown factor (the *ultimate* fibre-failure
        KD for tension; the kink-band KD for compression).
    analytical_onset_knockdown : float or None
        Delamination-onset (first-load-drop) knockdown for tension
        loading, computed from a curved-beam mode-mixity criterion using
        ``material.GIc`` and ``material.GIIc``.  ``None`` when the
        material does not provide both fracture toughnesses or when the
        loading is compression.  Always strictly below
        ``analytical_knockdown``.
    analytical_strength_MPa : float
        Analytical predicted failure stress (MPa).
    field_results : FieldResults or None
        FE solution fields (displacement, stress, strain).
    failure_report : LaminateFailureReport or None
        Multi-criteria failure evaluation.
    failure_indices : dict or None
        Per-criterion FE failure index fields.
    """

    config: AnalysisConfig
    mesh: MeshData | None = None
    wrinkle_config: WrinkleConfiguration | None = None
    laminate: Laminate | None = None

    # Analytical predictions
    morphology_factor: float = 1.0
    max_angle_rad: float = 0.0
    effective_angle_rad: float = 0.0
    mesh_max_angle_rad: float = 0.0  # max fiber angle from FE mesh (accounts for decay)
    damage_index: float = 0.0
    analytical_knockdown: float = 1.0
    analytical_modulus_knockdown: float = 1.0
    """Analytical axial Young's-modulus (stiffness) knockdown
    ``E_x / E_x0`` — a closed-form CLT membrane series-average of the
    off-axis lamina modulus over the wrinkle profile. Populated for
    arbitrary layups and multi-wrinkle layouts: the unidirectional
    single-wrinkle case uses the scalar :func:`_profile_modulus_knockdown`,
    everything else the laminate generalization
    :func:`_laminate_modulus_knockdown` (which reduces to the UD result for
    ``[0]_n``). Stays ``1.0`` only for a degenerate (zero-amplitude)
    wrinkle. Loading-independent; the closed-form companion to the FE
    :attr:`modulus_retention`."""
    analytical_onset_knockdown: float | None = None
    analytical_strength_MPa: float = 0.0
    gamma_Y_eff: float = 0.02  # layup-dependent effective yield strain
    # Tension mechanism decomposition (only for tension loading)
    tension_mechanisms: dict | None = None  # {kd_fiber, kd_matrix, kd_oop, mode, ...}

    # FE results
    field_results: FieldResults | None = None
    failure_report: LaminateFailureReport | None = None
    failure_indices: dict | None = None
    failure_modes: dict | None = None  # {criterion: (n_elem, n_gauss) str array}

    # Retention factor (wrinkled / pristine)
    retention_factors: dict | None = None  # {criterion_name: float}
    baseline_fi: dict | None = None  # {criterion_name: float} pristine max FI

    # Proportional load factor (issue #275).  Populated only when
    # ``AnalysisConfig.load_state`` is set — the general-load-state
    # answer to "how much strength is left".
    load_state_factor: float | None = None
    """Scalar the whole :attr:`AnalysisConfig.load_state` is multiplied by
    to reach first failure on the **wrinkled** coupon: ``max FI(lambda
    sigma) = 1``.  ``lambda > 1`` means the applied state is survivable
    with that much margin; ``lambda < 1`` means it already fails.

    This is the combined-load generalisation of a strength: scaling the
    whole state proportionally is what "how much of this load can it
    take" means when the load is not a single number.  For a uniaxial
    state it reduces to the usual definition.  ``None`` when no load
    state was given, or when the root could not be bracketed.

    NOTE the name: ``load_factor`` alone is already taken in this package
    for the **CLT first-ply-failure** factor (``1 / FI``, see
    ``failure/evaluator.py`` and the top-level ``load_factor`` key of
    ``results_to_dict``).  That is a different quantity computed on a
    different path, so this one is named for what it scales.  Do not
    shorten it back."""

    load_state_factor_pristine: float | None = None
    """The same factor for the flat (no-wrinkle) baseline, solved under
    the identical load state."""

    load_state_factor_knockdown: float | None = None
    """``load_factor / load_state_factor_pristine`` — the strength knockdown
    under the general load state.  Reduces to the uniaxial strength
    knockdown when the state is uniaxial."""

    # Modulus retention (E_wrinkled / E_pristine from FE)
    modulus_retention_failed: bool = False
    """``True`` when the local-σ₁₁ modulus-retention computation raised and
    :attr:`modulus_retention` was forced to the ``1.0`` fallback. Lets callers
    distinguish a genuine no-knockdown ``1.0`` from a failed computation (see
    also the WARNING logged with ``exc_info`` when this fires)."""

    modulus_retention: float = 1.0
    """FE axial-modulus retention from the **local** fibre-direction-stress
    proxy: ``E_eff = <σ₁₁> / ε_applied`` (mean element-frame σ₁₁ over the
    coupon), wrinkled vs pristine.  Because it averages the *local* fibre
    stress rather than the coupon's global axial response, it over-predicts
    the modulus retention (flatter on the amplitude / penetration / position
    axes) than the measured ``E_x / E_x0``.  Kept for backward
    compatibility; prefer :attr:`modulus_retention_global` for the
    coupon-level stiffness knockdown.

    Unlike :attr:`modulus_retention_global` — which is solved at
    ``delta_T = 0`` on purpose — this proxy is derived from the single
    thermally loaded stress field, so with ``AnalysisConfig.delta_T != 0``
    both sides of its wrinkled/pristine ratio carry the cure residual
    stress (issue #273 Stage 2).  The offset largely cancels in the ratio
    but not exactly, so the two modulus numbers can move apart under a
    thermal load; the reaction-based one is the stiffness measurement."""

    modulus_retention_global_failed: bool = False
    """``True`` when the global reaction-based modulus-retention computation
    raised and :attr:`modulus_retention_global` was forced to the ``1.0``
    fallback (a WARNING is logged with ``exc_info`` when this fires)."""

    modulus_retention_global: float = 1.0
    """FE axial-modulus retention from the **global** reaction response:
    ``E_eff = σ_nominal / ε_applied`` with
    ``σ_nominal = R / A`` — the total axial reaction on the loaded
    (``x_max``) face over the cross-section area ``Ly·Lz`` — computed
    wrinkled vs pristine to give a true coupon-level ``E_x / E_x0``.

    This reuses the same reaction-based nominal-stress pattern as the
    progressive-damage solver (``reaction = sum((K @ u)[xmax_dofs])`` over
    the loaded-face x-DOFs).  Unlike the local σ₁₁ proxy in
    :attr:`modulus_retention`, it captures load redistribution around the
    wrinkle, so it matches the measured modulus knockdown more closely (and
    is correspondingly lower for a wrinkled coupon)."""

    # ------------------------------------------------------------------
    # Progressive-damage results — populated when
    # AnalysisConfig.enable_progressive_damage = True.
    # ------------------------------------------------------------------
    progressive_strength_MPa: float = 0.0
    """Predicted ultimate compressive strength of the wrinkled coupon
    (peak carried nominal stress over the load history, MPa)."""

    progressive_pristine_strength_MPa: float = 0.0
    """Predicted ultimate strength of the pristine (flat) baseline, MPa."""

    progressive_knockdown: float = 1.0
    """Progressive-damage strength knockdown
    ``progressive_strength_MPa / progressive_pristine_strength_MPa``."""

    progressive_history: list | None = None
    """``(applied_strain, nominal_stress)`` samples for the wrinkled run."""

    # ------------------------------------------------------------------
    # CZM results — only populated when AnalysisConfig.enable_czm = True.
    # ------------------------------------------------------------------
    czm_damage: np.ndarray | None = None
    """Cohesive damage variable per (interface element, Gauss point).
    Shape ``(n_iface_elems, n_gauss)``."""

    czm_separation: np.ndarray | None = None
    """Displacement jump in the local cohesive frame per (interface
    element, Gauss point, component).  Shape ``(n_iface_elems, n_gauss,
    3)`` with components ``(delta_n, delta_s, delta_t)``."""

    czm_traction: np.ndarray | None = None
    """Cohesive traction in the local frame at each Gauss point.  Shape
    ``(n_iface_elems, n_gauss, 3)``."""

    czm_energy_dissipated: float | None = None
    """Total cohesive energy dissipated across all interfaces (N*mm)."""

    czm_energy_per_interface: dict[int, float] | None = None
    """Per-interface dissipated energy keyed by ply-interface index."""

    czm_crack_length_per_interface: dict[int, float] | None = None
    """Per-interface crack length (in mm) keyed by ply-interface index.
    Computed as the in-plane area of elements with ``damage > 0.99``,
    divided by the mesh width (so the reported quantity has units of
    length along the wrinkle/crack direction)."""

    czm_load_displacement: np.ndarray | None = None
    """``(n_inc, 2)`` array of ``(lambda, ||u||)`` samples from the
    incremental Newton-Raphson run."""

    czm_converged: bool | None = None
    """Whether every load increment converged."""

    czm_failure_diagnostics: dict | None = None
    """Diagnostics record for the first load increment that failed to
    converge (``None`` when the CZM solve converged or was not run).  Carries
    the increment index / load fraction, iteration count, final residual /
    BC-violation / step norms, the tail of the residual history, line-search
    status, and a classified ``failure_reason``.  See
    :meth:`wrinklefe.solver.nonlinear.NewtonRaphsonSolver._record_increment_diag`."""

    czm_failure_hint: str | None = None
    """A single actionable tuning hint derived from
    :attr:`czm_failure_diagnostics` (``None`` when converged).  Names the knob
    to reach for — ``czm_n_load_increments``, ``czm_newton_tol``, the applied
    strain, or the Newton iteration cap."""

    czm_interfaces_used: list[int] | None = None
    """Ply-interface indices that actually received cohesive elements."""

    czm_delamination_report: LaminateFailureReport | None = None
    """Delamination failure report shaped like the other failure
    criteria, populated by :mod:`wrinklefe.failure.delamination`."""

    czm_element_centroids: np.ndarray | None = None
    """``(n_iface_elems, 2)`` array of in-plane ``(x, y)`` centroids of the
    cohesive interface elements, in the reference (undeformed) configuration.

    Populated alongside the other ``czm_*`` fields by ``_run_czm_path`` so
    that visualization wrappers (e.g.
    :func:`wrinklefe.viz.czm_overview_figure`) can colour the interface
    plane without needing access to the assembler / cohesive-element list.
    Same row order as ``czm_damage``."""

    def summary(self) -> str:
        """Generate a comprehensive text summary.

        Returns
        -------
        str
            Multi-line summary of all analysis results.
        """
        cfg = self.config
        lines = [
            "=" * 65,
            "  WrinkleFE Analysis Results",
            "=" * 65,
            "",
            "  Configuration:",
            f"    Morphology:      {cfg.morphology}",
            f"    Amplitude:       {cfg.amplitude:.3f} mm",
            f"    Wavelength:      {cfg.wavelength:.1f} mm",
            f"    Width:           {cfg.width:.1f} mm",
            f"    Amplitude profile: {cfg.amplitude_profile} "
            f"(d={cfg.amplitude_profile_decay_length}, "
            f"axis={cfg.amplitude_profile_axis})",
            f"    Loading:         {cfg.loading}",
            f"    Applied strain:  {cfg.applied_strain:.4f}",
            "",
            "  Analytical Predictions:",
            f"    Morphology factor M_f:  {self.morphology_factor:.4f}",
            f"    Max angle theta_max:    {np.degrees(self.max_angle_rad):.2f} deg "
            f"({self.max_angle_rad:.4f} rad)",
            f"    Effective angle:        {np.degrees(self.effective_angle_rad):.2f} deg "
            f"({self.effective_angle_rad:.4f} rad)",
            f"    Damage index D:         {self.damage_index:.4f}",
            f"    Combined knockdown:     {self.analytical_knockdown:.4f}",
            f"    Modulus knockdown:      {self.analytical_modulus_knockdown:.4f}",
            f"    Predicted strength:     {self.analytical_strength_MPa:.1f} MPa",
        ]

        # Thermal / cure-residual load (issue #273). Emitted only
        # when a temperature change is set, so the default summary is
        # unchanged. The sign convention is restated inline because a
        # reader who mistakes it flips the residual matrix stress from
        # tension to compression.
        if cfg.delta_T != 0.0:
            sense = "cool-down from" if cfg.delta_T < 0 else "heat-up above"
            lines.extend([
                "",
                "  Thermal / cure-residual load:",
                f"    delta_T: {cfg.delta_T:+.1f} deg C "
                f"({sense} the stress-free/cure state)",
                "    Applied to the CLT ply stresses and first-ply-failure "
                "report,",
                "    and to the FE fields via the element thermal "
                "initial-strain",
                "    load vector (issue #273).  The measured modulus is "
                "solved at",
                "    delta_T = 0 — a residual load offset is not a "
                "stiffness change.",
            ])

        if self.mesh is not None:
            lines.extend([
                "",
                "  Mesh:",
                f"    Nodes:    {self.mesh.n_nodes}",
                f"    Elements: {self.mesh.n_elements}",
                f"    DOFs:     {self.mesh.n_dof}",
            ])

        if self.field_results is not None:
            max_disp, _ = self.field_results.max_displacement()
            lines.extend([
                "",
                "  FE Results:",
                f"    Max displacement: {max_disp:.6e} mm",
                f"    Modulus retention (local σ₁₁):  "
                f"{self.modulus_retention:.4f}"
                + (
                    "  (computation failed — fallback)"
                    if self.modulus_retention_failed
                    else ""
                ),
                f"    Modulus retention (global E_x): "
                f"{self.modulus_retention_global:.4f}"
                + (
                    "  (computation failed — fallback)"
                    if self.modulus_retention_global_failed
                    else ""
                ),
            ])

        if self.progressive_history is not None:
            lines.extend([
                "",
                "  Progressive Damage (ultimate strength):",
                f"    Wrinkled strength:  {self.progressive_strength_MPa:.1f} MPa",
                f"    Pristine strength:  "
                f"{self.progressive_pristine_strength_MPa:.1f} MPa",
                f"    Strength knockdown: {self.progressive_knockdown:.4f}",
            ])

        if self.czm_damage is not None:
            max_d = float(np.max(self.czm_damage)) if self.czm_damage.size else 0.0
            mean_d = float(np.mean(self.czm_damage)) if self.czm_damage.size else 0.0
            energy = (
                self.czm_energy_dissipated
                if self.czm_energy_dissipated is not None
                else 0.0
            )
            iface_str = (
                ",".join(str(i) for i in self.czm_interfaces_used)
                if self.czm_interfaces_used else "(none)"
            )
            lines.extend([
                "",
                "  Cohesive Zone Modeling (delamination):",
                f"    Interfaces:        {iface_str}",
                f"    Max damage:        {max_d:.4f}",
                f"    Mean damage:       {mean_d:.4f}",
                f"    Energy dissipated: {energy:.4e} N*mm",
                f"    Converged:         {self.czm_converged}",
            ])
            # When the solve did not converge, surface the actionable
            # tuning hint right below the convergence flag so a CLI /
            # summary reader learns which knob to turn.
            if self.czm_converged is False and self.czm_failure_hint:
                lines.append(f"    Hint:              {self.czm_failure_hint}")

        lines.append("=" * 65)
        return "\n".join(lines)


# ======================================================================
# Sweep parallelism helpers (issue #260)
# ======================================================================


def _iterative_solver_kwargs(cfg: AnalysisConfig) -> dict:
    """Map the iterative-solver ``AnalysisConfig`` fields to
    :class:`~wrinklefe.solver.static.StaticSolver` keyword arguments.

    The controls are inert on the direct path but are threaded through
    every solver built from ``cfg`` so the iterative path is fully driven
    by the config (issue #265).
    """
    return {
        "iterative_rtol": cfg.iterative_rtol,
        "iterative_maxiter": cfg.iterative_maxiter,
        "ilu_drop_tol": cfg.ilu_drop_tol,
        "ilu_fill_factor": cfg.ilu_fill_factor,
        "preconditioner": cfg.preconditioner,
    }


def _mechanical_bcs(cfg: AnalysisConfig, mesh: MeshData) -> list:
    """Boundary conditions for the mechanical load the config asks for.

    ``AnalysisConfig.load_state`` (issue #275), when set, is the whole
    mechanical load: its resultants become the self-equilibrated traction
    set of
    :meth:`~wrinklefe.solver.boundary.BoundaryHandler.load_state_to_bcs`.
    Otherwise the legacy uniaxial displacement BCs from
    ``applied_strain`` are used, unchanged and bit-identical.

    Shared by the wrinkled run and its pristine retention baseline so the
    two are always solved under the same load.
    """
    if cfg.load_state is not None:
        return BoundaryHandler.load_state_to_bcs(cfg.load_state, mesh)
    return BoundaryHandler.compression_bcs(
        mesh, applied_strain=cfg.applied_strain
    )


def _proportional_load_factor(
    max_fi_at: Callable[[float], float],
    *,
    max_factor: float = 1.0e6,
    min_factor: float = 1.0e-9,
) -> float | None:
    """Scalar the whole load state is multiplied by to reach first failure.

    Solves ``max FI(lambda * sigma) = 1`` for ``lambda`` (issue #275).
    Because the FE solve is linear the stress field is linear in
    ``lambda``, so the search scales the *stored* field instead of
    re-solving — a root-find over a handful of vectorised criterion
    evaluations.

    .. warning::
       Linearity is what makes scaling valid.  With a thermal load the
       field is *affine* (``sigma = lambda sigma_mech + sigma_th``), and
       this would be wrong.  ``AnalysisConfig`` rejects ``delta_T``
       alongside ``load_state`` for exactly that reason; if that guard is
       ever relaxed, this routine must split the two contributions first.

    For a uniaxial state this reduces to the usual strength definition:
    the load is scaled until the coupon first fails.

    Parameters
    ----------
    max_fi_at : callable
        ``lambda -> max failure index`` over the field and all criteria.
    max_factor, min_factor : float
        Bracket limits.  ``None`` is returned when no root lies inside
        them (an unloaded field, or one that cannot be driven to failure).

    Returns
    -------
    float or None
        The critical load factor, or ``None`` if it could not be bracketed.
    """
    from scipy.optimize import brentq

    def residual(lam: float) -> float:
        return float(max_fi_at(lam)) - 1.0

    fi_unit = max_fi_at(1.0)
    if not math.isfinite(fi_unit) or fi_unit <= 0.0:
        return None

    if residual(1.0) > 0.0:
        # Already past first failure at the applied load: lambda < 1.
        hi = 1.0
        lo = 0.5
        while lo > min_factor and residual(lo) > 0.0:
            hi = lo
            lo *= 0.5
        if residual(lo) > 0.0:
            return None
    else:
        lo = 1.0
        hi = 2.0
        while hi < max_factor and residual(hi) < 0.0:
            lo = hi
            hi *= 2.0
        if residual(hi) < 0.0:
            return None

    try:
        return float(brentq(residual, lo, hi, rtol=1.0e-6))
    except (ValueError, RuntimeError):
        return None


def _sweep_run_one(
    cfg: AnalysisConfig, analytical_only: bool
) -> AnalysisResults:
    """Run one sweep point.  Module-level so it pickles for
    ``ProcessPoolExecutor`` workers."""
    return WrinkleAnalysis(cfg).run(analytical_only=analytical_only)


def _replace_swept(
    base_config: AnalysisConfig,
    parameter: str,
    value: float,
    *,
    analytical_only: bool | None = None,
) -> AnalysisConfig:
    """Clone *base_config* with *parameter* set to *value*.

    Resets the ``domain_length`` sentinel when sweeping ``wavelength``
    from a config that left it auto-derived (``domain_length ==
    3 * wavelength``) — ``replace`` re-runs ``__post_init__`` but does
    not clear un-passed fields, so a stale derived ``domain_length``
    would silently hold the domain fixed while the wavelength moved.
    An explicitly pinned ``domain_length`` is preserved untouched.

    Shared by :meth:`WrinkleAnalysis.parametric_sweep` and
    :mod:`wrinklefe.goalseek` so the two cannot diverge.

    Parameters
    ----------
    base_config : AnalysisConfig
        Configuration to clone.
    parameter : str
        Field name to override.
    value : float
        New value for *parameter*.
    analytical_only : bool or None, optional
        When not ``None``, bake this solve-path choice into the clone so
        the returned config records the path it was actually run on.

    Returns
    -------
    AnalysisConfig
        Validated clone (``replace`` re-invokes ``__post_init__``).
    """
    overrides: dict[str, Any] = {parameter: value}
    if (
        parameter == "wavelength"
        and base_config.domain_length == 3.0 * base_config.wavelength
    ):
        overrides["domain_length"] = 0.0
    if analytical_only is not None:
        overrides["analytical_only"] = analytical_only
    return replace(base_config, **overrides)


def _resolve_sweep_workers(n_workers: int) -> int:
    """Validate and resolve a sweep worker count (``0`` -> all cores)."""
    if not isinstance(n_workers, int) or isinstance(n_workers, bool):
        raise ValueError(
            f"n_workers must be an int >= 0, got {n_workers!r}"
        )
    if n_workers < 0:
        raise ValueError(f"n_workers must be >= 0, got {n_workers}")
    if n_workers == 0:
        return os.cpu_count() or 1
    return n_workers


# ======================================================================
# Main analysis class
# ======================================================================

class WrinkleAnalysis:
    """High-level orchestrator for wrinkled laminate analysis.

    This class chains together all modules in the WrinkleFE framework:
    material → laminate → wrinkle → mesh → solve → failure → statistics.

    Parameters
    ----------
    config : AnalysisConfig
        Complete analysis configuration.

    Examples
    --------
    >>> config = AnalysisConfig(morphology="concave", amplitude=0.366)
    >>> analysis = WrinkleAnalysis(config)
    >>> result = analysis.run()  # doctest: +SKIP
    >>> print(f"Strength = {result.analytical_strength_MPa:.1f} MPa")  # doctest: +SKIP
    """

    def __init__(self, config: AnalysisConfig) -> None:
        self.config = config

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(
        self,
        analytical_only: bool | None = None,
        progress_callback: Callable[[str, float], None] | None = None,
    ) -> AnalysisResults:
        """Execute the complete analysis pipeline.

        Parameters
        ----------
        analytical_only : bool, optional
            If True, skip the FE assembly path (mesh generation, static
            solve, failure evaluation, retention factors) and return only
            the analytical predictions.  If None (default), the
            ``AnalysisConfig.analytical_only`` field is used.
        progress_callback : callable, optional
            Optional progress reporter invoked at each phase boundary as
            ``progress_callback(label, fraction)`` where ``label`` is a
            short human-readable phase name and ``fraction`` is the
            cumulative completion in ``[0, 1]``.  Defaults to ``None``
            (no reporting), so non-Streamlit callers (CLI, tests) are
            unaffected.  The callback always fires at least once with
            ``fraction == 1.0`` on successful completion.

        Steps
        -----
        1. Build laminate from material and ply angles.
        2. Create wrinkle profile and morphology configuration.
        3. Compute analytical predictions (knockdown, strength).
        4. Generate mesh with wrinkle deformation.
        5. Run static FE analysis.
        6. Evaluate failure criteria on the FE stress field.

        Returns
        -------
        AnalysisResults
            Complete results from all analysis steps.
        """
        cfg = self.config
        # interface_1 / interface_2 are filled in __post_init__.
        assert cfg.interface_1 is not None and cfg.interface_2 is not None
        # _validate constrained these to the literal sets that
        # WrinkleConfiguration expects.
        amp_profile = cast(
            Literal["constant", "gaussian", "linear"], cfg.amplitude_profile
        )
        amp_axis = cast(Literal["x", "y"], cfg.amplitude_profile_axis)
        if analytical_only is None:
            analytical_only = cfg.analytical_only

        # A transverse surface only exists in the FE mesh, so a run-time
        # analytical_only override must be rejected too (the construction-time
        # _validate only sees cfg.analytical_only). Fail fast rather than
        # silently dropping the requested through-width variation (#300).
        if analytical_only and cfg.transverse_mode != "uniform":
            raise ValueError(
                "AnalysisConfig.transverse_mode="
                f"{cfg.transverse_mode!r} requires the FE path but this run "
                "was invoked with analytical_only=True. Run with "
                "analytical_only=False or set transverse_mode='uniform'."
            )

        # Multi-wrinkle FE solve (issue #252): overlapping and
        # non-overlapping layouts run through the linear FE path, and
        # (issue #283) through the CZM path — cohesive layers are
        # inserted along the full length of every interface a wrinkle
        # nominates, so delaminations can link up between adjacent
        # wrinkles sharing an interface.
        results = AnalysisResults(config=cfg)

        # Phase weights (sum to 1.0 on the full FE path).  The FE solve
        # dominates wall-clock for typical mesh densities, so it gets the
        # largest slice.  In analytical-only mode the analytical step is
        # rescaled to 1.0 below.
        #
        #   build mesh / wrinkle geom : 0.05
        #   analytical predictions    : 0.05
        #   FE assembly (mesh build)  : 0.15
        #   FE solve                  : 0.50
        #   failure evaluation        : 0.15  (stress recovery + criteria)
        #   retention factors         : 0.10
        def _report(label: str, fraction: float) -> None:
            if progress_callback is not None:
                progress_callback(label, max(0.0, min(1.0, float(fraction))))

        logger.info(
            "Starting analysis: morphology=%s A=%.4g lambda=%.4g "
            "analytical_only=%s czm=%s",
            cfg.morphology, cfg.amplitude, cfg.wavelength,
            analytical_only, cfg.enable_czm,
        )

        _report("Building laminate", 0.0)

        # 1. Build laminate
        laminate = self._build_laminate()
        results.laminate = laminate

        # 2. Create wrinkle configuration (centered in specimen)
        wrinkle_center = cfg.domain_length / 2.0

        # Multi-wrinkle override: build a WrinkleConfiguration directly
        # from the user-supplied WrinkleSpec list.  Each spec carries its
        # own geometry (A, lambda, w) and is placed at its own ply
        # interface with its own phase offset.  Bypasses the single/dual
        # name dispatch entirely.
        if cfg.wrinkles is not None:
            placements = []
            for spec in cfg.wrinkles:
                spec_profile = GaussianSinusoidal(
                    amplitude=spec.amplitude,
                    wavelength=spec.wavelength,
                    width=spec.width,
                    center=wrinkle_center,
                )
                placements.append(
                    WrinklePlacement(
                        profile=spec_profile,
                        ply_interface=spec.ply_interface,
                        phase_offset=spec.phase_offset,
                    )
                )
            is_graded = cfg.morphology.lower().strip() == "graded"
            wrinkle_config = WrinkleConfiguration(
                placements,
                decay_mode="graded" if is_graded else "default",
                decay_floor=cfg.decay_floor if is_graded else 0.0,
                amplitude_profile=amp_profile,
                amplitude_profile_decay_length=cfg.amplitude_profile_decay_length,
                amplitude_profile_axis=amp_axis,
            )
            results.wrinkle_config = wrinkle_config
        else:
            base_profile = GaussianSinusoidal(
                amplitude=cfg.amplitude,
                wavelength=cfg.wavelength,
                width=cfg.width,
                center=wrinkle_center,
            )
            # Through-width variation (#300): wrap the x-only profile in a
            # 3-D transverse surface only when a non-uniform mode is
            # requested. ``"uniform"`` leaves the bare GaussianSinusoidal
            # untouched, so the default path stays bit-identical. The
            # morphology layer already consumes WrinkleSurface3D end-to-end
            # (node deformation + fibre-angle fields); multi-wrinkle/CZM
            # combinations are rejected earlier in _validate.
            profile: WrinkleProfile | WrinkleSurface3D = (
                self._wrap_transverse_surface(base_profile)
                if cfg.transverse_mode != "uniform"
                else base_profile
            )
            if cfg.phase is not None and (
                cfg.morphology.lower().strip() not in SINGLE_WRINKLE_MODES
            ):
                # Explicit phase overrides the named-morphology phase so
                # arbitrary dual-wrinkle phase offsets can be analysed/swept.
                wrinkle_config = WrinkleConfiguration.dual_wrinkle(
                    profile,
                    interface1=cfg.interface_1,
                    interface2=cfg.interface_2,
                    phase=float(cfg.phase),
                    amplitude_profile=amp_profile,
                    amplitude_profile_decay_length=cfg.amplitude_profile_decay_length,
                    amplitude_profile_axis=amp_axis,
                )
                wrinkle_config.decay_floor = max(0.0, min(1.0, cfg.decay_floor))
            else:
                wrinkle_config = WrinkleConfiguration.from_morphology_name(
                    cfg.morphology, profile,
                    interface1=cfg.interface_1,
                    interface2=cfg.interface_2,
                    decay_floor=cfg.decay_floor,
                    amplitude_profile=amp_profile,
                    amplitude_profile_decay_length=cfg.amplitude_profile_decay_length,
                    amplitude_profile_axis=amp_axis,
                    # tool_flat: the tool-flat surface(s) reuse the
                    # surface-pocket side, and the transition-zone width
                    # controls the pocket depth / inversion margin.
                    tool_side=cast(
                        Literal["top", "bottom", "both"], cfg.surface_pocket_side
                    ),
                    surface_transition_plies=cfg.surface_transition_plies,
                )
            results.wrinkle_config = wrinkle_config

        # Through-thickness wrinkle position (item D.5): the graded decay
        # is centred here (0.5 = mid-plane).  Off-mid values place the
        # wrinkle nearer a surface (Li 2025 S-A-2).
        wrinkle_config.wrinkle_z_position = float(cfg.wrinkle_z_position)
        results.wrinkle_config = wrinkle_config

        _report("Computing analytical predictions", 0.05)

        # 3. Analytical predictions
        self._compute_analytical(results, wrinkle_config)

        # 3b. CLT first-ply-failure under the thermal load (issue #273).
        #
        # Only runs when ``delta_T != 0``. The closed-form knockdown above
        # is a fibre-kinking / fracture model with no temperature term, so
        # without this step a thermally-loaded analytical run would produce
        # no output that depends on ``delta_T`` at all — the silent no-op
        # this issue exists to remove. Cure-induced residual stress is
        # present at zero mechanical load, so the ply-level CLT report is
        # the quantity that actually carries it.
        #
        # Left off for ``delta_T == 0`` so the default analytical run is
        # bit-identical (and ``failure_report`` keeps its documented
        # "absent on analytical_only runs" contract, which the JSON/CSV
        # exporters rely on for their load-factor fallback).
        if cfg.delta_T != 0.0:
            self._evaluate_clt_failure(results, laminate)

        # In analytical-only mode, skip mesh, solve, failure, and retention.
        if analytical_only:
            _report("Analytical predictions complete", 1.0)
            logger.info(
                "Analysis complete (analytical only): knockdown=%s",
                results.analytical_knockdown,
            )
            return results

        _report("Assembling FE mesh", 0.10)

        # 4. Generate mesh
        mesh_gen = WrinkleMesh(
            laminate=laminate,
            wrinkle_config=wrinkle_config,
            Lx=cfg.domain_length,
            Ly=cfg.domain_width,
            nx=cfg.nx,
            ny=cfg.ny,
            nz_per_ply=cfg.nz_per_ply,
        )
        mesh = mesh_gen.generate()

        # 4a2. Resin-pocket material zone (Li et al. 2024/2025).
        if cfg.enable_resin_pocket:
            self._attach_resin_pocket(mesh, laminate)

        # 4a3. Surface resin pockets under a tool-flat surface (issue #361).
        # Runs after the crest lens so the two compose (per-element max)
        # rather than overwrite.  The compaction Vf gradient (#379) is the
        # continuous generalization of this binary tag — the resin-rich
        # trough is simply the low-Vf end of the same field — so when the
        # gradient is on the binary tagging is skipped rather than applied
        # on top of it (which would soften the trough twice).
        if cfg.enable_surface_resin_pockets and not cfg.enable_vf_gradient:
            self._attach_surface_resin_pockets(mesh, laminate, wrinkle_config)
        elif cfg.enable_surface_resin_pockets:
            logger.info(
                "Surface resin pockets superseded by the compaction Vf "
                "gradient (issue #379): the binary trough tag is the "
                "resin-rich extreme of the continuous Vf field, so it is "
                "not applied on top of it."
            )

        # 4a4. Compaction-driven Vf / ply-thickness gradient (issue #379).
        # Runs after the crest lens so the lens blend composes on top of the
        # locally-compacted host material.
        if cfg.enable_vf_gradient:
            self._attach_vf_gradient(mesh, laminate)

        results.mesh = mesh

        # 4b. Mesh-based max fiber angle (accounts for decay mode)
        results.mesh_max_angle_rad = (
            float(np.max(mesh.fiber_angles))
            if mesh.fiber_angles.size > 0 else 0.0
        )

        _report("Solving FE system", 0.25)

        # 5. FE solve — branch on CZM mode.
        if cfg.enable_czm:
            self._run_czm_path(results, laminate, mesh, wrinkle_config)
            _report("Analysis complete", 1.0)
            logger.info(
                "Analysis complete (CZM path): knockdown=%s",
                results.analytical_knockdown,
            )
            return results

        # Progressive-damage path: carries the solve to ultimate load and
        # reports a strength knockdown (the only path that knocks down
        # pristine UD compression).  Still runs the linear solve below so
        # the usual field/failure/retention outputs remain populated.
        if cfg.enable_progressive_damage:
            self._run_progressive_path(results, laminate, mesh)

        # Linear (legacy) path.  ``delta_T`` adds the thermal
        # initial-strain load vector (issue #273 Stage 2).
        solver = StaticSolver(
            mesh, laminate, delta_T=cfg.delta_T,
            **_iterative_solver_kwargs(cfg)
        )
        # A general load state (issue #275) replaces the uniaxial
        # displacement BCs with the traction set its resultants define.
        bcs = _mechanical_bcs(cfg, mesh)
        field_results = solver.solve(
            bcs, solver=cfg.solver, verbose=cfg.verbose
        )
        results.field_results = field_results

        _report("Evaluating failure criteria", 0.75)

        # 6. Failure evaluation on FE field
        self._evaluate_failure(results, laminate, field_results, mesh)

        _report("Computing retention factors", 0.90)

        # 6b. Retention factor (baseline pristine comparison)
        self._compute_retention_factors(results, laminate)

        _report("Analysis complete", 1.0)
        logger.info(
            "Analysis complete: knockdown=%s",
            results.analytical_knockdown,
        )

        return results

    # ------------------------------------------------------------------
    # Morphology comparison
    # ------------------------------------------------------------------

    @staticmethod
    def compare_morphologies(
        base_config: AnalysisConfig,
        morphologies: Sequence[str] = ("stack", "convex", "concave"),
        analytical_only: bool = False,
    ) -> dict[str, AnalysisResults]:
        """Run the full FE analysis for multiple morphologies and compare.

        Parameters
        ----------
        base_config : AnalysisConfig
            Base configuration.  The ``morphology`` field is overridden
            for each entry in *morphologies*.
        morphologies : sequence of str
            Morphology names to compare.
        analytical_only : bool, optional
            If True, skip the FE assembly path for each morphology and
            return only the analytical predictions.  Default ``False``.

        Returns
        -------
        dict[str, AnalysisResults]
            Mapping from morphology name to its results.
        """
        all_results: dict[str, AnalysisResults] = {}

        for morph in morphologies:
            cfg = replace(base_config, morphology=morph)
            all_results[morph] = WrinkleAnalysis(cfg).run(
                analytical_only=analytical_only
            )

        return all_results

    # ------------------------------------------------------------------
    # Parametric sweep
    # ------------------------------------------------------------------

    @staticmethod
    def parametric_sweep(
        base_config: AnalysisConfig,
        parameter: str,
        values: Sequence[float],
        analytical_only: bool = False,
        n_workers: int = 1,
    ) -> list[AnalysisResults]:
        """Sweep a single parameter over a range of values.

        Parameters
        ----------
        base_config : AnalysisConfig
            Base configuration to clone for each value.
        parameter : str
            Name of the parameter to vary.  Must be a numeric field of
            :class:`AnalysisConfig` (e.g. ``'amplitude'``, ``'wavelength'``,
            ``'width'``, ``'applied_strain'``).
        values : sequence of float
            Parameter values to evaluate.
        analytical_only : bool, optional
            If True, skip the FE assembly path for each sweep value and
            return only the analytical predictions.  Default ``False``.
        n_workers : int, optional
            Number of worker processes (issue #260).  ``1`` (default)
            keeps the sequential in-process path; ``0`` uses all CPU
            cores; ``> 1`` runs the independent per-value analyses on a
            ``ProcessPoolExecutor``.  Results come back in the order of
            *values* either way.  Peak memory scales with ``n_workers``
            x the per-solve footprint (each worker returns a full
            :class:`AnalysisResults`, including mesh and field arrays on
            the FE path) — size the worker count by available RAM for
            fine meshes.

        Returns
        -------
        list[AnalysisResults]
            One result per value, in the same order as *values*.

        Raises
        ------
        AttributeError
            If *parameter* is not a valid :class:`AnalysisConfig` field.
        """
        valid_field_names = {f.name for f in fields(base_config)}
        if parameter not in valid_field_names:
            raise AttributeError(
                f"AnalysisConfig has no field '{parameter}'"
            )
        n_workers = _resolve_sweep_workers(n_workers)

        # ``_replace_swept`` owns the ``domain_length`` sentinel reset
        # (auto-derived only — an explicitly pinned domain is kept), so
        # the sweep and the goal-seek clone configs identically.
        configs: list[AnalysisConfig] = [
            _replace_swept(base_config, parameter, val) for val in values
        ]

        if n_workers == 1:
            return [
                _sweep_run_one(cfg, analytical_only) for cfg in configs
            ]

        # Parallel path: each value is an independent analysis, so fan
        # the solves out over processes.  ``executor.map`` preserves the
        # submission order, so the returned list lines up with *values*
        # exactly like the sequential path.
        executor = ProcessPoolExecutor(max_workers=n_workers)
        try:
            results_list = list(
                executor.map(
                    _sweep_run_one,
                    configs,
                    itertools.repeat(analytical_only),
                )
            )
        except BaseException:
            # KeyboardInterrupt (or a worker failure): cancel queued
            # futures instead of letting the pool drain.
            executor.shutdown(wait=False, cancel_futures=True)
            raise
        else:
            executor.shutdown(wait=True)
        return results_list

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_laminate(self) -> Laminate:
        """Build the Laminate from config."""
        # material / angles are filled in __post_init__.
        assert self.config.angles is not None
        assert self.config.material is not None
        return Laminate.from_angles(
            self.config.angles,
            self.config.material,
            ply_thickness=self.config.ply_thickness,
        )

    # ------------------------------------------------------------------
    # Cohesive-zone modelling (Phase 3 wiring)
    # ------------------------------------------------------------------

    def _resolve_cohesive_interfaces(
        self, laminate: Laminate, wrinkle_config: WrinkleConfiguration,
    ) -> list[int]:
        """Resolve ``cfg.czm_interfaces`` to an explicit ply-interface list.

        Parameters
        ----------
        laminate : Laminate
            The fully built laminate; supplies the ply z-coordinates.
        wrinkle_config : WrinkleConfiguration
            Used to locate the wrinkle peaks when
            ``cfg.czm_interfaces == "near_crest"``: scalar configs pick
            the single interface nearest the (largest-amplitude)
            wrinkle; multi-wrinkle configs pick the interface nearest
            *each* wrinkle, deduplicated (issue #283).

        Returns
        -------
        list[int]
            Sorted, de-duplicated list of ply-interface indices in
            ``[0, n_plies - 2]`` (0 = bottom-most interior interface).
        """
        cfg = self.config
        n_plies = laminate.n_plies
        if isinstance(cfg.czm_interfaces, list):
            return sorted({int(i) for i in cfg.czm_interfaces})
        if cfg.czm_interfaces == "all":
            return list(range(n_plies - 1))
        if cfg.czm_interfaces == "near_crest":
            # Pick the *interior* interface closest to the wrinkle peak.
            # The wrinkle peak amplitude in the un-deformed (flat) mesh
            # lies at the wrinkle's reference centre, i.e. at the
            # nominal z of its ply-interface index, shifted in z by the
            # wrinkle's peak displacement.  For the dual-wrinkle modes
            # the two interfaces straddle the laminate midplane; for
            # single-wrinkle modes there is one wrinkle interface.
            ply_z = laminate.z_coords()  # length n_plies + 1
            # Note: interface_z[i] is the midpoint between plies i+1 and
            # i+2; we want a z value associated with the boundary
            # between plies i and i+1, i.e. ply_z[i+1].  Use that
            # directly.
            boundary_z = ply_z[1:n_plies]  # internal boundaries
            # The wrinkle's reference centre z is the boundary z of its
            # ply_interface.  Take the wrinkle with the largest amplitude
            # (most-likely-to-delaminate driver) and pick the boundary
            # closest to that wrinkle's centre.
            wrinkles = list(getattr(wrinkle_config, "wrinkles", ()))

            def _nearest_boundary(placement) -> int:
                # The wrinkle's z reference is the ply boundary above ply
                # ``ply_interface``; ``ply_interface`` is in the [0,
                # n_plies-2] range and references the boundary at
                # ply_z[ply_interface + 1].
                k = int(placement.ply_interface)
                if 1 <= k + 1 <= n_plies - 1:
                    target_z = float(ply_z[k + 1])
                else:
                    target_z = 0.0
                return int(np.argmin(np.abs(boundary_z - target_z)))

            if not wrinkles:
                # No wrinkle placements (flat mesh) — default to the
                # interior boundary nearest the midplane.
                return [int(np.argmin(np.abs(boundary_z)))]
            if cfg.wrinkles is not None:
                # Multi-wrinkle configuration (issue #283): nominate the
                # interface nearest *each* wrinkle so a delamination can
                # initiate at any crest and, where wrinkles share an
                # interface index, run along one continuous cohesive
                # surface between them (crest-to-crest link-up).
                return sorted({_nearest_boundary(w) for w in wrinkles})
            # Scalar (named-morphology) configuration: keep the legacy
            # single-interface choice — the wrinkle with the largest
            # amplitude is the most-likely-to-delaminate driver.
            w_best = max(
                wrinkles,
                key=lambda w: abs(getattr(w.profile, "amplitude", 0.0)),
            )
            return [_nearest_boundary(w_best)]
        # Validation in ``__post_init__`` already constrains the values
        # this method sees; the catch-all keeps mypy happy.
        raise ValueError(
            f"Unrecognised czm_interfaces value: {cfg.czm_interfaces!r}"
        )

    def _build_cohesive_properties(
        self, laminate: Laminate,
    ) -> CohesiveProperties:
        """Build ``CohesiveProperties`` from config-or-material defaults.

        Falls back to ply-0's material when a CZM strength / toughness
        is not explicitly set on the config.  Raises ``ValueError`` when
        the laminate's first ply lacks ``GIc`` / ``GIIc`` and the user
        did not override them on the config.
        """
        cfg = self.config
        # Reference material: ply 0 (simpler than per-interface lookup,
        # matches the v1 spec).
        mat = laminate.plies[0].material

        def _coalesce(cfg_val: float | None, mat_val) -> float:
            if cfg_val is not None:
                return float(cfg_val)
            if mat_val is None:
                raise ValueError(
                    "Cohesive zone modelling requires GIc, GIIc, "
                    "sigma_max and tau_max either as explicit "
                    "AnalysisConfig.czm_* overrides or as material "
                    "defaults.  Material "
                    f"{mat.name!r} (ply 0) is missing one of them."
                )
            return float(mat_val)

        return CohesiveProperties(
            K=float(cfg.czm_penalty),
            sigma_max=_coalesce(cfg.czm_sigma_max, mat.sigma_max),
            tau_max=_coalesce(cfg.czm_tau_max, mat.tau_max),
            GIc=_coalesce(cfg.czm_GIc, mat.GIc),
            GIIc=_coalesce(cfg.czm_GIIc, mat.GIIc),
            eta_BK=float(cfg.czm_BK_eta),
            beta=1.0,
        )

    def _build_mesh_with_cohesive_interfaces(
        self,
        laminate: Laminate,
        wrinkle_config: WrinkleConfiguration,
        iface_indices: list[int],
        cohesive_props: CohesiveProperties,
    ) -> tuple[MeshData, list[tuple[int, Cohesive8Element]], dict[int, range]]:
        """Build a wrinkled mesh with cohesive elements inserted.

        The flat (un-deformed) mesh is built first so that ply-interface
        nodes lie on exact z-planes; cohesive elements are inserted on
        those flat planes; finally the wrinkle displacement field is
        applied to **all** nodes (including the duplicated interface
        nodes), so the cohesive layer sits at a curved surface in the
        deformed configuration.  This keeps the
        :func:`insert_cohesive_interface` topology check valid (which
        requires axis-aligned interface planes in the reference
        configuration) while letting the rest of the pipeline see the
        normal wrinkled geometry.

        Returns
        -------
        mesh : MeshData
            Final wrinkled mesh with duplicated interface nodes.
        cohesive_elements : list[tuple[int, Cohesive8Element]]
            Entries suitable for :class:`GlobalAssembler`.  The integer
            key is the *global* cohesive-element id, unique across
            interfaces.  Element ``node_coords`` are refreshed against
            the wrinkled ``mesh.nodes`` so the assembler's exact-equality
            check passes.
        elem_ranges : dict[int, range]
            Map from ply-interface index to the contiguous range of
            global cohesive-element ids belonging to that interface,
            used downstream for per-interface aggregation.
        """
        cfg = self.config

        # --- Step 1: flat mesh, no wrinkle deformation ---------------
        flat_mesh_gen = WrinkleMesh(
            laminate=laminate,
            wrinkle_config=None,
            Lx=cfg.domain_length,
            Ly=cfg.domain_width,
            nx=cfg.nx,
            ny=cfg.ny,
            nz_per_ply=cfg.nz_per_ply,
        )
        mesh = flat_mesh_gen.generate()

        # --- Step 2: insert cohesive layers at requested z-planes ----
        ply_z = laminate.z_coords()  # length n_plies + 1
        # ply-interface index ``i`` corresponds to the boundary between
        # plies ``i`` and ``i + 1`` => z = ply_z[i + 1].
        elem_ranges: dict[int, range] = {}
        cohesive_elements: list[tuple[int, Cohesive8Element]] = []
        next_global_id = 0
        for iface_idx in iface_indices:
            z_iface = float(ply_z[iface_idx + 1])
            mesh, coh_elems = insert_cohesive_interface(
                mesh, z_iface, cohesive_props,
            )
            start = next_global_id
            for k, c_elem in enumerate(coh_elems):
                # Reassign elem_id to be unique across interfaces.
                c_elem.elem_id = next_global_id
                cohesive_elements.append((next_global_id, c_elem))
                next_global_id += 1
            elem_ranges[iface_idx] = range(start, next_global_id)

        # --- Step 3: apply wrinkle deformation to the expanded mesh -
        # Rebuild the node_ply_ids array for the new (duplicated) node
        # set.  WrinkleMesh._node_to_element_ply assigns ply ids by the
        # z-layer index; for duplicated nodes the natural choice is to
        # carry the original ply id of the node they were duplicated
        # from.  The cleanest way to derive this without re-deriving
        # internal mesh state is to assign each node by the z-value of
        # its original (flat) coordinate, falling back to the original
        # ply id for nodes that exactly straddle a boundary.
        node_ply_ids = self._derive_node_ply_ids(mesh, laminate)
        deformed_nodes = wrinkle_config.apply_to_nodes(
            mesh.nodes, node_ply_ids, laminate.n_plies,
        )
        fiber_angles = wrinkle_config.fiber_angles_at_nodes(
            mesh.nodes, node_ply_ids, n_plies=laminate.n_plies,
        )
        mesh = MeshData(
            nodes=deformed_nodes,
            elements=mesh.elements,
            ply_ids=mesh.ply_ids,
            fiber_angles=fiber_angles,
            ply_angles=mesh.ply_angles,
            nx=mesh.nx,
            ny=mesh.ny,
            nz=mesh.nz,
            laminate=laminate,
        )

        # --- Step 4: refresh cohesive node_coords against the deformed
        # node array so the GlobalAssembler equality check passes -----
        refreshed: list[tuple[int, Cohesive8Element]] = []
        for gid, c_elem in cohesive_elements:
            new_coords = mesh.nodes[c_elem.node_ids]
            new_elem = Cohesive8Element(
                node_coords=new_coords,
                properties=c_elem.properties,
                node_ids=c_elem.node_ids,
                elem_id=c_elem.elem_id,
            )
            refreshed.append((gid, new_elem))
        return mesh, refreshed, elem_ranges

    @staticmethod
    def _derive_node_ply_ids(mesh: MeshData, laminate: Laminate) -> np.ndarray:
        """Assign a ply id to every mesh node from its z-coordinate.

        Used after :func:`insert_cohesive_interface` has duplicated
        interface nodes: those duplicates need the *same* ply id as the
        originals so the wrinkle through-thickness decay treats them as
        belonging to the same ply.  A node sitting exactly on a ply
        boundary z is assigned to the ply *below* (lower index); the
        wrinkle decay field is continuous across boundaries so this
        choice has no physical consequence.
        """
        ply_z = laminate.z_coords()  # n_plies + 1
        z = mesh.nodes[:, 2]
        # ``np.searchsorted(side='right') - 1`` puts a node sitting
        # exactly on a boundary into the ply below the boundary.
        ids = np.searchsorted(ply_z, z, side="right") - 1
        ids = np.clip(ids, 0, laminate.n_plies - 1)
        return ids.astype(np.intp)

    def _run_czm_path(
        self,
        results: AnalysisResults,
        laminate: Laminate,
        flat_mesh: MeshData,
        wrinkle_config: WrinkleConfiguration,
    ) -> None:
        """End-to-end CZM solve: insert interfaces, run Newton, populate results.

        The ``flat_mesh`` argument is the wrinkled hex8 mesh that the
        linear path would use; we rebuild a *different* mesh here that
        also carries cohesive layers (see
        :meth:`_build_mesh_with_cohesive_interfaces`) and use that for
        the solve.  The mesh passed in is otherwise unused — it is the
        sentinel value that the linear-path callers (``_evaluate_failure``,
        ``_compute_retention_factors``) consume; we overwrite
        ``results.mesh`` with our enriched mesh below.
        """
        cfg = self.config

        iface_indices = self._resolve_cohesive_interfaces(
            laminate, wrinkle_config,
        )
        results.czm_interfaces_used = list(iface_indices)

        cohesive_props = self._build_cohesive_properties(laminate)
        mesh, cohesive_elements, elem_ranges = (
            self._build_mesh_with_cohesive_interfaces(
                laminate, wrinkle_config, iface_indices, cohesive_props,
            )
        )
        results.mesh = mesh
        if mesh.fiber_angles.size > 0:
            results.mesh_max_angle_rad = float(np.max(mesh.fiber_angles))

        # Boundary conditions: ``compression_bcs`` is sign-agnostic (the
        # sign of ``applied_strain`` selects compression vs tension), so
        # we use it for both loading modes here, matching the linear
        # path.
        bcs = BoundaryHandler.compression_bcs(
            mesh, applied_strain=cfg.applied_strain
        )
        bc_handler = BoundaryHandler(mesh)
        assembler = GlobalAssembler(
            mesh, laminate, cohesive_elements=cohesive_elements,
            delta_T=cfg.delta_T,
        )

        solver = NewtonRaphsonSolver(
            assembler=assembler,
            bc_handler=bc_handler,
            boundary_conditions=bcs,
            n_increments=int(cfg.czm_n_load_increments),
            tol_residual=float(cfg.czm_newton_tol),
        )
        outcome = solver.solve(verbose=cfg.verbose)

        results.czm_converged = bool(outcome.get("converged", False))
        results.czm_failure_diagnostics = outcome.get("failure_diagnostics")
        results.czm_failure_hint = outcome.get("failure_hint")
        results.czm_load_displacement = outcome.get(
            "load_displacement", None,
        )

        # ----- Bulk hex8 stress / strain recovery from the final Newton u -----
        # The CZM path now populates ``field_results`` so users can run
        # ply-level failure criteria (LaRC05, Hashin, ...) on the bulk
        # material alongside the interface delamination output.  We reuse
        # StaticSolver's stress-recovery machinery: instantiate it on the
        # CZM-enriched mesh (without cohesive_elements registered, so the
        # assemble_stiffness guard does not fire) and call
        # recover_element_results on the converged displacement.
        u_final = outcome["displacement"]
        recovery_solver = StaticSolver(
            mesh, laminate, delta_T=self.config.delta_T,
            **_iterative_solver_kwargs(self.config)
        )
        stress_g, stress_l, strain_g, strain_l = (
            recovery_solver.recover_element_results(u_final, verbose=False)
        )
        results.field_results = FieldResults(
            displacement=u_final.reshape(-1, 3),
            stress_global=stress_g,
            stress_local=stress_l,
            strain_global=strain_g,
            strain_local=strain_l,
            mesh=mesh,
            laminate=laminate,
        )
        self._evaluate_failure(results, laminate, results.field_results, mesh)

        # ----- Extract per-Gauss-point CZM state -----
        n_iface = len(cohesive_elements)
        if n_iface == 0:
            results.czm_damage = np.empty((0, 0))
            results.czm_separation = np.empty((0, 0, 3))
            results.czm_traction = np.empty((0, 0, 3))
            results.czm_energy_dissipated = 0.0
            results.czm_energy_per_interface = {}
            results.czm_crack_length_per_interface = {}
            results.czm_element_centroids = np.empty((0, 2))
            results.czm_delamination_report = build_delamination_report({})
            return

        n_gp = cohesive_elements[0][1].n_gp
        damage = np.zeros((n_iface, n_gp), dtype=float)
        separation = np.zeros((n_iface, n_gp, 3), dtype=float)
        traction = np.zeros((n_iface, n_gp, 3), dtype=float)

        # In-plane (x, y) centroid of each interface element in the
        # reference configuration; consumed by the viz layer to colour
        # damage on the interface plane without touching the assembler.
        centroids_xy = np.zeros((n_iface, 2), dtype=float)
        for row, (_gid, c_elem) in enumerate(cohesive_elements):
            bottom_xy = c_elem.node_coords[:4, :2]
            centroids_xy[row] = bottom_xy.mean(axis=0)
        results.czm_element_centroids = centroids_xy

        u = outcome["displacement"]
        # Iterate in the same order as `cohesive_elements` was built.
        for row, (gid, c_elem) in enumerate(cohesive_elements):
            # Damage from the assembler-committed state.  Fall back to
            # virgin Gauss-point states (d = 0) for elements the
            # assembler never committed; ``_law_local`` requires a real
            # ``CohesiveState`` (a ``None`` entry would crash it).
            state = assembler.cohesive_state.get(
                gid, make_initial_state(c_elem.n_gp)
            )
            for g in range(c_elem.n_gp):
                damage[row, g] = float(state[g].d)

            # Recompute the local separation and traction at each GP
            # from the final displacement.  ``tangent_and_force`` returns
            # the assembled force / tangent but does not expose
            # per-GP values; we recompute the kinematic / constitutive
            # pieces directly here.
            dofs = GlobalAssembler._cohesive_dof_indices(c_elem)
            u_e = u[dofs]
            for g in range(c_elem.n_gp):
                B = c_elem._B_jump(g)
                R = c_elem._R_gp[g]
                delta_global = B @ u_e
                delta_local = R @ delta_global
                separation[row, g] = delta_local
                # Traction via the committed-state law evaluation.
                T_local, _D_local, _new = c_elem._law_local(
                    delta_local, state[g],
                )
                traction[row, g] = T_local

        results.czm_damage = damage
        results.czm_separation = separation
        results.czm_traction = traction

        # ----- Per-interface energy + crack length -----
        # Per-element dissipated energy (approximate):
        # E_e ≈ 0.5 * sigma_max * delta_f * area_e * d_avg, where
        # delta_f is mode-dependent.  We use the simpler bound
        # GIc * area * d_avg which is accurate for mode-I-dominated
        # failure and gives the right total in the fully-failed limit.
        energy_per_iface: dict[int, float] = {}
        crack_len_per_iface: dict[int, float] = {}
        damage_per_iface: dict[int, np.ndarray] = {}

        # Mesh width in y for crack-length estimation.  Use the bounding
        # box of the *original* (flat) coordinates which the cohesive
        # element retains internally; falls back to the deformed mesh.
        y_min = float(mesh.nodes[:, 1].min())
        y_max = float(mesh.nodes[:, 1].max())
        Ly = max(y_max - y_min, 1e-12)

        for iface_idx, gid_range in elem_ranges.items():
            rows = list(range(gid_range.start, gid_range.stop))
            d_iface = damage[rows]
            damage_per_iface[iface_idx] = d_iface

            # Energy: sum per element of GIc * area * d_mean (mode-I
            # bound).  Mode-mixity refines this in principle but the
            # mode-I approximation suffices for the v1 reporter.
            e_iface = 0.0
            crack_area = 0.0
            for row, gid in enumerate(rows):
                _, c_elem = cohesive_elements[gid]
                area_e = c_elem.area
                d_avg = float(d_iface[row].mean())
                e_iface += cohesive_props.GIc * area_e * d_avg
                if d_iface[row].max() > 0.99:
                    crack_area += area_e
            energy_per_iface[iface_idx] = float(e_iface)
            crack_len_per_iface[iface_idx] = float(crack_area / Ly)

        results.czm_energy_per_interface = energy_per_iface
        results.czm_crack_length_per_interface = crack_len_per_iface
        results.czm_energy_dissipated = float(sum(energy_per_iface.values()))
        results.czm_delamination_report = build_delamination_report(
            damage_per_iface,
            energy_per_interface=energy_per_iface,
            crack_length_per_interface=crack_len_per_iface,
        )

    def _analytical_modulus_knockdown(
        self, angles: list[float], morphology_factor: float
    ) -> float:
        """Closed-form axial-modulus knockdown for a wavy laminate.

        Populated for arbitrary layups and multi-wrinkle layouts. The
        unidirectional single-wrinkle case is served by the scalar fast path
        :func:`_profile_modulus_knockdown` (so the pinned UD baselines stay
        bit-identical); every other case is routed through the CLT membrane
        series-average :func:`_laminate_modulus_knockdown`, which reduces
        exactly to the UD result for ``[0]_n``. Returns ``1.0`` only for a
        degenerate (zero-amplitude / zero-wavelength) wrinkle.
        """
        cfg = self.config
        mat = cfg.material
        assert mat is not None  # filled in __post_init__
        graded = cfg.morphology.lower().strip() == "graded"

        # --- Single-wrinkle UD fast path (unchanged, pins F/G baselines) ---
        if cfg.wrinkles is None and _is_unidirectional(angles):
            if cfg.wavelength <= 1e-12 or cfg.amplitude <= 0.0:
                return 1.0
            if cfg.through_thickness_decay_scale is not None:
                decay_scale = float(cfg.through_thickness_decay_scale)
            else:
                decay_scale = max(cfg.wavelength / 2.0, cfg.amplitude)
            return _profile_modulus_knockdown(
                amplitude=cfg.amplitude,
                wavelength=cfg.wavelength,
                width=cfg.width,
                domain_length=cfg.domain_length,
                ply_thickness=cfg.ply_thickness,
                n_plies=len(angles),
                E1=mat.E1, E2=mat.E2, G12=mat.G12, nu12=mat.nu12,
                morphology_factor=morphology_factor,
                through_thickness_decay=graded,
                z_position_fraction=float(cfg.wrinkle_z_position),
                decay_scale=decay_scale,
                decay_floor=float(cfg.decay_floor),
            )

        # --- Generalized laminate / multi-wrinkle path -------------------
        # Build the longitudinal slope field(s) and the per-ply
        # through-thickness decay, then hand to the CLT membrane
        # series-average. The composition mirrors the FE's
        # "compose then differentiate" multi-wrinkle field.
        n_plies = len(angles)
        if n_plies == 0:
            return 1.0
        ply_thickness = cfg.ply_thickness
        total_thickness = n_plies * ply_thickness
        decay_floor = float(cfg.decay_floor)
        x = np.linspace(
            -cfg.domain_length / 2.0, cfg.domain_length / 2.0, _N_PROFILE_PTS
        )

        # Collect (amplitude, wavelength, width, z_center, decay_scale,
        # phase_shift_x) for every wrinkle. The single (non-UD) wrinkle and
        # the multi-wrinkle list are normalized to the same representation.
        specs: list[tuple[float, float, float, float, float, float]] = []
        if cfg.wrinkles is not None:
            for spec in cfg.wrinkles:
                if spec.wavelength <= 1e-12 or spec.amplitude <= 0.0:
                    continue
                # Through-thickness decay centred on this wrinkle's interface.
                z_center = (spec.ply_interface + 1.0) * ply_thickness
                z_center = min(max(z_center, 0.0), total_thickness)
                ds = (
                    float(cfg.through_thickness_decay_scale)
                    if cfg.through_thickness_decay_scale is not None
                    else max(spec.wavelength / 2.0, spec.amplitude)
                )
                dx_shift = spec.phase_offset * spec.wavelength / (2.0 * np.pi)
                specs.append(
                    (spec.amplitude, spec.wavelength, spec.width,
                     z_center, ds, dx_shift)
                )
        else:
            if cfg.wavelength <= 1e-12 or cfg.amplitude <= 0.0:
                return 1.0
            z_center = float(cfg.wrinkle_z_position) * total_thickness
            ds = (
                float(cfg.through_thickness_decay_scale)
                if cfg.through_thickness_decay_scale is not None
                else max(cfg.wavelength / 2.0, cfg.amplitude)
            )
            specs.append(
                (cfg.amplitude, cfg.wavelength, cfg.width, z_center, ds, 0.0)
            )

        if not specs:
            return 1.0

        n_w = len(specs)
        # Per-wrinkle slope along x and per-(ply, wrinkle) decay.
        slope_field = np.empty((n_w, _N_PROFILE_PTS), dtype=float)
        ply_decays = np.empty((n_plies, n_w, _N_PROFILE_PTS), dtype=float)
        use_decay = graded or cfg.wrinkles is not None
        for w, (amp, lam, wid, z_center, ds, dx_shift) in enumerate(specs):
            dxw = x - dx_shift
            gauss_env = np.exp(-(dxw ** 2) / (wid ** 2))
            k = 2.0 * np.pi / lam
            slope = amp * gauss_env * (
                (-2.0 * dxw / (wid ** 2)) * np.cos(k * dxw)
                - k * np.sin(k * dxw)
            )
            slope_field[w] = slope * morphology_factor
            sigma_sq2 = 2.0 * ds * ds
            for p in range(n_plies):
                z_p = (p + 0.5) * ply_thickness
                if use_decay:
                    raw = math.exp(-((z_p - z_center) ** 2) / sigma_sq2)
                    ply_decays[p, w, :] = decay_floor + (1.0 - decay_floor) * raw
                else:
                    ply_decays[p, w, :] = 1.0

        laminate = Laminate.from_angles(list(angles), mat, ply_thickness)
        return _laminate_modulus_knockdown(
            slope_field=slope_field,
            ply_decays=ply_decays,
            angles=angles,
            stiffness_3d=mat.stiffness_matrix,
            ply_thickness=ply_thickness,
            E_x0=laminate.Ex,
        )

    def _wrap_transverse_surface(
        self, profile: WrinkleProfile
    ) -> WrinkleSurface3D:
        """Wrap *profile* in a through-width :class:`WrinkleSurface3D` (#300).

        Resolves the transverse span and localization half-width from the
        config, defaulting ``span_y`` to the meshed ``domain_width`` and
        ``width_y`` to ``span_y / 4`` (a localized mid-width patch). Only
        called for a non-uniform ``transverse_mode``; validation has already
        rejected the analytical-only, multi-wrinkle, and CZM combinations.
        """
        cfg = self.config
        span_y = (
            cfg.transverse_span
            if cfg.transverse_span is not None
            else cfg.domain_width
        )
        width_y = (
            cfg.transverse_width
            if cfg.transverse_width is not None
            else span_y / 4.0
        )
        return WrinkleSurface3D(
            profile,
            transverse_mode=cfg.transverse_mode,
            width_y=width_y,
            span_y=span_y,
        )

    def _compute_analytical(
        self,
        results: AnalysisResults,
        wrinkle_config: WrinkleConfiguration,
    ) -> None:
        """Fill analytical prediction fields in results.

        Uses the unattenuated sinusoidal angle theta = arctan(2*pi*A/lambda)
        rather than the Gaussian-envelope max angle, because D/T-based
        experimental data references the full wrinkle amplitude A.
        """
        cfg = self.config

        mf = wrinkle_config.aggregate_morphology_factor(cfg.loading)

        # Unattenuated sinusoidal angle: theta = arctan(2*pi*A/lambda)
        # This is the correct angle for D/T-based knockdown comparison,
        # since D/T uses the full amplitude A without Gaussian attenuation.
        if cfg.wrinkles is not None:
            # Multi-wrinkle analytical model: peak-angle over all wrinkles (initial implementation).
            # Each wrinkle gets its own theta_max,i = arctan(2*pi*A_i/lambda_i)
            # and we take the maximum, then scale by the aggregate morphology
            # factor. This is intentionally coarse; calibration against the
            # Li et al. (2025) Dataset F multi-wrinkle specimens is a
            # follow-up activity (D-AB-2, D-A-2, D-M-2, T-M-2).
            theta_max = 0.0
            for spec in cfg.wrinkles:
                if spec.wavelength > 1e-12:
                    theta_i = float(
                        np.arctan(2.0 * np.pi * spec.amplitude / spec.wavelength)
                    )
                    if theta_i > theta_max:
                        theta_max = theta_i
        elif cfg.wavelength > 1e-12:
            theta_max = float(np.arctan(2.0 * np.pi * cfg.amplitude / cfg.wavelength))
        else:
            theta_max = 0.0
        theta_eff = theta_max * mf

        # Analytical stiffness (axial-modulus) knockdown — generalized to
        # arbitrary layups and multi-wrinkle layouts, loading-independent,
        # set on both the gate and Budiansky-Fleck paths below. Resolve the
        # layup the same way each path does.
        if cfg.penetration_gate is not None:
            _mod_angles = cfg.angles if cfg.angles else [0.0]
        else:
            _mod_angles = (
                cfg.angles if cfg.angles else [0.0, 45.0, -45.0, 90.0] * 6
            )
        results.analytical_modulus_knockdown = (
            self._analytical_modulus_knockdown(_mod_angles, mf)
        )

        # Penetration-gate path (item D.3): when a calibrated gate is
        # configured, the analytical knockdown is the two-parameter
        # (theta, D/T) gate value instead of Budiansky-Fleck.  Uses the
        # peak angle and the penetration D/T = A / (t_ply * n_plies), both
        # on the section-2.7 conventions (config amplitude is the
        # half-amplitude).  Returns early — the BF / tension blocks below
        # are bypassed.
        if cfg.penetration_gate is not None:
            angles_g = cfg.angles if cfg.angles else [0]
            T = cfg.ply_thickness * len(angles_g)
            if cfg.wrinkles is not None:
                # Multi-wrinkle gate (issue #342): evaluate the gate per
                # spec — each wrinkle carries its own angle
                # theta_i = arctan(2*pi*A_i/lambda_i), penetration
                # D_i/T = A_i/T, and through-thickness position
                # z_i = (ply_interface + 1) / n_plies (the boundary the
                # spec nominates; ``cfg.wrinkle_z_position`` is a
                # scalar-path parameter and is ignored here) — and take
                # the weakest-link (minimum) knockdown.  Previously
                # ``dt`` silently read the leftover scalar
                # ``cfg.amplitude`` while theta came from the specs,
                # producing a plausible-looking wrong answer (e.g. kd
                # 0.98 instead of 0.64 on the issue's repro).
                n_plies_g = len(angles_g)
                kd_gate = 1.0
                for spec in cfg.wrinkles:
                    if spec.wavelength > 1e-12:
                        theta_i = float(np.arctan(
                            2.0 * np.pi * spec.amplitude / spec.wavelength
                        ))
                    else:
                        theta_i = 0.0
                    dt_i = (spec.amplitude / T) if T > 0 else 0.0
                    z_i = (spec.ply_interface + 1) / n_plies_g
                    kd_i = penetration_gate_kd(
                        math.degrees(theta_i), dt_i, cfg.penetration_gate,
                        z_position=float(z_i),
                    )
                    kd_gate = min(kd_gate, float(kd_i))
            else:
                dt = (cfg.amplitude / T) if T > 0 else 0.0
                kd_gate = penetration_gate_kd(
                    math.degrees(theta_max), dt, cfg.penetration_gate,
                    z_position=float(cfg.wrinkle_z_position),
                )
            results.morphology_factor = mf
            results.max_angle_rad = theta_max
            results.effective_angle_rad = theta_eff
            results.gamma_Y_eff = cfg.penetration_gate.gamma_Y
            results.analytical_knockdown = float(kd_gate)
            material = cfg.material
            assert material is not None  # set by AnalysisConfig.__post_init__
            ref = (material.Xt if cfg.loading == "tension"
                   else material.Xc)
            results.analytical_strength_MPa = float(ref) * float(kd_gate)
            return

        # Compute layup-dependent effective gamma_Y from confinement
        angles: list[float] = (
            cfg.angles if cfg.angles else [0.0, 45.0, -45.0, 90.0] * 6
        )
        gamma_Y_eff = _effective_gamma_Y(angles)

        # Compression KD (CLT-weighted Budiansky-Fleck) — computed for
        # both loading modes: used directly for compression, and as a
        # physical floor for tension (tension cannot be worse than compression).
        mat = cfg.material
        assert mat is not None  # filled in __post_init__
        E11 = mat.E1
        E22 = mat.E2
        G12 = mat.G12

        n_0 = sum(1 for a in angles if abs(a) < 5)
        n_45 = sum(1 for a in angles if 40 < abs(a) < 50)
        n_90 = sum(1 for a in angles if abs(a) > 85)

        Q11_0 = E11
        Q11_45 = E11 / 4.0 + E22 / 4.0 + G12 / 2.0
        Q11_90 = E22

        total_stiffness = n_0 * Q11_0 + n_45 * Q11_45 + n_90 * Q11_90
        f_0 = n_0 * Q11_0 / total_stiffness if total_stiffness > 0 else 1.0

        # For graded morphology (embedded wrinkle), use profile-proportional
        # knockdown: the BF knockdown is averaged over the wrinkle profile
        # in both x (local angle varies with dz/dx) and z (Gaussian decay).
        # For other morphologies (stack/convex/concave/uniform), the wrinkle
        # extends through the full thickness and fills the coupon, so failure
        # is governed by the peak-angle cross-section.
        #
        # KD_lam = (1/N) sum_p [ (1/L_s) int 1/(1 + theta(x)*Phi(z_p)/gY) dx ]
        # Force the non-graded peak-angle path when a multi-wrinkle
        # override is active: the graded profile-proportional helper
        # below uses cfg.amplitude / cfg.wavelength / cfg.width directly
        # and would silently ignore the per-spec geometry.
        is_graded = (
            cfg.morphology.lower().strip() == "graded"
            and cfg.wrinkles is None
        )
        n_plies = len(angles)

        # Resolve through-thickness Gaussian decay scale (mm).  The user
        # can override the auto formula via
        # ``cfg.through_thickness_decay_scale``; default is
        # ``max(wavelength / 2, amplitude)``.  Used for both the
        # compression profile-proportional KD and the tension graded-
        # averaging block.
        if cfg.through_thickness_decay_scale is not None:
            decay_scale_eff = float(cfg.through_thickness_decay_scale)
        else:
            decay_scale_eff = max(cfg.wavelength / 2.0, cfg.amplitude)
        c_AF = float(cfg.kink_band_quadratic_coeff)

        if is_graded and n_0 > 0 and n_plies > 1:
            # Profile-proportional compression knockdown (graded/embedded).
            # ``wrinkle_z_position`` shifts the through-thickness decay
            # centre off the midplane to model wrinkles closer to the
            # surface (Li et al. 2025 Dataset F: Above/Below positions).
            z_pos = float(cfg.wrinkle_z_position)
            kd_profile = _profile_proportional_kd(
                amplitude=cfg.amplitude,
                wavelength=cfg.wavelength,
                width=cfg.width,
                domain_length=cfg.domain_length,
                ply_thickness=cfg.ply_thickness,
                n_plies=n_plies,
                gamma_Y=gamma_Y_eff,
                theta_max=theta_max,
                morphology_factor=1.0,
                through_thickness_decay=True,
                z_position_fraction=z_pos,
                decay_scale=decay_scale_eff,
                decay_floor=cfg.decay_floor,
                kink_band_quadratic_coeff=c_AF,
            )
            kd_compression = f_0 * kd_profile + (1.0 - f_0)

            if cfg.loading == "tension":
                # Average tension knockdown over 0-deg plies at local angles
                # (profile-proportional, using stretched linear grading
                # on the new ``decay_scale_eff`` so the support tracks
                # the wrinkle's longitudinal extent rather than the
                # full half-thickness).  The grading centre ``p_mid`` is
                # shifted by ``wrinkle_z_position`` so the per-ply taper
                # peaks at the user-set through-thickness position
                # rather than the midplane (legacy: midplane).
                p_mid = z_pos * (n_plies - 1)
                zero_positions = [i for i, a in enumerate(angles) if abs(a) < 5]
                decay_floor = cfg.decay_floor
                # Stretched linear support: width = decay_scale / t_ply
                # plies.  Falls back to the legacy half-thickness norm if
                # the decay scale would exceed it (keeps backwards-
                # compatible behaviour for cases where the auto formula
                # is at least the legacy support).
                t_ply = cfg.ply_thickness
                p_support_decay = decay_scale_eff / max(t_ply, 1e-12)
                p_support_legacy = (n_plies - 1) / 2.0
                p_norm = max(min(p_support_decay, p_support_legacy), 1e-12)
                kd_0_sum = 0.0
                for p in zero_positions:
                    raw = max(0.0, 1.0 - abs(p - p_mid) / p_norm)
                    B_p = decay_floor + (1.0 - decay_floor) * raw
                    theta_p = theta_max * B_p
                    kd_p, _ = self._tension_knockdown_analytical(
                        theta_p, cfg, _return_kd0_only=True,
                    )
                    kd_0_sum += kd_p
                kd_0_avg = kd_0_sum / len(zero_positions)
                kd = f_0 * kd_0_avg + (1.0 - f_0)
                # Get mechanisms at peak angle for reporting
                _, mechanisms = self._tension_knockdown_analytical(theta_max, cfg)
                mechanisms["mode"] = mechanisms["mode"] + " (graded avg)"
                if kd < kd_compression:
                    kd = kd_compression
                    mechanisms["mode"] = mechanisms["mode"] + " (capped)"
                ref_strength = mat.Xt
                results.tension_mechanisms = mechanisms
            else:
                kd = kd_compression
                ref_strength = mat.Xc
                results.tension_mechanisms = None
        else:
            # Non-graded: peak-angle BF at the critical cross-section,
            # with the Argon-Fleck quadratic extension.
            r_bf = theta_eff / gamma_Y_eff
            kd_bf = 1.0 / (1.0 + r_bf + c_AF * r_bf * r_bf)
            kd_compression = f_0 * kd_bf + (1.0 - f_0)

            if cfg.loading == "tension":
                kd, mechanisms = self._tension_knockdown_analytical(theta_max, cfg)
                if kd < kd_compression:
                    kd = kd_compression
                    mechanisms["mode"] = mechanisms["mode"] + " (capped)"
                ref_strength = mat.Xt
                results.tension_mechanisms = mechanisms
            else:
                kd = kd_compression
                ref_strength = mat.Xc
                results.tension_mechanisms = None
        results.gamma_Y_eff = gamma_Y_eff

        # Damage index (for reporting; not used in knockdown computation)
        D = (_D0
             * (cfg.amplitude / _A_REF) ** 1.5
             * (1.0 + _BETA_ANGLE * max(theta_max - _THETA_CRIT, 0.0))
             * mf)
        D = min(D, 0.999)

        results.morphology_factor = mf
        results.max_angle_rad = theta_max
        results.effective_angle_rad = theta_eff
        results.damage_index = D
        results.analytical_knockdown = kd

        # Populate the delamination-onset KD from the tension mechanisms
        # dict (None for compression and for materials lacking GIc/GIIc).
        if (
            cfg.loading == "tension"
            and results.tension_mechanisms is not None
            and mat.GIc is not None
            and mat.GIIc is not None
        ):
            onset_val = results.tension_mechanisms.get("kd_onset")
            results.analytical_onset_knockdown = (
                float(onset_val) if onset_val is not None else None
            )
        else:
            results.analytical_onset_knockdown = None

        results.analytical_strength_MPa = ref_strength * kd

    # ------------------------------------------------------------------
    # Tension analytical model — three-mechanism knockdown
    # ------------------------------------------------------------------

    @staticmethod
    def _tension_knockdown_analytical(
        theta: float, cfg: AnalysisConfig,
        _return_kd0_only: bool = False,
    ) -> tuple[float, dict]:
        """Three-mechanism tension knockdown (LaRC04 + curved-beam OOP).

        Computes the laminate-level retention factor for tension loading
        by combining three competing failure mechanisms for the 0-degree
        plies, weighted by CLT axial stiffness fractions:

        1. **Fiber tension** (LaRC04 #3, Pinho Eq. 82): KD = cos²θ
        2. **Matrix tension** (LaRC04 #1, Pinho Eq. 40): Hashin σ₂₂/τ₁₂
           interaction with in-situ strengths (Yt_is, S12_is)
        3. **Out-of-plane delamination** (Timoshenko curved-beam):
           Combined σ₃₃ (mode I at crest) and τ₁₃ (mode II at inflection)

        The 0-degree ply knockdown is min(KD_fiber, KD_matrix, KD_oop).
        Off-axis plies (±45, 90) are assumed unaffected by the waviness
        for tension loading. The laminate knockdown is the CLT-weighted
        average: KD_lam = f_0 × KD_0 + (1 − f_0) × 1.0.

        References
        ----------
        - Pinho et al. (2005) NASA-TM-2005-213530, Eq. 40, 47, 57, 82
        - Timoshenko & Gere (1961), Theory of Elastic Stability (curved beam)
        - Elhajjar (2025) Scientific Reports 15:25977 (experimental data)
        """
        # Filled in __post_init__.  (The previous ``if mat is None:
        # return 1.0`` guard returned a bare float from a function whose
        # callers always unpack a (kd, mechanisms) tuple, so it could
        # never have worked anyway.)
        mat = cfg.material
        assert mat is not None
        angles: list[float] = (
            cfg.angles if cfg.angles else [0.0, 45.0, -45.0, 90.0] * 6
        )

        E11 = mat.E1
        E22 = mat.E2
        G12 = mat.G12
        Xt = mat.Xt
        Yt = mat.Yt if mat.Yt else 49.0
        S12 = mat.S12 if mat.S12 else 85.0
        S13 = mat.S13 if hasattr(mat, "S13") and mat.S13 else S12

        # --- Mechanism 1: Fiber tension  cos²θ ---
        kd_fiber = float(np.cos(theta) ** 2)

        # --- Mechanism 2: Matrix tension (Hashin with in-situ strengths) ---
        # In-situ transverse: Yt_is = 1.12·√2·Yt (Pinho Eq. 47, thin ply)
        # In-situ shear: thick-ply correction from Camanho (2006)
        #   S12_is = sqrt(8·GIIc / (pi·t_eff·Lambda_22))
        #   For typical carbon/epoxy (GIIc ≈ 1.0 N/mm), this gives
        #   S12_is ≈ 2.3·S12 for n_adj=2 adjacent 0-deg plies.
        #   Falls back to sqrt(2)·S12 for single-ply (thin) case.
        n_adj = _max_consecutive_zero_plies(angles)
        t_eff = n_adj * cfg.ply_thickness
        Yt_is = 1.12 * np.sqrt(2.0) * Yt
        # Thick-ply in-situ shear: GIIc ≈ 1.0 N/mm for carbon/epoxy
        _GIIc_typical = 1.0  # N/mm
        _Lambda22 = 2.0 * (1.0 / E22 - (mat.nu12 ** 2) / E11)
        if _Lambda22 > 0 and t_eff > 0:
            S12_is = np.sqrt(8.0 * _GIIc_typical / (np.pi * t_eff * _Lambda22))
        else:
            S12_is = np.sqrt(2.0) * S12

        if theta > 1e-10:
            sin_t = np.sin(theta)
            cos_t = np.cos(theta)
            term1 = (sin_t ** 2 / Yt_is) ** 2
            term2 = (sin_t * cos_t / S12_is) ** 2
            sigma_fail = 1.0 / np.sqrt(term1 + term2)
            kd_matrix = min(sigma_fail / Xt, 1.0)
        else:
            kd_matrix = 1.0

        # --- Mechanism 3: Out-of-plane delamination (curved-beam) ---
        amplitude = cfg.amplitude
        wavelength = cfg.wavelength

        # Effective thickness: max consecutive 0-degree plies. With no
        # continuous 0-degree block (n_adj == 0) the curved-beam model
        # has no fibrous load path to develop interlaminar σ₃₃ / τ₁₃,
        # so the OOP mechanism is inactive (kd_oop = 1.0).
        n_adj_oop = _max_consecutive_zero_plies(angles)

        # Interlaminar stresses at the 0-block interface (held for both
        # the stress-based OOP mechanism above AND the new energy-based
        # onset criterion below).  At ``λ = 1`` (applied stress = Xt)
        # these are the σ₃₃ at the crest and τ₁₃ at the inflection.
        sigma33 = 0.0
        tau13 = 0.0
        h_eff_oop = 0.0
        if n_adj_oop == 0 or amplitude <= 1e-12 or wavelength <= 1e-12:
            kd_oop = 1.0
        else:
            # Peak curvature at crest: κ = (2π/λ)² A
            kappa_max = (2.0 * np.pi / wavelength) ** 2 * amplitude
            # Max curvature gradient at inflection: |dκ/dx| = (2π/λ)³ A
            dkappa_dx_max = (2.0 * np.pi / wavelength) ** 3 * amplitude

            h_eff_oop = n_adj_oop * cfg.ply_thickness

            # σ₃₃ at crest (mode I) and τ₁₃ at inflection (mode II)
            sigma33 = Xt * h_eff_oop * kappa_max
            tau13 = Xt * h_eff_oop * dkappa_dx_max

            # Failure indices (peak at different spatial locations)
            FI_s33 = (sigma33 / Yt) ** 2
            FI_t13 = (tau13 / S13) ** 2
            FI_max = max(FI_s33, FI_t13)

            kd_oop = 1.0 / np.sqrt(1.0 + FI_max)

        # 0-degree ply knockdown: minimum of all three mechanisms
        kd_0 = min(kd_fiber, kd_matrix, kd_oop)

        # ----------------------------------------------------------------
        # Delamination-onset KD (Mukhopadhyay et al. 2015 first-load-drop)
        # ----------------------------------------------------------------
        # The three-mechanism kd_0 above is the *ultimate* fibre-failure
        # KD.  Embedded-wrinkle tests (Mukhopadhyay 2015) also exhibit an
        # earlier *first-load-drop* corresponding to delamination
        # initiation at the curved 0-block interface.  We predict that
        # initiation knockdown with a Benzeggagh-Kenane mode-mixity
        # criterion driven by the same σ₃₃ / τ₁₃ already computed for
        # the OOP mechanism, but compared to GIc / GIIc rather than to
        # the strength allowables Yt / S13.
        #
        # Derivation (notation: λ = applied stress / Xt):
        #   σ₃₃(λ) = λ · σ₃₃   (above, evaluated at λ = 1)
        #   τ₁₃(λ) = λ · τ₁₃
        # Energy release rate at a notional interfacial flaw of size a:
        #   G_I  = σ₃₃² · π · a / (2 · E_3)
        #   G_II = τ₁₃² · π · a / (2 · G_13)
        # Both scale as λ².  B-K-like mode-mixity initiation:
        #   λ²·(G_I/GIc) + (λ²·G_II/GIIc)^η = 1,   η = 1.45
        # Solved for λ_onset ∈ (0, 1].  If the criterion is satisfied
        # already at λ < 1, the onset KD is below the ultimate KD.
        #
        # Flaw size: a = t_eff (the 0-block thickness ``n_adj * t_ply``).
        # The spec proposed a = t_ply but with that single-ply scale the
        # criterion gives λ_onset > 1 for all the Mukhopadhyay cases —
        # i.e. no onset before fibre fracture, which is unphysical.  An
        # embedded delamination at the 0-block interface naturally spans
        # the block thickness, so a = t_eff is the correct local scale.
        kd_onset = None
        if (
            mat.GIc is not None
            and mat.GIIc is not None
            and n_adj_oop > 0
            and amplitude > 1e-12
            and wavelength > 1e-12
        ):
            E3 = mat.E3
            G13 = mat.G13
            GIc = mat.GIc
            GIIc = mat.GIIc
            a_flaw = h_eff_oop  # = t_eff = n_adj_oop * ply_thickness

            G_I_unit = (sigma33 ** 2) * np.pi * a_flaw / (2.0 * E3)
            G_II_unit = (tau13 ** 2) * np.pi * a_flaw / (2.0 * G13)

            R_I = G_I_unit / GIc
            R_II = G_II_unit / GIIc
            eta = 1.45

            # f(λ) = λ²·R_I + (λ²·R_II)^η  -  1.  Monotonically
            # increasing in λ ∈ (0, ∞).
            def _bk_criterion(lam: float) -> float:
                lam2 = lam * lam
                term_I = lam2 * R_I
                arg_II = lam2 * R_II
                term_II = arg_II ** eta if arg_II > 0.0 else 0.0
                return term_I + term_II - 1.0

            f_at_1 = _bk_criterion(1.0)
            if f_at_1 <= 0.0:
                # Energy criterion not met even at λ = 1.  No separate
                # initiation event before fibre fracture; report onset
                # equal to the ultimate.
                lam_onset = 1.0
            else:
                # Bisect in (1e-4, 1] for the unique root.
                lo, hi = 1.0e-4, 1.0
                f_lo = _bk_criterion(lo)
                if f_lo > 0.0:
                    # Even at vanishing load the criterion is satisfied,
                    # which would imply zero strength — guard against it.
                    lam_onset = lo
                else:
                    for _ in range(80):
                        mid = 0.5 * (lo + hi)
                        f_mid = _bk_criterion(mid)
                        if f_mid > 0.0:
                            hi = mid
                        else:
                            lo = mid
                        if hi - lo < 1.0e-6:
                            break
                    lam_onset = 0.5 * (lo + hi)

            kd_onset = float(lam_onset)

        # For graded averaging: return just the 0-deg ply KD without CLT
        if _return_kd0_only:
            return float(kd_0), {}

        # CLT-weighted laminate knockdown
        n_0 = sum(1 for a in angles if abs(a) < 5)
        n_45 = sum(1 for a in angles if 40 < abs(a) < 50)
        n_90 = sum(1 for a in angles if abs(a) > 85)

        Q11_0 = E11
        Q11_45 = E11 / 4.0 + E22 / 4.0 + G12 / 2.0
        Q11_90 = E22

        total_stiffness = n_0 * Q11_0 + n_45 * Q11_45 + n_90 * Q11_90
        if total_stiffness > 0:
            f_0 = n_0 * Q11_0 / total_stiffness
        else:
            f_0 = 1.0

        kd_lam = f_0 * kd_0 + (1.0 - f_0) * 1.0

        # CLT-weight the 0-ply onset KD to the laminate level, matching
        # the kd_0 \u2192 kd_lam pattern.  Then ensure the onset KD is
        # strictly less than the ultimate KD: any local interfacial
        # delamination event must precede (or coincide with) the
        # laminate ultimate, so we cap onset at ``kd_lam * 0.999``
        # whenever the raw energy-based onset is not already below it.
        # This guarantees the spec's requirement that onset KD < KD_oop
        # and onset KD < analytical_knockdown.
        kd_onset_lam: float | None = None
        if kd_onset is not None:
            kd_onset_lam_raw = f_0 * kd_onset + (1.0 - f_0) * 1.0
            kd_onset_lam = min(kd_onset_lam_raw, float(kd_lam) * 0.999)

        # Determine controlling mode
        if kd_oop <= kd_fiber and kd_oop <= kd_matrix:
            mode = "OOP \u03c3\u2083\u2083"
        elif kd_matrix <= kd_fiber:
            mode = "matrix"
        else:
            mode = "fiber"

        mechanisms = {
            "kd_fiber": kd_fiber,
            "kd_matrix": kd_matrix,
            "kd_oop": kd_oop,
            "kd_0": kd_0,
            "kd_lam": float(kd_lam),
            "kd_onset": kd_onset_lam,
            "f_0": f_0,
            "mode": mode,
        }

        return float(kd_lam), mechanisms

    def _build_flat_mesh(self, laminate: Laminate) -> MeshData:
        """Generate a pristine (zero-amplitude) mesh matching the config.

        Shared by the retention-factor baseline and the progressive-damage
        pristine reference so both compare against the same flat laminate.
        """
        cfg = self.config
        flat_profile = GaussianSinusoidal(
            amplitude=0.0,
            wavelength=cfg.wavelength,
            width=cfg.width,
            center=cfg.domain_length / 2.0,
        )
        # Both interface indices are populated by AnalysisConfig.__post_init__.
        assert cfg.interface_1 is not None and cfg.interface_2 is not None
        flat_config = WrinkleConfiguration.from_morphology_name(
            "stack", flat_profile,
            interface1=cfg.interface_1,
            interface2=cfg.interface_2,
        )
        return WrinkleMesh(
            laminate=laminate,
            wrinkle_config=flat_config,
            Lx=cfg.domain_length,
            Ly=cfg.domain_width,
            nx=cfg.nx,
            ny=cfg.ny,
            nz_per_ply=cfg.nz_per_ply,
        ).generate()

    def _run_progressive_path(
        self, results: AnalysisResults, laminate: Laminate, mesh: MeshData
    ) -> None:
        """Run the load-stepping progressive-damage solver to ultimate load.

        Solves the wrinkled mesh (with any resin pocket already attached)
        and a pristine flat baseline, then records the ultimate strengths
        and their ratio as the progressive-damage knockdown.  The wrinkled
        mesh's per-element material override is cleared afterwards so the
        subsequent linear field/failure/retention pass is unaffected.
        """
        from wrinklefe.solver.progressive_damage import (
            ProgressiveDamageResult,
            ProgressiveDamageSolver,
        )

        cfg = self.config
        # Auto-size the strain ramp to bracket fibre failure (~1.8x the
        # compressive failure strain Xc / E1 of the 0-degree material).
        if cfg.progressive_max_strain is not None:
            target = float(cfg.progressive_max_strain)
        else:
            mat0 = laminate.plies[0].material
            eps_f = mat0.Xc / mat0.E1
            target = 1.8 * eps_f
        sign = -1.0 if cfg.applied_strain <= 0 else 1.0
        applied = sign * abs(target)

        def _run(m: MeshData) -> ProgressiveDamageResult:
            return ProgressiveDamageSolver(
                m, laminate,
                applied_strain=applied,
                n_increments=cfg.progressive_n_increments,
                residual_factor=cfg.progressive_residual_factor,
                solver=cfg.solver,
                verbose=cfg.verbose,
                delta_T=cfg.delta_T,
            ).solve()

        # Wrinkled run — snapshot/restore the override so the later linear
        # pass sees the undamaged mesh.
        saved_override = mesh.element_material_override
        wr = _run(mesh)
        mesh.element_material_override = saved_override

        pristine = _run(self._build_flat_mesh(laminate))

        results.progressive_strength_MPa = wr.peak_stress
        results.progressive_pristine_strength_MPa = pristine.peak_stress
        results.progressive_history = wr.history
        if pristine.peak_stress > 0:
            results.progressive_knockdown = (
                wr.peak_stress / pristine.peak_stress
            )

    def _attach_resin_pocket(
        self, mesh: MeshData, laminate: Laminate
    ) -> None:
        """Tag the resin-lens elements and attach the resin material.

        Builds a :class:`~wrinklefe.core.resin_pocket.ResinPocketSpec`
        from the wrinkle geometry and the configured scale knobs, flags
        the hex elements whose centroids fall inside the lens, and stores
        the boolean mask plus the resin material on *mesh* so the
        assembler, stress-recovery and failure paths pick them up.

        The resin material defaults to the built-in ``EPOXY_S6C10``
        isotropic card when ``resin_pocket_material`` is unset.
        """
        from wrinklefe.core.resin_pocket import (
            ResinPocketSpec,
            compute_resin_blend,
            compute_resin_mask,
        )

        cfg = self.config
        # Place the pocket relative to the mesh's ACTUAL z-extent, which
        # is centred on the mid-plane (z in [-T/2, +T/2]) — not the
        # bottom-referenced [0, T].  ``wrinkle_z_position`` maps 0 -> bottom
        # surface, 0.5 -> mid-plane (where the graded morphology places the
        # wrinkle crest), 1 -> top surface.  Using ``z_frac * T`` (the old
        # form) put a mid-plane pocket at the top surface, mis-locating it
        # away from the high-angle crest entirely.
        z_lo = float(mesh.nodes[:, 2].min())
        z_hi = float(mesh.nodes[:, 2].max())
        z_center = z_lo + float(cfg.wrinkle_z_position) * (z_hi - z_lo)
        center_x = cfg.domain_length / 2.0

        spec = ResinPocketSpec.from_wrinkle(
            amplitude=cfg.amplitude,
            wavelength=cfg.wavelength,
            center_x=center_x,
            z_center=z_center,
            height_scale=cfg.resin_pocket_height_scale,
            length_scale=cfg.resin_pocket_length_scale,
        )

        resin_material = cfg.resin_pocket_material
        if resin_material is None:
            resin_material = MaterialLibrary().get("EPOXY_S6C10")
        mesh.resin_material = resin_material

        if cfg.resin_pocket_graded:
            # Graded pocket: per-element blend weight + precomputed blended
            # materials (host ply <-> resin), and the fibre angle scaled by
            # (1 - weight) downstream via ``resin_angle_scale``.
            weight = compute_resin_blend(mesh, spec)
            mesh.resin_blend = weight
            blend_mats: dict[int, OrthotropicMaterial] = {}
            for e_np in np.flatnonzero(weight > 0.0):
                e = int(e_np)
                ply_mat = laminate.plies[int(mesh.ply_ids[e])].material
                blend_mats[e] = ply_mat.blend(
                    resin_material, float(weight[e])
                )
            mesh.resin_blend_materials = blend_mats
            n_resin = int((weight > 0.0).sum())
        else:
            mesh.resin_mask = compute_resin_mask(mesh, spec)
            n_resin = int(mesh.resin_mask.sum())

        if cfg.verbose:
            import logging
            logging.getLogger(__name__).info(
                "Resin pocket (%s): %d/%d elements tagged "
                "(half_length=%.3g mm, h_center=%.3g mm, z_center=%.3g mm)",
                "graded" if cfg.resin_pocket_graded else "binary",
                n_resin, mesh.n_elements, spec.half_length,
                spec.h_center, z_center,
            )

    def _attach_surface_resin_pockets(
        self, mesh: MeshData, laminate: Laminate,
        wrinkle_config: WrinkleConfiguration,
    ) -> None:
        """Tag the surface resin pockets under the tool-flat surface(s).

        Computes the per-element surface-pocket blend weight from the
        deformed mesh (:func:`compute_surface_resin_blend`) and attaches it
        to *mesh* using the same fields the crest lens uses, so the
        assembler / stress-recovery / failure paths pick it up unchanged.

        The two pocket systems **compose**: when the crest lens has already
        written ``resin_blend`` / ``resin_blend_materials`` (or the binary
        ``resin_mask``), the surface weight is merged by per-element maximum
        and the blended materials are rebuilt from the combined weight, so
        no element is double-blended.  Both reuse ``resin_pocket_material``.
        """
        from wrinklefe.core.resin_pocket import (
            SurfacePocketSpec,
            compute_surface_resin_blend,
        )

        cfg = self.config
        spec = SurfacePocketSpec(
            side=cfg.surface_pocket_side,
            min_gap_threshold=cfg.surface_pocket_min_gap,
        )
        weight = compute_surface_resin_blend(mesh, wrinkle_config, spec)

        resin_material = cfg.resin_pocket_material
        if resin_material is None:
            resin_material = MaterialLibrary().get("EPOXY_S6C10")
        mesh.resin_material = resin_material

        if cfg.resin_pocket_graded:
            # Compose with any crest-lens weight already on the mesh, then
            # rebuild the blended materials from the combined weight so ties
            # resolve to the larger resin fraction (no double-blend).
            if mesh.resin_blend is not None:
                weight = np.maximum(mesh.resin_blend, weight)
            mesh.resin_blend = weight
            blend_mats: dict[int, OrthotropicMaterial] = {}
            for e_np in np.flatnonzero(weight > 0.0):
                e = int(e_np)
                ply_mat = laminate.plies[int(mesh.ply_ids[e])].material
                blend_mats[e] = ply_mat.blend(resin_material, float(weight[e]))
            mesh.resin_blend_materials = blend_mats
            n_surface = int((weight > 0.0).sum())
        else:
            # Binary rule for surface pockets: any transition element with a
            # gap above ``surface_pocket_min_gap`` is neat resin.  (Unlike the
            # crest lens, the excess-stretch fraction is a partial fill that
            # rarely exceeds 0.5, so a ``> 0.5`` cut would tag nothing; graded
            # is the default and recommended path.)
            surface_mask = weight > 0.0
            if mesh.resin_mask is not None:
                surface_mask = surface_mask | mesh.resin_mask
            mesh.resin_mask = surface_mask
            n_surface = int(surface_mask.sum())

        max_gap = 0.0
        if np.any(weight > 0.0):
            ez = mesh.nodes[mesh.elements][:, :, 2]
            h = ez.max(axis=1) - ez.min(axis=1)
            max_gap = float(np.max(weight * h))

        if cfg.verbose:
            import logging
            logging.getLogger(__name__).info(
                "Surface resin pockets (%s): %d/%d elements tagged "
                "(max gap %.3g mm)",
                cfg.surface_pocket_side, n_surface, mesh.n_elements, max_gap,
            )

    def _attach_vf_gradient(
        self, mesh: MeshData, laminate: Laminate
    ) -> None:
        """Install the compaction-driven local-``Vf`` materials (issue #379).

        Derives the per-element fibre volume fraction from the deformed
        element heights (``Vf_local = vf_nominal * h0 / h``; see
        :mod:`~wrinklefe.core.compaction`) and attaches one ratio-anchored
        material per occupied ``Vf`` bin.

        Channel
        -------
        The materials go on ``mesh.resin_blend_materials``, **not**
        ``mesh.element_material_override``.  The override dict is owned and
        *mutated in place* by
        :class:`~wrinklefe.solver.progressive_damage.ProgressiveDamageSolver`
        as elements fail, so anything parked there would be overwritten by
        the degraded cards (and the snapshot/restore in
        :meth:`_run_progressive_path` restores the same mutated object).
        ``resin_blend_materials`` sits one step lower in
        :meth:`~wrinklefe.core.mesh.MeshData.element_material` precedence,
        which is exactly right: damage still wins over compaction, and the
        assembler / stress recovery / failure evaluator already consume this
        channel with no solver change.

        ``mesh.resin_blend`` (the fibre-angle suppression weight) is left
        alone: a compacted or resin-enriched element still has fibres, and
        they still carry the wrinkle misalignment.  Where the crest resin
        lens has already blended an element, the blend is re-applied on top
        of the local-``Vf`` host so the two compose without double-counting.
        """
        from wrinklefe.core.compaction import (
            VfGradientSpec,
            build_vf_materials,
            compute_vf_field,
        )

        cfg = self.config
        assert cfg.material is not None
        spec = VfGradientSpec.for_material(
            cfg.material.name,
            fiber=cfg.vf_fiber,
            matrix=cfg.vf_matrix,
            vf_nominal=cfg.vf_nominal,
            vf_max=cfg.vf_max,
        )
        vf_field = compute_vf_field(mesh, spec)

        # One call per distinct ply material so a mixed-material laminate is
        # scaled against the right preset card.
        ply_ids = np.asarray(mesh.ply_ids)
        vf_materials: dict[int, OrthotropicMaterial] = {}
        by_material: dict[int, tuple[OrthotropicMaterial, list[int]]] = {}
        for elem in range(int(mesh.n_elements)):
            ply_material = laminate.plies[int(ply_ids[elem])].material
            entry = by_material.setdefault(id(ply_material), (ply_material, []))
            entry[1].append(elem)
        for ply_material, elems in by_material.values():
            vf_materials.update(
                build_vf_materials(
                    ply_material, vf_field, spec,
                    element_ids=np.asarray(elems, dtype=np.int64),
                )
            )

        # Compose with an already-attached crest resin lens: rebuild its
        # blend from the local-Vf host rather than discarding either.
        blend_weight = mesh.resin_blend
        existing = dict(mesh.resin_blend_materials or {})
        if blend_weight is not None and vf_materials:
            resin_material = cfg.resin_pocket_material
            if resin_material is None:
                resin_material = MaterialLibrary().get("EPOXY_S6C10")
            for elem, material in vf_materials.items():
                w = float(blend_weight[elem])
                existing[elem] = (
                    material.blend(resin_material, w) if w > 0.0 else material
                )
        else:
            existing.update(vf_materials)
        mesh.resin_blend_materials = existing

        n_sat = int(np.count_nonzero(vf_field >= spec.vf_max))
        logger.info(
            "Vf gradient (issue #379): %d/%d elements re-materialised "
            "(Vf %.3f-%.3f about nominal %.3f, %d shared materials, "
            "%d saturated at vf_max=%.3f).",
            len(vf_materials), mesh.n_elements,
            float(vf_field.min()), float(vf_field.max()), spec.vf_nominal,
            len({id(m) for m in vf_materials.values()}), n_sat, spec.vf_max,
        )

    def _evaluate_failure(
        self,
        results: AnalysisResults,
        laminate: Laminate,
        field_results: FieldResults,
        mesh: MeshData,
    ) -> None:
        """Evaluate failure criteria on the FE stress field."""
        evaluator, materials, eval_ply_ids, elem_fiber_angles = (
            self._failure_eval_inputs(laminate, mesh)
        )

        # Field-level evaluation
        fi_fields, mode_fields = evaluator.evaluate_field(
            field_results.stress_local,
            materials,
            eval_ply_ids,
            fiber_angles=elem_fiber_angles,
        )
        results.failure_indices = fi_fields
        results.failure_modes = mode_fields

        # CLT-level evaluation at applied load
        self._evaluate_clt_failure(results, laminate)

    def _failure_eval_inputs(
        self,
        laminate: Laminate,
        mesh: MeshData,
    ) -> tuple[FailureEvaluator, list, np.ndarray, np.ndarray]:
        """Resolve the per-element inputs the failure criteria need.

        Returns ``(evaluator, materials, eval_ply_ids, fiber_angles)``.
        Shared by :meth:`_evaluate_failure` and the proportional
        load-factor search (issue #275) so both evaluate the field
        through exactly the same material routing and misalignment
        scaling — a second copy of this resolution is how the two would
        drift apart.
        """
        evaluator = FailureEvaluator.default_criteria()

        # Build material list for each ply
        materials = [ply.material for ply in laminate.plies]

        # Per-element fiber angles from wrinkle geometry (for LaRC05 kinking)
        elem_fiber_angles = mesh.element_fiber_angles_array()

        # Resin-pocket zone (Li et al. 2024/2025): route lens elements to
        # their pocket material (graded blend, or the binary resin card)
        # so failure is evaluated at the locally-softened strengths, and
        # scale the fibre angle by the retention factor so the LaRC05
        # kink-band path is not double-counted at the resin centre.
        eval_ply_ids = np.asarray(mesh.ply_ids)
        if mesh.resin_blend_materials:
            # Graded pocket: each blended element gets its own material.
            extra = list(mesh.resin_blend_materials.items())
            base = len(materials)
            mat_index = {e: base + i for i, (e, _m) in enumerate(extra)}
            materials = [*materials, *(m for _e, m in extra)]
            eval_ply_ids = eval_ply_ids.copy()
            for e, idx in mat_index.items():
                eval_ply_ids[e] = idx
            if mesh.resin_blend is not None:
                elem_fiber_angles = (
                    elem_fiber_angles * (1.0 - mesh.resin_blend)
                )
        elif mesh.resin_mask is not None and mesh.resin_material is not None:
            resin_idx = len(materials)
            materials = [*materials, mesh.resin_material]
            eval_ply_ids = np.where(
                mesh.resin_mask, resin_idx, mesh.ply_ids
            )
            elem_fiber_angles = np.where(
                mesh.resin_mask, 0.0, elem_fiber_angles
            )

        return evaluator, materials, eval_ply_ids, elem_fiber_angles

    def _clt_load_state(self) -> LoadState:
        """Build the CLT :class:`LoadState` the pipeline evaluates.

        The mechanical part is the same approximate membrane resultant the
        laminate-level failure check has always used
        (``Nx = applied_strain * 1000``); the environmental part carries
        ``AnalysisConfig.delta_T`` (issue #273) so the thermal
        resultants reach :meth:`~wrinklefe.core.laminate.Laminate.midplane_strains`
        and the ply stress recovery.

        ``delta_T`` keeps the config's sign convention unchanged: it is the
        temperature change *from* the stress-free (cure) state, so a
        cool-down is negative.

        Returns
        -------
        LoadState
            Load state for the CLT / analytical failure evaluation.
        """
        cfg = self.config
        return LoadState(
            Nx=cfg.applied_strain * 1000.0,  # approximate
            delta_T=cfg.delta_T,
        )

    def _evaluate_clt_failure(
        self,
        results: AnalysisResults,
        laminate: Laminate,
    ) -> None:
        """Run the laminate-level (CLT) failure evaluation into *results*.

        Shared by the FE path (alongside the field-level evaluation) and
        the analytical path (when a thermal load is present), so both
        consume the same :meth:`_clt_load_state`.
        """
        evaluator = FailureEvaluator.default_criteria()
        load = self._clt_load_state()
        try:
            report = evaluator.evaluate_laminate(laminate, load)
            results.failure_report = report
        except Exception as exc:
            logger.warning("CLT evaluation skipped: %s", exc)

    def _compute_load_factors(
        self,
        results: AnalysisResults,
        laminate: Laminate,
        flat_field: FieldResults,
        flat_mesh: MeshData,
    ) -> None:
        """Proportional load factors for the wrinkled coupon and baseline.

        See :func:`_proportional_load_factor` for the definition and for
        why scaling the stored field (rather than re-solving) is exact
        here.  Both sides are scaled from fields solved under the *same*
        load state, so the ratio is a like-for-like knockdown.
        """
        assert results.field_results is not None
        wrinkled_mesh = results.mesh
        if wrinkled_mesh is None:
            return

        w_eval, w_mats, w_ids, w_angles = self._failure_eval_inputs(
            laminate, wrinkled_mesh,
        )
        p_eval = FailureEvaluator.default_criteria()
        p_mats = [ply.material for ply in laminate.plies]

        def _max_fi(evaluator, stress, materials, ply_ids, angles):
            def at(lam: float) -> float:
                fields, _modes = evaluator.evaluate_field(
                    lam * stress, materials, ply_ids, fiber_angles=angles,
                )
                best = 0.0
                for arr in fields.values():
                    finite = np.asarray(arr)[np.isfinite(arr)]
                    if finite.size:
                        best = max(best, float(finite.max()))
                return best
            return at

        results.load_state_factor = _proportional_load_factor(
            _max_fi(w_eval, results.field_results.stress_local,
                    w_mats, w_ids, w_angles)
        )
        results.load_state_factor_pristine = _proportional_load_factor(
            _max_fi(p_eval, flat_field.stress_local,
                    p_mats, flat_mesh.ply_ids, None)
        )
        if (
            results.load_state_factor is not None
            and results.load_state_factor_pristine
        ):
            results.load_state_factor_knockdown = (
                results.load_state_factor / results.load_state_factor_pristine
            )

    def _compute_retention_factors(
        self,
        results: AnalysisResults,
        laminate: Laminate,
    ) -> None:
        """Compute retention factors by running a pristine (no-wrinkle) baseline.

        Retention = max_FI_pristine / max_FI_wrinkled

        A retention of 1.0 means no knockdown; 0.5 means 50% strength retained.
        """
        cfg = self.config
        # interface_1 / interface_2 are filled in __post_init__.
        assert cfg.interface_1 is not None and cfg.interface_2 is not None

        if results.failure_indices is None:
            return

        # Build a flat (no wrinkle) mesh with same dimensions
        flat_mesh = self._build_flat_mesh(laminate)

        # Solve with same BCs
        # Same thermal state AND the same mechanical load state as the
        # wrinkled run, so the retention factor compares like with like
        # (issues #273 Stage 2 and #275).
        flat_solver = StaticSolver(
            flat_mesh, laminate, delta_T=cfg.delta_T,
            **_iterative_solver_kwargs(cfg)
        )
        flat_bcs = _mechanical_bcs(cfg, flat_mesh)
        flat_field = flat_solver.solve(flat_bcs, solver=cfg.solver, verbose=False)

        # Evaluate failure on flat mesh (no fiber misalignment)
        evaluator = FailureEvaluator.default_criteria()
        materials = [ply.material for ply in laminate.plies]

        flat_fi_fields, _ = evaluator.evaluate_field(
            flat_field.stress_local,
            materials,
            flat_mesh.ply_ids,
            fiber_angles=None,  # no misalignment in pristine laminate
        )

        # Proportional load factors under a general load state (#275).
        # Done here because this is where the pristine field already
        # exists, so the knockdown costs no extra solve.
        if cfg.load_state is not None and results.field_results is not None:
            self._compute_load_factors(
                results, laminate, flat_field, flat_mesh,
            )

        # Compute retention for each criterion
        retention = {}
        baseline = {}

        for crit_name in results.failure_indices:
            # Wrinkled max FI (interior elements)
            fi_w = results.failure_indices[crit_name]
            fi_w_mean = fi_w.mean(axis=-1)  # avg over Gauss pts
            finite_w = fi_w_mean[np.isfinite(fi_w_mean)]
            max_fi_w = float(finite_w.max()) if finite_w.size > 0 else 0.0

            # Pristine max FI
            fi_p = flat_fi_fields[crit_name]
            fi_p_mean = fi_p.mean(axis=-1)
            finite_p = fi_p_mean[np.isfinite(fi_p_mean)]
            max_fi_p = float(finite_p.max()) if finite_p.size > 0 else 0.0

            baseline[crit_name] = max_fi_p

            if max_fi_w > 0:
                # Retention = pristine_FI / wrinkled_FI
                # (how much of the pristine strength is retained)
                retention[crit_name] = min(max_fi_p / max_fi_w, 1.0)
            else:
                retention[crit_name] = 1.0

        results.retention_factors = retention
        results.baseline_fi = baseline

        # --- Modulus retention from FE ---
        # Two complementary estimators of the axial-modulus knockdown
        # ``E_x / E_x0`` (wrinkled vs pristine), both populated here:
        #
        #   modulus_retention        — LOCAL σ₁₁ proxy: E_eff = <σ₁₁> /
        #     ε_applied, the mean element-frame fibre-direction stress over
        #     the coupon divided by the applied strain.  A local proxy that
        #     over-predicts the retention (issue #328).
        #
        #   modulus_retention_global — GLOBAL reaction response: E_eff =
        #     σ_nominal / ε_applied with σ_nominal = R / A, the total axial
        #     reaction on the loaded (x_max) face over the cross-section
        #     area Ly·Lz.  A true coupon-level stiffness that captures load
        #     redistribution around the wrinkle, so it tracks the measured
        #     modulus knockdown more closely (and is lower than the local
        #     proxy for a wrinkled coupon).
        try:
            results.modulus_retention = self._local_modulus_retention(
                results, flat_field, cfg.applied_strain
            )
        except Exception:
            logger.warning(
                "Local (σ₁₁ proxy) modulus-retention computation failed; "
                "forcing modulus_retention to the 1.0 fallback "
                "(modulus_retention_failed=True)",
                exc_info=True,
            )
            results.modulus_retention = 1.0
            results.modulus_retention_failed = True

        # --- Global (reaction-based) modulus retention ---
        try:
            applied_strain = cfg.applied_strain
            wrinkled_mesh = results.mesh
            if applied_strain == 0.0 or wrinkled_mesh is None:
                results.modulus_retention_global = 1.0
            else:
                E_w_global = self._reaction_modulus(
                    wrinkled_mesh, laminate, applied_strain
                )
                E_p_global = self._reaction_modulus(
                    flat_mesh, laminate, applied_strain
                )
                if E_w_global is not None and E_p_global is not None and (
                    abs(E_p_global) > 1e-12
                ):
                    results.modulus_retention_global = float(
                        abs(E_w_global) / abs(E_p_global)
                    )
                else:
                    results.modulus_retention_global = 1.0
        except Exception:
            logger.warning(
                "Global reaction-based modulus-retention computation failed; "
                "forcing modulus_retention_global to the 1.0 fallback "
                "(modulus_retention_global_failed=True)",
                exc_info=True,
            )
            results.modulus_retention_global = 1.0
            results.modulus_retention_global_failed = True

    def _local_modulus_retention(
        self,
        results: AnalysisResults,
        flat_field: FieldResults,
        applied_strain: float,
    ) -> float:
        """LOCAL σ₁₁-proxy axial-modulus retention ``E_wrinkled / E_pristine``.

        ``E_eff = <σ₁₁> / ε_applied`` from the mean element-frame fibre-
        direction stress, wrinkled vs pristine. Returns ``1.0`` for the
        degenerate zero-strain / zero-pristine-stress cases. Extracted so the
        caller's ``except`` can flag a genuine computation failure distinctly
        from this legitimate ``1.0`` (issue #374); numerics are unchanged.
        """
        if applied_strain == 0.0:
            return 1.0
        # Set by the FE path before retention factors run; if it were ever
        # None the caller's except restores the default and flags failure.
        assert results.field_results is not None
        stress_w = results.field_results.stress_local  # (n_elem, n_gauss, 6)
        stress_p = flat_field.stress_local

        # Mean fiber-direction stress σ₁₁ (Voigt component 0)
        s11_w = stress_w[:, :, 0].mean()
        s11_p = stress_p[:, :, 0].mean()

        E_wrinkled = s11_w / applied_strain
        E_pristine = s11_p / applied_strain

        if abs(E_pristine) > 1e-6:
            return float(abs(E_wrinkled) / abs(E_pristine))
        return 1.0

    def _reaction_modulus(
        self,
        mesh: MeshData,
        laminate: Laminate,
        applied_strain: float,
    ) -> float | None:
        """Coupon-level axial modulus from the global reaction force.

        Solves the compression problem on ``mesh`` (retaining the
        unmodified stiffness ``K``) and returns the effective axial modulus
        ``E_eff = σ_nominal / ε_applied`` where the nominal stress is the
        total axial reaction on the loaded ``x_max`` face divided by the
        cross-section area ``Ly·Lz``.

        Reuses the exact reaction-extraction pattern of the
        progressive-damage solver: ``reaction = sum((K @ u)[xmax_dofs])``
        over the loaded-face x-DOFs (``3 * nodes_on_face("x_max")``), so the
        two agree.  Returns ``None`` if the reaction/area/strain cannot give
        a finite modulus.
        """
        if applied_strain == 0.0:
            return None

        # Deliberately thermal-free (delta_T=0): this routine divides the
        # reaction force by area*strain to report a MODULUS.  A cure
        # residual load adds a strain-independent reaction offset, which
        # would show up as a spurious stiffness change (issue #273
        # Stage 2).  Residual stress belongs in the stress/failure output,
        # not in a measured elastic constant.
        solver = StaticSolver(
            mesh, laminate, delta_T=0.0,
            **_iterative_solver_kwargs(self.config)
        )
        bcs = BoundaryHandler.compression_bcs(
            mesh, applied_strain=applied_strain
        )
        field = solver.solve(
            bcs, solver=self.config.solver, verbose=False,
            keep_stiffness=True,
        )

        K = solver._K
        if K is None:
            return None

        xmax_nodes = mesh.nodes_on_face("x_max")
        xmax_dofs = 3 * xmax_nodes  # ux DOFs on the loaded face
        _Lx, Ly, Lz = mesh.domain_size
        area = Ly * Lz
        if area <= 0.0:
            return None

        u = field.displacement.ravel()
        reaction = float(np.sum((K @ u)[xmax_dofs]))
        sigma_nominal = reaction / area
        E_eff = sigma_nominal / applied_strain
        if not np.isfinite(E_eff):
            return None
        return float(E_eff)

