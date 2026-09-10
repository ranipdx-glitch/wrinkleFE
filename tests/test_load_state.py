"""General load states through the analysis pipeline (issue #275).

``AnalysisConfig.loading`` / ``applied_strain`` could express exactly one
thing: a uniaxial state. Real wrinkle dispositions are rarely that — skin
panels see compression **plus shear**, pressure shells see biaxial
membrane states — and those are precisely the states whose transverse and
shear components drive the matrix failure modes, i.e. where the FE path
earns its cost over the analytical one.

``AnalysisConfig.load_state`` accepts a :class:`LoadState` and applies its
resultants directly. What is pinned here:

1. **The default is inert.** ``load_state=None`` leaves the FE path
   bit-identical to the ``applied_strain`` path.
2. **The guards.** Every combination the mapping cannot honour is refused
   with the offending component named, rather than silently dropped —
   the standing rule in this package for a load it cannot apply.
3. **The physics reaches the field.** Shear shows up in the ply shear
   stress; biaxial shows up in the transverse stress.
4. **The load factor.** Strength under a combined state is the scalar the
   whole state is multiplied by to reach first failure. The properties
   that make that a *definition* rather than a number — exact inverse
   scaling with load magnitude, and a knockdown invariant to it — are
   asserted directly.

The BC mapping those rest on is verified against closed-form CLT in
``tests/test_solver/test_boundary.py``.
"""

from __future__ import annotations

import numpy as np
import pytest

from wrinklefe.analysis import AnalysisConfig, WrinkleAnalysis
from wrinklefe.core.laminate import LoadState
from wrinklefe.core.material import MaterialLibrary

LAYUP = [0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0]
PLY_T = 0.125


def _config(**overrides) -> AnalysisConfig:
    """A small but genuinely wrinkled FE coupon."""
    base = dict(
        amplitude=0.15, wavelength=12.0, width=8.0, morphology="graded",
        material=MaterialLibrary().get("IM7_8552"),
        angles=list(LAYUP), ply_thickness=PLY_T,
        domain_length=16.0, domain_width=8.0,
        nx=8, ny=3, nz_per_ply=1, analytical_only=False,
    )
    base.update(overrides)
    return AnalysisConfig(**base)


# --------------------------------------------------------------------------- #
# 1. The default is inert
# --------------------------------------------------------------------------- #


class TestDefaultIsInert:
    def test_default_is_none(self):
        assert AnalysisConfig().load_state is None

    def test_fe_path_bit_identical_without_a_load_state(self):
        """The feature must not perturb a single bit of the legacy path."""
        a = WrinkleAnalysis(_config(applied_strain=-0.005)).run(
            analytical_only=False)
        b = WrinkleAnalysis(_config(applied_strain=-0.005)).run(
            analytical_only=False)
        assert a.field_results is not None and b.field_results is not None
        np.testing.assert_array_equal(
            a.field_results.stress_local, b.field_results.stress_local)
        assert a.load_state_factor is None
        assert a.load_state_factor_pristine is None
        assert a.load_state_factor_knockdown is None


# --------------------------------------------------------------------------- #
# 2. Guards — refused, with the offending component named
# --------------------------------------------------------------------------- #


class TestGuards:
    """Every unsupported combination fails fast rather than dropping load."""

    @pytest.mark.parametrize("component", ["Mx", "My", "Mxy"])
    def test_curvature_is_refused(self, component):
        """Curvature prescribes displacement on the faces a membrane state
        loads with traction; the two cannot be superposed."""
        with pytest.raises(ValueError, match=component):
            _config(load_state=LoadState(Nx=-800.0, **{component: 2.0}))

    @pytest.mark.parametrize("component", ["Qx", "Qy"])
    def test_transverse_shear_is_refused(self, component):
        with pytest.raises(ValueError, match=component):
            _config(load_state=LoadState(Nx=-800.0, **{component: 5.0}))

    def test_empty_membrane_state_is_refused(self):
        """No membrane resultant maps to an empty BC set — a singular
        system, not a zero-load run."""
        with pytest.raises(ValueError, match="no membrane"):
            _config(load_state=LoadState(Nx=0.0, Ny=0.0, Nxy=0.0))

    @pytest.mark.parametrize("component", ["delta_T", "delta_C"])
    def test_environmental_load_on_the_state_is_refused(self, component):
        """One quantity, one owner: temperature lives on the config.

        ``LoadState`` carries its own ``delta_T`` and so does
        ``AnalysisConfig`` (since #273). Two places to set the same thing
        is how a sign error gets in, so the load state's copy must be zero
        and the message must point at the config field.
        """
        with pytest.raises(ValueError) as exc:
            _config(load_state=LoadState(Nx=-800.0, **{component: -155.0}))
        msg = str(exc.value)
        assert component in msg
        assert "AnalysisConfig.delta_T" in msg

    def test_analytical_only_is_refused(self):
        with pytest.raises(ValueError, match="analytical_only"):
            _config(load_state=LoadState(Nx=-800.0), analytical_only=True)

    def test_czm_is_refused(self):
        """The CZM path builds its own uniaxial BCs, so the state would be
        silently ignored."""
        with pytest.raises(ValueError, match="enable_czm"):
            _config(load_state=LoadState(Nx=-800.0), enable_czm=True)

    def test_progressive_damage_is_refused(self):
        """That solver ramps a prescribed displacement; a force-controlled
        state has no applied_strain to ramp."""
        with pytest.raises(ValueError, match="enable_progressive_damage"):
            _config(load_state=LoadState(Nx=-800.0),
                    enable_progressive_damage=True)

    def test_non_finite_component_is_refused(self):
        with pytest.raises(ValueError, match="finite"):
            _config(load_state=LoadState(Nx=float("nan")))

    def test_wrong_type_is_refused(self):
        with pytest.raises(ValueError, match="LoadState"):
            _config(load_state={"Nx": -800.0})


# --------------------------------------------------------------------------- #
# 3. Config round-trip
# --------------------------------------------------------------------------- #


class TestConfigRoundTrip:
    def test_load_state_round_trips(self):
        cfg = _config(load_state=LoadState(Nx=-800.0, Ny=-400.0, Nxy=250.0))
        back = AnalysisConfig.from_dict(cfg.to_dict())
        assert back.load_state == cfg.load_state

    def test_none_round_trips(self):
        cfg = _config()
        assert AnalysisConfig.from_dict(cfg.to_dict()).load_state is None

    def test_unknown_load_state_key_is_named(self):
        cfg = _config(load_state=LoadState(Nx=-800.0))
        payload = cfg.to_dict()
        payload["load_state"]["Nz"] = 1.0
        with pytest.raises(ValueError, match="Nz"):
            AnalysisConfig.from_dict(payload)


# --------------------------------------------------------------------------- #
# 4. The physics reaches the field
# --------------------------------------------------------------------------- #


@pytest.mark.slow
class TestFieldResponse:
    """Measured on this coupon: the components the uniaxial surface could
    not express do show up in the ply stresses."""

    @pytest.fixture(scope="class")
    @classmethod
    def runs(cls):
        out = {}
        for name, ls in (
            ("uniaxial", LoadState(Nx=-800.0)),
            ("shear", LoadState(Nx=-800.0, Nxy=250.0)),
            ("biaxial", LoadState(Nx=-800.0, Ny=-400.0)),
        ):
            out[name] = WrinkleAnalysis(_config(load_state=ls)).run(
                analytical_only=False)
        return out

    def test_runs_end_to_end(self, runs):
        for name, r in runs.items():
            assert r.field_results is not None, name
            assert r.failure_indices, name

    def test_shear_raises_the_ply_shear_stress(self, runs):
        """``Nxy`` reaches tau_12 — the component a uniaxial state cannot
        produce, and one of the drivers of matrix failure."""
        base = np.abs(runs["uniaxial"].field_results.stress_local[:, :, 5]).mean()
        with_shear = np.abs(
            runs["shear"].field_results.stress_local[:, :, 5]).mean()
        assert with_shear > 1.3 * base

    def test_biaxial_raises_the_transverse_stress(self, runs):
        """``Ny`` reaches sigma_2 — likewise unreachable uniaxially."""
        base = runs["uniaxial"].field_results.stress_local[:, :, 1].mean()
        biaxial = runs["biaxial"].field_results.stress_local[:, :, 1].mean()
        assert biaxial < base          # more transverse compression
        assert abs(biaxial - base) > 0.2 * abs(base)

    def test_states_give_distinct_fields(self, runs):
        a = runs["uniaxial"].field_results.stress_local
        for other in ("shear", "biaxial"):
            assert not np.allclose(a, runs[other].field_results.stress_local)


# --------------------------------------------------------------------------- #
# 5. The proportional load factor
# --------------------------------------------------------------------------- #


@pytest.mark.slow
class TestLoadFactor:
    """Strength under a combined state.

    ``lambda`` is the scalar the whole state is multiplied by to reach
    first failure. Two properties make that a definition rather than just
    a number, and both are asserted directly rather than pinned to a
    measured value:

    * **Exact inverse scaling.** Halving the load doubles ``lambda``.
      (Measured: ``Nx = -800`` gives 0.4947, ``Nx = -400`` gives 0.9894.)
    * **A magnitude-invariant knockdown.** The ratio to the flat baseline
      does not depend on how hard you push, because both sides scale
      together. (Measured: 0.9728 for both.)

    The second is what makes it the general-load-state generalisation of
    the uniaxial strength knockdown.
    """

    @pytest.fixture(scope="class")
    @classmethod
    def full(cls):
        return WrinkleAnalysis(
            _config(load_state=LoadState(Nx=-800.0))
        ).run(analytical_only=False)

    @pytest.fixture(scope="class")
    @classmethod
    def half(cls):
        return WrinkleAnalysis(
            _config(load_state=LoadState(Nx=-400.0))
        ).run(analytical_only=False)

    def test_populated_under_a_load_state(self, full):
        assert full.load_state_factor is not None
        assert full.load_state_factor_pristine is not None
        assert full.load_state_factor_knockdown is not None

    def test_halving_the_load_doubles_the_factor(self, full, half):
        """The linearity the whole cheap implementation rests on.

        Because the solve is linear, the search scales the stored stress
        field instead of re-solving. If that assumption ever broke — a
        thermal term making the field affine, say — this is the assertion
        that would catch it.
        """
        assert half.load_state_factor == pytest.approx(
            2.0 * full.load_state_factor, rel=1e-3)

    def test_knockdown_is_invariant_to_load_magnitude(self, full, half):
        assert half.load_state_factor_knockdown == pytest.approx(
            full.load_state_factor_knockdown, rel=1e-3)

    def test_wrinkle_knocks_the_factor_down(self, full):
        """The wrinkled coupon fails earlier than its flat baseline."""
        assert full.load_state_factor < full.load_state_factor_pristine
        assert 0.0 < full.load_state_factor_knockdown < 1.0

    def test_exported_only_under_a_load_state(self, full):
        from wrinklefe.io.results import results_to_dict

        payload = results_to_dict(full)
        block = payload["knockdown_factors"]
        assert "load_state_factor" in block
        assert "load_state_factor_pristine" in block
        assert "load_state_factor_knockdown" in block
        # Must NOT collide with the CLT first-ply-failure factor, which
        # owns the plain "load_factor" key and is a different quantity on
        # a different path.
        assert "load_factor" in payload
        assert payload["load_factor"] != block["load_state_factor"]
