"""Tabular / structured export of :class:`AnalysisResults` to CSV and JSON.

This module complements :mod:`wrinklefe.io.export` (which already provides
the legacy ``export_results_json`` and the mesh-oriented Abaqus/VTK
writers).  The exporters here are designed for downstream comparison and
plotting in Excel, Pandas or Jupyter:

- :func:`export_results_json` writes a *structured*, schema-versioned
  JSON document.  Numeric scalars and short 1D arrays are serialised in
  full; larger arrays (e.g. per-Gauss-point stress fields) are reduced
  to summary statistics so the file stays a few-KB artefact rather than
  a multi-MB dump.
- :func:`export_results_csv` writes a per-ply, Pandas-friendly CSV table
  with one row per ply.
- :func:`results_to_dict` returns the same JSON-shaped dict in-memory,
  so callers (CLI, Streamlit, tests) can post-process the structure
  without round-tripping through a file.

The JSON output is deterministic (``sort_keys=True``) and all numeric
values are plain Python ``float`` / ``int`` so the file round-trips
cleanly through :func:`json.load`.

A schema version field is embedded so future evolutions don't silently
break consumers — bump :data:`SCHEMA_VERSION` whenever the public shape
changes.
"""

from __future__ import annotations

import csv
import json
from dataclasses import fields, is_dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from wrinklefe.analysis import AnalysisConfig, AnalysisResults


# ----------------------------------------------------------------------
# Schema
# ----------------------------------------------------------------------

#: Bump this when the public JSON layout changes in a non-additive way.
SCHEMA_VERSION = "1.1"

#: 1D arrays at or below this length are serialised in full; longer
#: arrays are reduced to {min, max, mean, p95}.  A few-hundred-element
#: ceiling keeps the JSON small without truncating per-ply data for
#: realistic layups.
_ARRAY_INLINE_LIMIT = 256


# ----------------------------------------------------------------------
# JSON-friendly conversion helpers
# ----------------------------------------------------------------------

def _to_jsonable(obj: Any) -> Any:
    """Coerce numpy / dataclass / set values into JSON-native types.

    Used as ``json.dumps(..., default=_to_jsonable)`` for anything the
    inline conversion missed.  Raises ``TypeError`` for genuinely
    unsupported types so silent data loss is impossible.
    """
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (set, frozenset, tuple)):
        return list(obj)
    if isinstance(obj, Path):
        return str(obj)
    if is_dataclass(obj):
        return {f.name: _to_jsonable(getattr(obj, f.name)) for f in fields(obj)}
    raise TypeError(
        f"Object of type {type(obj).__name__!r} is not JSON-serialisable"
    )


def _f(value: Any) -> float:
    """Convert a numeric scalar to a plain Python float."""
    return float(value)


def _array_summary(arr: np.ndarray) -> dict:
    """Reduce a numeric array to (min, max, mean, p95) summary stats.

    The exporter inlines short 1D arrays in full; long arrays
    (per-Gauss-point fields) are summarised here so the JSON stays small.
    """
    a = np.asarray(arr, dtype=np.float64).ravel()
    if a.size == 0:
        return {"min": None, "max": None, "mean": None, "p95": None, "n": 0}
    return {
        "min": _f(np.min(a)),
        "max": _f(np.max(a)),
        "mean": _f(np.mean(a)),
        "p95": _f(np.percentile(a, 95)),
        "n": int(a.size),
    }


def _array_or_summary(arr: Any) -> Any:
    """Inline short 1D arrays; summarise larger / multi-D arrays."""
    if arr is None:
        return None
    a = np.asarray(arr)
    if a.ndim == 1 and a.size <= _ARRAY_INLINE_LIMIT:
        return [_f(x) for x in a]
    return _array_summary(a)


# ----------------------------------------------------------------------
# Config + results -> dict
# ----------------------------------------------------------------------

def _config_to_dict(cfg: AnalysisConfig) -> dict:
    """Serialise an :class:`AnalysisConfig` to a plain JSON-able dict.

    The material is reduced to its name (full property tables are not
    needed in a results export — callers who need them already have the
    material library).
    """
    out: dict = {}
    for f in fields(cfg):
        value = getattr(cfg, f.name)
        if f.name == "material":
            out[f.name] = value.name if value is not None else None
        elif isinstance(value, np.ndarray):
            out[f.name] = value.tolist()
        elif isinstance(value, (list, tuple)):
            out[f.name] = [
                _f(v) if isinstance(v, (np.integer, np.floating)) else v
                for v in value
            ]
        elif isinstance(value, (np.integer,)):
            out[f.name] = int(value)
        elif isinstance(value, (np.floating,)):
            out[f.name] = _f(value)
        else:
            out[f.name] = value
    return out


def _per_ply_rows(results: AnalysisResults) -> list[dict]:
    """Build the per-ply summary table.

    One row per ply with:
        index, angle_deg, max_FI, min_RF, critical_mode, critical_criterion

    If no failure_report is attached (e.g. ``analytical_only=True``) only
    ``index`` and ``angle_deg`` are populated — the FI / mode columns are
    written as ``None`` (JSON ``null``) so the schema is stable.
    """
    cfg = results.config
    angles = list(cfg.angles or [])
    n_plies = len(angles)

    rows: list[dict] = []
    rep = results.failure_report

    # Pre-collect per-criterion FI arrays so we can find the worst across
    # all criteria for each ply.
    per_ply_fi: dict[str, np.ndarray] = {}
    if rep is not None and getattr(rep, "ply_failure_indices", None):
        for crit_name, arr in rep.ply_failure_indices.items():
            per_ply_fi[crit_name] = np.asarray(arr, dtype=np.float64).ravel()

    for k in range(n_plies):
        row: dict = {
            "index": int(k),
            "angle_deg": _f(angles[k]),
            "max_FI": None,
            "min_RF": None,
            "critical_mode": None,
            "critical_criterion": None,
        }
        if per_ply_fi:
            # Worst (largest) FI across criteria for this ply.
            crit_best = None
            fi_best = -np.inf
            for crit_name, fi_arr in per_ply_fi.items():
                if k < fi_arr.size and fi_arr[k] > fi_best:
                    fi_best = float(fi_arr[k])
                    crit_best = crit_name
            if crit_best is not None and np.isfinite(fi_best):
                row["max_FI"] = fi_best
                # Reserve factor = 1 / FI for linear criteria; fall back
                # to None when FI is zero (infinite RF) to keep the CSV
                # round-trippable through csv.DictReader.
                row["min_RF"] = (1.0 / fi_best) if fi_best > 0 else None
                row["critical_criterion"] = crit_best
                # Mode: pull from the FPF/LPF dict if available; the
                # LaminateFailureReport stores per-criterion FPF/LPF
                # which is keyed by criterion, not by ply, so we only
                # populate the mode when the ply happens to be the FPF
                # ply for that criterion.
                fpf = (getattr(rep, "fpf", {}) or {}).get(crit_best, {})
                if fpf.get("ply") == k:
                    row["critical_mode"] = fpf.get("mode") or None
                else:
                    lpf = (getattr(rep, "lpf", {}) or {}).get(crit_best, {})
                    if lpf.get("ply") == k:
                        row["critical_mode"] = lpf.get("mode") or None
        rows.append(row)
    return rows


def _first_ply_failure(results: AnalysisResults) -> dict | None:
    """First-ply failure summary across all criteria (lowest load factor)."""
    rep = results.failure_report
    if rep is None or not getattr(rep, "fpf", None):
        return None
    crit = rep.critical_criterion or next(iter(rep.fpf))
    data = rep.fpf.get(crit, {})
    return {
        "criterion": crit,
        "ply_index": int(data.get("ply", rep.critical_ply or 0)),
        "mode": data.get("mode") or rep.critical_mode or None,
        "load_factor": (
            _f(data.get("load_factor"))
            if data.get("load_factor") is not None
            else None
        ),
    }


def _knockdown_factors(results: AnalysisResults) -> dict:
    """Bundle the analytical + FE knockdown / retention factors."""
    out: dict = {
        "analytical": _f(results.analytical_knockdown),
        "analytical_modulus": _f(results.analytical_modulus_knockdown),
    }
    # FE-derived strength retention factors are keyed by criterion in
    # results.retention_factors; we surface a single ``fe`` scalar as the
    # minimum (most conservative) retention across criteria, plus the
    # per-criterion breakdown for callers who want detail.
    retention = results.retention_factors or {}
    if retention:
        rf_values = [float(v) for v in retention.values() if v is not None]
        if rf_values:
            out["fe"] = min(rf_values)
        out["fe_per_criterion"] = {k: float(v) for k, v in retention.items()}
    if results.modulus_retention is not None:
        out["modulus_retention"] = _f(results.modulus_retention)
    if results.modulus_retention_global is not None:
        out["modulus_retention_global"] = _f(results.modulus_retention_global)
    # Only emitted on failure so machine consumers can tell a fallback 1.0
    # from a genuinely computed 1.0; absent (and byte-identical) for valid
    # runs, preserving ledger zero-drift.
    if results.modulus_retention_failed:
        out["modulus_retention_failed"] = True
    if results.modulus_retention_global_failed:
        out["modulus_retention_global_failed"] = True
    # Proportional load factor under a general load state (issue #275).
    # Emitted only when AnalysisConfig.load_state was set — absent (and
    # byte-identical) for every applied_strain run, preserving ledger
    # zero-drift.
    if results.load_state_factor is not None:
        out["load_state_factor"] = _f(results.load_state_factor)
    if results.load_state_factor_pristine is not None:
        out["load_state_factor_pristine"] = _f(results.load_state_factor_pristine)
    if results.load_state_factor_knockdown is not None:
        out["load_state_factor_knockdown"] = _f(results.load_state_factor_knockdown)
    # CZM/Newton convergence-failure diagnostics — present only on a
    # non-converged nonlinear solve (both None when converged), so they are
    # absent and byte-identical for valid runs, preserving ledger zero-drift.
    if results.czm_failure_diagnostics is not None:
        out["czm_failure_diagnostics"] = results.czm_failure_diagnostics
    if results.czm_failure_hint is not None:
        out["czm_failure_hint"] = results.czm_failure_hint
    return out


def _load_factor(results: AnalysisResults) -> float:
    """Top-level "load factor" surfaced in the export.

    For a wrinkled laminate the most useful single number is the FPF
    load factor from the failure report.  When no FE/failure report is
    available (analytical_only runs) we fall back to the analytical
    knockdown so the field is always populated.
    """
    rep = results.failure_report
    if rep is not None and getattr(rep, "fpf", None):
        crit = rep.critical_criterion or next(iter(rep.fpf))
        lf = rep.fpf.get(crit, {}).get("load_factor")
        if lf is not None and np.isfinite(lf):
            return _f(lf)
    return _f(results.analytical_knockdown)


def _stress_field_summary(results: AnalysisResults) -> dict | None:
    """Reduce the (n_elem, n_gauss, 6) stress tensor to summary stats."""
    fr = results.field_results
    if fr is None or fr.stress_global is None:
        return None
    stress = np.asarray(fr.stress_global, dtype=np.float64)
    if stress.size == 0:
        return None
    component_names = (
        "stress_11", "stress_22", "stress_33",
        "stress_23", "stress_13", "stress_12",
    )
    out: dict = {}
    # Voigt index along last axis.
    flat = stress.reshape(-1, stress.shape[-1])
    for i, name in enumerate(component_names):
        out[name] = _array_summary(flat[:, i])
    return out


def results_to_dict(results: AnalysisResults) -> dict:
    """Build the JSON-shaped dict for :func:`export_results_json`.

    Exposed so callers (CLI, Streamlit, tests) can introspect the same
    payload that would be written to disk.
    """
    cfg = results.config

    payload: dict = {
        "schema_version": SCHEMA_VERSION,
        "config": _config_to_dict(cfg),
        "load_factor": _load_factor(results),
        "analytical": {
            "morphology_factor": _f(results.morphology_factor),
            "max_angle_rad": _f(results.max_angle_rad),
            "max_angle_deg": _f(np.degrees(results.max_angle_rad)),
            "effective_angle_rad": _f(results.effective_angle_rad),
            "effective_angle_deg": _f(np.degrees(results.effective_angle_rad)),
            "damage_index": _f(results.damage_index),
            "analytical_knockdown": _f(results.analytical_knockdown),
            "analytical_onset_knockdown": (
                _f(results.analytical_onset_knockdown)
                if results.analytical_onset_knockdown is not None
                else None
            ),
            "analytical_modulus_knockdown": _f(
                results.analytical_modulus_knockdown
            ),
            "analytical_strength_MPa": _f(results.analytical_strength_MPa),
            "gamma_Y_eff": _f(results.gamma_Y_eff),
        },
        "first_ply_failure": _first_ply_failure(results),
        "per_ply": _per_ply_rows(results),
        "knockdown_factors": _knockdown_factors(results),
    }

    # Optional mesh / FE summaries — only included when present so
    # analytical-only runs produce a compact file.
    if results.mesh is not None:
        payload["mesh"] = {
            "n_nodes": int(results.mesh.n_nodes),
            "n_elements": int(results.mesh.n_elements),
            "n_dof": int(results.mesh.n_dof),
        }

    if results.field_results is not None:
        max_disp, max_disp_node = results.field_results.max_displacement()
        payload["fe"] = {
            "max_displacement_mm": _f(max_disp),
            "max_displacement_node": int(max_disp_node),
            "stress_field_summary": _stress_field_summary(results),
        }

    if results.tension_mechanisms:
        payload["tension_mechanisms"] = {
            k: (_f(v) if isinstance(v, (int, float, np.floating)) else v)
            for k, v in results.tension_mechanisms.items()
        }

    # Progressive-damage block — only when a progressive run happened.
    # ``progressive_history`` is the reliable sentinel (the scalar fields
    # default to 0.0 / 1.0); it holds ``(applied_strain, nominal_stress)``
    # samples from the wrinkled coupon's load history.
    if results.progressive_history is not None:
        payload["progressive"] = {
            "strength_MPa": _f(results.progressive_strength_MPa),
            "pristine_strength_MPa": _f(
                results.progressive_pristine_strength_MPa
            ),
            "knockdown": _f(results.progressive_knockdown),
            "n_increments": len(results.progressive_history),
            "history": [
                [_f(strain), _f(stress)]
                for strain, stress in results.progressive_history
            ],
        }

    return payload


# ----------------------------------------------------------------------
# Public file writers
# ----------------------------------------------------------------------

def export_results_json(
    results: AnalysisResults,
    path: str | Path,
) -> None:
    """Write an :class:`AnalysisResults` as a schema-versioned JSON file.

    The output is deterministic (``sort_keys=True``) and uses plain
    Python ``int`` / ``float`` everywhere so it round-trips cleanly
    through :func:`json.load`.

    Parameters
    ----------
    results : AnalysisResults
        Output of :meth:`WrinkleAnalysis.run`.
    path : str or Path
        Output file path; parent directories are created if needed.

    Notes
    -----
    Numpy arrays longer than ``_ARRAY_INLINE_LIMIT`` (or multi-
    dimensional) are reduced to ``{min, max, mean, p95, n}`` summary
    stats so a stress field at every Gauss point does not bloat the
    file.

    See Also
    --------
    export_results_csv : Pandas-friendly per-ply tabular export.
    """
    payload = results_to_dict(results)
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(
        json.dumps(
            payload, sort_keys=True, indent=2, default=_to_jsonable
        ),
        encoding="utf-8",
    )


# CSV column order — kept stable so downstream consumers can pin to it.
_CSV_COLUMNS = (
    "ply_index",
    "angle_deg",
    "max_FI",
    "min_RF",
    "critical_mode",
    "critical_criterion",
)


def export_results_csv(
    results: AnalysisResults,
    path: str | Path,
) -> None:
    """Write the per-ply summary table as a Pandas-friendly CSV.

    Columns
    -------
    ``ply_index, angle_deg, max_FI, min_RF, critical_mode, critical_criterion``

    The row count equals ``len(results.config.angles)``.  When no
    failure report is attached (e.g. ``analytical_only=True``) the FI /
    mode columns are written as empty strings, but the table shape is
    preserved.

    Parameters
    ----------
    results : AnalysisResults
        Output of :meth:`WrinkleAnalysis.run`.
    path : str or Path
        Output file path; parent directories are created if needed.

    See Also
    --------
    export_results_json : Structured JSON export with all numeric scalars.
    """
    rows = _per_ply_rows(results)
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=_CSV_COLUMNS)
        writer.writeheader()
        for r in rows:
            writer.writerow({
                "ply_index": r["index"],
                "angle_deg": r["angle_deg"],
                "max_FI": (
                    f"{r['max_FI']:.10g}" if r["max_FI"] is not None else ""
                ),
                "min_RF": (
                    f"{r['min_RF']:.10g}" if r["min_RF"] is not None else ""
                ),
                "critical_mode": r["critical_mode"] or "",
                "critical_criterion": r["critical_criterion"] or "",
            })
