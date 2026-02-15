#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import statistics
import subprocess
import time
from pathlib import Path


BASE_INPUT = Path("inputs/current/bath_calibration.i")
OVERRIDE_COMMON = Path("inputs/solver_overrides/quick_test_common.i")

RUN_MATRIX = [
    ("baseline", Path("inputs/solver_overrides/solve_baseline_smp_lu.i")),
    ("fsp+amg", Path("inputs/solver_overrides/solve_fsp3_amg.i")),
    (
        "fsp+amg+globalize+lag",
        Path("inputs/solver_overrides/solve_fsp3_amg_globalize_lag.i"),
    ),
]

FALLBACK_OVERRIDE = Path("inputs/solver_overrides/solve_fsp3_gamg.i")

SOLVER_CSV = Path("outputs/bench_solver_watch.csv")
MESH_CSV = Path("outputs/bench_mesh_watch.csv")
SUMMARY_JSON = Path("outputs/bench_summary.json")
LOG_DIR = Path("outputs/bench_logs")


def _is_hypre_unavailable(log: str) -> bool:
    text = log.lower()
    if "hypre" not in text:
        return False
    needles = (
        "unknown type hypre",
        "pc type hypre is invalid",
        "not configured with hypre",
        "hypre support is disabled",
        "unable to locate hypre",
    )
    return any(n in text for n in needles)


def _run_case(name: str, override: Path) -> dict:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    SOLVER_CSV.unlink(missing_ok=True)
    MESH_CSV.unlink(missing_ok=True)

    cmd = ["./collie-opt", "-i", str(BASE_INPUT), str(OVERRIDE_COMMON), str(override)]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    wall = time.perf_counter() - t0

    log_path = LOG_DIR / f"{name.replace('+', '_')}.log"
    log_path.write_text(proc.stdout)

    return {
        "name": name,
        "override": str(override),
        "cmd": cmd,
        "returncode": proc.returncode,
        "wall_time_s": wall,
        "log_path": str(log_path),
        "log_text": proc.stdout,
    }


def _read_csv_rows(path: Path) -> list[dict]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def _col_values(rows: list[dict], key: str) -> list[float]:
    vals = []
    for row in rows:
        if key not in row:
            continue
        raw = row[key].strip()
        if not raw:
            continue
        try:
            vals.append(float(raw))
        except ValueError:
            continue
    return vals


def _final_value(rows: list[dict], key: str) -> float:
    if not rows:
        return math.nan
    raw = rows[-1].get(key, "").strip()
    try:
        return float(raw)
    except ValueError:
        return math.nan


def _metric_mean(vals: list[float]) -> float:
    return statistics.fmean(vals) if vals else math.nan


def _metric_max(vals: list[float]) -> float:
    return max(vals) if vals else math.nan


def _metric_min(vals: list[float]) -> float:
    return min(vals) if vals else math.nan


def _fmt(x: float, prec: int = 6) -> str:
    if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
        return "nan"
    return f"{x:.{prec}g}"


def _collect_metrics(run_info: dict) -> dict:
    if run_info["returncode"] != 0:
        return {
            "status": "failed",
            "returncode": run_info["returncode"],
            "wall_time_s": run_info["wall_time_s"],
            "log_path": run_info["log_path"],
        }

    if not SOLVER_CSV.exists() or not MESH_CSV.exists():
        return {
            "status": "failed_missing_csv",
            "returncode": run_info["returncode"],
            "wall_time_s": run_info["wall_time_s"],
            "log_path": run_info["log_path"],
        }

    solver_rows = _read_csv_rows(SOLVER_CSV)
    mesh_rows = _read_csv_rows(MESH_CSV)

    nonlinear_its = _col_values(solver_rows, "nonlinear_its")
    linear_its = _col_values(solver_rows, "linear_its")
    mesh_dt = _col_values(mesh_rows, "dt")
    solver_dt = _col_values(solver_rows, "dt")
    dt_vals = mesh_dt if mesh_dt else solver_dt
    mesh_warns = _col_values(solver_rows, "mesh_distortion_warning")

    metrics = {
        "status": "ok",
        "returncode": run_info["returncode"],
        "wall_time_s": run_info["wall_time_s"],
        "timesteps": len(solver_rows),
        "nonlinear_its_mean": _metric_mean(nonlinear_its),
        "nonlinear_its_max": _metric_max(nonlinear_its),
        "linear_its_mean": _metric_mean(linear_its),
        "linear_its_max": _metric_max(linear_its),
        "dt_min": _metric_min(dt_vals),
        "dt_max": _metric_max(dt_vals),
        "bath_value_final": _final_value(solver_rows, "bath_value"),
        "flux_out_final": _final_value(solver_rows, "flux_out"),
        "n_min_final": _final_value(solver_rows, "n_min"),
        "n_max_final": _final_value(solver_rows, "n_max"),
        "mesh_distortion_warning_count": int(sum(1 for x in mesh_warns if abs(x - 1.0) < 1e-12)),
        "log_path": run_info["log_path"],
    }
    return metrics


def _print_commands() -> None:
    print("Commands:")
    print(
        "./collie-opt -i inputs/current/bath_calibration.i "
        "inputs/solver_overrides/quick_test_common.i "
        "inputs/solver_overrides/solve_baseline_smp_lu.i"
    )
    print(
        "./collie-opt -i inputs/current/bath_calibration.i "
        "inputs/solver_overrides/quick_test_common.i "
        "inputs/solver_overrides/solve_fsp3_amg.i"
    )
    print(
        "./collie-opt -i inputs/current/bath_calibration.i "
        "inputs/solver_overrides/quick_test_common.i "
        "inputs/solver_overrides/solve_fsp3_amg_globalize_lag.i"
    )
    print(
        "./collie-opt -i inputs/current/bath_calibration.i "
        "inputs/solver_overrides/quick_test_common.i "
        "inputs/solver_overrides/solve_fsp3_gamg.i"
    )


def main() -> int:
    if not BASE_INPUT.exists():
        raise FileNotFoundError(f"Missing base input: {BASE_INPUT}")

    results = {}
    metadata = {}

    for name, override in RUN_MATRIX:
        run_info = _run_case(name, override)
        used_override = override
        fallback_used = False

        if name == "fsp+amg" and run_info["returncode"] != 0 and _is_hypre_unavailable(run_info["log_text"]):
            run_info = _run_case("fsp+gamg_fallback", FALLBACK_OVERRIDE)
            used_override = FALLBACK_OVERRIDE
            fallback_used = True

        metrics = _collect_metrics(run_info)
        results[name] = metrics
        metadata[name] = {
            "override_requested": str(override),
            "override_used": str(used_override),
            "fallback_used": fallback_used,
            "command": run_info["cmd"],
        }

    print("\nQuick Solver Benchmark")
    header = (
        f"{'config':<28} {'nl_mean':>9} {'nl_max':>8} {'lin_mean':>9} {'lin_max':>8} "
        f"{'dt_min':>9} {'dt_max':>9} {'bath_f':>11} {'flux_f':>11} {'wall_s':>9}"
    )
    print(header)
    print("-" * len(header))
    for name in ("baseline", "fsp+amg", "fsp+amg+globalize+lag"):
        m = results.get(name, {})
        if m.get("status") != "ok":
            print(f"{name:<28} FAILED (status={m.get('status')}, rc={m.get('returncode')})")
            continue
        print(
            f"{name:<28} "
            f"{_fmt(m['nonlinear_its_mean']):>9} {_fmt(m['nonlinear_its_max']):>8} "
            f"{_fmt(m['linear_its_mean']):>9} {_fmt(m['linear_its_max']):>8} "
            f"{_fmt(m['dt_min']):>9} {_fmt(m['dt_max']):>9} "
            f"{_fmt(m['bath_value_final']):>11} {_fmt(m['flux_out_final']):>11} "
            f"{_fmt(m['wall_time_s']):>9}"
        )

    invariance = {"warnings": [], "ok": True}
    base = results.get("baseline", {})
    if base.get("status") == "ok":
        bath_base = base["bath_value_final"]
        flux_base = base["flux_out_final"]
        flux_tol = 1e-6 * max(1.0, abs(flux_base))
        for name in ("fsp+amg", "fsp+amg+globalize+lag"):
            m = results.get(name, {})
            if m.get("status") != "ok":
                invariance["ok"] = False
                invariance["warnings"].append(
                    f"{name}: no invariance check (status={m.get('status')})"
                )
                continue

            bath_diff = abs(m["bath_value_final"] - bath_base)
            flux_diff = abs(m["flux_out_final"] - flux_base)
            if bath_diff >= 1e-6 or flux_diff >= flux_tol:
                invariance["ok"] = False
                severity = "LARGE" if (bath_diff > 1e-3 or flux_diff > 1e-3 * max(1.0, abs(flux_base))) else "WARN"
                invariance["warnings"].append(
                    f"{severity}: {name} differs from baseline "
                    f"(|dbath|={bath_diff:.6g}, |dflux|={flux_diff:.6g}, flux_tol={flux_tol:.6g})"
                )

    if invariance["warnings"]:
        print("\nInvariance checks:")
        for msg in invariance["warnings"]:
            print(f"  - {msg}")
    else:
        print("\nInvariance checks: OK")

    _print_commands()

    summary = {
        "base_input": str(BASE_INPUT),
        "common_override": str(OVERRIDE_COMMON),
        "results": results,
        "metadata": metadata,
        "invariance": invariance,
    }
    SUMMARY_JSON.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2))
    print(f"\nWrote: {SUMMARY_JSON}")

    failed = [k for k, v in results.items() if v.get("status") != "ok"]
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
