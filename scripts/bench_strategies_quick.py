#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import os
import shutil
import statistics
import subprocess
import tempfile
import time
from pathlib import Path


SUMMARY_JSON = Path("outputs/bench_strategies_summary.json")
LOG_DIR = Path("outputs/bench_logs")
DEFAULT_TIMEOUT_S = int(os.environ.get("BENCH_CASE_TIMEOUT", "300"))


def detect_exe() -> str:
    env = os.environ.get("COLLIE_EXE")
    if env:
        return env
    for cand in ("./collie-opt", "./collie-openmp-opt", "./collie-dbg"):
        if Path(cand).exists():
            return cand
    raise FileNotFoundError("No executable found. Set COLLIE_EXE or build collie-opt.")


def detect_base_input() -> Path:
    preferred = Path("inputs/bath_calibration.i")
    if preferred.exists():
        return preferred
    fallback = Path("inputs/current/bath_calibration.i")
    if fallback.exists():
        return fallback
    raise FileNotFoundError("Missing base input: inputs/bath_calibration.i or inputs/current/bath_calibration.i")


def _read_csv_rows(path: Path) -> list[dict]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def _col_values(rows: list[dict], key: str) -> list[float]:
    vals = []
    for row in rows:
        raw = row.get(key, "")
        if raw is None:
            continue
        raw = raw.strip()
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
    raw = rows[-1].get(key, "")
    if raw is None:
        return math.nan
    raw = raw.strip()
    if not raw:
        return math.nan
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


def _write_case_override(tag: str) -> Path:
    Path("outputs").mkdir(parents=True, exist_ok=True)
    tmp = tempfile.NamedTemporaryFile("w", suffix=".i", prefix=f"{tag}_case_", delete=False)
    path = Path(tmp.name)
    tmp.write(
        "\n".join(
            [
                "Outputs/exodus := false",
                f"Outputs/solver_watch/file_base := outputs/{tag}_solver_watch",
                f"Outputs/mesh_watch/file_base := outputs/{tag}_mesh_watch",
                # Use case-local checkpoint roots to avoid cross-case interference.
                f"Outputs/chk/file_base := outputs/{tag}_chk",
                "",
            ]
        )
    )
    tmp.close()
    return path


def _clean_tag_outputs(tag: str) -> None:
    out = Path("outputs")
    if not out.exists():
        return
    for p in out.glob(f"{tag}_*"):
        try:
            if p.is_dir():
                shutil.rmtree(p, ignore_errors=True)
            else:
                p.unlink(missing_ok=True)
        except OSError:
            pass


def _run_and_log(name: str, cmd: list[str], timeout_s: int = DEFAULT_TIMEOUT_S) -> tuple[int, float, Path, bool]:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"{name}.log"
    t0 = time.perf_counter()
    timed_out = False
    out = ""
    rc = 0
    try:
        proc = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=timeout_s
        )
        out = proc.stdout
        rc = proc.returncode
    except subprocess.TimeoutExpired as e:
        timed_out = True
        out = e.stdout or ""
        if isinstance(out, bytes):
            out = out.decode(errors="replace")
        rc = 124
    wall = time.perf_counter() - t0
    log_path.write_text(out)
    return rc, wall, log_path, timed_out


def _collect_metrics(
    tag: str, rc: int, wall: float, log_path: Path, command: list[str], timed_out: bool = False
) -> dict:
    solver_csv = Path(f"outputs/{tag}_solver_watch.csv")
    mesh_csv = Path(f"outputs/{tag}_mesh_watch.csv")

    if timed_out:
        return {
            "status": "timed_out",
            "returncode": rc,
            "wall_time_s": wall,
            "log_path": str(log_path),
            "command": command,
        }
    if rc != 0:
        return {
            "status": "failed",
            "returncode": rc,
            "wall_time_s": wall,
            "log_path": str(log_path),
            "command": command,
        }
    if not solver_csv.exists():
        return {
            "status": "failed_missing_solver_csv",
            "returncode": rc,
            "wall_time_s": wall,
            "log_path": str(log_path),
            "command": command,
        }

    solver_rows = _read_csv_rows(solver_csv)
    mesh_rows = _read_csv_rows(mesh_csv) if mesh_csv.exists() else []

    nonlinear_its = _col_values(solver_rows, "nonlinear_its")
    linear_its = _col_values(solver_rows, "linear_its")
    dt_solver = _col_values(solver_rows, "dt")
    dt_mesh = _col_values(mesh_rows, "dt")
    dt_vals = dt_solver if dt_solver else dt_mesh

    metrics = {
        "status": "ok",
        "returncode": rc,
        "wall_time_s": wall,
        "timesteps": len(solver_rows),
        "nonlinear_its_mean": _metric_mean(nonlinear_its),
        "nonlinear_its_max": _metric_max(nonlinear_its),
        "linear_its_mean": _metric_mean(linear_its),
        "linear_its_max": _metric_max(linear_its),
        "dt_min": _metric_min(dt_vals),
        "dt_max": _metric_max(dt_vals),
        "bath_value_final": _final_value(solver_rows, "bath_value"),
        "flux_full_final": _final_value(solver_rows, "flux_full"),
        "n_min_final": _final_value(solver_rows, "n_min"),
        "n_max_final": _final_value(solver_rows, "n_max"),
        "log_path": str(log_path),
        "solver_csv": str(solver_csv),
        "mesh_csv": str(mesh_csv),
        "command": command,
    }
    return metrics


def run_direct_case(
    name: str,
    tag: str,
    exe: str,
    base: Path,
    common: Path,
    overrides: list[Path],
) -> dict:
    _clean_tag_outputs(tag)
    case_override = _write_case_override(tag)
    try:
        cmd = [exe, "-w", "-i", str(base), str(common)] + [str(p) for p in overrides] + [str(case_override)]
        rc, wall, log_path, timed_out = _run_and_log(name, cmd)
        return _collect_metrics(tag, rc, wall, log_path, cmd, timed_out=timed_out)
    finally:
        case_override.unlink(missing_ok=True)


def run_driver_case(
    name: str,
    tag: str,
    script: Path,
    exe: str,
    base: Path,
    common: Path,
    extra_args: list[str] | None = None,
) -> dict:
    extra_args = extra_args or []
    _clean_tag_outputs(tag)
    cmd = [
        "python3",
        str(script),
        "--exe",
        exe,
        "--base",
        str(base),
        "--common",
        str(common),
        "--tag",
        tag,
    ] + extra_args
    rc, wall, log_path, timed_out = _run_and_log(name, cmd)
    return _collect_metrics(tag, rc, wall, log_path, cmd, timed_out=timed_out)


def recommendation(results: dict) -> str:
    candidates = []
    for name, m in results.items():
        if m.get("status") != "ok":
            continue
        if m.get("n_min_final", 0.0) < -1e-10:
            continue
        candidates.append((name, m))
    if not candidates:
        return "No successful positive-strict strategy; keep baseline and inspect logs."
    best_name, _ = min(
        candidates,
        key=lambda kv: (
            kv[1].get("nonlinear_its_max", math.inf),
            kv[1].get("linear_its_mean", math.inf),
            kv[1].get("dt_min", -math.inf) * -1.0,
        ),
    )
    return f"Recommended: {best_name} (best nonlinear/linear iteration profile among successful positive runs)."


def main() -> int:
    exe = detect_exe()
    base = detect_base_input()
    common = Path("inputs/solver_overrides/common_quick.i")

    paths = {
        "A_diff": Path("inputs/solver_overrides/A_diff.i"),
        "A_rxn": Path("inputs/solver_overrides/A_rxn.i"),
        "B_bounds_vi": Path("inputs/solver_overrides/B_bounds_vi.i"),
        "B_artdiff": Path("inputs/solver_overrides/B_artdiff.i"),
        "C_vlc": Path("inputs/solver_overrides/C_vlc.i"),
        "C_quad8": Path("inputs/solver_overrides/C_quad8.i"),
        "D_ptc": Path("inputs/solver_overrides/D_ptc_ramp_lag.i"),
        "imex_driver": Path("scripts/imex_split_driver.py"),
        "ptc_driver": Path("scripts/ptc_ramp_driver.py"),
    }

    results: dict[str, dict] = {}
    metadata: dict[str, dict] = {}

    # baseline
    results["baseline"] = run_direct_case(
        "baseline",
        "bench_baseline",
        exe,
        base,
        common,
        [],
    )
    metadata["baseline"] = {"fallback_used": False}

    # A via IMEX split driver
    results["A_imex_split"] = run_driver_case(
        "A_imex_split",
        "bench_A",
        paths["imex_driver"],
        exe,
        base,
        common,
        # Keep IMEX subprocess splitting cheap in quick diagnostics.
        extra_args=["--macro-dt", "0.05", "--end-time", "0.05"],
    )
    metadata["A_imex_split"] = {"fallback_used": "imex log indicates fallback if split failed"}

    # B: VI bounds
    results["B_bounds_vi"] = run_direct_case(
        "B_bounds_vi",
        "bench_B_bounds_vi",
        exe,
        base,
        common,
        [paths["B_bounds_vi"]],
    )
    metadata["B_bounds_vi"] = {"fallback_used": False}
    b_vi = results["B_bounds_vi"]
    if b_vi.get("status") != "ok" or b_vi.get("n_min_final", 0.0) < -1e-10:
        results["B_bounds_vi_artdiff"] = run_direct_case(
            "B_bounds_vi_artdiff",
            "bench_B_bounds_vi_artdiff",
            exe,
            base,
            common,
            [paths["B_bounds_vi"], paths["B_artdiff"]],
        )
        metadata["B_bounds_vi_artdiff"] = {"fallback_used": True}

    # C: locking correction, fallback quad8 if unsupported or no effect
    results["C_vlc"] = run_direct_case(
        "C_vlc",
        "bench_C_vlc",
        exe,
        base,
        common,
        [paths["C_vlc"]],
    )
    metadata["C_vlc"] = {"fallback_used": False}

    c_need_fallback = results["C_vlc"].get("status") != "ok"
    if c_need_fallback:
        results["C_quad8"] = run_direct_case(
            "C_quad8",
            "bench_C_quad8",
            exe,
            base,
            common,
            [paths["C_quad8"]],
        )
        metadata["C_quad8"] = {"fallback_used": True}

    # D: run dedicated driver (direct control first, then staged fallback)
    results["D_ptc_ramp_lag"] = run_driver_case(
        "D_ptc_ramp_lag",
        "bench_D",
        paths["ptc_driver"],
        exe,
        base,
        common,
    )
    metadata["D_ptc_ramp_lag"] = {"fallback_used": "ptc driver log indicates fallback if direct control failed"}

    # Invariance checks against baseline
    invariance = {"warnings": [], "ok": True}
    base_metrics = results.get("baseline", {})
    if base_metrics.get("status") == "ok":
        bath_base = base_metrics.get("bath_value_final", math.nan)
        flux_base = base_metrics.get("flux_full_final", math.nan)
        flux_tol = 1e-6 * max(1.0, abs(flux_base)) if not math.isnan(flux_base) else math.nan
        for name, m in results.items():
            if name == "baseline" or m.get("status") != "ok":
                continue
            bath = m.get("bath_value_final", math.nan)
            flux = m.get("flux_full_final", math.nan)
            if math.isnan(bath) or math.isnan(flux) or math.isnan(bath_base) or math.isnan(flux_base):
                continue
            db = abs(bath - bath_base)
            df = abs(flux - flux_base)
            if db >= 1e-6 or df >= flux_tol:
                invariance["ok"] = False
                invariance["warnings"].append(
                    f"{name}: |dbath|={db:.6g}, |dflux_full|={df:.6g}, flux_tol={flux_tol:.6g}"
                )

    # Print concise summary table
    order = [
        "baseline",
        "A_imex_split",
        "B_bounds_vi",
        "B_bounds_vi_artdiff",
        "C_vlc",
        "C_quad8",
        "D_ptc_ramp_lag",
    ]
    print("Quick Strategy Benchmark (end_time <= 0.25)")
    header = (
        f"{'config':<22} {'nl_mean':>8} {'nl_max':>8} {'lin_mean':>9} {'lin_max':>8} "
        f"{'dt_min':>9} {'bath_f':>10} {'flux_full_f':>12} {'n_min_f':>10} {'wall_s':>9}"
    )
    print(header)
    print("-" * len(header))
    for name in order:
        if name not in results:
            continue
        m = results[name]
        if m.get("status") != "ok":
            print(f"{name:<22} FAILED ({m.get('status')}, rc={m.get('returncode')})")
            continue
        print(
            f"{name:<22} {_fmt(m['nonlinear_its_mean']):>8} {_fmt(m['nonlinear_its_max']):>8} "
            f"{_fmt(m['linear_its_mean']):>9} {_fmt(m['linear_its_max']):>8} "
            f"{_fmt(m['dt_min']):>9} {_fmt(m['bath_value_final']):>10} "
            f"{_fmt(m['flux_full_final']):>12} {_fmt(m['n_min_final']):>10} "
            f"{_fmt(m['wall_time_s']):>9}"
        )

    if invariance["warnings"]:
        print("\nInvariance warnings:")
        for w in invariance["warnings"]:
            print(f"  - {w}")
    else:
        print("\nInvariance checks: OK")

    rec = recommendation(results)
    print(f"\n{rec}")

    summary = {
        "exe": exe,
        "base_input": str(base),
        "common_override": str(common),
        "results": results,
        "metadata": metadata,
        "invariance": invariance,
        "recommendation": rec,
    }
    SUMMARY_JSON.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2))
    print(f"Wrote {SUMMARY_JSON}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
