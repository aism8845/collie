#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def _to_float(raw: str) -> float:
    try:
        s = str(raw).strip()
        if not s:
            return math.nan
        return float(s)
    except Exception:
        return math.nan


def _read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as f:
        text = f.read()
    if not text.strip():
        return []

    # First try normal CSV with header.
    rows = list(csv.DictReader(text.splitlines()))
    if rows and rows[0].get("time") is not None:
        return rows

    # Fallback: headerless restart segment from CSV output.
    # Observed order for the active solver_watch in this repo:
    # time, J_min_1, NumFixedPointIterations, bath_value, dt, flux_full,
    # flux_right, flux_top, linear_its, n_max, n_min, nonlinear_its, bath
    fallback_names = [
        "time",
        "J_min_1",
        "NumFixedPointIterations",
        "bath_value",
        "dt",
        "flux_full",
        "flux_right",
        "flux_top",
        "linear_its",
        "n_max",
        "n_min",
        "nonlinear_its",
        "bath",
    ]
    parsed = []
    for rec in csv.reader(text.splitlines()):
        if not rec:
            continue
        if len(rec) != len(fallback_names):
            continue
        parsed.append({k: v for k, v in zip(fallback_names, rec)})
    return parsed


def _col(rows: list[dict], key: str) -> list[float]:
    vals = []
    for r in rows:
        v = _to_float(r.get(key, ""))
        if not math.isnan(v):
            vals.append(v)
    return vals


def _closest_at(rows: list[dict], t_target: float, key: str) -> float:
    best = None
    best_dt = None
    for r in rows:
        t = _to_float(r.get("time", ""))
        v = _to_float(r.get(key, ""))
        if math.isnan(t) or math.isnan(v):
            continue
        d = abs(t - t_target)
        if best_dt is None or d < best_dt:
            best_dt = d
            best = v
    return best if best is not None else math.nan


def _find_bath_col(rows: list[dict]) -> str | None:
    if not rows:
        return None
    if "bath_value" in rows[0]:
        return "bath_value"
    if "bath" in rows[0]:
        return "bath"
    return None


def _find_flux_col(rows: list[dict]) -> str | None:
    if not rows:
        return None
    for key in ("flux_full", "flux_out", "flux_pp"):
        if key in rows[0]:
            return key
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description="Audit one solver_watch CSV.")
    ap.add_argument("solver_watch", type=Path)
    ap.add_argument("--mesh-watch", type=Path, default=None)
    ap.add_argument("--refresh-start", type=float, default=24.0)
    ap.add_argument("--refresh-duration", type=float, default=0.10)
    ap.add_argument("--refresh-tol", type=float, default=1e-8)
    ap.add_argument("--json-out", type=Path, default=None)
    args = ap.parse_args()

    rows = _read_csv(args.solver_watch)
    if not rows:
        print(json.dumps({"status": "error", "reason": "empty_solver_watch"}, indent=2))
        return 2

    bath_key = _find_bath_col(rows)
    flux_key = _find_flux_col(rows)
    if bath_key is None:
        print(json.dumps({"status": "error", "reason": "no_bath_column"}, indent=2))
        return 2

    time = _col(rows, "time")
    dt = _col(rows, "dt")
    nli = _col(rows, "nonlinear_its")
    li = _col(rows, "linear_its")
    n_min = _col(rows, "n_min")
    bath = _col(rows, bath_key)
    flux = _col(rows, flux_key) if flux_key else []

    end_time = time[-1] if time else math.nan
    end_bath = bath[-1] if bath else math.nan
    end_flux = flux[-1] if flux else math.nan
    bath_depletion = 1.0 - end_bath if not math.isnan(end_bath) else math.nan

    # Pre-refresh monotonicity (bath should not increase for t < refresh_start).
    pre = [(t, _to_float(r.get(bath_key, ""))) for r in rows for t in [_to_float(r.get("time", ""))] if not math.isnan(t)]
    pre = [(t, b) for t, b in pre if t < args.refresh_start and not math.isnan(b)]
    pre_mono = True
    if len(pre) >= 2:
        pre_mono = all(pre[i + 1][1] <= pre[i][1] + 1e-12 for i in range(len(pre) - 1))

    # Refresh window increase check (if that window is in data).
    r0 = args.refresh_start
    r1 = args.refresh_start + args.refresh_duration
    win = [(t, _to_float(r.get(bath_key, ""))) for r in rows for t in [_to_float(r.get("time", ""))]
           if not math.isnan(t) and r0 <= t <= r1]
    win = [(t, b) for t, b in win if not math.isnan(b)]
    refresh_check = "skipped"
    refresh_increase = math.nan
    refresh_ok = None
    if len(win) >= 2:
        start_b = win[0][1]
        max_b = max(b for _, b in win)
        refresh_increase = max_b - start_b
        refresh_ok = refresh_increase > args.refresh_tol
        refresh_check = "pass" if refresh_ok else "fail"

    n_min_min = min(n_min) if n_min else math.nan
    n_nonnegative = (not math.isnan(n_min_min)) and (n_min_min >= -1e-12)

    mesh_summary = {
        "mesh_watch": str(args.mesh_watch) if args.mesh_watch else None,
        "min_J_min_1": math.nan,
        "mesh_jmin_gt_0p2": None,
    }
    if args.mesh_watch and args.mesh_watch.exists():
        mrows = _read_csv(args.mesh_watch)
        mJ = _col(mrows, "J_min_1")
        if mJ:
            mesh_summary["min_J_min_1"] = min(mJ)
            mesh_summary["mesh_jmin_gt_0p2"] = mesh_summary["min_J_min_1"] > 0.2

    checks = {
        "pre_refresh_bath_monotone_nonincreasing": pre_mono,
        "refresh_window_bath_increase": refresh_check,
        "nutrient_nonnegative": n_nonnegative,
        "mesh_jmin_gt_0p2": mesh_summary["mesh_jmin_gt_0p2"],
    }

    overall = True
    for key, value in checks.items():
        if value is None or value == "skipped":
            continue
        if value is False:
            overall = False

    summary = {
        "status": "ok" if overall else "check_failed",
        "solver_watch": str(args.solver_watch),
        "rows": len(rows),
        "end_time": end_time,
        "end_bath": end_bath,
        "end_flux": end_flux,
        "bath_depletion": bath_depletion,
        "max_nonlinear_its": max(nli) if nli else math.nan,
        "max_linear_its": max(li) if li else math.nan,
        "dt_min": min(dt) if dt else math.nan,
        "dt_max": max(dt) if dt else math.nan,
        "n_min_min": n_min_min,
        "bath_at_24": _closest_at(rows, 24.0, bath_key),
        "bath_at_36": _closest_at(rows, 36.0, bath_key),
        "bath_at_48": _closest_at(rows, 48.0, bath_key),
        "refresh_window": {
            "start": r0,
            "end": r1,
            "bath_increase": refresh_increase,
            "check": refresh_check,
        },
        "checks": checks,
        "mesh": mesh_summary,
    }

    text = json.dumps(summary, indent=2)
    print(text)
    if args.json_out:
        args.json_out.write_text(text + "\n")

    return 0 if overall else 2


if __name__ == "__main__":
    raise SystemExit(main())
