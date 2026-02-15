#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path


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


def write_case_override(path: Path, tag: str) -> None:
    path.write_text(
        "\n".join(
            [
                "Outputs/exodus := false",
                f"Outputs/solver_watch/file_base := outputs/{tag}_solver_watch",
                f"Outputs/mesh_watch/file_base := outputs/{tag}_mesh_watch",
                f"Outputs/chk/file_base := outputs/{tag}_chk",
                "",
            ]
        )
    )


def clean_tag_outputs(tag: str) -> None:
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


def run_cmd(cmd: list[str], log_path: Path) -> tuple[int, float]:
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    wall = time.perf_counter() - t0
    with log_path.open("a") as f:
        f.write("\n$ " + " ".join(cmd) + "\n")
        f.write(proc.stdout)
        if not proc.stdout.endswith("\n"):
            f.write("\n")
    return proc.returncode, wall


def latest_checkpoint_prefix(chk_base: Path) -> str | None:
    cp_dir = Path(str(chk_base) + "_cp")
    if not cp_dir.is_dir():
        return None

    best = None
    best_key = None
    re_step = re.compile(r"^(\d+)-restart-(\d+)\.rd$")
    for p in cp_dir.glob("*-restart-*.rd"):
        if not p.is_dir():
            continue
        if not (p / "data").exists() or not (p / "header").exists():
            continue
        m = re_step.match(p.name)
        if not m:
            continue
        step = int(m.group(1))
        part = int(m.group(2))
        key = (step, part, p.stat().st_mtime)
        if best_key is None or key > best_key:
            best_key = key
            best = p

    if best is None:
        return None

    # Recover expects the prefix without "-restart-<rank>.rd".
    return str(best).rsplit("-restart-", 1)[0]


def main() -> int:
    parser = argparse.ArgumentParser(description="Run Strang-style IMEX split (A_diff/A_rxn) with checkpoint recover.")
    parser.add_argument("--exe", default=detect_exe())
    parser.add_argument("--base", type=Path, default=detect_base_input())
    parser.add_argument("--common", type=Path, default=Path("inputs/solver_overrides/common_quick.i"))
    parser.add_argument("--diff", type=Path, default=Path("inputs/solver_overrides/A_diff.i"))
    parser.add_argument("--rxn", type=Path, default=Path("inputs/solver_overrides/A_rxn.i"))
    parser.add_argument("--macro-dt", type=float, default=1e-2)
    parser.add_argument("--end-time", type=float, default=0.25)
    parser.add_argument("--tag", default="bench_A")
    args = parser.parse_args()

    clean_tag_outputs(args.tag)
    Path("outputs").mkdir(parents=True, exist_ok=True)
    log_path = Path(f"outputs/{args.tag}_imex.log")
    if log_path.exists():
        log_path.unlink()

    chk_base = Path(f"outputs/{args.tag}_chk")
    with tempfile.NamedTemporaryFile("w", suffix=".i", prefix=f"{args.tag}_case_", delete=False) as tf:
        case_override = Path(tf.name)
    write_case_override(case_override, args.tag)

    total_wall = 0.0
    t_cur = 0.0
    recover = False
    stage_id = 0
    failed = False
    fail_reason = ""

    try:
        while t_cur < args.end_time - 1e-14:
            h = min(args.macro_dt, args.end_time - t_cur)
            macro_start = t_cur
            stage_targets = (
                ("diff_half_1", args.diff, macro_start + 0.5 * h),
                ("rxn_full", args.rxn, macro_start + h),
                ("diff_half_2", args.diff, macro_start + h),
            )
            for stage_name, stage_file, target in stage_targets:
                if target <= t_cur + 1e-16:
                    continue
                stage_id += 1
                cmd = [args.exe, "-w"]
                if recover:
                    rec = latest_checkpoint_prefix(chk_base)
                    if not rec:
                        failed = True
                        fail_reason = f"stage={stage_name} id={stage_id} no_checkpoint_found"
                        break
                    cmd.append(f"--recover={rec}")
                cmd += [
                    "-i",
                    str(args.base),
                    str(args.common),
                    str(stage_file),
                    str(case_override),
                    f"Executioner/end_time={target:.15g}",
                ]
                rc, wall = run_cmd(cmd, log_path)
                total_wall += wall
                if rc != 0:
                    failed = True
                    fail_reason = f"stage={stage_name} id={stage_id} rc={rc}"
                    break
                t_cur = target
                recover = True
            if failed:
                break

        if failed:
            # Fallback: short monolithic run using common quick settings.
            with log_path.open("a") as f:
                f.write(f"\nIMEX recover/split failed ({fail_reason}). Falling back to monolithic quick run.\n")
            cmd = [
                args.exe,
                "-w",
                "-i",
                str(args.base),
                str(args.common),
                str(case_override),
                f"Executioner/end_time={args.end_time:.15g}",
                "Executioner/petsc_options=-snes_converged_reason -ksp_converged_reason",
                "Executioner/petsc_options_iname=-snes_type -snes_linesearch_type -ksp_type -ksp_rtol -pc_type -pc_factor_mat_solver_type -snes_rtol -snes_atol",
                "Executioner/petsc_options_value=newtonls bt fgmres 1e-6 lu mumps 1e-4 1e-7",
            ]
            rc, wall = run_cmd(cmd, log_path)
            total_wall += wall
            if rc != 0:
                print(f"IMEX fallback failed (rc={rc}). See {log_path}")
                return rc
            print(
                f"IMEX fallback completed (monolithic). tag={args.tag} end_time={args.end_time:.6g} wall_s={total_wall:.3f}"
            )
            return 0

        print(f"IMEX split completed. tag={args.tag} end_time={t_cur:.6g} wall_s={total_wall:.3f}")
        return 0
    finally:
        try:
            case_override.unlink(missing_ok=True)
        except OSError:
            pass


if __name__ == "__main__":
    sys.exit(main())
