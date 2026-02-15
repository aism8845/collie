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

    return str(best).rsplit("-restart-", 1)[0]


def main() -> int:
    parser = argparse.ArgumentParser(description="Run D strategy with control-ramp, fallback to 2-stage continuation.")
    parser.add_argument("--exe", default=detect_exe())
    parser.add_argument("--base", type=Path, default=detect_base_input())
    parser.add_argument("--common", type=Path, default=Path("inputs/solver_overrides/common_quick.i"))
    parser.add_argument("--direct", type=Path, default=Path("inputs/solver_overrides/D_ptc_ramp_lag.i"))
    parser.add_argument("--stage1", type=Path, default=Path("inputs/solver_overrides/D_stage1.i"))
    parser.add_argument("--stage2", type=Path, default=Path("inputs/solver_overrides/D_stage2.i"))
    parser.add_argument("--tag", default="bench_D")
    args = parser.parse_args()

    clean_tag_outputs(args.tag)
    Path("outputs").mkdir(parents=True, exist_ok=True)
    log_path = Path(f"outputs/{args.tag}_ptc.log")
    if log_path.exists():
        log_path.unlink()

    chk_base = Path(f"outputs/{args.tag}_chk")
    with tempfile.NamedTemporaryFile("w", suffix=".i", prefix=f"{args.tag}_case_", delete=False) as tf:
        case_override = Path(tf.name)
    write_case_override(case_override, args.tag)

    total_wall = 0.0
    try:
        # Try direct ramp+lag strategy first.
        direct_cmd = [args.exe, "-w", "-i", str(args.base), str(args.common), str(args.direct), str(case_override)]
        rc, wall = run_cmd(direct_cmd, log_path)
        total_wall += wall
        if rc == 0:
            print(f"D strategy direct run completed. tag={args.tag} wall_s={total_wall:.3f}")
            return 0

        with log_path.open("a") as f:
            f.write("\nDirect D strategy failed; trying staged continuation fallback.\n")

        stage1_cmd = [args.exe, "-w", "-i", str(args.base), str(args.common), str(args.stage1), str(case_override)]
        rc1, wall1 = run_cmd(stage1_cmd, log_path)
        total_wall += wall1
        if rc1 != 0:
            print(f"D stage1 failed (rc={rc1}). See {log_path}")
            return rc1

        rec = latest_checkpoint_prefix(chk_base)
        if not rec:
            print(f"D stage2 skipped: no checkpoint found under {chk_base}_cp. See {log_path}")
            return 1

        stage2_cmd = [
            args.exe,
            "-w",
            f"--recover={rec}",
            "-i",
            str(args.base),
            str(args.common),
            str(args.stage2),
            str(case_override),
        ]
        rc2, wall2 = run_cmd(stage2_cmd, log_path)
        total_wall += wall2
        if rc2 != 0:
            print(f"D stage2 failed (rc={rc2}). See {log_path}")
            return rc2

        print(f"D strategy fallback (staged continuation) completed. tag={args.tag} wall_s={total_wall:.3f}")
        return 0
    finally:
        try:
            case_override.unlink(missing_ok=True)
        except OSError:
            pass


if __name__ == "__main__":
    sys.exit(main())
