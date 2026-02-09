#!/bin/bash
#SBATCH --job-name=yp3d5
#SBATCH --account=ucb735_peak1
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=12:00:00
#SBATCH --mem=900G
#SBATCH --output=/projects/%u/projects/collie/logs/yp3d5-%j.out
#SBATCH --error=/projects/%u/projects/collie/logs/yp3d5-%j.err

# Strict for script logic; keep conda activation safe
set -eo pipefail

BASE=/projects/$USER/projects/collie
RUN=/scratch/alpine/$USER/collie_runs/YP3D5/$SLURM_JOB_ID

mkdir -p "$RUN"
cd "$RUN"

echo "JobID: ${SLURM_JOB_ID}"
echo "Host : $(hostname)"
echo "Run  : $RUN"

# --- Conda ---
source "$HOME/miniforge/etc/profile.d/conda.sh"
conda activate "/projects/$USER/conda_envs/moose"
set -euo pipefail

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# --- Stage executable + inputs ---
cp -p "$BASE/collie-opt" .
chmod +x ./collie-opt

mkdir -p inputs/3D
cp -p "$BASE/inputs/3D/YP3D5.i" inputs/3D/

# Mesh/field file expected in CWD by your FileMeshGenerator + SolutionUserObject
cp -p "$BASE/inputs/3D/phi_ref_filter3d_debugmesh2_out.e" .

# (Optional) keep a copy in inputs too, if you want
# cp -p "$BASE/inputs/3D/phi_ref_filter3d_debugmesh2_out.e" inputs/3D/

# MUMPS out-of-core scratch
mkdir -p "$RUN/mumps_ooc"

echo "conda env: $CONDA_PREFIX"
which mpirun || true
mpirun --version || true

# IMPORTANT:
#  - Add PETSc/MUMPS out-of-core to reduce FACTOR_OUTMEMORY risk.
#  - Force errors if it truly cannot converge (better than silent garbage).
mpirun -np "$SLURM_NTASKS" ./collie-opt -i inputs/3D/YP3D5.i -w --n-threads=1 \
  -mat_mumps_icntl_22 1 \
  -mat_mumps_ooc_tmpdir "$RUN/mumps_ooc" \
  -snes_error_if_not_converged 1 \
  -ksp_error_if_not_converged 1 \
  -ksp_converged_reason

echo "Run dir: $RUN"

# --- Copy outputs back to /projects ---
OUTDIR="$BASE/results/YP3D5/$SLURM_JOB_ID"
mkdir -p "$OUTDIR"

# Nemesis produces multiple files; grab them all
cp -p ./*.e*  "$OUTDIR"/ 2>/dev/null || true
cp -p ./*.csv "$OUTDIR"/ 2>/dev/null || true
cp -p ./*.log "$OUTDIR"/ 2>/dev/null || true
cp -p inputs/3D/YP3D5.i "$OUTDIR"/ 2>/dev/null || true

echo "Copied outputs to: $OUTDIR"
