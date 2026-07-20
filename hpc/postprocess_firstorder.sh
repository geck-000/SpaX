#!/bin/bash -l
#SBATCH --job-name=fo_post
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --account=project_XXXXXX
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
# Standard first-order extraction (E_x, E_z, E_z/E_x ...) for one study CSV.
# Env: WORKDIR, CSV (rve_*.csv), RESULTS (results_*.csv).
unset SLURM_GTIDS
export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
module load abaqus/2026
WORKDIR=${WORKDIR:?set WORKDIR}
CSV=${CSV:?set CSV}
RESULTS=${RESULTS:?set RESULTS}
cd "$WORKDIR" || exit 1
abaqus python SpaX_PostProcess.py "$CSV" "$WORKDIR" "$RESULTS"
echo "===== wrote $RESULTS ====="
