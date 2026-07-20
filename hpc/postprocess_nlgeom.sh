#!/bin/bash -l
#SBATCH --job-name=nlg_post
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --account=project_XXXXXX
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
# Large-deformation extraction (nominal sigma-eps path from the RP reaction) for
# one nlgeom study CSV. Env: WORKDIR, CSV (rve_nlgeom_*.csv), SUMM, CURVES.
unset SLURM_GTIDS
export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
module load abaqus/2026
WORKDIR=${WORKDIR:?set WORKDIR}
CSV=${CSV:?set CSV}
SUMM=${SUMM:?set SUMM}
CURVES=${CURVES:?set CURVES}
cd "$WORKDIR" || exit 1
abaqus python nlgeom_extract.py "$CSV" "$WORKDIR" "$SUMM" "$CURVES"
echo "===== wrote $SUMM / $CURVES ====="
