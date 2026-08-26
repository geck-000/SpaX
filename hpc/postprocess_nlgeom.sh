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
# Common post-processing contract (see hpc/README.md): every postprocess_*.sh
# takes WORKDIR, CSV, RESULTS and OUTDIR, and uses whichever of them it needs.
# Anything extra a script wants is derived from RESULTS rather than demanded
# separately, so one caller can drive them all. SUMM/CURVES may still be set
# explicitly to override.
WORKDIR=${WORKDIR:?set WORKDIR}
CSV=${CSV:?set CSV}
SUMM=${SUMM:-${RESULTS:?set RESULTS (or SUMM)}}
CURVES=${CURVES:-curves_${SUMM#results_}}
cd "$WORKDIR" || exit 1
abaqus python nlgeom_extract.py "$CSV" "$WORKDIR" "$SUMM" "$CURVES"
echo "===== wrote $SUMM / $CURVES ====="
