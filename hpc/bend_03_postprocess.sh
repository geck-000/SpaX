#!/bin/bash -l
#SBATCH --job-name=bnd_post
#SBATCH --output=logs/bnd_post_%A_%a.out
#SBATCH --error=logs/bnd_post_%A_%a.err
#SBATCH --account=project_XXXXXX
#SBATCH --cpus-per-task=1
#SBATCH --mem=120G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --partition=small
# Reads a multi-million-element bending ODB (120 GB on a standard 192 GB node is
# enough; bump --mem if the largest OOMs). Extracts D_rve per RVE; E_eff/G_eff
# come out MISSING because uniaxial/shear were never solved here (expected).
# Uses THIS study's CSV + its own post_parts_bend dir.
unset SLURM_GTIDS
module load abaqus/2025
WORKDIR=${WORKDIR:-/scratch/project_XXXXXX/test_rve}
cd "$WORKDIR" || exit 1
export PYTHONUNBUFFERED=1
CSV=rve_bending_extended.csv

if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    mkdir -p post_parts_bend
    abaqus python Spatium_PostProcess.py "$CSV" "$WORKDIR" \
        "post_parts_bend/row_${SLURM_ARRAY_TASK_ID}.csv" "${SLURM_ARRAY_TASK_ID}"
else
    abaqus python Spatium_PostProcess.py "$CSV" "$WORKDIR" results_bending.csv
fi
