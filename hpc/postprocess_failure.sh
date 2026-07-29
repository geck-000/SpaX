#!/bin/bash -l
#SBATCH --job-name=fail_post
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --account=project_XXXXXX
#SBATCH --partition=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
# Run failure_extract.py over every column-slice utx ODB -> results_failure.csv.
# Column-slice decks (ICE_z*) all use box size L = 0.50. Mohr-Coulomb friction
# angle phi defaults to 30 deg (override with SPAX_MC_PHI_DEG).
unset SLURM_GTIDS
module load abaqus/2025
WORKDIR=${WORKDIR:-/scratch/project_XXXXXX/test_rve}
cd "$WORKDIR"
export PYTHONUNBUFFERED=1
export SPAX_MC_PHI_DEG=${SPAX_MC_PHI_DEG:-30}
L=0.50
RES=results_failure.csv
rm -f "$RES"

while read -r JOB; do
    [ -z "$JOB" ] && continue
    id=${JOB#Job-}; id=${id%-utx}
    if [ ! -f "${JOB}.odb" ]; then echo "MISSING odb: ${JOB}"; continue; fi
    abaqus python failure_extract.py "${JOB}.odb" "$L" "$id" "$RES"
done < GlobalJobList_failure
echo "wrote $RES (phi=${SPAX_MC_PHI_DEG} deg)"
