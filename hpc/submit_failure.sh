#!/bin/bash -l
# Submit the strength / failure-onset uniaxial solves for the 10 column slices
# (ICE_z05..ICE_z95). Linear C3D4, ~50k elements each -> tiny/cheap. One uniaxial
# tension step per deck; failure_extract.py then gives per-slice SCF + Mohr-Coulomb
# stress-concentration percentiles. Run on the Puhti login node: bash submit_failure.sh
set -e
WORKDIR=${WORKDIR:-/scratch/project_XXXXXX/test_rve}
cd "$WORKDIR"
mkdir -p logs

ls Job-ICE_z*-utx.inp 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_failure
N=$(wc -l < GlobalJobList_failure)
echo "failure (column-slice utx) jobs: $N"
[ "$N" -ge 1 ] || { echo "ERROR: no Job-ICE_z*-utx.inp in $WORKDIR"; exit 1; }

# Trivially small: 4 cores, 12G, 30 min is generous (each solves in ~1 min).
SOLVE=$(sbatch --parsable \
  --partition=small --cpus-per-task=4 --mem=12G --time=00:30:00 \
  --array=1-${N}%10 \
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_failure \
  csc_solve_array.sh)
echo "solve array: $SOLVE"

POST=$(sbatch --parsable --dependency=afterany:${SOLVE} \
  --export=ALL,WORKDIR=$WORKDIR postprocess_failure.sh)
echo "postprocess: $POST  -> results_failure.csv"
