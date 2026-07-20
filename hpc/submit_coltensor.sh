#!/bin/bash -l
# Submit the full-6x6-tensor solves for ALL 10 column slices (CTEN_z05..z95),
# full_tensor=Yes -> 6 load cases each (utx,uty,utz,ss12,ss13,ss23) = 60 decks.
# Linear C3D4H, ~46k elements each -> small/cheap. extract_elasticity_tensor then
# assembles each slice's C_ij. Run on a Roihu CPU login node: bash submit_coltensor.sh
set -e
WORKDIR=/scratch/project_XXXXXX/test_rve
mkdir -p "$WORKDIR"
cd "$WORKDIR"
mkdir -p logs

ls Job-CTEN_z*-*.inp 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_coltensor
N=$(wc -l < GlobalJobList_coltensor)
echo "coltensor (full-tensor load-case) jobs: $N   (expect 60 = 10 slices x 6)"
[ "$N" -ge 1 ] || { echo "ERROR: no Job-CTEN_z*-*.inp in $WORKDIR"; exit 1; }

# Small: 4 cores, 8G (~4 cores on Roihu's 2 GB/core), 30 min generous per linear solve.
SOLVE=$(sbatch --parsable \
  --partition=small --cpus-per-task=4 --mem=8G --time=00:30:00 \
  --array=1-${N}%15 \
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_coltensor \
  csc_solve_array.sh)
echo "solve array: $SOLVE"

POST=$(sbatch --parsable --dependency=afterany:${SOLVE} \
  --export=ALL,WORKDIR=$WORKDIR postprocess_coltensor.sh)
echo "postprocess: $POST  -> elasticity_tensor_CTEN_z*.csv"
