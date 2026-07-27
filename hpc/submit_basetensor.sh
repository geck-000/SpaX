#!/bin/bash -l
# Submit the full-6x6-tensor solves for the 5 warm-base replicates
# (BTEN_z95_s1..s5), full_tensor=Yes -> 6 load cases each = 30 decks.
#
# Settles whether the in-plane split seen in the single base cell (E_x 4.85 vs
# E_y 5.02 GPa) is a realisation effect or a genuine in-plane texture: with five
# packings, E_y/E_x becomes an ensemble statement instead of one draw.
#
# Same shape as submit_coltensor.sh but ~250k elements per deck (the base slice
# is the densest in the column), so the walltime is raised to 1 h.
# Run on a Roihu CPU login node: bash submit_basetensor.sh
set -e
WORKDIR=/scratch/project_XXXXXX/test_rve
mkdir -p "$WORKDIR"
cd "$WORKDIR"
mkdir -p logs

ls Job-BTEN_z95_s*-*.inp 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_basetensor
N=$(wc -l < GlobalJobList_basetensor)
echo "basetensor (full-tensor load-case) jobs: $N   (expect 30 = 5 packings x 6)"
[ "$N" -ge 1 ] || { echo "ERROR: no Job-BTEN_z95_s*-*.inp in $WORKDIR"; exit 1; }

SOLVE=$(sbatch --parsable \
  --partition=small --cpus-per-task=4 --mem=12G --time=01:00:00 \
  --array=1-${N}%15 \
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_basetensor \
  csc_solve_array.sh)
echo "solve array: $SOLVE"

POST=$(sbatch --parsable --dependency=afterany:${SOLVE} \
  --export=ALL,WORKDIR=$WORKDIR postprocess_basetensor.sh)
echo "postprocess: $POST  -> post_basetensor/elasticity_tensor_BTEN_z95_s*.csv"
