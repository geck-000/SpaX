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
#
# PREFIX selects the slice. The base, BTEN_z95, is the default and the one the
# paper reports; PREFIX=BTEN_z85 solves the slice above it, which carries the
# other half of the anisotropy claim -- Table 4's replicate scatter rises in the
# bottom *two* slices, and only the lowest of them has an ensemble tensor.
# Decks are generated on the workstation and rsynced in; this script solves and
# extracts only.
set -e
WORKDIR=${WORKDIR:-/scratch/project_2019020/test_rve}
PREFIX=${PREFIX:-BTEN_z95}
TENSOR_DIR=${TENSOR_DIR:-tensors/basetensor_seeds}
TAG=$(echo "$PREFIX" | tr 'A-Z' 'a-z')
mkdir -p "$WORKDIR"
cd "$WORKDIR"
mkdir -p logs

ls Job-${PREFIX}_s*-*.inp 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_${TAG}
N=$(wc -l < GlobalJobList_${TAG})
echo "basetensor (full-tensor load-case) jobs for ${PREFIX}: $N   (expect 30 = 5 packings x 6)"
[ "$N" -ge 1 ] || { echo "ERROR: no Job-${PREFIX}_s*-*.inp in $WORKDIR"; exit 1; }

ACCT=${SPAX_ACCT:-project_2019020}
SOLVE=$(sbatch --parsable --account=$ACCT \
  --partition=small --cpus-per-task=4 --mem=12G --time=01:00:00 \
  --array=1-${N}%15 \
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_${TAG} \
  csc_solve_array.sh)
echo "solve array: $SOLVE"

POST=$(sbatch --parsable --account=$ACCT --dependency=afterany:${SOLVE} \
  --export=ALL,WORKDIR=$WORKDIR,PREFIX=$PREFIX,TENSOR_DIR=$TENSOR_DIR \
  postprocess_basetensor.sh)
echo "postprocess: $POST  -> ${TENSOR_DIR}/elasticity_tensor_${PREFIX}_s*.csv"
