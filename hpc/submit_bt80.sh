#!/bin/bash -l
# Submit the full-6x6-tensor solves for the 5 warm-base replicates at L=0.80
# (BT80_z95_s1..s5), full_tensor=Yes -> 6 load cases each = 30 decks.
#
# Settles the in-plane isotropy question the L=0.50 ensemble could not. At the
# base volume fraction a 0.50 cell holds only 3-5 channels, so the two in-plane
# directions are not equivalent within one realisation and E_y/E_x scatters by
# ~1% per cell. check_channel_isotropy.py showed the generator itself is
# unbiased, so the fix is a larger cell, not more packings: L=0.80 is 2.56x the
# cross-section and holds proportionally more channels at the same L_mesh.
#
# Much heavier than the L=0.50 run: ~7.8e5 elements per deck against 1.9e5
# (cubic in L at fixed mesh size), so cpus, memory and walltime are all raised.
# Run on a Roihu CPU login node: bash submit_bt80.sh
set -e
WORKDIR=${WORKDIR:-/scratch/project_XXXXXX/test_rve}
mkdir -p "$WORKDIR"
cd "$WORKDIR"
mkdir -p logs

ls Job-BT80_z95_s*-*.inp 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_bt80
N=$(wc -l < GlobalJobList_bt80)
echo "bt80 (full-tensor load-case) jobs: $N   (expect 30 = 5 packings x 6)"
[ "$N" -ge 1 ] || { echo "ERROR: no Job-BT80_z95_s*-*.inp in $WORKDIR"; exit 1; }

# %10 rather than %15: these are large enough that packing many onto one node
# starves them of memory bandwidth (the L=0.50 run slowed ~7x when 9 landed
# together on a full node).
SOLVE=$(sbatch --parsable \
  --partition=small --cpus-per-task=8 --mem=48G --time=04:00:00 \
  --array=1-${N}%10 \
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_bt80 \
  csc_solve_array.sh)
echo "solve array: $SOLVE"

POST=$(sbatch --parsable --dependency=afterany:${SOLVE} \
  --export=ALL,WORKDIR=$WORKDIR postprocess_bt80.sh)
echo "postprocess: $POST  -> post_bt80/elasticity_tensor_BT80_z95_s*.csv"
