#!/bin/bash -l
# Solve the seeded first-year C-shape column (statistical replicate campaign,
# rve_colseeds.csv): 10 depth slices x 5 packings x 2 load cases (utx, utz),
# linear C3D4H. One solve array over all Job-CSEED_* decks, then a first-order
# postprocess (E_x, E_z, E_z/E_x per replicate) -> results_colseeds.csv, from
# which analysis/make_rev_figs.py draws the per-slice mean +/- spread envelopes.
# Run on a Roihu login node from the staging dir. Mirrors submit_firstorder.sh.
set -e
WORKDIR=${WORKDIR:-/scratch/project_XXXXXX/test_rve}
cd "$WORKDIR"; mkdir -p logs

ls Job-CSEED_*.inp 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_cseeds
N=$(wc -l < GlobalJobList_cseeds)
echo "colseeds decks: $N   (expect 100 = 10 slices x 5 seeds x 2 loads)"
[ "$N" -ge 1 ] || { echo "ERROR: no Job-CSEED_*.inp in $WORKDIR"; exit 1; }

# test partition: tiny linear solves (~30 s) fit the 15-min cap; cheap & fast to
# schedule. If the array size is rejected, rerun with --partition=small.
SOLVE=$(sbatch --parsable \
  --partition=test --cpus-per-task=4 --mem=8G --time=00:15:00 \
  --array=1-${N}%20 \
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_cseeds \
  csc_solve_array.sh)
echo "solve array: $SOLVE"

POST=$(sbatch --parsable --dependency=afterany:${SOLVE} \
  --export=ALL,WORKDIR=$WORKDIR,CSV=rve_colseeds.csv,RESULTS=results_colseeds.csv \
  postprocess_firstorder.sh)
echo "postprocess: $POST  -> results_colseeds.csv"
