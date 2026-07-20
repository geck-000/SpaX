#!/bin/bash -l
# Submit the two first-order campaigns together (#3 sizechan + #4 salfamily):
# utx + utz per RVE, linear C3D4H, small/cheap. One solve array over all decks,
# then a per-study postprocess (E_x, E_z, E_z/E_x). Run on a Roihu login node.
set -e
WORKDIR=/scratch/project_XXXXXX/test_rve
mkdir -p "$WORKDIR"; cd "$WORKDIR"; mkdir -p logs

ls Job-SZCH_*.inp Job-SAL_*.inp 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_fo
N=$(wc -l < GlobalJobList_fo)
echo "first-order (sizechan + salfamily) decks: $N   (expect 80 = 40 RVEs x 2)"
[ "$N" -ge 1 ] || { echo "ERROR: no Job-SZCH_*/Job-SAL_* in $WORKDIR"; exit 1; }

# test partition: tiny linear solves (~30 s) fit the 15-min cap; faster to schedule
# and cheaper. If test rejects the array size, rerun with --partition=small.
SOLVE=$(sbatch --parsable \
  --partition=test --cpus-per-task=4 --mem=8G --time=00:15:00 \
  --array=1-${N}%20 \
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_fo \
  csc_solve_array.sh)
echo "solve array: $SOLVE"

P1=$(sbatch --parsable --dependency=afterany:${SOLVE} \
  --export=ALL,WORKDIR=$WORKDIR,CSV=rve_sizechan.csv,RESULTS=results_sizechan.csv \
  postprocess_firstorder.sh)
P2=$(sbatch --parsable --dependency=afterany:${SOLVE} \
  --export=ALL,WORKDIR=$WORKDIR,CSV=rve_salfamily.csv,RESULTS=results_salfamily.csv \
  postprocess_firstorder.sh)
echo "postprocess sizechan: $P1  -> results_sizechan.csv"
echo "postprocess salfamily: $P2 -> results_salfamily.csv"
