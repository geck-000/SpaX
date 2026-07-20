#!/bin/bash -l
# Study #8 (nlgeom): 3 slices x {linear ref, nlgeom tension, nlgeom compression},
# single-axis utx. One solve array over all 9 decks, then a per-case extraction of
# the reaction-based nominal sigma-eps path. Run on a Roihu login node.
set -e
WORKDIR=/scratch/project_XXXXXX/test_rve
mkdir -p "$WORKDIR"; cd "$WORKDIR"; mkdir -p logs

ls Job-NLG*-utx.inp 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_nlg
N=$(wc -l < GlobalJobList_nlg)
echo "nlgeom decks: $N   (expect 9 = 3 slices x 3 cases)"
[ "$N" -ge 1 ] || { echo "ERROR: no Job-NLG*-utx.inp in $WORKDIR"; exit 1; }

# nlgeom ramps take more increments than the linear solves; give 30 min on small.
SOLVE=$(sbatch --parsable \
  --partition=small --cpus-per-task=4 --mem=8G --time=00:30:00 \
  --array=1-${N}%9 \
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_nlg \
  csc_solve_array.sh)
echo "solve array: $SOLVE"

for case in lin ten cmp; do
  P=$(sbatch --parsable --dependency=afterany:${SOLVE} \
    --export=ALL,WORKDIR=$WORKDIR,CSV=rve_nlgeom_${case}.csv,SUMM=results_nlgeom_${case}.csv,CURVES=curves_nlgeom_${case}.csv \
    postprocess_nlgeom.sh)
  echo "postprocess $case: $P -> results_nlgeom_${case}.csv + curves_nlgeom_${case}.csv"
done
