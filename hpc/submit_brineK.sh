#!/bin/bash -l
# Study #5 (brine K(T)): paired first-order column, Kconst vs Ktemp.
# utx + utz per RVE, linear C3D4H hybrid, tiny/cheap. One solve array over both
# deck sets, then a per-study postprocess (E_x, E_z, E_z/E_x).
# Run on a Roihu login node from the scratch dir.
set -e
WORKDIR=/scratch/project_XXXXXX/test_rve
mkdir -p "$WORKDIR"; cd "$WORKDIR"; mkdir -p logs

ls Job-BKC_*.inp Job-BKT_*.inp 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_bk
N=$(wc -l < GlobalJobList_bk)
echo "brine-K decks: $N   (expect 40 = 20 RVEs x 2 modes)"
[ "$N" -ge 1 ] || { echo "ERROR: no Job-BKC_*/Job-BKT_* in $WORKDIR"; exit 1; }

# small partition: 40-task array exceeds the test-partition submit limit
# (AssocMaxSubmitJobLimit); small takes it at identical BU cost. Tiny linear
# solves (~30 s) still fit well under the walltime.
SOLVE=$(sbatch --parsable \
  --partition=small --cpus-per-task=4 --mem=8G --time=00:15:00 \
  --array=1-${N}%20 \
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_bk \
  csc_solve_array.sh)
echo "solve array: $SOLVE"

P1=$(sbatch --parsable --dependency=afterany:${SOLVE} \
  --export=ALL,WORKDIR=$WORKDIR,CSV=rve_brineKconst.csv,RESULTS=results_brineKconst.csv \
  postprocess_firstorder.sh)
P2=$(sbatch --parsable --dependency=afterany:${SOLVE} \
  --export=ALL,WORKDIR=$WORKDIR,CSV=rve_brineKtemp.csv,RESULTS=results_brineKtemp.csv \
  postprocess_firstorder.sh)
echo "postprocess Kconst: $P1  -> results_brineKconst.csv"
echo "postprocess Ktemp:  $P2  -> results_brineKtemp.csv"
