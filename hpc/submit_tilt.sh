#!/bin/bash -l
# Channel-tilt study (#6): straight (0) vs wavy-inclined (15, 30 deg) channels.
# First-order utx+utz per RVE, 4 seeds per tilt. One solve array over all decks,
# then a per-tilt postprocess (E_x, E_z, E_z/E_x). ODBs are KEPT (not cleaned) so
# a couple can be exported to VTK for ParaView figures. Run on a Roihu login node.
set -e
WORKDIR=/scratch/project_XXXXXX/test_rve
mkdir -p "$WORKDIR"; cd "$WORKDIR"; mkdir -p logs

ls Job-TLT00_*.inp Job-TLT15_*.inp Job-TLT30_*.inp 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_tilt
N=$(wc -l < GlobalJobList_tilt)
echo "tilt decks: $N   (expect 24 = 3 tilts x 4 seeds x 2 modes)"
[ "$N" -ge 1 ] || { echo "ERROR: no Job-TLT*_*.inp in $WORKDIR"; exit 1; }

SOLVE=$(sbatch --parsable \
  --partition=small --cpus-per-task=4 --mem=8G --time=00:20:00 \
  --array=1-${N}%24 \
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_tilt \
  csc_solve_array.sh)
echo "solve array: $SOLVE"

for tag in 00 15 30; do
  P=$(sbatch --parsable --dependency=afterany:${SOLVE} \
    --export=ALL,WORKDIR=$WORKDIR,CSV=rve_tilt${tag}.csv,RESULTS=results_tilt${tag}.csv \
    postprocess_firstorder.sh)
  echo "postprocess tilt${tag}: $P -> results_tilt${tag}.csv"
done
