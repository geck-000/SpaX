#!/bin/bash -l
# Base-microstructure size sweep: does the warm-base modulus converge with cell
# size, and where?
#
# The full-tensor ensembles gave E_x = 4.79 GPa at L=0.50 and 2.94 GPa at
# L=0.80 -- a 39% drop with the soft-phase fraction matched to ~2 points, so it
# is not a composition difference. At ~30% soft phase the base sits near the
# continuum percolation threshold, where a small cell suppresses connectivity
# and holds the stiffness up. The published box-size sweep shows E_x flat over
# the same L range but at a third of this soft fraction, so it does not settle
# the question here.
#
# This fills in L=0.65 and L=1.00 (5 packings each, first-order utx+utz only --
# E_x and E_z are all the convergence curve needs). Combined with the existing
# L=0.50 and L=0.80 ensembles that gives four points.
#
# The two sizes are submitted as separate arrays because they are very
# different jobs: ~3.7e5 elements at L=0.65 against ~1.35e6 at L=1.00.
# Run on a Roihu CPU login node: bash submit_basesweep.sh
set -e
WORKDIR=/scratch/project_XXXXXX/test_rve
mkdir -p "$WORKDIR"
cd "$WORKDIR"
mkdir -p logs

ls Job-BS65_z95_s*-*.inp 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_bs65
ls Job-BS10_z95_s*-*.inp 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_bs10
N65=$(wc -l < GlobalJobList_bs65)
N10=$(wc -l < GlobalJobList_bs10)
echo "L=0.65 jobs: $N65   L=1.00 jobs: $N10   (expect 10 each = 5 packings x 2)"
[ "$N65" -ge 1 ] || { echo "ERROR: no Job-BS65_*.inp in $WORKDIR"; exit 1; }
[ "$N10" -ge 1 ] || { echo "ERROR: no Job-BS10_*.inp in $WORKDIR"; exit 1; }

S65=$(sbatch --parsable \
  --partition=small --cpus-per-task=8 --mem=24G --time=01:00:00 \
  --array=1-${N65}%10 \
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_bs65 \
  csc_solve_array.sh)
echo "L=0.65 solve array: $S65"

# L=1.00 is ~2x the L=0.80 job, which ran comfortably on 8 cpus / 48G; double
# both and give it 6 h, since a timeout here would leave a truncated ODB that
# csc_solve_array.sh would then SKIP on resubmit.
S10=$(sbatch --parsable \
  --partition=small --cpus-per-task=16 --mem=96G --time=06:00:00 \
  --array=1-${N10}%5 \
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_bs10 \
  csc_solve_array.sh)
echo "L=1.00 solve array: $S10"

POST=$(sbatch --parsable --dependency=afterany:${S65}:${S10} \
  --export=ALL,WORKDIR=$WORKDIR postprocess_basesweep.sh)
echo "postprocess: $POST  -> results_basesweep_L065.csv, results_basesweep_L100.csv"
