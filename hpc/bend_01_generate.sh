#!/bin/bash -l
#SBATCH --job-name=bnd_gen
#SBATCH --output=logs/bnd_gen_%j.out
#SBATCH --error=logs/bnd_gen_%j.err
#SBATCH --account=project_XXXXXX
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=18:00:00
#SBATCH --partition=small
# Runs in the SHARED test_rve dir but scopes EVERY file op to this study's
# BND_L* run_ids and to its own job-list / post / results names, so it never
# clobbers the other runs already living in test_rve.
unset SLURM_GTIDS
module load python-data
export PYTHONUSERBASE=${PYTHONUSERBASE:-/projappl/project_XXXXXX/my-python-env}
export OMP_NUM_THREADS=1; export OPENBLAS_NUM_THREADS=1; export MKL_NUM_THREADS=1
export SPAX_RESUME=1
export SPAX_GEN_WORKERS=4        # few workers: each meshes a 2-5 M-element RVE
export SPAX_MAX_RETRIES=8
export SPAX_MESH_TIMEOUT=5400

WORKDIR=${WORKDIR:-/scratch/project_XXXXXX/test_rve}
cd "$WORKDIR" || exit 1
CSV=rve_bending_extended.csv
echo "CSV: $CSV   WORKDIR: $WORKDIR"

python3 Spatium_Standalone.py "$CSV" "$WORKDIR"

# Bending job list SCOPED to this study (BND_L*-ben). The -utx decks generated
# alongside are left on disk (unsolved) so SPAX_RESUME stays correct on resubmit.
ls Job-BND_L*-ben.inp 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_bend
N=$(wc -l < GlobalJobList_bend)
echo "Bending jobs: $N -> GlobalJobList_bend"
[ "$N" -lt 1 ] && { echo "ERROR: no Job-BND_L*-ben.inp generated"; exit 1; }

# Solve array — LOW concurrency (%3): big/hugemem nodes are scarce.
ARRAY_JID=$(sbatch --parsable --array=1-${N}%3 bend_02_solve.sh)
echo "Submitted solve array: $ARRAY_JID (1-$N)"
[ -z "$ARRAY_JID" ] && { echo "ERROR: solve submit failed"; exit 1; }

# One post task per RVE (CSV rows) -> partials in post_parts_bend -> merge.
N_RVE=$(python3 -c "import csv,sys;print(sum(1 for _ in csv.DictReader(open(sys.argv[1]))))" "$CSV")
mkdir -p post_parts_bend; rm -f post_parts_bend/row_*.csv 2>/dev/null
POST_JID=$(sbatch --parsable --dependency=afterany:${ARRAY_JID} --array=1-${N_RVE}%10 bend_03_postprocess.sh)
echo "Submitted postprocess array: $POST_JID (1-$N_RVE)"
[ -z "$POST_JID" ] && { echo "ERROR: post submit failed"; exit 1; }

MERGE_JID=$(sbatch --parsable --dependency=afterany:${POST_JID} bend_04_merge.sh)
echo "Submitted merge: $MERGE_JID  -> results_bending.csv"
