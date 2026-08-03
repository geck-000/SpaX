#!/bin/bash -l
#SBATCH --job-name=spax_gen
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --account=project_XXXXXX
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
# Generate decks on the cluster, one array task per CSV row.
#
# Deck generation had been done on the workstation, on the grounds that it needs
# only Python and Gmsh and no Abaqus licence. That is true but it is not a
# reason to keep it there: generation is the slow half of this campaign, and on
# a laptop it is serial and memory-bound. Quadratic cells at L/d=8-10 need
# ~2.6 GB each, so a 16 GB machine runs one or two at a time and the large
# Eringen cells took hours; the solves they feed then finish in minutes.
#
# Here each row becomes an independent array task, so the wall time is that of
# the single slowest cell rather than the sum over all of them, and each task
# gets as much memory as it needs. numpy and gmsh are installed under
# PYTHONUSERBASE below (python3 -m pip install --user numpy gmsh).
#
# Env: WORKDIR, CSV (deck), OUTDIR. Optional: SPAX_MESH_ORDER (2 for bending),
# SPAX_SEED, PYTHONUSERBASE. Set the account with sbatch --account=<acct>, which
# overrides the placeholder directive above; see hpc/README.md.
#
#   sbatch --array=1-$(($(wc -l < params/rve_x.csv)-1)) \
#          --export=ALL,WORKDIR=$W,CSV=rve_x.csv,OUTDIR=out_x generate_array.sh
#
# SpaX_Standalone.py resolves its own seed per row, so one task per row
# reproduces exactly what a single serial run over the whole deck would build.
set -e
unset SLURM_GTIDS
export PYTHONUSERBASE=${PYTHONUSERBASE:-/projappl/project_XXXXXX/spax_py}
export PATH="$PYTHONUSERBASE/bin:$PATH"

WORKDIR=${WORKDIR:?set WORKDIR}
CSV=${CSV:?set CSV}
OUTDIR=${OUTDIR:?set OUTDIR}
cd "$WORKDIR" || exit 1
mkdir -p "$OUTDIR" logs

# Gmsh and the BLAS underneath numpy both thread by default; left alone they
# oversubscribe the 4 cores this task asked for and slow each other down.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SPAX_GEN_WORKERS=1          # one row per task; the array is the parallelism

ROW=${SLURM_ARRAY_TASK_ID}
RUN_ID=$(awk -F, -v r=$((ROW+1)) 'NR==1{for(i=1;i<=NF;i++) if($i=="run_id") c=i} NR==r{print $c}' "$CSV")
[ -z "$RUN_ID" ] && { echo "ERROR: no row $ROW in $CSV"; exit 1; }

# Each task generates into a private directory and the finished decks are moved
# up afterwards. Sharing one output directory across concurrent array tasks does
# not work: the generator writes intermediate artefacts there (the Gmsh mesh and
# the <run_id>_periodic_pairs.csv that the .inp writer consumes) and tidies them
# between rows, so parallel tasks delete each other's files and the writer fails
# with FileNotFoundError on a pairs file that existed moments earlier. The
# failures are scattered rather than reproducible, which is the signature.
TASKDIR="$OUTDIR/.task_${SLURM_ARRAY_JOB_ID}_${ROW}"
mkdir -p "$TASKDIR"

# One-row slice of the deck, so the task builds exactly its own RVE.
SLICE="$TASKDIR/row.csv"
head -1 "$CSV" > "$SLICE"
awk -v r=$((ROW+1)) 'NR==r' "$CSV" >> "$SLICE"

echo "===== row ${ROW}: ${RUN_ID}  start $(date) ====="
python3 SpaX_Standalone.py "$SLICE" "$TASKDIR/"
RC=$?

N=$(ls "$TASKDIR"/Job-*.inp 2>/dev/null | wc -l)
if [ "$N" -gt 0 ]; then
    mv -f "$TASKDIR"/Job-*.inp "$OUTDIR"/
fi
rm -rf "$TASKDIR"

echo "===== row ${ROW}: ${RUN_ID}  done $(date)  rc=${RC} ====="
echo "decks: ${N}"
[ "$N" -eq 0 ] && exit 1
exit $RC
