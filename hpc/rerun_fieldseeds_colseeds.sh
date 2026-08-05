#!/bin/bash -l
#SBATCH --job-name=rerun_ctl
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=small
#SBATCH --ntasks=1 --cpus-per-task=1 --mem=2G --time=00:20:00
# Stage $STAGE of the re-run. Kept as a chain of small controllers because the
# small partition caps submitted jobs at 200, so the 300-deck solve arrays
# cannot be queued alongside the generation arrays that produce them.
set -e
# Placeholders as in the rest of hpc/ (see hpc/README.md): override per run with
# WORKDIR / SPAX_ACCOUNT / PYTHONUSERBASE rather than editing this file.
W=${WORKDIR:-/scratch/project_XXXXXX/test_rve}
cd "$W"
A=--account=${SPAX_ACCOUNT:-project_XXXXXX}
PUB=${PYTHONUSERBASE:-/projappl/project_XXXXXX/spax_py}
# Each stage re-submits this same controller, so it must be staged in $W.
SELF=rerun_fieldseeds_colseeds.sh

collect () {  # $1 = outdir, $2 = deck glob, $3 = joblist name
  find "$1" -name 'Job-*.inp' -exec mv -t "$W" {} +
  find "$1" -name '*_periodic_pairs.csv' -exec mv -t "$W" {} + 2>/dev/null || true
  ls $2 2>/dev/null | sed 's/\.inp$//' | sort > "$3"
  wc -l < "$3"
}

case "$STAGE" in
  fs_solve)
    N=$(collect out_fs2 'Job-MSEED_*.inp Job-NSEED_*.inp Job-N2SEED_*.inp' GJ_fs2)
    echo "fieldseeds decks: $N (expect 300)"
    S1=$(sbatch --parsable $A --partition=small --cpus-per-task=4 --mem=16G \
         --time=00:30:00 --array=1-150%30 \
         --export=ALL,WORKDIR=$W,JOBLIST=GJ_fs2 csc_solve_array.sh)
    echo "fs wave1: $S1"
    sbatch $A --dependency=afterany:$S1 --export=ALL,STAGE=fs_wave2 "$SELF"
    ;;
  fs_wave2)
    N=$(wc -l < GJ_fs2)
    S2=$(sbatch --parsable $A --partition=small --cpus-per-task=4 --mem=16G \
         --time=00:30:00 --array=151-${N}%30 \
         --export=ALL,WORKDIR=$W,JOBLIST=GJ_fs2 csc_solve_array.sh)
    echo "fs wave2: $S2"
    P=$(sbatch --parsable $A --dependency=afterany:$S2 \
        --export=ALL,WORKDIR=$W,CSV=rve_fieldseeds.csv,RESULTS=results_fieldseeds.csv \
        postprocess_firstorder.sh)
    echo "fs post: $P"
    G=$(sbatch --parsable $A --dependency=afterany:$P --array=1-50%40 \
        --export=ALL,WORKDIR=$W,CSV=params/rve_colseeds_extra.csv,OUTDIR=out_cse,SPAX_SEED=20260723,PYTHONUSERBASE=$PUB \
        generate_array.sh)
    echo "cs gen: $G"
    sbatch $A --dependency=afterany:$G --export=ALL,STAGE=cs_solve "$SELF"
    ;;
  cs_solve)
    N=$(collect out_cse 'Job-CSEED_*.inp' GJ_cse)
    echo "colseeds-extra decks: $N (expect 100)"
    S=$(sbatch --parsable $A --partition=small --cpus-per-task=4 --mem=16G \
        --time=00:30:00 --array=1-${N}%30 \
        --export=ALL,WORKDIR=$W,JOBLIST=GJ_cse csc_solve_array.sh)
    echo "cs solve: $S"
    P=$(sbatch --parsable $A --dependency=afterany:$S \
        --export=ALL,WORKDIR=$W,CSV=rve_colseeds_extra.csv,RESULTS=results_colseeds_extra.csv \
        postprocess_firstorder.sh)
    echo "cs post: $P -> results_colseeds_extra.csv"
    ;;
esac
