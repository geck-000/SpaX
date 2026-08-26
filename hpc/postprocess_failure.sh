#!/bin/bash -l
#SBATCH --job-name=fail_post
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
# Run failure_extract.py over every column-slice utx ODB -> results_failure.csv.
# Mohr-Coulomb friction angle phi defaults to 30 deg (override SPAX_MC_PHI_DEG).
#
# Follows the common post-processing contract (see hpc/README.md): WORKDIR, CSV,
# RESULTS, OUTDIR. It previously read a hand-built GlobalJobList_failure and
# hard-coded its output name, which meant the staged re-run controller could not
# drive it like the other postprocess_*.sh scripts.
set -e
unset SLURM_GTIDS
WORKDIR=${WORKDIR:-/scratch/project_XXXXXX/test_rve}
CSV=${CSV:-parametric_sea_ice_column.csv}
RES=${RESULTS:-results_failure.csv}
cd "$WORKDIR"
mkdir -p logs

# Lmod does not work in batch on this cluster, so source the snapshot of the
# working interactive environment, exactly as csc_solve_array.sh does.
source "$HOME/abaqus_env.sh"

export PYTHONUNBUFFERED=1
export SPAX_MC_PHI_DEG=${SPAX_MC_PHI_DEG:-30}

rm -f "$RES"
n=0
for odb in Job-ICE_z*-utx.odb; do
    [ -e "$odb" ] || { echo "no Job-ICE_z*-utx.odb in $WORKDIR"; exit 1; }
    rid=$(basename "$odb" -utx.odb); rid=${rid#Job-}
    # Cell edge is a per-deck property; read it back from the deck rather than
    # assuming, so a deck at another box size still extracts correctly.
    Lr=$(awk -F, -v id="$rid" 'NR==1{for(i=1;i<=NF;i++){if($i=="run_id")c=i; if($i=="L")l=i}; next} $c==id{print $l; exit}' "params/$CSV" 2>/dev/null || true)
    Lr=${Lr:-0.50}
    echo "== $rid  (L=$Lr)"
    abaqus python failure_extract.py "$odb" "$Lr" "$rid" "$RES"
    n=$((n+1))
done
echo "wrote $RES from $n ODBs (phi=${SPAX_MC_PHI_DEG} deg)"
