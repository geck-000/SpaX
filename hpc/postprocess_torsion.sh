#!/bin/bash -l
#SBATCH --job-name=tors_post
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
# Extract the torsional rigidity from every solved -tor ODB matching GLOB.
# Follows the common contract (see hpc/README.md): WORKDIR, CSV, RESULTS, GLOB.
set -e
unset SLURM_GTIDS
WORKDIR=${WORKDIR:-/scratch/project_XXXXXX/test_rve}
CSV=${CSV:-rve_torsion.csv}
RES=${RESULTS:-results_torsion_K.csv}
GLOB=${GLOB:-Job-TOR*-tor.odb}
cd "$WORKDIR"
mkdir -p logs

# Lmod does not work in batch here; source the snapshot of the working
# interactive environment, exactly as csc_solve_array.sh does.
source "$HOME/abaqus_env.sh"
export PYTHONUNBUFFERED=1

rm -f "$RES"
n=0
bad=0
for odb in $GLOB; do
    [ -e "$odb" ] || { echo "no ODB matching $GLOB in $WORKDIR"; exit 1; }
    rid=$(basename "$odb" -tor.odb); rid=${rid#Job-}
    # cell edge is a per-deck property; read it back rather than assuming
    Lr=$(awk -F, -v id="$rid" 'NR==1{for(i=1;i<=NF;i++){if($i=="run_id")c=i; if($i=="L")l=i}; next} $c==id{print $l; exit}' "params/$CSV" 2>/dev/null || true)
    Lr=${Lr:-0.50}
  # A truncated ODB -- an Abaqus killed mid-write, most often because the
  # filesystem filled -- opens as corrupt and takes the whole extraction down
  # with it under set -e, losing the cells that did solve. Skip it and carry on;
  # the short count in the summary is what flags the cells needing a re-solve.
    if abaqus python torsion_extract.py "$odb" "$Lr" "$rid" "$RES"; then
      n=$((n+1))
    else
      echo "  SKIPPED $rid: unreadable ODB"; bad=$((bad+1))
    fi
done
echo "extracted $n torsion ODBs ($bad unreadable) -> $RES"
