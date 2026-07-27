#!/bin/bash -l
#SBATCH --job-name=basetensor_post
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --account=project_XXXXXX
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
# Assemble the full 6x6 effective stiffness C_ij for each of the 5 base-slice
# replicates from its 6 solved load-case ODBs. extract_elasticity_tensor writes
# one elasticity_tensor_<run_id>.csv per packing into post_basetensor/. L=0.50.
unset SLURM_GTIDS
# CSC software stack v2026_03 broke `module load abaqus` in batch (the Lmod dance
# leaves the Tcl env-modules, which cannot parse Abaqus's .lua -> "Magic cookie
# missing"). Source the snapshot of the working interactive env instead, exactly
# as csc_solve_array.sh does.
source "$HOME/abaqus_env.sh"
WORKDIR=${WORKDIR:?set WORKDIR}
cd "$WORKDIR" || exit 1
mkdir -p post_basetensor

for s in 1 2 3 4 5; do
    RID="BTEN_z95_s${s}"
    n=$(ls Job-${RID}-*.odb 2>/dev/null | wc -l)
    echo "===== ${RID}: ${n}/6 ODBs present ====="
    [ "$n" -ge 1 ] || { echo "  SKIP ${RID}: no ODBs"; continue; }
    abaqus python SpaX_PostProcess.py elasticity "$WORKDIR" post_basetensor 0.50 "$RID"
done

echo "===== tensors written ====="
ls -1 post_basetensor/elasticity_tensor_BTEN_z95_s*.csv 2>/dev/null
