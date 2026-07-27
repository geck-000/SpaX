#!/bin/bash -l
#SBATCH --job-name=basesweep_post
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --account=project_XXXXXX
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
# First-order extraction (E_x, E_z, E_z/E_x) for the base-microstructure size
# sweep: one results CSV per cell size. The decks are 2 load cases (utx, utz),
# so this is the standard first-order path rather than the tensor assembler.
unset SLURM_GTIDS
# CSC v2026_03 broke `module load abaqus` in batch (Lmod leaves the Tcl
# env-modules, which cannot parse Abaqus's .lua -> "Magic cookie missing").
# Source the snapshot of the working interactive env, as csc_solve_array.sh does.
source "$HOME/abaqus_env.sh"
WORKDIR=${WORKDIR:?set WORKDIR}
cd "$WORKDIR" || exit 1

for spec in "rve_basetensor_bs65.csv results_basesweep_L065.csv" \
            "rve_basetensor_bs10.csv results_basesweep_L100.csv"; do
    set -- $spec
    CSV=$1; RES=$2
    if [ ! -f "$CSV" ]; then
        echo "SKIP $CSV: not present"
        continue
    fi
    echo "===== $CSV -> $RES ====="
    abaqus python SpaX_PostProcess.py "$CSV" "$WORKDIR" "$RES"
done

echo "===== results written ====="
ls -1 results_basesweep_L*.csv 2>/dev/null
