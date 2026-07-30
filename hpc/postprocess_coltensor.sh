#!/bin/bash -l
#SBATCH --job-name=coltensor_post
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --account=project_XXXXXX
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
# Assemble the full 6x6 effective stiffness C_ij for each of the 10 column slices
# from its 6 solved load-case ODBs. extract_elasticity_tensor writes one
# elasticity_tensor_<run_id>.csv per slice into tensors/column/. L=0.50.
unset SLURM_GTIDS
# See csc_solve_array.sh: force CSC Lmod init in non-interactive/batch context.
export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
module load abaqus/2026
WORKDIR=${WORKDIR:?set WORKDIR}
cd "$WORKDIR" || exit 1
mkdir -p tensors/column

for z in 05 15 25 35 45 55 65 75 85 95; do
    RID="CTEN_z${z}"
    # need all 6 ODBs for a full tensor; extractor skips missing columns but warn here
    n=$(ls Job-${RID}-*.odb 2>/dev/null | wc -l)
    echo "===== ${RID}: ${n}/6 ODBs present ====="
    [ "$n" -ge 1 ] || { echo "  SKIP ${RID}: no ODBs"; continue; }
    abaqus python SpaX_PostProcess.py elasticity "$WORKDIR" tensors/column 0.50 "$RID"
done

echo "===== tensors written ====="
ls -1 tensors/column/elasticity_tensor_CTEN_z*.csv 2>/dev/null
