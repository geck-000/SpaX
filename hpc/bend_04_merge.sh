#!/bin/bash -l
#SBATCH --job-name=bnd_merge
#SBATCH --output=logs/bnd_merge_%j.out
#SBATCH --error=logs/bnd_merge_%j.err
#SBATCH --account=project_XXXXXX
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=00:15:00
#SBATCH --partition=small
# Union this study's per-RVE partials into results_bending.csv. Pure Python.
unset SLURM_GTIDS
module load python-data
export PYTHONUSERBASE=${PYTHONUSERBASE:-/projappl/project_XXXXXX/my-python-env}
WORKDIR=${WORKDIR:-/scratch/project_XXXXXX/test_rve}
cd "$WORKDIR" || exit 1
export PYTHONUNBUFFERED=1
python3 Spatium_PostProcess.py --merge post_parts_bend results_bending.csv
echo "Merged -> $WORKDIR/results_bending.csv"
echo "Next (with the small-RVE first-order results.csv):"
echo "  python3 Spatium_PostProcess.py analyze hybrid results.csv results_bending.csv"
