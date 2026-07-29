#!/bin/bash -l
#SBATCH --job-name=si2nd_post
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --account=project_XXXXXX
#SBATCH --partition=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
# Extract D_rve / E_bending / E_x etc. for the 4 L=0.40 RVEs only (rows 9-12 of
# rve_seaice_2nd.csv = SI2_L400_s1..s4). Writes one partial per row, then appends
# the data rows to results_si2nd.csv IF that file (with the L240/L320 rows) was
# copied here. Eq.19 (3-size MCST slope) is best run locally on the full CSV.
unset SLURM_GTIDS
module load abaqus/2025
WORKDIR=${WORKDIR:-/scratch/project_XXXXXX/test_rve}
cd "$WORKDIR"
export PYTHONUNBUFFERED=1
CSV=rve_seaice_2nd.csv
RES=results_si2nd.csv
mkdir -p post_parts_si2nd

for i in 9 10 11 12; do
    abaqus python Spatium_PostProcess.py "$CSV" "$WORKDIR" "post_parts_si2nd/row_${i}.csv" "$i"
done

# Append L400 data rows to results_si2nd.csv (keep its header; create if absent).
if [ ! -f "$RES" ] && [ -f post_parts_si2nd/row_9.csv ]; then
    head -1 post_parts_si2nd/row_9.csv > "$RES"
fi
for i in 9 10 11 12; do
    [ -f "post_parts_si2nd/row_${i}.csv" ] && tail -n +2 "post_parts_si2nd/row_${i}.csv" >> "$RES"
done
echo "appended L400 rows -> $RES"

# Try the 3-size Eq.19 MCST slope (pandas-only; harmless if it can't run here).
python3 Spatium_PostProcess.py analyze eq19 "$RES" 2>/dev/null \
  || abaqus python Spatium_PostProcess.py analyze eq19 "$RES" \
  || echo "run Eq.19 locally: python3 Spatium_PostProcess.py analyze eq19 results_si2nd.csv"
