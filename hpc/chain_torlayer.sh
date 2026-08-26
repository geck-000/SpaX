#!/bin/bash -l
#SBATCH --job-name=torl_collect
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=small
#SBATCH --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:30:00
# Collect the layered torsion cells and their matched phi=0 controls, solve
# them together, then extract with torsion_extract.py.
#
# The extraction is the part that matters. The first torsion campaign was
# postprocessed with the generic first-order route, which returns E and nu and
# leaves G_eff MISSING -- so the sweep produced no shear modulus at all and the
# odbs were reclaimed before anyone noticed. postprocess_torsion.sh is the one
# that runs torsion_extract.py, and it is what both decks get here.
set -e
WORKDIR=/scratch/project_2019020/test_rve
cd "$WORKDIR"

for d in out_torlayer out_torlayerh; do
  find $d -name 'Job-*.inp' -exec mv -t "$WORKDIR" {} + 2>/dev/null || true
  find $d -name '*_periodic_pairs.csv' -exec mv -t "$WORKDIR" {} + 2>/dev/null || true
done

ls Job-TORL_*.inp Job-TORLH_*.inp 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_torlayer
N=$(wc -l < GlobalJobList_torlayer)
echo "layered torsion decks: $N (expect 20 = 15 layered + 5 control)"
[ "$N" -ge 1 ] || { echo "ERROR: no TORL decks collected"; exit 1; }

S=$(sbatch --parsable --account=project_2019020 \
  --partition=small --cpus-per-task=4 --mem=32G --time=01:30:00 \
  --array=1-${N}%20 \
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_torlayer \
  csc_solve_array.sh)
echo "solve: $S"

P1=$(sbatch --parsable --account=project_2019020 --dependency=afterany:${S} \
  --export=ALL,WORKDIR=$WORKDIR,CSV=rve_torsion_layer.csv,RESULTS=results_torsion_layer_K.csv,GLOB='Job-TORL_*.odb' \
  postprocess_torsion.sh)
P2=$(sbatch --parsable --account=project_2019020 --dependency=afterany:${S} \
  --export=ALL,WORKDIR=$WORKDIR,CSV=rve_torsion_layer_homog.csv,RESULTS=results_torsion_layer_homog_K.csv,GLOB='Job-TORLH_*.odb' \
  postprocess_torsion.sh)
echo "torsion extraction: $P1 (layered) $P2 (control)"
