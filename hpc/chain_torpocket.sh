#!/bin/bash -l
#SBATCH --job-name=torp_collect
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=small
#SBATCH --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:30:00
# The pocket-and-channel torsion sweep, regenerated so its torsional rigidity
# can actually be extracted.
#
# The decks were always right -- Bending_Plane=torsion, Kappa=0.11, six sizes,
# 33 cells plus 6 matched phi=0 controls. What went wrong the first time was
# downstream: the campaign was postprocessed through the generic first-order
# route, which returns E and nu and leaves G_eff empty, and the odbs were
# reclaimed before the omission surfaced. Nothing survives to re-extract from,
# so the meshes are rebuilt from the same decks and put through
# postprocess_torsion.sh, which is the route that calls torsion_extract.py.
set -e
WORKDIR=/scratch/project_2019020/test_rve
cd "$WORKDIR"

for d in out_torpocket out_torpocketh; do
  find $d -name 'Job-*.inp' -exec mv -t "$WORKDIR" {} + 2>/dev/null || true
  find $d -name '*_periodic_pairs.csv' -exec mv -t "$WORKDIR" {} + 2>/dev/null || true
done

ls Job-TOR_*.inp Job-TORH_*.inp 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_torpocket
N=$(wc -l < GlobalJobList_torpocket)
echo "pocket torsion decks: $N (expect 39 = 33 heterogeneous + 6 control)"
[ "$N" -ge 1 ] || { echo "ERROR: no TOR decks collected"; exit 1; }

# L=0.80 at quadratic order is the expensive member; size the whole array for it.
S=$(sbatch --parsable --account=project_2019020 \
  --partition=small --cpus-per-task=4 --mem=32G --time=01:30:00 \
  --array=1-${N}%20 \
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_torpocket \
  csc_solve_array.sh)
echo "solve: $S"

P1=$(sbatch --parsable --account=project_2019020 --dependency=afterany:${S} \
  --export=ALL,WORKDIR=$WORKDIR,CSV=rve_torsion.csv,RESULTS=results_torsion_K.csv,GLOB='Job-TOR_*-tor.odb' \
  postprocess_torsion.sh)
P2=$(sbatch --parsable --account=project_2019020 --dependency=afterany:${S} \
  --export=ALL,WORKDIR=$WORKDIR,CSV=rve_torsion_homog.csv,RESULTS=results_torsion_homog_K.csv,GLOB='Job-TORH_*-tor.odb' \
  postprocess_torsion.sh)
echo "torsion extraction: $P1 (heterogeneous) $P2 (control)"
