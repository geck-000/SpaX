#!/bin/bash
# Generate, solve and read the (b, t) grid on a workstation.
#
# The grid asks whether the closure's bridge factor is a function of b alone.
# LAYERB said no, but its sweep moved t by 29% while it moved b, because
# phi_slab = t(1-b) ties them; this deck samples the plane instead of walking one
# line across it. See hpc/make_bt_grid.py.
#
# Local because the cells are small -- 72k to 96k elements to the edge against
# 140k-350k for RAMP -- and because nothing here needs to be compared against a
# cell solved elsewhere. Every point on the grid carries the same size bias and
# the question is the shape of the surface, not its level.
#
# One increment. Safe now that extract_first_order has the single-point branch,
# and verified on the cluster: a one-increment ODB there extracted E = 1.713 GPa
# where the unpatched extractor returned exactly 0.0. Ten increments would cost
# roughly 2.2x the solve for a number the analytic argument says is identical --
# for an exactly linear response the least-squares slope over the 10-40% window
# equals the secant through any single point.
#
# Only the x load case is solved. n is a drained transverse quantity; nothing
# here reads E_z, and skipping it halves the solve.
#
#   ./run_btgrid_local.sh [workdir]
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=${1:-${CLAUDE_JOB_DIR:-/tmp}/btgrid}
mkdir -p "$WORK"

ABQ=${ABAQUS_CMD:-/c/SIMULIA/Commands/abaqus.bat}
CPUS=${ABQ_CPUS:-4}
# Three generator workers rather than one: these are order-1 cells at ~1.4M tets,
# far lighter than the order-2 cells that made a low worker count necessary.
WORKERS=${SPAX_GEN_WORKERS:-3}

cd "$WORK"
cp -f "$HERE/params/rve_btgrid.csv" .

if [ "$(ls Job-BT_*.inp 2>/dev/null | wc -l)" -lt 18 ]; then
  echo "=== generating (18 cells, $WORKERS workers) ==="
  OMP_NUM_THREADS=1 SPAX_LINEAR_ONE_STEP=1 SPAX_MESH_TIMEOUT=5400 \
    SPAX_GEN_WORKERS=$WORKERS \
    python "$HERE/SpaX_Standalone.py" rve_btgrid.csv . > gen.log 2>&1 || true
fi
echo "generated: $(ls Job-BT_*-utx.inp 2>/dev/null | wc -l) x-decks of 18"

echo "=== solving (utx only) ==="
for f in Job-BT_*-utx.inp; do
  [ -e "$f" ] || continue
  j="${f%.inp}"
  [ -e "${j}.odb" ] && { echo "  skip $j"; continue; }
  s=$(date +%s)
  "$ABQ" job="$j" input="$f" cpus=$CPUS interactive > "${j}.solvelog" 2>&1 || true
  e=$(date +%s)
  echo "  $j  $((e-s))s  odb=$([ -e ${j}.odb ] && echo yes || echo NO)"
done

echo "=== reading ==="
"$ABQ" python "$HERE/analysis/bt_surface.py" "$WORK"
