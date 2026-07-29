#!/bin/bash -l
# Unattended full quadratic porous length-scale run:
#   regen (all 3 load cases) -> solve 90 jobs (utx+ss13+ben) -> post -> extract.
# Quadratic C3D10H config with the cap-rejection + minSICN validity gate.
set -u
cd "$(dirname "$0")"
ROOT="$(pwd)"
OUT=porous_q3
CSV=rve_porous_q3.csv

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SPAX_MESH_ORDER=2 SPAX_MESH_TIMEOUT=1800 SPAX_OPT_PASSES=2
export SPAX_SLIVER_START=1 SPAX_MAX_RETRIES=12 SPAX_SLIVER_MULT_Q=1.0
export SPAX_GEN_WORKERS=5
unset SPAX_LC_FINE_MULT SPAX_MESH_ALGO2D 2>/dev/null || true

echo "==== [1/4] GENERATE (all 3 load cases, quadratic) $(date) ===="
rm -rf "$OUT"; mkdir -p "$OUT/logs"
SPAX_RESUME=0 python -u Spatium_Standalone.py "$CSV" "$OUT" > "$OUT/logs/gen.log" 2>&1
nben=$(ls "$OUT"/Job-*-ben.inp 2>/dev/null | wc -l)
nall=$(ls "$OUT"/Job-*.inp 2>/dev/null | wc -l)
echo "  decks: $nall total, $nben bending  (expect 90 / 30)"

echo "==== [2/4] SOLVE all $nall jobs (utx+ss13+ben) $(date) ===="
# solve_local solves every Job-*.inp present; bending-only stripping NOT done.
AUTO_POST=0 SOLVE_JOBS=2 SOLVE_CPUS=3 bash solve_local.sh "$OUT" > "$OUT/logs/solve.log" 2>&1
echo "  odbs: $(ls "$OUT"/Job-*.odb 2>/dev/null | wc -l)/$nall"

echo "==== [3/4] POST (fixed multi-IP processor) $(date) ===="
bash post_local.sh "$CSV" "$OUT" "$ROOT/results_porous_q3.csv" > "$OUT/logs/post.log" 2>&1
echo "  -> results_porous_q3.csv"

echo "==== [4/4] EXTRACT length scale $(date) ===="
python Spatium_PostProcess.py analyze lengthscale results_porous_q3.csv | tee "$ROOT/porous_q_extract.txt"
echo "==== DONE $(date) ===="
