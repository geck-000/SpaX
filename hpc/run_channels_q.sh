#!/bin/bash -l
set -u
cd "$(dirname "$0")/.."; ROOT="$(pwd)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SPAX_MESH_ORDER=2 SPAX_MESH_TIMEOUT=1800 SPAX_OPT_PASSES=2
export SPAX_SLIVER_START=1 SPAX_MAX_RETRIES=12 SPAX_SLIVER_MULT_Q=1.0 SPAX_GEN_WORKERS=5
echo "==== [1/5] GEN homogeneous xy baseline $(date) ===="
rm -rf homog_qxy; mkdir -p homog_qxy/logs
SPAX_RESUME=0 python -u Spatium_Standalone.py rve_homog_qxy.csv homog_qxy > homog_qxy/logs/gen.log 2>&1
echo "==== [2/5] GEN channels (xy, quadratic, all 3 load cases) $(date) ===="
rm -rf channels_q; mkdir -p channels_q/logs
SPAX_RESUME=0 python -u Spatium_Standalone.py rve_channels_q.csv channels_q > channels_q/logs/gen.log 2>&1
echo "  channel decks: $(ls channels_q/Job-*-ben.inp 2>/dev/null|wc -l)/18 ; homog: $(ls homog_qxy/Job-*-ben.inp 2>/dev/null|wc -l)/3"
echo "==== [3/5] SOLVE homog baseline (9 jobs) $(date) ===="
AUTO_POST=0 SOLVE_JOBS=2 SOLVE_CPUS=3 bash solve_local.sh homog_qxy > homog_qxy/logs/solve.log 2>&1
echo "==== [4/5] SOLVE channels (54 jobs) $(date) ===="
AUTO_POST=0 SOLVE_JOBS=2 SOLVE_CPUS=3 bash solve_local.sh channels_q > channels_q/logs/solve.log 2>&1
echo "  odbs: homog=$(ls homog_qxy/Job-*.odb 2>/dev/null|wc -l)/9 channels=$(ls channels_q/Job-*.odb 2>/dev/null|wc -l)/54"
echo "==== [5/5] POST + EXTRACT $(date) ===="
bash post_local.sh rve_homog_qxy.csv homog_qxy "$ROOT/results_homog_qxy.csv" > homog_qxy/logs/post.log 2>&1
bash post_local.sh rve_channels_q.csv channels_q "$ROOT/results_channels_q.csv" > channels_q/logs/post.log 2>&1
echo "--- f_quad(xy) baseline vs xz (should match for isotropic cube) ---" | tee "$ROOT/channels_q_extract.txt"
python Spatium_PostProcess.py analyze lengthscale results_channels_q.csv results_homog_qxy.csv 2>&1 | tee -a "$ROOT/channels_q_extract.txt"
echo "==== DONE $(date) ===="
