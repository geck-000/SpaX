#!/bin/bash -l
# Submit the L=0.40 second-order (quadratic C3D10H) sea-ice bending jobs on Puhti.
# Small footprint: 12 decks (4 seeds x {utx,ss13,ben}), ~78k quadratic tets each.
# These are NOT the multi-million-element extended-bending study -> NO hugemem.
# Run from the SHARED test_rve dir on the Puhti login node:  bash submit_si2nd_l400.sh
set -e
WORKDIR=${WORKDIR:-/scratch/project_XXXXXX/test_rve}
cd "$WORKDIR"
mkdir -p logs

ls Job-SI2_L400_s*-utx.inp Job-SI2_L400_s*-ss13.inp Job-SI2_L400_s*-ben.inp 2>/dev/null \
  | sed 's/\.inp$//' | sort > GlobalJobList_si2nd
N=$(wc -l < GlobalJobList_si2nd)
echo "si2nd L400 jobs: $N"
[ "$N" -eq 12 ] || echo "WARN: expected 12 decks, found $N (check the .inp were copied)"
[ "$N" -ge 1 ] || { echo "ERROR: no Job-SI2_L400_*.inp in $WORKDIR"; exit 1; }

# 8 cores keeps billing == cores (Puhti standard node ~4.8 GB/core, 32G < 8*4.8).
# time=2h caps worst-case BU; quadratic ben usually finishes well under that.
SOLVE=$(sbatch --parsable \
  --partition=small --cpus-per-task=8 --mem=32G --time=02:00:00 \
  --gres=nvme:100 --array=1-${N}%6 \
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_si2nd \
  csc_solve_array.sh)
echo "solve array: $SOLVE"

POST=$(sbatch --parsable --dependency=afterany:${SOLVE} \
  --export=ALL,WORKDIR=$WORKDIR postprocess_si2nd_l400.sh)
echo "postprocess: $POST  -> post_parts_si2nd/row_{9..12}.csv (+ append if results_si2nd.csv present)"
echo "When done, pull the partials (or results_si2nd.csv) back and run Eq.19 locally:"
echo "  python3 Spatium_PostProcess.py analyze eq19 results_si2nd.csv"
