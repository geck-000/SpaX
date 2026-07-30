#!/bin/bash -l
# Skeletal basal layer campaign (studies/make_skeletal.py):
#
#   rve_skeletal.csv          27 RVEs x 2 loads  E(phi_b) into the high-porosity
#                                                regime, 2 morphologies
#   rve_skeletal_laminae.csv  12 RVEs x 2 loads  the four sub-laminae of the
#                                                bottom 5%, utx + ss12 for CLT
#   rve_steep_column.csv      42 RVEs x 2 loads  steeply monotonic column carried
#                                                down to a resolved skeletal base
#
# The last of these is the discriminating run. analysis/skeletal_clt.py shows
# that inserting a skeletal lamina under the existing column moves the neutral
# plane only from z0/H=0.466 to 0.452, against a measured 0.37-0.39: a thin
# compliant layer at the plate surface cannot shift the centroid of E(z). What
# controls the neutral plane is the profile SHAPE, so the test is a steeper
# interior, not merely a softer base.
#
# Generate first, from the staging dir (linear elements are exact for the
# volume-averaged constant-strain measure, so order 1 is right here):
#
#   SPAX_SEED=20260730 OMP_NUM_THREADS=1 python3 SpaX_Standalone.py \
#       params/rve_skeletal.csv out_skeletal/
#   ... likewise for the other two decks ...
#
# High brine fractions are the hard case for the packer: phi_b=0.50 with
# channel-dominated morphology puts ~17 lamellae at 42% areal fraction in the
# cell. SPAX_GAP_REFINE and the off-axis sliver repair are on by default and
# are what keep the ligaments meshable; expect some rows to need retries.
set -e
WORKDIR=${WORKDIR:-/scratch/project_XXXXXX/test_rve}
cd "$WORKDIR"; mkdir -p logs

for spec in "SKEL:rve_skeletal.csv:results_skeletal.csv" \
            "SKLM:rve_skeletal_laminae.csv:results_skeletal_laminae.csv" \
            "STEEP:rve_steep_column.csv:results_steep_column.csv"; do
  TAG=${spec%%:*}; rest=${spec#*:}; CSV=${rest%%:*}; RES=${rest#*:}

  ls Job-${TAG}_*.inp 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_${TAG}
  N=$(wc -l < GlobalJobList_${TAG})
  echo "${TAG}: $N decks   (deck ${CSV})"
  [ "$N" -ge 1 ] || { echo "  WARNING: no Job-${TAG}_*.inp -- generate first"; continue; }

  # High-porosity cells carry more elements than the column slices; give them
  # the small partition rather than the 15-minute test cap.
  SOLVE=$(sbatch --parsable \
    --partition=small --cpus-per-task=4 --mem=12G --time=01:00:00 \
    --array=1-${N}%20 \
    --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_${TAG} \
    csc_solve_array.sh)
  echo "  solve array: $SOLVE"

  POST=$(sbatch --parsable --dependency=afterany:${SOLVE} \
    --export=ALL,WORKDIR=$WORKDIR,CSV=${CSV},RESULTS=${RES} \
    postprocess_firstorder.sh)
  echo "  postprocess: $POST  -> ${RES}"
done

cat <<'EOF'

Then, offline (no Abaqus licence needed):

  cd results
  python3 ../analysis/skeletal_clt.py results_column.csv \
          --skeletal results_skeletal_laminae.csv
  python3 ../analysis/skeletal_clt.py results_steep_column.csv --probe

The question each answers: does resolving the base, or steepening the profile,
bring the neutral plane into the 0.37-0.39 band measured by Kujala et al. (1990)?
EOF
