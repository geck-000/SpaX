#!/bin/bash -l
# The in-plane percolation transition, where the closure rests on one cell.
#
# Eq. (5)'s ramp runs from phi_c = 0.09 to phi_sat = 0.104 and the endpoint is
# read off a single layered cell at phi = 0.099 whose exponent came out at 0.66
# against a plateau of 0.98. Across that window E(phi) moves by a factor of
# three, and seven of the twelve Gogolaze slices sit in or beside it. Three
# decks, generated together, solved in dependency order because the first one
# decides whether the other two can be pooled with the existing LCOL cells.
#
#   rve_layerb.csv  48 cells  b swept at fixed phi -- is n independent of b?
#                             Written by make_layer_decks.py long ago and never
#                             solved. It is the prerequisite: every n in the
#                             paper comes from a cell whose b was Assur's value
#                             at the SLAB fraction while phi was read back as
#                             the realised total, ~0.019 higher, so the LCOL
#                             cells sit off the Assur curve by construction.
#                             31 of the 48 are ALREADY MESHED in out_lb from an
#                             abandoned pass and are reused, not rebuilt; run
#                             reuse_layerb.sh first. The other 17 are
#                             rve_layerb_missing.csv, and they failed for two
#                             different reasons -- see make_layerb_missing.py.
#   rve_rampn.csv   24 cells  phi = 0.092, 0.096, 0.104, 0.110 -- n(phi) through
#                             the ramp, five points with the existing 0.099
#   rve_subc.csv    18 cells  phi = 0.075, 0.082, 0.088 -- below phi_c, where
#                             the closure asserts w = 0 and nothing measures it.
#                             Includes LCOL p060's condition, which would not
#                             mesh at min_distance 0.002.
#
# The control for the one-increment solve runs on a workstation, not here; see
# analysis/local_control.md and the note in collect_ramp.sh below.
#
# Build the last two first, on the workstation or here:
#
#   python3 hpc/make_ramp_decks.py params/
#
# Every cell in the two new decks runs 0.005-0.0061 element size on a 0.5 cube,
# so 0.55-1.0 M cells to the edge count and several million tets. That is two to
# four times the LCOL cells and is why the solve asks for 32G rather than the
# 12G those used. The generation array is the slow half regardless.
#
# The wider min_distance is the fix for the sliver that killed LCOL p060: at
# these layer thicknesses a pocket straddling a plane leaves a facet Gmsh
# reports as an overlapping boundary and cannot repair. Expect a few rows to
# need retries anyway; verify_solves.sh lists them.
set -e
WORKDIR=${WORKDIR:-/scratch/project_2019020/test_rve}
cd "$WORKDIR"; mkdir -p logs

ACCT=${SPAX_ACCT:-project_2019020}
MAXARR=${SPAX_MAX_ARRAY:-25}
# SpaX_Standalone.py:2402 gives the mesh subprocess 900s and calls the overrun
# "degenerate geometry". For five of the LAYERB cells that diagnosis was wrong:
# they are the b = 0.10 cells, where mesh_for hits its 0.005 floor and the cell
# carries a million elements to the edge, and the mesher simply wanted longer.
# The RAMP and SUBC cells are the same size -- the pilot spent the better part
# of an hour meshing one -- so the whole campaign gets the longer fuse. It costs
# nothing on a cell that does not need it, and the generation array's own
# --time=04:00:00 still bounds the task.
MESH_TIMEOUT=${SPAX_MESH_TIMEOUT:-5400}
export PYTHONUSERBASE=/projappl/project_2019020/spax_py

# Linear elements: the volume-averaged constant-strain measure is exact for
# them, and this is a first-order homogenisation exactly as LCOL was. Anything
# else would make the new exponents incomparable with the four they extend.
#
# SPAX_LINEAR_ONE_STEP IS NOT USED HERE, and the reason is worth recording
# because it cost a full campaign of solves to learn.
#
# The flag collapses the ten 0.1 increments to one, on the argument that the
# cells are linear and the extractor reads only the last frame. The first half
# is true -- a workstation control put one increment against ten on identical
# meshes and the reaction forces agreed to every digit printed. The second half
# is false. extract_first_order, which is what produces E_x and E_z, does not
# read the last frame: it walks the frame series and fits sigma against epsilon
# by polyfit over the 10-40% window of peak strain. With one increment there is
# a single point, at 100% of peak, the window is empty, and E comes back exactly
# zero. (SpaX_PostProcess.py:1381, last_frame_only=True, is extract_principals
# -- a different route, used for the full-tensor and bending paths.)
#
# Ten increments is not arbitrary either: at ten, strains land at 0.1..1.0 of
# peak and the 10-40% window captures four points. Four increments would capture
# one, which is no better than one.
#
# The lesson for the control: it compared reaction forces with a purpose-written
# reader and so verified the SOLVE, never the extraction. A control for a
# pipeline change has to run the pipeline.
gen () {   # $1 deck  $2 outdir
  local N=$(($(wc -l < "$1") - 1))
  sbatch --parsable --account=$ACCT \
    --array=1-${N}%${MAXARR} \
    --export=ALL,WORKDIR=$WORKDIR,CSV=$1,OUTDIR=$2,SPAX_MESH_ORDER=1,SPAX_SEED=20260819,SPAX_MESH_TIMEOUT=$MESH_TIMEOUT,PYTHONUSERBASE=$PYTHONUSERBASE \
    generate_array.sh
}

# out_lb's 31 cells are moved into place and collapsed to one increment before
# anything is submitted, so the generation array below only covers what is
# genuinely absent.
./reuse_layerb.sh "$WORKDIR"

GB=$(gen rve_layerb_missing.csv out_layerb)
GR=$(gen rve_rampn.csv  out_rampn)
GS=$(gen rve_subc.csv   out_subc)
echo "gen: $GB (7 layerb, 17 drained reused) $GR (12 ramp) $GS (9 subc)"

cat > collect_ramp.sh <<'EOS'
#!/bin/bash -l
#SBATCH --job-name=ramp_collect
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=small
#SBATCH --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:30:00
set -e
WORKDIR=${WORKDIR:-/scratch/project_2019020/test_rve}
ACCT=${SPAX_ACCT:-project_2019020}
cd "$WORKDIR"
for d in out_layerb out_rampn out_subc; do
  find $d -name 'Job-*.inp' -exec mv -t "$WORKDIR" {} + 2>/dev/null || true
  find $d -name '*_periodic_pairs.csv' -exec mv -t "$WORKDIR" {} + 2>/dev/null || true
done

# The control for the single-increment solve is NOT run here. The claim being
# tested -- that one increment and ten return the same moduli -- is a property
# of the equations and not of the discretisation, so a coarse copy of the same
# cell settles it as well as a production one and settles it on a workstation
# for nothing. Run control_local.sh there and read it with
# analysis/local_control_compare.py. Running it at full size here would have
# bought the same yes/no answer at six cells of billing.

sub () {   # $1 tag  $2 deck  $3 results  $4 mem  $5 time  [$6 glob]
  # $6 lets a tag solve a subset of what it generated. LAYERB was meshed in both
  # drainage states before the campaign narrowed to drained, so its undrained
  # decks stay on disk -- they cost nothing to keep and would cost a regeneration
  # to recover -- but they are not solved, because nothing reads them.
  ls ${6:-Job-${1}_*.inp} 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_${1}
  local N=$(wc -l < GlobalJobList_${1})
  echo "${1}: $N decks"
  [ "$N" -ge 1 ] || { echo "  WARNING: none generated -- check logs/spax_gen_*"; return; }
  local S=$(sbatch --parsable --account=$ACCT \
    --partition=small --cpus-per-task=4 --mem=$4 --time=$5 \
    --array=1-${N}%20 \
    --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_${1} \
    csc_solve_array.sh)
  echo "  solve: $S"
  local P=$(sbatch --parsable --account=$ACCT --dependency=afterany:${S} \
    --export=ALL,WORKDIR=$WORKDIR,CSV=$2,RESULTS=$3 \
    postprocess_firstorder.sh)
  echo "  post:  $P -> $3"
}

# layerb spans 0.005 to 0.012 element size, so its largest members are the same
# size as the ramp cells and it gets the same allocation.
sub LB    rve_layerb.csv  results_layerb.csv  240G 04:00:00 'Job-LB_*_drn_*.inp'
sub RAMP  rve_rampn.csv   results_rampn.csv   240G 04:00:00
sub SUBC  rve_subc.csv    results_subc.csv    240G 04:00:00
EOS
chmod +x collect_ramp.sh

C=$(sbatch --parsable --account=$ACCT \
  --dependency=afterany:${GB}:${GR}:${GS} \
  --export=ALL,WORKDIR=$WORKDIR,SPAX_ACCT=$ACCT collect_ramp.sh)
echo "collect+solve+post chained as $C"

cat <<'EOF'

Then, offline:

  python3 analysis/ramp_exponent.py

which pools results_rampn.csv and results_subc.csv with the existing
results_layercol.csv, refits n(phi) through the window, and reports what the
ramp width and phi_sat become. It checks results_layerb.csv first and refuses
to pool if n drifts with b, since that would invalidate the pooling rather than
merely widen the band.
EOF
