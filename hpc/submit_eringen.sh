#!/bin/bash -l
# Nonlocal length-scale campaign (studies/make_eringen.py):
#
#   rve_eringen.csv        33 RVEs   channelled base at L/d = 3,4,5,6,8,10
#   rve_eringen_homog.csv   6 RVEs   the same cells with VoF=0 (geometric baseline)
#
# Each RVE yields a utx deck (first-order reference E_eff) and, through
# Kappa>0, a bending deck, so the apparent-modulus ratio is formed per cell
# rather than against one global plate modulus.
#
# MUST be generated with quadratic elements. Linear tetrahedra lock in bending
# and the whole measurement is a few-percent modulus difference:
#
#   SPAX_MESH_ORDER=2 SPAX_SEED=20260730 OMP_NUM_THREADS=1 \
#       python3 SpaX_Standalone.py params/rve_eringen.csv out_eringen/
#   SPAX_MESH_ORDER=2 SPAX_SEED=20260730 \
#       python3 SpaX_Standalone.py params/rve_eringen_homog.csv out_eringen/
#
# The homogeneous cells are the control the existing calibration lacks: they
# have no microstructure and therefore no length scale of any kind, so whatever
# size dependence they show is cube-versus-plate kinematics and discretisation.
# Matched sizes let that be subtracted point-by-point instead of interpolated.
#
# L=0.80 is the practical ceiling: the L=0.96 cells of both earlier bending
# studies failed to solve (results_bending.csv and results_lscale.csv carry
# zeros there), which is why L/d=10 runs at reduced seed count and anything
# beyond it is not attempted.
set -e
WORKDIR=${WORKDIR:-/scratch/project_XXXXXX/test_rve}
cd "$WORKDIR"; mkdir -p logs

ls Job-ERG_*.inp Job-ERGH_*.inp 2>/dev/null | sed 's/\.inp$//' | sort > GlobalJobList_erg
N=$(wc -l < GlobalJobList_erg)
echo "eringen decks: $N   (expect 78 = 39 RVEs x {utx, ben})"
[ "$N" -ge 1 ] || { echo "ERROR: no Job-ERG*.inp in $WORKDIR -- generate first"; exit 1; }

# Quadratic tets on cells up to L/d=10 are the most expensive solves in the
# campaign; %10 keeps the memory footprint of concurrent C3D10H jobs sane.
SOLVE=$(sbatch --parsable \
  --partition=small --cpus-per-task=8 --mem=24G --time=04:00:00 \
  --array=1-${N}%10 \
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_erg \
  csc_solve_array.sh)
echo "solve array: $SOLVE"

POST=$(sbatch --parsable --dependency=afterany:${SOLVE} \
  --export=ALL,WORKDIR=$WORKDIR,CSV=rve_eringen.csv,RESULTS=results_eringen.csv \
  postprocess_firstorder.sh)
echo "postprocess: $POST  -> results_eringen.csv"

POSTH=$(sbatch --parsable --dependency=afterany:${SOLVE} \
  --export=ALL,WORKDIR=$WORKDIR,CSV=rve_eringen_homog.csv,RESULTS=results_eringen_homog.csv \
  postprocess_firstorder.sh)
echo "postprocess: $POSTH -> results_eringen_homog.csv"

cat <<'EOF'

Then, offline:

  cd results
  python3 ../analysis/fit_nonlocal.py results_eringen.csv results_eringen_homog.csv

which fits both nonclassical families to the same sweep:

  III.A  gradient / couple stress   E_app/E_inf = 1 + 12 l^2 / L^2     (stiffening)
  III.B  Eringen integral nonlocal  E_inf/E_app = 1 + (e0a)^2 / L^2    (softening)

On the published three-size sweep these give l^2 < 0 (family ruled out) and
e0a = 1.12 d with an intercept of 0.966 where the model requires 1.000. The
intercept is the large-cell asymptote, so it is the added sizes at L/d = 6, 8, 10
that decide whether the fitted nonlocal length is real or is absorbing the
ordinary first-order dilution.
EOF
