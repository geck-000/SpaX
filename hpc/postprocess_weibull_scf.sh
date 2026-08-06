#!/bin/bash -l
#SBATCH --job-name=wbl_scf
#SBATCH --partition=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/wbl_scf_%j.out
# Extract the matrix stress-concentration field from every solved WBL utx ODB,
# writing BOTH the summary row (results_weibull_scf.csv) and the per-element
# (scf, volume) dump that the offline Weibull m-sweep needs.
#
# The split matters: opening an ODB needs an Abaqus licence and is the whole
# cost here, while the m-sweep itself is arithmetic. Dumping the arrays once
# means analysis/weibull_sensitivity.py can be re-run for any Weibull modulus,
# on any machine, without touching the cluster again.
set -e
# Common post-processing contract (see hpc/README.md): WORKDIR, CSV, RESULTS,
# OUTDIR, each used where relevant, so one caller can drive every
# postprocess_*.sh the same way.
WORKDIR=${WORKDIR:-/scratch/project_XXXXXX/test_rve}
CSV=${CSV:-rve_weibull.csv}
DUMPDIR=${DUMPDIR:-$WORKDIR/weibull_dumps}
L=${L:-0.50}
cd "$WORKDIR"
mkdir -p "$DUMPDIR" logs

# Lmod does not work in batch on Roihu (CSC's zz-csc-env.sh bails out for
# non-interactive shells, leaving Tcl env-modules that cannot parse Abaqus's
# .lua modulefile). Source the snapshot of the working interactive environment
# instead, exactly as csc_solve_array.sh does.
source "$HOME/abaqus_env.sh"

OUT=${RESULTS:-results_weibull_scf.csv}
n=0
for odb in Job-WBL_*-utx.odb; do
  [ -e "$odb" ] || { echo "no Job-WBL_*-utx.odb in $WORKDIR"; exit 1; }
  rid=$(basename "$odb" -utx.odb); rid=${rid#Job-}
  # cell edge is a per-deck property; read it back from the deck rather than
  # assuming, since the six WBL cases inherit different L from their sources
  Lr=$(awk -F, -v id="$rid" 'NR==1{for(i=1;i<=NF;i++){if($i=="run_id")c=i; if($i=="L")l=i}; next} $c==id{print $l; exit}' "$CSV" 2>/dev/null || true)
  Lr=${Lr:-$L}
  echo "== $rid  (L=$Lr)"
  abaqus python scf_extract.py "$odb" "$Lr" "$rid" "$OUT" "$DUMPDIR/$rid.npz"
  n=$((n+1))
done
echo "extracted $n ODBs -> $OUT  and  $DUMPDIR/*.npz"

cat <<EOF

Then, offline (plain python3, no Abaqus):

  python3 weibull_sensitivity.py $DUMPDIR results_weibull

which reports the effective concentration factor as a volume-weighted m-norm,

  SCF_eff(m) = [ (1/V) INT SCF^m dV ]^(1/m)

over m = 1..50, together with the ranking of the six cases at each m. The P99
quoted in the manuscript equals SCF_eff at one particular m; this makes that
choice explicit and shows which conclusions survive it. Note that the ranking
is already known to move with the statistic: gas voids exceed brine channels at
P50 and P90 but fall below them at P99 and in the maximum, so connectivity
populates the extreme tail without shifting the bulk.
EOF
