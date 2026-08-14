#!/bin/bash -l
# Submit a campaign end to end: generate, collect, solve, extract.
#
# This exists because of a failure mode that has now cost two campaigns. Every
# postprocessor takes the same environment contract, so any of them can be
# pointed at any campaign and will run happily to completion. The generic
# first-order route returns E and nu for anything; asked for a torsional
# rigidity or a stress concentration factor it returns nothing, writes its CSV,
# and exits zero. The reclaim stage then deletes the ODBs, and by the time the
# empty column is noticed there is nothing left to re-extract from -- the
# campaign has to be solved again from the mesh up.
#
# The fix is to stop choosing the postprocessor at submit time. A campaign's
# extraction is a property of the campaign, so it is recorded here, next to the
# deck name, and cannot be omitted by forgetting a flag.
#
#   ./spax_submit.sh <campaign> [extra sbatch args for the solve]
#
# Add a campaign by adding one line to the table below. A campaign whose name
# is not in the table is refused rather than run with a default, because the
# default is exactly what caused the problem.
set -e
WORKDIR=${WORKDIR:-/scratch/project_2019020/test_rve}
ACCT=${SPAX_ACCT:-project_2019020}
export PYTHONUSERBASE=${PYTHONUSERBASE:-/projappl/project_2019020/spax_py}
cd "$WORKDIR"; mkdir -p logs

# --- preflight: is there room to write the answer? ----------------------------
# A full /scratch does not stop a solve, it corrupts it. Abaqus is killed
# mid-write and leaves a truncated ODB that opens as "database file is corrupt",
# so the compute is spent and the result unrecoverable without a re-solve -- and
# it takes down unrelated campaigns writing at the same time. Cheaper to refuse
# to start than to discover it from a stack trace two campaigns later.
QUOTA_LINE=$(lfs quota -q -p 602019020 /scratch 2>/dev/null | head -1)
if [ -n "$QUOTA_LINE" ]; then
  USED_KB=$(echo "$QUOTA_LINE" | awk '{gsub(/\*/,"",$2); print $2}')
  LIM_KB=$(echo "$QUOTA_LINE" | awk '{print $4}')
  if [ -n "$USED_KB" ] && [ -n "$LIM_KB" ] && [ "$LIM_KB" -gt 0 ] 2>/dev/null; then
    PCT=$(( 100 * USED_KB / LIM_KB ))
    USED_G=$(( USED_KB / 1024 / 1024 )); LIM_G=$(( LIM_KB / 1024 / 1024 ))
    echo "scratch   : ${USED_G}G of ${LIM_G}G (${PCT}%)"
    if [ "$PCT" -ge "${SPAX_QUOTA_HARD:-90}" ]; then
      echo "ERROR: /scratch is ${PCT}% full. Solves would write truncated ODBs."
      echo "       Clear extracted ODBs first, or set SPAX_QUOTA_HARD to override."
      exit 3
    elif [ "$PCT" -ge "${SPAX_QUOTA_WARN:-75}" ]; then
      echo "WARNING: /scratch is ${PCT}% full -- clear extracted ODBs soon."
    fi
  fi
fi

CAMPAIGN=${1:?usage: spax_submit.sh <campaign> [solve sbatch args]}
shift || true
SOLVE_EXTRA="$@"

# campaign | deck | job prefix | postprocessor | results file | mesh order | solve mem
read -r -d '' TABLE <<'EOT' || true
torsion|rve_torsion.csv|TOR|postprocess_torsion.sh|results_torsion_K.csv|2|64G
torsion_big|rve_torsion_big.csv|TORB|postprocess_torsion.sh|results_torsion_big_K.csv|2|180G
torsion_homog|rve_torsion_homog.csv|TORH|postprocess_torsion.sh|results_torsion_homog_K.csv|2|64G
torsion_layer|rve_torsion_layer.csv|TORL|postprocess_torsion.sh|results_torsion_layer_K.csv|2|64G
torsion_layer_homog|rve_torsion_layer_homog.csv|TORLH|postprocess_torsion.sh|results_torsion_layer_homog_K.csv|2|64G
weibull|rve_weibull.csv|WBL|postprocess_weibull_scf.sh|results_weibull_scf.csv|1|32G
weibull_layer|rve_weibull_layer.csv|WBLL|postprocess_weibull_scf.sh|results_weibull_layer_scf.csv|1|180G
nlgeom_layer|rve_nlgeom_layer.csv|NLGL|postprocess_nlgeom.sh|results_nlgeom_layer.csv|1|32G
eringen_layer|rve_eringen_layer.csv|ERGL|postprocess_firstorder.sh|results_eringen_layer.csv|2|64G
eringen_layer_homog|rve_eringen_layer_homog.csv|ERGLH|postprocess_firstorder.sh|results_eringen_layer_homog.csv|2|64G
layercol_p060|rve_layercol_p060.csv|LCOL_p060|postprocess_firstorder.sh|results_layercol_p060.csv|1|32G
EOT

ROW=$(printf '%s\n' "$TABLE" | awk -F'|' -v c="$CAMPAIGN" '$1==c{print; exit}')
if [ -z "$ROW" ]; then
  echo "ERROR: '$CAMPAIGN' is not a known campaign. Known:"
  printf '%s\n' "$TABLE" | awk -F'|' '{printf "  %-22s -> %s\n", $1, $4}'
  exit 2
fi
DECK=$(echo "$ROW" | cut -d'|' -f2)
PREFIX=$(echo "$ROW" | cut -d'|' -f3)
POST=$(echo "$ROW" | cut -d'|' -f4)
RESULTS=$(echo "$ROW" | cut -d'|' -f5)
ORDER=$(echo "$ROW" | cut -d'|' -f6)
MEM=$(echo "$ROW" | cut -d'|' -f7)
case "$POST" in
  postprocess_torsion.sh)      SUFFIX='-tor' ;;
  postprocess_weibull_scf.sh)  SUFFIX='-utx' ;;
  postprocess_nlgeom.sh)       SUFFIX='-utx' ;;
  *)                           SUFFIX='' ;;
esac

[ -f "$DECK" ] || { echo "ERROR: deck $DECK not found in $WORKDIR"; exit 1; }
# postprocess_torsion.sh reads the deck from params/; keep the two in step.
mkdir -p params && cp -f "$DECK" params/

N=$(($(wc -l < "$DECK") - 1))
OUTDIR="out_${CAMPAIGN}"
echo "campaign  : $CAMPAIGN"
echo "deck      : $DECK ($N cells, mesh order $ORDER)"
echo "extraction: $POST -> $RESULTS"

GEN=$(sbatch --parsable --account=$ACCT --array=1-${N}%${SPAX_MAX_ARRAY:-25} --mem=$MEM \
  --export=ALL,WORKDIR=$WORKDIR,CSV=$DECK,OUTDIR=$OUTDIR,SPAX_MESH_ORDER=$ORDER,PYTHONUSERBASE=$PYTHONUSERBASE \
  generate_array.sh)
echo "generate  : $GEN"

CHAIN=$(mktemp "chain_${CAMPAIGN}_XXXX.sh")
cat > "$CHAIN" <<EOS
#!/bin/bash -l
#SBATCH --job-name=${CAMPAIGN}_chain
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=small --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:30:00
set -e
cd "$WORKDIR"
find $OUTDIR -name 'Job-*.inp' -exec mv -t "$WORKDIR" {} + 2>/dev/null || true
find $OUTDIR -name '*_periodic_pairs.csv' -exec mv -t "$WORKDIR" {} + 2>/dev/null || true
ls Job-${PREFIX}_*.inp 2>/dev/null | sed 's/\.inp\$//' | sort > GlobalJobList_${CAMPAIGN}
M=\$(wc -l < GlobalJobList_${CAMPAIGN})
echo "collected \$M of $N decks"
# A short collection is a generation failure, not something to solve around.
[ "\$M" -ge 1 ] || { echo "ERROR: nothing collected for $CAMPAIGN"; exit 1; }
[ "\$M" -eq "$N" ] || echo "WARNING: $((N)) expected, \$M collected -- check logs/spax_gen_${GEN}_*.err"
S=\$(sbatch --parsable --account=$ACCT --partition=small --cpus-per-task=4 \\
  --mem=$MEM --time=02:00:00 --array=1-\${M}%20 $SOLVE_EXTRA \\
  --export=ALL,WORKDIR=$WORKDIR,JOBLIST=GlobalJobList_${CAMPAIGN} csc_solve_array.sh)
echo "solve: \$S"
P=\$(sbatch --parsable --account=$ACCT --dependency=afterany:\${S} \\
  --export=ALL,WORKDIR=$WORKDIR,CSV=$DECK,RESULTS=$RESULTS,SUMM=$RESULTS,CURVES=curves_${CAMPAIGN}.csv,GLOB='Job-${PREFIX}_*${SUFFIX}.odb' \\
  $POST)
echo "extract: \$P -> $RESULTS   (via $POST)"
EOS
chmod +x "$CHAIN"
C=$(sbatch --parsable --account=$ACCT --dependency=afterany:${GEN} "$CHAIN")
echo "chain     : $C"
