#!/bin/bash -l
#SBATCH --job-name=rerun2
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#
# Staged re-run controller, v2.
#
# Three things v1 got wrong, each of which cost a chain:
#
#   1. SUBMIT CAP. v1 chunked the solve array and dependency-chained the chunks,
#      but submitted every chunk at once. The cap counts SUBMITTED tasks
#      (pending+running), not running ones, so a campaign larger than the cap
#      was rejected outright however it was throttled. v2 submits ONE chunk and
#      re-enters itself for the next, so at most CHUNK tasks are ever submitted.
#
#   2. SCRATCH. ODBs are ~200 MB each and the project quota is 250 GB; the last
#      full re-run peaked at 230 GB and had to be cleaned by hand. v2 deletes a
#      campaign's ODBs and decks as soon as its results file is written and
#      verified, and refuses to start a solve stage without headroom.
#
#   3. TORSION. The stock post-processor looks for a step named 'Bending' and
#      silently mis-extracts a torsion ODB. v2 lets a campaign name several
#      post-processors, '+'-separated, so torsion runs the first-order pass for
#      the moduli and torsion_extract.py for the rigidity.
#
# USAGE, from $WORKDIR:
#   sbatch --account=<acct> --export=ALL,WORKDIR=$PWD,STAGE=gen,IDX=1,\
#          MANIFEST=<file>,SPAX_ACCOUNT=<acct>,PYTHONUSERBASE=<pub> rerun_paper_v2.sh
#
set -e
unset SLURM_GTIDS

W=${WORKDIR:?set WORKDIR}
cd "$W"
A="--account=${SPAX_ACCOUNT:?set SPAX_ACCOUNT}"
PUB=${PYTHONUSERBASE:?set PYTHONUSERBASE}
SELF=rerun_paper_v2.sh
MAN=${MANIFEST:-rerun_paper_manifest.tsv}
SEED=${SPAX_SEED:-20260810}
TOL=${VOLUME_TOL:-1.15}
IDX=${IDX:?set IDX}
STAGE=${STAGE:?set STAGE}
CHUNK=${SPAX_MAX_ARRAY:-150}        # tasks submitted at once; cap is 200
LO=${LO:-1}                          # first task of the chunk being submitted
KEEP_ODB=${SPAX_KEEP_ODB:-0}         # 1 disables the cleanup, for debugging
MIN_FREE_GB=${SPAX_MIN_FREE_GB:-60}  # refuse to solve with less headroom

export PYTHONUSERBASE=$PUB
export PATH="$PUB/bin:$PATH"

read_row () { grep -v '^#' "$MAN" | grep -v '^[[:space:]]*$' | sed -n "${1}p"; }
ROW=$(read_row "$IDX")
if [ -z "$ROW" ]; then
  echo "=== manifest exhausted at index $IDX -- re-run complete ==="
  exit 0
fi
NAME=$(echo   "$ROW" | cut -f1)
DECK=$(echo   "$ROW" | cut -f2)
OUTDIR=$(echo "$ROW" | cut -f3)
GLOB=$(echo   "$ROW" | cut -f4)
RESULTS=$(echo "$ROW" | cut -f5)
POST=$(echo   "$ROW" | cut -f6)
EXTRA=$(echo  "$ROW" | cut -f7)
[ "$EXTRA" = "-" ] && EXTRA=""

echo "=== [$IDX] $NAME : $DECK, stage $STAGE${LO:+ (from task $LO)} ==="

# --- scratch accounting -------------------------------------------------------
free_gb () {
  local used
  used=$(du -s --block-size=1G "$W" 2>/dev/null | cut -f1)
  echo $(( ${SPAX_QUOTA_GB:-250} - ${used:-0} ))
}

# Remove a finished campaign's heavy artefacts. Only ever called once its
# results file exists and is non-trivial, so extraction has already happened.
reclaim () {
  local deck="$1" res="$2"
  [ "$KEEP_ODB" = "1" ] && { echo "reclaim: disabled"; return 0; }
  if [ ! -s "$res" ] || [ "$(wc -l < "$res")" -lt 2 ]; then
    echo "reclaim: $res missing or empty -- keeping ODBs for diagnosis"
    return 0
  fi
  local n=0
  python3 - "params/$deck" > ".reclaim_$NAME" <<'PYEOF'
import csv, sys
seen = set()
for row in csv.DictReader(open(sys.argv[1], encoding='utf8', errors='replace')):
    rid = (row.get('run_id') or '').strip()
    if rid and rid not in seen:
        seen.add(rid); print(rid)
PYEOF
  while read -r rid; do
    for f in Job-${rid}-*.odb Job-${rid}-*.inp Job-${rid}-*.sta Job-${rid}-*.msg \
             Job-${rid}-*.dat Job-${rid}-*.prt Job-${rid}-*.com Job-${rid}-*.sim \
             Job-${rid}-*.stt Job-${rid}-*.mdl Job-${rid}-*.lck; do
      [ -e "$f" ] && { rm -f "$f"; n=$((n+1)); }
    done
    rm -rf Job-${rid}-*.simdir 2>/dev/null || true
  done < ".reclaim_$NAME"
  rm -f ".reclaim_$NAME"
  echo "reclaim: removed $n files for $NAME; free now ~$(free_gb) GB"
}

case "$STAGE" in

  gen)
    mkdir -p "$OUTDIR" logs
    NR=$(python3 -c "
import csv,sys
print(sum(1 for _ in csv.DictReader(open('params/$DECK',encoding='utf8',errors='replace'))))")
    echo "$NR RVEs"
    G=$(sbatch --parsable $A --array=1-${NR}%40 \
        --export=ALL,WORKDIR=$W,CSV=params/$DECK,OUTDIR=$OUTDIR,SPAX_SEED=$SEED,PYTHONUSERBASE=$PUB${EXTRA:+,$EXTRA} \
        generate_array.sh)
    echo "gen: $G"
    sbatch $A --dependency=afterany:$G \
      --export=ALL,WORKDIR=$W,STAGE=solve,IDX=$IDX,LO=1,MANIFEST=$MAN,SPAX_ACCOUNT=${SPAX_ACCOUNT},PYTHONUSERBASE=$PUB "$SELF"
    ;;

  solve)
    if [ "$LO" -eq 1 ]; then
      # audit and stage the decks once, on the first chunk only
      python3 audit_volume.py "params/$DECK" "$OUTDIR" "$TOL"; rc=$?
      if [ "$rc" -eq 2 ]; then
        echo "!!! VOLUME AUDIT FAILED for $NAME -- chain stopped."; exit 1
      elif [ "$rc" -ne 0 ]; then
        echo "!!! audit_volume.py could not run (exit $rc) -- environment fault."; exit 1
      fi
      find "$OUTDIR" -name 'Job-*.inp' -exec mv -t "$W" {} + 2>/dev/null || true
      find "$OUTDIR" -name '*_periodic_pairs.csv' -exec mv -t "$W" {} + 2>/dev/null || true
      python3 - "params/$DECK" > "GJ_$NAME" <<'PYEOF'
import csv, glob, os, sys
seen = set()
for row in csv.DictReader(open(sys.argv[1], encoding='utf8', errors='replace')):
    rid = (row.get('run_id') or '').strip()
    if not rid or rid in seen:
        continue
    seen.add(rid)
    for f in sorted(glob.glob('Job-%s-*.inp' % rid)):
        print(os.path.basename(f)[:-4])
PYEOF
    fi

    N=$(wc -l < "GJ_$NAME")
    [ "$N" -ge 1 ] || { echo "!!! no decks for $NAME"; exit 1; }
    [ "$LO" -eq 1 ] && echo "decks collected: $N"

    FREE=$(free_gb)
    if [ "$FREE" -lt "$MIN_FREE_GB" ]; then
      echo "!!! only ${FREE} GB free against a ${SPAX_QUOTA_GB:-250} GB quota."
      echo "!!! Refusing to solve; earlier campaigns should have reclaimed."
      exit 1
    fi

    LMAX=$(python3 - "params/$DECK" <<'PYEOF'
import csv, sys
v = []
for r in csv.DictReader(open(sys.argv[1], encoding='utf8', errors='replace')):
    try: v.append(float(r.get('L') or 0))
    except ValueError: pass
print('%.3f' % (max(v) if v else 0.5))
PYEOF
)
    BEND=$(grep -cE -- '-(ben|tor)$' "GJ_$NAME" 2>/dev/null || echo 0)
    read SOLVE_MEM SOLVE_TIME <<EOF
$(python3 -c "
L=float('$LMAX') or 0.5
s=(L/0.5)**3
heavy = $BEND > 0
mem=max(16,min(360, 16*s*(3 if heavy else 1)))
hrs=max(0.67,min(8.0,(40/60.0)*s*(2 if heavy else 1)))
h=int(hrs); m=int(round((hrs-h)*60))
print('%dG %02d:%02d:00' % (int(mem), h, m))
")
EOF
    HI=$(( LO + CHUNK - 1 )); [ "$HI" -gt "$N" ] && HI=$N
    echo "cell edge $LMAX -> mem $SOLVE_MEM, walltime $SOLVE_TIME; free ${FREE} GB"
    S=$(sbatch --parsable $A --partition=small --cpus-per-task=4 \
        --mem=$SOLVE_MEM --time=$SOLVE_TIME --array=${LO}-${HI}%30 \
        --export=ALL,WORKDIR=$W,JOBLIST=GJ_$NAME csc_solve_array.sh)
    echo "solve ${LO}-${HI} of $N: $S"

    if [ "$HI" -lt "$N" ]; then
      # Re-enter for the next chunk rather than submitting it now: the cap
      # counts submitted tasks, so only one chunk may be in flight.
      sbatch $A --dependency=afterany:$S \
        --export=ALL,WORKDIR=$W,STAGE=solve,IDX=$IDX,LO=$(( HI + 1 )),MANIFEST=$MAN,SPAX_ACCOUNT=${SPAX_ACCOUNT},PYTHONUSERBASE=$PUB "$SELF"
    else
      sbatch $A --dependency=afterany:$S \
        --export=ALL,WORKDIR=$W,STAGE=post,IDX=$IDX,MANIFEST=$MAN,SPAX_ACCOUNT=${SPAX_ACCOUNT},PYTHONUSERBASE=$PUB "$SELF"
    fi
    ;;

  post)
    # A campaign may name several post-processors, '+'-separated. Torsion needs
    # two: the first-order pass for E_eff/nu_eff, and torsion_extract.py for the
    # rigidity, because the first-order pass cannot read a torsion step.
    DEP=""; LASTP=""
    echo "$POST" | tr '+' '\n' | while read -r p; do
      [ -z "$p" ] && continue
      echo "  post script: $p"
    done
    for p in $(echo "$POST" | tr '+' ' '); do
      RES_P="$RESULTS"
      case "$p" in
        postprocess_torsion.sh) RES_P="${RESULTS%.csv}_K.csv" ;;
      esac
      P=$(sbatch --parsable $A $DEP \
          --export=ALL,WORKDIR=$W,CSV=params/$DECK,RESULTS=$RES_P,OUTDIR=$OUTDIR,GLOB=$GLOB \
          "$p")
      echo "post: $P -> $RES_P  ($p)"
      DEP="--dependency=afterany:$P"
      LASTP=$P
    done
    # reclaim once every extractor for this campaign has run
    sbatch $A --dependency=afterany:$LASTP \
      --export=ALL,WORKDIR=$W,STAGE=reclaim,IDX=$IDX,MANIFEST=$MAN,SPAX_ACCOUNT=${SPAX_ACCOUNT},PYTHONUSERBASE=$PUB "$SELF"
    ;;

  reclaim)
    reclaim "$DECK" "$RESULTS"
    NEXT=$(( IDX + 1 ))
    sbatch $A --export=ALL,WORKDIR=$W,STAGE=gen,IDX=$NEXT,MANIFEST=$MAN,SPAX_ACCOUNT=${SPAX_ACCOUNT},PYTHONUSERBASE=$PUB "$SELF"
    ;;

  *)
    echo "unknown STAGE=$STAGE"; exit 1;;
esac
