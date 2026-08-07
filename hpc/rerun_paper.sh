#!/bin/bash -l
#SBATCH --job-name=rerun_paper
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=small
#SBATCH --ntasks=1 --cpus-per-task=1 --mem=2G --time=00:20:00
#
# Re-run every campaign the paper rests on, against the fixed mesher.
#
# WHY: SpaX_GmshPeriodic.py built each inclusion as its bounding sphere rather
# than the ellipsoid the packer placed. Every solved RVE therefore carried
# 1/sphericity^2 too much inclusion volume (1.6x cold, 2.2x warm) and had
# spherical pockets instead of the depth-varying elongated ones the setup
# specifies. The moduli were right for the cells built; the cells were wrong.
#
# STRUCTURE: one campaign at a time, three chained jobs each --
#   gen   : array over the deck rows, one task per RVE
#   solve : collect decks, build a job list, array over it
#   post  : extraction to results_<campaign>.csv, then chain to the next campaign
# Chained rather than submitted at once because the small partition caps
# submitted jobs at 200 and the whole re-run is ~1100 solves.
#
# GATING: the first campaign is followed by a volume audit (deckvol.py against
# the deck's own target). If the meshed inclusion fraction still exceeds the
# analytic one by more than TOL, the chain stops rather than spending the rest.
#
# USAGE, from $WORKDIR:
#   sbatch --account=<acct> --export=ALL,WORKDIR=$PWD,STAGE=gen,IDX=1 rerun_paper.sh
#
set -e
unset SLURM_GTIDS

W=${WORKDIR:?set WORKDIR}
cd "$W"
A="--account=${SPAX_ACCOUNT:-project_XXXXXX}"
PUB=${PYTHONUSERBASE:-/projappl/project_XXXXXX/spax_py}
SELF=rerun_paper.sh
MAN=${MANIFEST:-rerun_paper_manifest.tsv}
SEED=${SPAX_SEED:-20260806}
TOL=${VOLUME_TOL:-1.15}       # meshed/analytic must fall below this
IDX=${IDX:?set IDX}
STAGE=${STAGE:?set STAGE}

# ---- read row IDX of the manifest (comments and blanks skipped) -------------
read_row () {
  grep -v '^#' "$MAN" | grep -v '^[[:space:]]*$' | sed -n "${1}p"
}
ROW=$(read_row "$IDX")
if [ -z "$ROW" ]; then
  echo "=== manifest exhausted at index $IDX -- re-run complete ==="
  exit 0
fi
NAME=$(echo "$ROW"  | cut -f1)
DECK=$(echo "$ROW"  | cut -f2)
OUTDIR=$(echo "$ROW"| cut -f3)
GLOB=$(echo "$ROW"  | cut -f4)
RESULTS=$(echo "$ROW"| cut -f5)
POST=$(echo "$ROW"  | cut -f6)
EXTRA=$(echo "$ROW" | cut -f7)
[ "$EXTRA" = "-" ] && EXTRA=""

NROW=$(( $(wc -l < "params/$DECK") - 1 ))
echo "=== [$IDX] $NAME : $DECK ($NROW RVEs), stage $STAGE ==="

case "$STAGE" in

  gen)
    mkdir -p logs "$OUTDIR"
    G=$(sbatch --parsable $A --array=1-${NROW}%40 \
        --export=ALL,WORKDIR=$W,CSV=params/$DECK,OUTDIR=$OUTDIR,SPAX_SEED=$SEED,PYTHONUSERBASE=$PUB${EXTRA:+,$EXTRA} \
        generate_array.sh)
    echo "gen: $G"
    sbatch $A --dependency=afterany:$G \
      --export=ALL,WORKDIR=$W,STAGE=solve,IDX=$IDX,MANIFEST=$MAN "$SELF"
    ;;

  solve)
    # Audit the geometry before spending solve time on it.
    export PYTHONUSERBASE=$PUB; export PATH="$PUB/bin:$PATH"
    if [ ! -f audit_volume.py ]; then
      echo "!!! audit_volume.py is missing -- refusing to run ungated."
      exit 1
    fi
    echo "--- volume audit ---"
    python3 audit_volume.py "params/$DECK" "$OUTDIR" "$TOL" || {
      echo "!!! VOLUME AUDIT FAILED for $NAME -- chain stopped."
      echo "!!! The meshed inclusion fraction still exceeds what the deck asks for."
      exit 1
    }
    find "$OUTDIR" -name 'Job-*.inp' -exec mv -t "$W" {} + 2>/dev/null || true
    find "$OUTDIR" -name '*_periodic_pairs.csv' -exec mv -t "$W" {} + 2>/dev/null || true
    # Build the job list from the deck's own run_ids, not from a filename glob.
    # A glob cannot always separate two campaigns: colseeds and colseeds_extra
    # are both Job-CSEED_*, so the second would collect the first's decks as
    # well -- harmless while they run in sequence, because the solve array skips
    # decks whose ODB exists and the post-processor filters by the deck's
    # run_ids, but wasteful, and wrong outright if the two ever run at the same
    # time. The deck is the authority on which RVEs belong to the campaign.
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
    N=$(wc -l < "GJ_$NAME")
    echo "decks collected: $N"
    if [ "$N" -lt 1 ]; then
      echo "!!! no decks for $NAME -- chain stopped."; exit 1
    fi
    # Bending is memory-bound in a way the uniaxial cases are not: the mesh
    # grows as (L/d)^3 and the direct solver dominates. Four eringen bending
    # solves were killed at 16.3-16.7 GB against a 16 GB allocation, and an
    # Abaqus killed for memory still leaves an ODB behind -- which the solve
    # array reads as success, so the loss is silent until the results table
    # comes up short. Give any campaign carrying bending decks real headroom
    # rather than discovering the ceiling one cell at a time.
    # Size the solve from the cell edge, not from the load case. Element count
    # goes as L^3, and so does the direct solver's appetite, so a campaign at
    # L=1.00 needs roughly eight times what one at L=0.50 does. Keying this on
    # bending alone was wrong: the uniaxial L=1.00 campaigns ran out of memory
    # at 16 GB (MaxRSS 16.7 GB) and out of walltime at forty minutes, while the
    # bending flag never applied to them. Bending on top costs about another
    # factor of three. Nodes carry 762 GB, so being generous costs nothing but
    # queue position.
    LMAX=$(python3 - "params/$DECK" <<'PYEOF'
import csv, sys
vals = []
for r in csv.DictReader(open(sys.argv[1], encoding='utf8', errors='replace')):
    try:
        vals.append(float(r.get('L') or 0))
    except ValueError:
        pass
print('%.3f' % (max(vals) if vals else 0.5))
PYEOF
)
    read SOLVE_MEM SOLVE_TIME <<EOF
$(python3 -c "
import sys
L=float('$LMAX') or 0.5
bend=$(grep -qc -- '-ben\$' "GJ_$NAME" 2>/dev/null && echo 1 || echo 0)
s=(L/0.5)**3
mem=16*s*(3 if bend else 1)
mem=max(16,min(360,mem))
hrs=(40/60.0)*s*(2 if bend else 1)
hrs=max(0.67,min(8.0,hrs))
h=int(hrs); m=int(round((hrs-h)*60))
print('%dG %02d:%02d:00' % (int(mem), h, m))
")
EOF
    echo "cell edge $LMAX -> mem $SOLVE_MEM, walltime $SOLVE_TIME"
    S=$(sbatch --parsable $A --partition=small --cpus-per-task=4 --mem=$SOLVE_MEM \
        --time=$SOLVE_TIME --array=1-${N}%30 \
        --export=ALL,WORKDIR=$W,JOBLIST=GJ_$NAME csc_solve_array.sh)
    echo "solve: $S"
    sbatch $A --dependency=afterany:$S \
      --export=ALL,WORKDIR=$W,STAGE=post,IDX=$IDX,MANIFEST=$MAN "$SELF"
    ;;

  post)
    P=$(sbatch --parsable $A \
        --export=ALL,WORKDIR=$W,CSV=params/$DECK,RESULTS=$RESULTS,OUTDIR=$OUTDIR \
        "$POST")
    echo "post: $P -> $RESULTS"
    NEXT=$(( IDX + 1 ))
    sbatch $A --dependency=afterany:$P \
      --export=ALL,WORKDIR=$W,STAGE=gen,IDX=$NEXT,MANIFEST=$MAN "$SELF"
    ;;

  *)
    echo "unknown STAGE=$STAGE"; exit 1;;
esac
