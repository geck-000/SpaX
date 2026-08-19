#!/bin/bash -l
# Bring the 31 already-meshed LAYERB cells into the campaign instead of
# rebuilding them.
#
# An abandoned pass on 2026-08-13 left 31 of the 48 LAYERB cells in out_lb as
# complete utx+utz pairs, 19 GB of them, none truncated. Generation is the slow
# half of this campaign -- the pilot cell spent the better part of an hour in
# the mesher -- so rebuilding them would cost roughly 31 cells of it and change
# nothing about the answer.
#
# Two things are done here and nothing else.
#
# 1. The decks are moved to WORKDIR, where csc_solve_array.sh looks for them.
#
# 2. Their step card is collapsed from ten increments to one, by editing the
#    line rather than regenerating. These decks predate SPAX_LINEAR_ONE_STEP, so
#    they still ask for ten 0.1 increments and would write nine field frames the
#    first-order extractor never opens. Editing keeps the mesh, the packing and
#    the periodic equations exactly as generated -- which is the point of
#    reusing them at all -- and touches only the increment.
#
#    A copy of the original card is kept in reuse_layerb.manifest so the edit
#    can be checked, and the loop refuses to edit a deck that does not carry the
#    ten-increment card, rather than assuming it did.
#
# What this script does NOT do is decide whether the collapse is sound. That is
# settled by control_local.sh on a workstation, and the answer should be in hand
# before the solves are launched. Run with SPAX_KEEP_TEN=1 to reuse the decks
# as they are, at ten increments, if you would rather not depend on it.
#
#   ./reuse_layerb.sh [workdir]
set -e
W=${1:-${WORKDIR:-/scratch/project_2019020/test_rve}}
cd "$W"
# Idempotent by design: submit_ramp.sh calls this under `set -e`, and after the
# first run out_lb no longer exists. A second run must therefore be a no-op that
# still refreshes have_lb.txt and the missing deck, not a failure that aborts the
# submission.
if [ -d out_lb ]; then
  n=$(ls out_lb/Job-LB_*.inp 2>/dev/null | wc -l)
  echo "out_lb holds $n decks ($((n/2)) cells x 2 load cases)"
  if [ "$n" -ge 2 ]; then
    mv -f out_lb/Job-LB_*.inp . 2>/dev/null || true
    find out_lb -name '*_periodic_pairs.csv' -exec mv -t . {} + 2>/dev/null || true
  fi
  rmdir out_lb 2>/dev/null || true
else
  echo "no out_lb: already reused, or nothing was left there"
fi

: > reuse_layerb.manifest
edited=0; kept=0; odd=0
for f in Job-LB_*.inp; do
  [ -e "$f" ] || continue
  if [ -n "${SPAX_KEEP_TEN:-}" ]; then kept=$((kept+1)); continue; fi
  if grep -q '^0\.1, 1\., 1e-10, 0\.1' "$f"; then
    kept=$((kept+1))          # ten increments, which is what the extractor needs
  elif grep -q '^1\., 1\., 1e-10, 1\.' "$f"; then
    # A deck left at one increment by an earlier pass. Put it back: E_x comes
    # from a polyfit over the frame series, so a two-frame ODB extracts as zero.
    sed -i 's/^1\., 1\., 1e-10, 1\./0.1, 1., 1e-10, 0.1/' "$f"
    echo "$f  1.,1.,1e-10,1. -> 0.1,1.,1e-10,0.1 (restored)" >> reuse_layerb.manifest
    edited=$((edited+1))
  else
    echo "$f  UNRECOGNISED STEP CARD -- left alone" >> reuse_layerb.manifest
    odd=$((odd+1))
  fi
done

echo "restored to ten increments: $edited"
echo "already ten increments:      $kept"
[ "$odd" -gt 0 ] && echo "UNRECOGNISED, left as found:  $odd  (see reuse_layerb.manifest)"
echo "manifest: $W/reuse_layerb.manifest"

# The deck for what out_lb does not hold, derived from the filesystem so a
# partial rerun narrows it automatically rather than being hand-maintained.
ls Job-LB_*-utx.inp 2>/dev/null | sed 's|.*Job-||; s|-utx.inp||' | sort > have_lb.txt
echo "have: $(wc -l < have_lb.txt) cells"
python3 hpc/make_layerb_missing.py have_lb.txt . 2>/dev/null \
  || python3 make_layerb_missing.py have_lb.txt . \
  || echo "  (run make_layerb_missing.py by hand: it needs rve_layerb.csv beside it)"
