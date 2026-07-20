#!/bin/bash
# Study #8 deck builder. Mesh generation is non-deterministic (inclusions are
# randomly placed), so we CANNOT get identical geometry from separate generation
# runs. Instead we generate ONE mesh per slice (the nlgeom tension deck, disp
# +0.02) and stamp all three load cases onto that SAME deck, editing only:
#   * the *Step nlgeom flag         (OFF for the linear reference)
#   * the driving RP-1 displacement (magnitude EPS_DISP, sign per case)
#   * the *Static increment line    (finer stepping for nlgeom convergence + curve)
# All three load cases of a slice therefore share byte-identical geometry, mesh,
# material and PBC.
#
# Target |strain| = EPS_DISP / L = 0.010 / 0.50 = 2%  (4% tension de-stabilises the
# porous channelled base; 2% converges for all slices and still exposes the
# geometric nonlinearity, which scales ~linearly with strain).
#
# In:  out_nlgeom_ten/Job-NLGTEN_z{25,65,95}-utx.inp   (nlgeom=YES, disp=+0.02)
# Out: out_nlgeom/  per slice: LIN (OFF,+2%), TEN (ON,+2%), CMP (ON,-2%)
set -e
SRC=out_nlgeom_ten
OUT=out_nlgeom
EPS_DISP=0.010
STATIC_OLD='0.1, 1., 1e-10, 0.1'
STATIC_NEW='0.05, 1., 1e-08, 0.1'      # initial 0.1% strain steps, max 0.2%
rm -rf "$OUT"; mkdir -p "$OUT"
for base in "$SRC"/Job-NLGTEN_z*-utx.inp; do
  slice=$(basename "$base" | sed -E 's/Job-NLGTEN_(z[0-9]+)-utx\.inp/\1/')
  [ "$(grep -c 'RP-1, 1, 1, 0.02$' "$base")" = "1" ] || { echo "ERR: BC line not unique in $base"; exit 1; }
  [ "$(grep -c 'nlgeom=YES' "$base")" = "1" ]        || { echo "ERR: nlgeom flag not unique in $base"; exit 1; }
  # common: retarget displacement magnitude and refine the increment line
  common_sed="s/RP-1, 1, 1, 0.02\$/RP-1, 1, 1, ${EPS_DISP}/; s/${STATIC_OLD}/${STATIC_NEW}/"
  sed "$common_sed; s/nlgeom=YES/nlgeom=NO/"                         "$base" > "$OUT/Job-NLGLIN_${slice}-utx.inp"
  sed "$common_sed"                                                  "$base" > "$OUT/Job-NLGTEN_${slice}-utx.inp"
  sed "$common_sed; s/RP-1, 1, 1, ${EPS_DISP}\$/RP-1, 1, 1, -${EPS_DISP}/" "$base" > "$OUT/Job-NLGCMP_${slice}-utx.inp"
  echo "$slice: TEN/LIN/CMP stamped (|disp|=${EPS_DISP}) from the same mesh"
done
echo "built $(ls "$OUT"/*.inp | wc -l) decks in $OUT/ (expect 9)"
# verify LIN differs from TEN by exactly the nlgeom flag line, CMP by the disp sign
for slice in z25 z65 z95; do
  nLIN=$(diff "$OUT/Job-NLGTEN_${slice}-utx.inp" "$OUT/Job-NLGLIN_${slice}-utx.inp" | grep -c '^<')
  nCMP=$(diff "$OUT/Job-NLGTEN_${slice}-utx.inp" "$OUT/Job-NLGCMP_${slice}-utx.inp" | grep -c '^<')
  { [ "$nLIN" = "1" ] && [ "$nCMP" = "1" ]; } \
    && echo "$slice: geometry identical across all 3 cases OK (1-line delta each)" \
    || echo "$slice: UNEXPECTED delta (LIN=$nLIN CMP=$nCMP lines)"
done