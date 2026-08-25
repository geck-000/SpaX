#!/bin/bash -l
# The direct C3D4 vs C3D4H measurement, without an Abaqus licence.
#
# Everything before this inferred the hybrid element's worth from an
# order-1-versus-order-2 comparison inside CalculiX. That can only say the two
# element orders differ. The question is what ABAQUS gets from C3D4H that
# CalculiX cannot, and the tree already holds Abaqus's answer: the layered
# campaigns stored drained AND undrained E_x for hundreds of cells.
#
# The campaign's generation seed is not recorded, so a fresh run cannot
# reproduce its packing and E_x cannot be compared cell for cell
# (validate_gas.sh hit the same wall). But the decks come in DRAINED/UNDRAINED
# PAIRS, and
#
#     R = E_x(undrained) / E_x(drained)
#
# is the quantity the hybrid element acts on. Both codes mesh the drained cell
# (nu = 0.406) with the plain C3D4; only the undrained cell (nu = 0.49993)
# differs, C3D4H in Abaqus against C3D4 here. Geometry enters R through the
# ratio, where it largely cancels -- and the stored tables carry two or three
# SEEDS per condition, so ABAQUS'S OWN SEED SPREAD IN R sets the noise floor
# the CalculiX-vs-Abaqus difference has to beat.
#
# Locking makes a displacement element too STIFF, and only the undrained cell
# can lock, so it inflates R:
#
#     R_ccx >> R_abq   ->  C3D4 locks; C3D4H is buying Abaqus a compliance
#                          CalculiX cannot reproduce
#     R_ccx ~= R_abq   ->  nothing is lost by not having it
#
# One asymmetry to keep in view: here the drained twin is the SAME MESH with
# one elastic card rewritten, so R carries no packing noise at all, while the
# campaign re-packed for each drainage state (the stored porosity and
# phi_inclusion differ by 1-10% within a seed label). That makes the Abaqus R
# the noisier of the two, which is the safe direction -- it inflates the noise
# floor rather than the signal.
#
# Usage: one or more CASES, each the run_id stem the deck uses before
# _und_<seed> / _drn_<seed>.
#
#   CASES="LMESH_m0p0240 LMESH_m0p0120" \
#   DECK=params/rve_layermesh.csv REF=results/results_layermesh.csv \
#       calculix/layered_abaqus_ratio.sh
#
#   CASES="BRKB_b020 BRKB_b280" SEEDS="s1 s2 s3" \
#   DECK=params/rve_bracket_bridge.csv REF=results/results_bracket_bridge.csv \
#       calculix/layered_abaqus_ratio.sh
set -eu
cd "$(dirname "$0")/.."

PY=${PY:-/home/giacomo/venvs/sci/bin/python}
CCX=${CCX:-ccx_spax}
ROOT=${ROOT:-out_layerabq}
DECK=${DECK:-params/rve_layermesh.csv}
REF=${REF:-results/results_layermesh.csv}
CASES=${CASES:-"LMESH_m0p0240 LMESH_m0p0120"}
SEEDS=${SEEDS:-"s1 s2"}
GENSEED=${GENSEED:-s1}          # which stored seed's deck row to copy
CPUS=${CPUS:-8}

export SPAX_CCX="$CCX"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SPAX_MESH_TIMEOUT=14400 SPAX_MAX_RETRIES=12 SPAX_SEED=20260821
export CCX_ITER_TOL=${CCX_ITER_TOL:-1e-5}
# The layered campaigns ran on linear elements (5196ff1, the face-constraint
# limit), and C3D4 is the element under test.
export SPAX_MESH_ORDER=${SPAX_MESH_ORDER:-1}

mkdir -p "$ROOT"

# Pull the campaign's own undrained row verbatim -- every deck field is the
# campaign's, not a hand-written approximation of it.
extract_row () {   # $1 = case stem, $2 = seed, $3 = out csv
    $PY - "$DECK" "$1" "$2" "$3" <<'PYEOF'
import csv, os, sys
deck, stem, seed, out = sys.argv[1:5]
want = '%s_und_%s' % (stem, seed)
rows = list(csv.DictReader(open(deck)))
hit = [r for r in rows if r['run_id'] == want]
if not hit:
    raise SystemExit('no row %s in %s' % (want, deck))
r = dict(hit[0])
r['run_id'] = 'LAY'
# LMESH overrides the campaign's own element size.  Only for probing an
# element's behaviour at a resolution the machine can actually solve -- the
# Abaqus reference is at the campaign's L_mesh, so R is no longer comparable
# to it once this is set.
_lm = os.environ.get('LMESH', '').strip()
if _lm:
    r['L_mesh'] = _lm
w = csv.DictWriter(open(out, 'w', newline=''), fieldnames=list(rows[0].keys()))
w.writeheader()
w.writerow(r)
print('  deck row %s -> L=%s L_mesh=%s n_slabs=%s slab_vof=%s '
      'bridge_fraction=%s K=%s'
      % (want, r['L'], r['L_mesh'], r['n_slabs'], r['slab_vof'],
         r.get('bridge_fraction', '-'), r['K_inclusion']))
PYEOF
}

drain_deck () {          # rewrite the inclusion elastic card in place
    $PY - "$1" <<'PYEOF'
import sys
K, G = 2.2e6, 440029.33528897085
E = 9.0 * K * G / (3.0 * K + G)
nu = (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G))
p = sys.argv[1]
lines = open(p).readlines()
out, i, done = [], 0, False
while i < len(lines):
    out.append(lines[i])
    if lines[i].strip().lower().startswith('*material, name=mat_inclusion'):
        out.append(lines[i + 1])                 # *Elastic
        out.append('%r, %r\n' % (E, nu))         # replaces the E, nu line
        i += 3
        done = True
        continue
    i += 1
if not done:
    raise SystemExit('no Mat_Inclusion card in ' + p)
open(p, 'w').writelines(out)
PYEOF
}

for stem in $CASES; do
    echo "======== $stem ========"
    extract_row "$stem" "$GENSEED" "$ROOT/$stem.und.csv"
    # The drained twin, by column name -- E_sphere_inclusion is also 2.2e9,
    # so a positional sed rewrites the wrong field.
    $PY - "$ROOT/$stem.und.csv" "$ROOT/$stem.drn.csv" <<'PYEOF'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
for r in rows:
    r['K_inclusion'] = '2.2e+06'
w = csv.DictWriter(open(sys.argv[2], 'w', newline=''), fieldnames=list(rows[0].keys()))
w.writeheader()
w.writerows(rows)
PYEOF

    gen="$ROOT/${stem}_gen"
    if [ ! -f "$gen/Job-LAY-utx.inp" ]; then
        "$PY" -u SpaX_Standalone.py "$ROOT/$stem.und.csv" "$gen" \
            > "$ROOT/$stem.gen.log" 2>&1 || { echo "  generation FAILED"; continue; }
    fi
    grep -m1 "Element types" "$ROOT/$stem.gen.log" | sed 's/^ */  /'

    u="$ROOT/${stem}_und"; rm -rf "$u"; mkdir -p "$u"; cp "$gen"/Job-LAY-ut*.inp "$u/"
    d="$ROOT/${stem}_drn"; rm -rf "$d"; mkdir -p "$d"; cp "$gen"/Job-LAY-ut*.inp "$d/"
    for f in "$d"/Job-LAY-ut*.inp; do drain_deck "$f"; done

    for st in und drn; do
        w="$ROOT/${stem}_$st"
        echo "  -- $st"
        python3 SpaX_CalculiX.py convert "$w" > /dev/null
        python3 SpaX_CalculiX.py solve "$w" --cpus "$CPUS" --jobs "${JOBS:-2}" | sed 's/^/     /'
        python3 SpaX_PostProcess.py "$ROOT/$stem.$st.csv" "$w" "$ROOT/${stem}_$st.out.csv" \
            > "$ROOT/${stem}_$st.post.log" 2>&1

        # ---- U5+U6 arm: the same converted decks, with the inclusion
        # retyped to the deviatoric tet plus its nodal B-bar patches.  Same
        # mesh, same equations, same drainage state -- the ONLY difference
        # is the element, which is what R is being asked to isolate.
        # U6_ARMS: which patch coverages to run.
        #   u6     -- patches over the SOFT PHASE ONLY (how the campaign has
        #             been run so far)
        #   u6all  -- patches over BOTH phases.  nodalbbar's own comment warns
        #             that soft-phase-only leaves interface nodes with
        #             ONE-SIDED patches, and in a slab a few elements thick
        #             most soft-phase nodes ARE interface nodes.  On the
        #             frozen-geometry convergence cell u6all beat u6 at every
        #             one of five meshes, so it is carried here as its own arm
        #             rather than assumed.
        for arm in ${U6_ARMS:-u6 u6all}; do
            [ "${WITH_U6:-1}" = 1 ] || break
            case "$arm" in
                u6)    ES="--elset Sphere_Only" ;;
                u6all) ES="--elset Sphere_Only --elset Matrix_Only" ;;
                *) echo "  unknown U6 arm $arm"; continue ;;
            esac
            wu="${w}_${arm}"; rm -rf "$wu"; mkdir -p "$wu"
            ok=1
            for f in "$w"/*-ccx.inp; do
                SPAX_BBAR_SOLVER="${BBAR_SOLVER:-PARDISO}" "$PY" \
                    elements_ccx/nodalbbar.py "$f" "$wu/$(basename "$f")" $ES \
                    >> "$ROOT/${stem}_${st}.${arm}gen.log" 2>&1 || ok=0
            done
            if [ "$ok" = 1 ]; then
                echo "  -- $st ($arm)"
                python3 SpaX_CalculiX.py solve "$wu" --cpus "$CPUS" --jobs "${JOBS:-2}" \
                    | sed 's/^/     /'
                python3 SpaX_PostProcess.py "$ROOT/$stem.$st.csv" "$wu" \
                    "$ROOT/${stem}_${st}_${arm}.out.csv" \
                    > "$ROOT/${stem}_${st}_${arm}.post.log" 2>&1
            else
                echo "  -- $st ($arm) deck generation FAILED"
            fi
        done
    done
done

echo
$PY calculix/report_abaqus_ratio.py "$ROOT" "$REF" "$SEEDS" $CASES
