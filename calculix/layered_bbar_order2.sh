#!/bin/bash -l
# Can C3D10 + B-bar reach what Abaqus gets from C3D4H?
#
# layered_abaqus_ratio.sh established the problem: at L_mesh=0.0120 CalculiX's
# C3D4 undrained cell reads +12.28% above Abaqus's C3D4H, and refining makes it
# worse (-1.06% of movement against Abaqus's -10.55%). The displacement element
# cannot relieve the volumetric constraint by refining.
#
# B-bar is a no-op on C3D4 -- one integration point, so the element mean of the
# divergence IS the pointwise value. C3D10 has four, so it is not. This runs
# the same campaign cell three ways on ONE tetrahedralisation:
#
#   C3D10            does order 2 alone close the gap?
#   C3D10 + B-bar    does the mean-dilatation volumetric term close it?
#   C3D10 drained    the twin, for the ratio R = E_x(und)/E_x(drn)
#
# CCX_BBAR_NU is a Poisson-ratio threshold, so setting it to 0.45 turns B-bar
# on for the brine (nu = 0.49993) and leaves the ice (0.330) and the drained
# brine (0.406) on the stock assembly -- the same criterion SPAX_HYBRID_NU uses
# to ask Abaqus for a hybrid element. The drained twin is therefore unaffected
# by the switch and only needs solving once.
#
# Targets, from results/results_layermesh.csv:
#   Abaqus C3D4H at the same mesh (0.0120):  E_x 5.5656e9, R 2.3786
#   Abaqus C3D4H converged      (0.0060):    E_x 4.6102e9, R 1.9897
#
# WARNING: with B-bar active the equilibrium_gap check is invalid. Only the
# stiffness is patched; resultsmech.f still recovers internal forces with the
# standard B, so the volume-averaged stress stays consistent (the element mean
# of lambda*div(u) IS the B-bar pressure) but the reaction cross-check does
# not. Expect ~1.0 instead of ~1e-7 on the B-bar runs and do not read it as a
# failed solve. The unpatched runs beside it keep the check.
set -eu
cd "$(dirname "$0")/.."

PY=${PY:-/home/giacomo/venvs/sci/bin/python}
ROOT=${ROOT:-out_bbar2}
DECK=${DECK:-params/rve_layermesh.csv}
REF=${REF:-results/results_layermesh.csv}
STEM=${STEM:-LMESH_m0p0120}
GENSEED=${GENSEED:-s1}
BBNU=${BBNU:-0.45}
CPUS=${CPUS:-8}
JOBS=${JOBS:-1}                 # ~4.9M equations at order 2; keep peak memory down

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SPAX_MESH_TIMEOUT=21600 SPAX_MAX_RETRIES=12 SPAX_SEED=20260821
export CCX_ITER_TOL=${CCX_ITER_TOL:-1e-5}
export SPAX_MESH_ORDER=2

mkdir -p "$ROOT"

$PY - "$DECK" "$STEM" "$GENSEED" "$ROOT/$STEM.und.csv" <<'PYEOF'
import csv, sys
deck, stem, seed, out = sys.argv[1:5]
want = '%s_und_%s' % (stem, seed)
rows = list(csv.DictReader(open(deck)))
hit = [r for r in rows if r['run_id'] == want]
if not hit:
    raise SystemExit('no row %s in %s' % (want, deck))
r = dict(hit[0])
r['run_id'] = 'LAY'
w = csv.DictWriter(open(out, 'w', newline=''), fieldnames=list(rows[0].keys()))
w.writeheader()
w.writerow(r)
print('  %s -> L=%s L_mesh=%s n_slabs=%s slab_vof=%s'
      % (want, r['L'], r['L_mesh'], r['n_slabs'], r['slab_vof']))
PYEOF
$PY - "$ROOT/$STEM.und.csv" "$ROOT/$STEM.drn.csv" <<'PYEOF'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
for r in rows:
    r['K_inclusion'] = '2.2e+06'
w = csv.DictWriter(open(sys.argv[2], 'w', newline=''), fieldnames=list(rows[0].keys()))
w.writeheader()
w.writerows(rows)
PYEOF

drain_deck () {
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
        out.append(lines[i + 1])
        out.append('%r, %r\n' % (E, nu))
        i += 3
        done = True
        continue
    i += 1
if not done:
    raise SystemExit('no Mat_Inclusion card in ' + p)
open(p, 'w').writelines(out)
PYEOF
}

gen="$ROOT/${STEM}_gen"
echo "==== generating order 2 ===="
if [ ! -f "$gen/Job-LAY-utx.inp" ]; then
    "$PY" -u SpaX_Standalone.py "$ROOT/$STEM.und.csv" "$gen" \
        > "$ROOT/$STEM.gen.log" 2>&1 || { echo "  generation FAILED"; exit 1; }
fi
grep -m1 "Element types" "$ROOT/$STEM.gen.log" | sed 's/^ */  /'

# und_o2   : C3D10, stock assembly
# und_bbar : C3D10, B-bar on the brine only
# drn_o2   : C3D10 drained twin, identical mesh, one elastic card rewritten
for w in und_o2 und_bbar drn_o2; do
    d="$ROOT/${STEM}_$w"; rm -rf "$d"; mkdir -p "$d"
    cp "$gen"/Job-LAY-ut*.inp "$d/"
done
for f in "$ROOT/${STEM}_drn_o2"/Job-LAY-ut*.inp; do drain_deck "$f"; done

run_one () {   # $1 = subdir tag, $2 = binary, $3 = bbar threshold or empty
    w="$ROOT/${STEM}_$1"
    echo "  -- $1  ($2${3:+, CCX_BBAR_NU=$3})"
    SPAX_CCX="$2" python3 SpaX_CalculiX.py convert "$w" > /dev/null
    if [ -n "${3:-}" ]; then
        SPAX_CCX="$2" CCX_BBAR_NU="$3" python3 SpaX_CalculiX.py solve "$w" \
            --cpus "$CPUS" --jobs "$JOBS" | sed 's/^/     /'
    else
        SPAX_CCX="$2" python3 SpaX_CalculiX.py solve "$w" \
            --cpus "$CPUS" --jobs "$JOBS" | sed 's/^/     /'
    fi
    src=und; case "$1" in drn*) src=drn;; esac
    python3 SpaX_PostProcess.py "$ROOT/$STEM.$src.csv" "$w" "$ROOT/${STEM}_$1.out.csv" \
        > "$ROOT/${STEM}_$1.post.log" 2>&1
}

run_one und_o2   ccx_spax
run_one drn_o2   ccx_spax
run_one und_bbar ccx_bbar "$BBNU"

echo
$PY calculix/report_bbar_order2.py "$ROOT" "$REF" "$STEM"
