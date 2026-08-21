#!/bin/bash -l
# Does the missing hybrid element bite in the LAYERED cells?
#
# hybrid_locking_test.sh answered "no" for spherical brine inclusions, and gave
# the reason: the brine is an isolated soft pocket that neither carries load nor
# constrains the matrix, so the element's isochoric behaviour inside it barely
# reaches the effective stiffness. It also named the conditions under which that
# stops holding -- a soft phase that PERCOLATES and carries load, or geometric
# CONFINEMENT.
#
# The undrained layered cell has both, and is the sharpest case in the tree:
#
#   drained    K=2.2 MPa  G=0.44 MPa  ->  nu = 0.406   (Abaqus: C3D4, no hybrid)
#   undrained  K=2.2 GPa  G=0.44 MPa  ->  nu = 0.49993 (Abaqus: C3D4H, HYBRID)
#
# The brine is a cell-spanning slab confined between ice plates, with the ice
# bridges as constrictions -- and pressure transmission through that confined
# brine is the mechanism the layered closure rests on. Undrained, K=2.2 GPa is
# only 4x below the ice's 9.25 GPa, so unlike the spherical case the brine IS a
# load path. This is exactly where a displacement element without the mixed
# formulation is expected to lock.
#
# Design: ONE geometry per element order. The drained twin is made by rewriting
# the single inclusion elastic card on the IDENTICAL mesh (the paired design of
# analysis/make_drained_twin.py), so drainage is the only difference and no
# re-packing scatter enters.
#
#   gap_und = E_eff(order 1) vs E_eff(order 2), undrained
#   gap_drn = the same, drained
#
# Locking shows as gap_und >> gap_drn. If the two gaps match, the order effect
# is ordinary linear-tet stiffness, as it was for the spheres.
set -eu
cd "$(dirname "$0")/.."

PY=${PY:-/home/giacomo/venvs/sci/bin/python}
CCX=${CCX:-ccx_spax}
ROOT=${ROOT:-out_layerinc}
L=${L:-0.30}
LMESH=${LMESH:-0.015}
CPUS=${CPUS:-8}

export SPAX_CCX="$CCX"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SPAX_MESH_TIMEOUT=7200 SPAX_MAX_RETRIES=12 SPAX_SEED=20260821
export CCX_ITER_TOL=${CCX_ITER_TOL:-1e-5}

mkdir -p "$ROOT"

# Morphology copied from params/rve_layermesh.csv (4 slabs, 2 bridges, x-normal),
# at a cell size that fits this machine. Mode X loads ACROSS the layers, which is
# the constriction-dominated direction; Mode2 Z is in-plane.
DISP=$($PY -c "print(0.01 * $L)")
HDR='run_id,L,L_mesh,Is_Porous,E_matrix,nu_matrix,VoF_sphere,r_avg,r_std,Mode,Disp,Mode2,Disp2,VoF_void_sphere,VoF_incl_sphere,E_sphere_inclusion,nu_sphere_inclusion,sphericity_avg,sphericity_std,min_distance,max_iterations,nlgeom_flag,PBC_Method,Kappa,Bending_Plane,Bending_PBC_Type,generate_channels,channel_vof_target,r_channel_avg,r_channel_std,Growth_Direction,Growth_Concentration,Inclusion_Type,K_inclusion,G_inclusion,n_slabs,slab_vof,bridge_fraction,n_bridges,slab_axis,bridge_correlation'
printf '%s\nLAY,%s,%s,Composite,9.43e9,0.33,0.0325,0.035,0.008,Uniaxial Tension X,%s,Uniaxial Tension Z,%s,0.0100,0.0225,2.2e9,0.48,0.80,0.1,0.002,200000,OFF,Gmsh,0,xz,Lesicar,No,0,0.02,0.005,Z,0.5,Liquid,2.2e9,440029.33528897085,4,0.1000,0.2929,2,x,0.0\n' \
    "$HDR" "$L" "$LMESH" "$DISP" "$DISP" > "$ROOT/und.csv"

# The drained twin: same everything, K lowered. Used only for the post-processor;
# the DECK is produced by rewriting the elastic card, not by regenerating.
sed 's/,Liquid,2.2e9,440029/,Liquid,2.2e6,440029/' "$ROOT/und.csv" > "$ROOT/drn.csv"

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

for order in 1 2; do
    gen="$ROOT/o$order"
    echo "======== order $order ========"
    if [ ! -f "$gen/Job-LAY-utx.inp" ]; then
        SPAX_MESH_ORDER=$order "$PY" -u SpaX_Standalone.py "$ROOT/und.csv" "$gen" \
            > "$ROOT/o$order.gen.log" 2>&1 || { echo "  generation FAILED"; continue; }
    fi
    grep -m1 "Done: LAY" "$ROOT/o$order.gen.log" | sed 's/^ */  /'
    grep -m1 "Element types" "$ROOT/o$order.gen.log" | sed 's/^ */  /'

    # undrained: solve the generated deck as-is
    u="$ROOT/o${order}_und"; rm -rf "$u"; mkdir -p "$u"
    cp "$gen"/Job-LAY-ut*.inp "$u/"
    # drained: identical mesh, one elastic card rewritten
    d="$ROOT/o${order}_drn"; rm -rf "$d"; mkdir -p "$d"
    cp "$gen"/Job-LAY-ut*.inp "$d/"
    for f in "$d"/Job-LAY-ut*.inp; do drain_deck "$f"; done

    for tag in und drn; do
        w="$ROOT/o${order}_$tag"
        echo "  -- $tag"
        python3 SpaX_CalculiX.py convert "$w" > /dev/null
        python3 SpaX_CalculiX.py solve "$w" --cpus "$CPUS" --jobs 2 | sed 's/^/     /'
        python3 SpaX_PostProcess.py "$ROOT/$tag.csv" "$w" "$ROOT/o${order}_$tag.csv" \
            > "$ROOT/o${order}_$tag.post.log" 2>&1
    done
done

echo
echo "==== confined layered brine: does order 1 lose more when undrained? ===="
python3 - "$ROOT" <<'PYEOF'
import csv, os, sys
root = sys.argv[1]


def get(name, col):
    p = os.path.join(root, name + '.csv')
    if not os.path.isfile(p):
        return float('nan')
    for r in csv.DictReader(open(p)):
        try:
            return float(r[col])
        except (KeyError, ValueError, TypeError):
            return float('nan')
    return float('nan')


print('%-22s %14s %14s %12s' % ('case', 'E_x (across)', 'E_z (in-plane)', 'gap'))
for n in ('o1_und', 'o2_und', 'o1_drn', 'o2_drn'):
    print('%-22s %14.6e %14.6e %12.2e' % (
        n, get(n, 'E_x'), get(n, 'E_z'), get(n, 'equilibrium_gap')))


def d(a, b, col):
    x, y = get(a, col), get(b, col)
    return 100.0 * (x - y) / y if y else float('nan')


print()
print('order 1 vs order 2, same mesh:')
print('  undrained  nu=0.49993   E_x %+7.3f %%   E_z %+7.3f %%'
      % (d('o1_und', 'o2_und', 'E_x'), d('o1_und', 'o2_und', 'E_z')))
print('  drained    nu=0.406     E_x %+7.3f %%   E_z %+7.3f %%'
      % (d('o1_drn', 'o2_drn', 'E_x'), d('o1_drn', 'o2_drn', 'E_z')))
print()
print('(spherical-inclusion reference, hybrid_locking_test.sh:')
print('  brine nu=0.490 +3.62 %, twin nu=0.300 +4.03 % -- indistinguishable)')
PYEOF
