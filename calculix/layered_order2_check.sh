#!/bin/bash -l
# Is order 2 actually converged on the confined undrained layered cell?
#
# layered_incompressible.sh measured order 1 losing +2.1 points more than its
# drained twin, and I concluded order 2 "removes it by construction". That is
# circular: it was an order-1-versus-order-2 comparison, so it can only say the
# two differ, not which one is right. If C3D10 locks too, both are wrong and
# the gap between them understates the damage.
#
# The non-circular test is mesh refinement AT ORDER 2 on a fixed geometry.
# Volumetric locking is a discretisation error, so it must shrink as the mesh
# refines. And it needs the same paired control as before: the drained twin on
# the identical mesh sets how much drift is ordinary discretisation, so only
# the EXCESS drift of the undrained cell is locking.
#
#   drift_und >> drift_drn   ->  C3D10 is still locking; B-bar is warranted
#   drift_und ~= drift_drn   ->  order 2 is converged and is the fix
#
# Geometry is frozen with SPAX_SAVE_PACKING/SPAX_LOAD_PACKING so only L_mesh
# changes; the slabs and bridges are deterministic from the deck.
set -eu
cd "$(dirname "$0")/.."

PY=${PY:-/home/giacomo/venvs/sci/bin/python}
CCX=${CCX:-ccx_spax}
ROOT=${ROOT:-out_layero2}
L=${L:-0.30}
MESHES=${MESHES:-"0.020 0.015 0.011"}
CPUS=${CPUS:-8}

export SPAX_CCX="$CCX"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SPAX_MESH_TIMEOUT=7200 SPAX_MAX_RETRIES=12 SPAX_SEED=20260821
export SPAX_MESH_ORDER=2
export CCX_ITER_TOL=${CCX_ITER_TOL:-1e-5}

mkdir -p "$ROOT/pack"
DISP=$($PY -c "print(0.01 * $L)")
HDR='run_id,L,L_mesh,Is_Porous,E_matrix,nu_matrix,VoF_sphere,r_avg,r_std,Mode,Disp,Mode2,Disp2,VoF_void_sphere,VoF_incl_sphere,E_sphere_inclusion,nu_sphere_inclusion,sphericity_avg,sphericity_std,min_distance,max_iterations,nlgeom_flag,PBC_Method,Kappa,Bending_Plane,Bending_PBC_Type,generate_channels,channel_vof_target,r_channel_avg,r_channel_std,Growth_Direction,Growth_Concentration,Inclusion_Type,K_inclusion,G_inclusion,n_slabs,slab_vof,bridge_fraction,n_bridges,slab_axis,bridge_correlation'
mkrow () {   # L_mesh -> csv on stdout path $2
    printf '%s\nLAY,%s,%s,Composite,9.43e9,0.33,0.0325,0.035,0.008,Uniaxial Tension X,%s,Uniaxial Tension Z,%s,0.0100,0.0225,2.2e9,0.48,0.80,0.1,0.002,200000,OFF,Gmsh,0,xz,Lesicar,No,0,0.02,0.005,Z,0.5,Liquid,2.2e9,440029.33528897085,4,0.1000,0.2929,2,x,0.0\n' \
        "$HDR" "$L" "$1" "$DISP" "$DISP" > "$2"
}

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

FIRST=$(echo $MESHES | awk '{print $1}')
mkrow "$FIRST" "$ROOT/seed.csv"
echo "==== freeze the geometry ===="
if [ ! -f "$ROOT/pack/LAY.npy" ]; then
    SPAX_SAVE_PACKING="$ROOT/pack" "$PY" -u SpaX_Standalone.py \
        "$ROOT/seed.csv" "$ROOT/seedgen" > "$ROOT/seed.log" 2>&1
fi
ls "$ROOT/pack"

for lm in $MESHES; do
    tag="m${lm/./}"
    echo "==== L_mesh = $lm (order 2) ===="
    mkrow "$lm" "$ROOT/$tag.csv"
    sed 's/,Liquid,2.2e9,440029/,Liquid,2.2e6,440029/' "$ROOT/$tag.csv" > "$ROOT/${tag}_drn.csv"
    if [ ! -f "$ROOT/$tag/Job-LAY-utx.inp" ]; then
        SPAX_LOAD_PACKING="$ROOT/pack" "$PY" -u SpaX_Standalone.py \
            "$ROOT/$tag.csv" "$ROOT/$tag" > "$ROOT/$tag.gen.log" 2>&1 \
            || { echo "  generation FAILED"; continue; }
    fi
    grep -m1 "Done: LAY" "$ROOT/$tag.gen.log" | sed 's/^ */  /'

    for var in und drn; do
        w="$ROOT/${tag}_$var"; rm -rf "$w"; mkdir -p "$w"
        cp "$ROOT/$tag"/Job-LAY-utx.inp "$w/"
        if [ "$var" = drn ]; then drain_deck "$w/Job-LAY-utx.inp"; fi
        python3 SpaX_CalculiX.py convert "$w" > /dev/null
        python3 SpaX_CalculiX.py solve "$w" --cpus "$CPUS" --jobs 1 | sed 's/^/    /'
        csv="$ROOT/$tag.csv"; [ "$var" = drn ] && csv="$ROOT/${tag}_drn.csv"
        python3 SpaX_PostProcess.py "$csv" "$w" "$ROOT/${tag}_$var.out.csv" \
            > "$ROOT/${tag}_$var.post.log" 2>&1
    done
done

echo
echo "==== order-2 mesh convergence, confined layered cell ===="
python3 - "$ROOT" <<'PYEOF'
import csv, glob, os, re, sys
root = sys.argv[1]


def val(tag, var, col='E_x'):
    p = os.path.join(root, '%s_%s.out.csv' % (tag, var))
    if not os.path.isfile(p):
        return float('nan')
    for r in csv.DictReader(open(p)):
        try:
            return float(r[col])
        except (KeyError, ValueError, TypeError):
            return float('nan')
    return float('nan')


tags = []
for f in sorted(glob.glob(os.path.join(root, 'm*_und.out.csv'))):
    tag = os.path.basename(f).replace('_und.out.csv', '')
    p = os.path.join(root, tag + '.csv')
    lm = float(list(csv.DictReader(open(p)))[0]['L_mesh'])
    gen = os.path.join(root, tag + '.gen.log')
    nel = '?'
    if os.path.isfile(gen):
        m = re.search(r'Done: LAY \(\d+ nodes, (\d+) elements\)', open(gen).read())
        if m:
            nel = m.group(1)
    tags.append((lm, tag, nel))
tags.sort(reverse=True)

print('%-8s %10s %16s %16s' % ('L_mesh', 'elements', 'E_x undrained', 'E_x drained'))
for lm, tag, nel in tags:
    print('%-8.3f %10s %16.6e %16.6e' % (lm, nel, val(tag, 'und'), val(tag, 'drn')))

if len(tags) >= 2:
    (lc, tc, _), (lf, tf, _) = tags[0], tags[-1]
    def drift(var):
        a, b = val(tc, var), val(tf, var)
        return 100.0 * (a - b) / b if b else float('nan')
    print()
    print('coarsest -> finest (%.3f -> %.3f):' % (lc, lf))
    print('  undrained nu=0.49993   E_x drift %+7.3f %%' % drift('und'))
    print('  drained   nu=0.406     E_x drift %+7.3f %%' % drift('drn'))
    print('  excess attributable to incompressibility: %+7.3f points'
          % (drift('und') - drift('drn')))
    print()
    print('order-1 reference (layered_incompressible.sh): +2.11 points')
    print('spherical order-2 refinement was +0.22 %% (converged)')
PYEOF
