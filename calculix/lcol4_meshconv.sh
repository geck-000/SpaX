#!/bin/bash -l
# Mesh convergence for the 4-bridge closure cell (LCOL4: b=0.2929, slab=0.10).
#
# WHY THIS STUDY EXISTS. The paper's Appendix A.6 caveat -- drained modulus
# 8.7% high, undrained 35% high at 0.7 elements across the brine layer -- was
# measured on the TWO-bridge rve_layermesh gate. But the closure's n(b) is
# calibrated on FOUR-bridge cells (N4_LCOL_p100 et al.) solved at mesh_for
# resolution (2.5 elements/layer), which the gate already put in the converged
# regime. This study re-runs the gate geometry at four bridges so the caveat
# that belongs in the paper is about the cells the closure actually uses.
#
#   geometry  : L=0.50, n_slabs=4, slab_vof=0.10, bridge_fraction=0.2929,
#               n_bridges=4   (identical layer thickness to the 2-bridge gate)
#   sweep     : L_mesh = 0.0240 0.0180 0.0120 0.0090 0.0060
#   drained   : C3D4  (nu ~ 0.406, no locking)          -> ccx_spax
#   undrained : F-barES-FEM-T4 c=1 (nu ~ 0.4999, locks) -> ccx_fbar
#
# Geometry is frozen with SPAX_SAVE_PACKING/SPAX_LOAD_PACKING so only L_mesh
# changes between points; slabs and bridges are deterministic from the deck.
# run_id is held at LCOL4 across the whole sweep so the frozen packing loads.
set -eu
cd "${SPAX_ROOT:-/home/giacomo/SpaX}"

PY=${PY:-/home/giacomo/venvs/sci/bin/python}
CCX=${CCX:-ccx_spax}
CCX_FBAR=${CCX_FBAR:-ccx_fbar}
ROOT=${ROOT:-out_lcol4_meshconv}
CSV=params/rve_lcol4_meshconv.csv
# Undrained (F-bar c=1) is limited by the mastruct insertion wall (2^31):
# 0.0240 -> 3.6e8, 0.0120 -> 1.6e9, and 0.0080 would be ~1e10 (refused).
# Drained (plain C3D4, nu ~ 0.406, no locking) has no such wall but 0.0060 is
# ~22M elements and out of reach for this workstation.
MESHES=${MESHES:-"0.0240 0.0120 0.0080"}
FBAR_MESHES=${FBAR_MESHES:-"0.0240 0.0120"}
CPUS=${CPUS:-8}

export SPAX_CCX="$CCX"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SPAX_MESH_TIMEOUT=7200 SPAX_MAX_RETRIES=12 SPAX_SEED=20260828
export SPAX_MESH_ORDER=1
export CCX_ITER_TOL=${CCX_ITER_TOL:-1e-5}

mkdir -p "$ROOT/pack"

# --- freeze the geometry at the coarsest mesh --------------------------------
FIRST=$(echo $MESHES | awk '{print $1}')
$PY - "$CSV" "$ROOT/seed.csv" "$FIRST" <<'PYEOF'
import csv, sys
src, dst, lm = sys.argv[1], sys.argv[2], sys.argv[3]
rows = list(csv.DictReader(open(src)))
fn = list(rows[0].keys())
for r in rows:
    if r['L_mesh'] == lm and r['K_inclusion'] == '2.2e+09':
        w = csv.DictWriter(open(dst, 'w', newline=''), fieldnames=fn)
        w.writeheader()
        w.writerow(r)
        break
else:
    raise SystemExit('no undrained row at L_mesh=%s' % lm)
PYEOF
echo "==== freeze geometry (L_mesh=$FIRST, undrained s1) ===="
if [ ! -f "$ROOT/pack/LCOL4.npy" ]; then
    SPAX_SAVE_PACKING="$ROOT/pack" "$PY" -u SpaX_Standalone.py \
        "$ROOT/seed.csv" "$ROOT/seedgen" > "$ROOT/seed.log" 2>&1
fi
ls "$ROOT/pack"

# --- sweep L_mesh -----------------------------------------------------------
for lm in $MESHES; do
    tag="m${lm/./p}"
    echo "==== L_mesh = $lm ===="
    # One row per mesh (run_id is constant, so generating more would collide).
    $PY - "$CSV" "$ROOT/$tag.csv" "$lm" <<'PYEOF'
import csv, sys
src, dst, lm = sys.argv[1], sys.argv[2], sys.argv[3]
rows = list(csv.DictReader(open(src)))
fn = list(rows[0].keys())
for r in rows:
    if r['L_mesh'] == lm and r['K_inclusion'] == '2.2e+09':
        w = csv.DictWriter(open(dst, 'w', newline=''), fieldnames=fn)
        w.writeheader()
        w.writerow(r)
        break
else:
    raise SystemExit('no undrained row at L_mesh=%s' % lm)
PYEOF

    if [ ! -f "$ROOT/$tag/Job-LCOL4-utx.inp" ]; then
        SPAX_LOAD_PACKING="$ROOT/pack" "$PY" -u SpaX_Standalone.py \
            "$ROOT/$tag.csv" "$ROOT/$tag" > "$ROOT/$tag.gen.log" 2>&1 \
            || { echo "  generation FAILED"; continue; }
    fi
    grep -m1 "Done: LCOL4" "$ROOT/$tag.gen.log" | sed 's/^ */  /'

    # undrained: convert + F-bar, solve with ccx_fbar (only where c=1 fits)
    if echo "$FBAR_MESHES" | grep -qw "$lm"; then
      w="$ROOT/${tag}_und"; rm -rf "$w"; mkdir -p "$w"
      cp "$ROOT/$tag"/Job-LCOL4-utx.inp "$w/"
      cp "$ROOT/$tag"/Job-LCOL4-utz.inp "$w/" 2>/dev/null || true
      "$PY" SpaX_CalculiX.py convert "$w" > /dev/null
      for inp in "$w"/Job-LCOL4-*-ccx.inp; do
          [ -f "$inp" ] || continue
          base="${inp%-ccx.inp}"
          SPAX_CCX="$CCX_FBAR" "$PY" elements_ccx/fbares.py "$inp" \
              "${base}-fbar.inp" --elset Sphere_Only --cycles 1 > /dev/null
          mv "${base}-fbar.inp" "$inp"
      done
      # F-barES-FEM-T4 c=1 on the unsymmetric path stores both triangles of
      # the factor, which OOMs this 30 GB machine past ~0.0120 (in-core factor
      # ~44 GB).  Spill the factor to disk: CCX_PARDISO_OOC is the in-core
      # budget in MB (MKL_PARDISO_OOC_MAX_CORE_SIZE); 4096 and 8192 return
      # PARDISO error -9 ("not enough memory" for the core), 16384 works with
      # a 22.6 GB peak RSS.  MKL_PARDISO_OOC_PATH must be real disk (not /tmp
      # tmpfs), so point it at the workdir itself, which lives under /home.
      mkdir -p "$w/ooc_temp"
      SPAX_CCX="$CCX_FBAR" CCX_FBAR_C=1 \
          CCX_PARDISO_OOC=${CCX_PARDISO_OOC:-16384} \
          MKL_PARDISO_OOC_PATH="ooc_temp" \
          "$PY" -u SpaX_CalculiX.py solve "$w" --cpus "$CPUS" --jobs 1 \
          | sed 's/^/    /'
      SPAX_SOLVER=calculix SPAX_CCX_SIGMA_FROM_RF=1 \
          "$PY" SpaX_PostProcess.py "$ROOT/$tag.csv" "$w" "$ROOT/${tag}_und.results.csv" \
          > "$ROOT/${tag}_und.post.log" 2>&1 || { echo "  und post FAILED"; }
    else
      echo "  undrained skipped: L_mesh=$lm not in FBAR_MESHES"
    fi

    # drained: same deck with the brine bulk modulus released
    w="$ROOT/${tag}_drn"; rm -rf "$w"; mkdir -p "$w"
    cp "$ROOT/$tag"/Job-LCOL4-utx.inp "$w/"
    cp "$ROOT/$tag"/Job-LCOL4-utz.inp "$w/" 2>/dev/null || true
    for inp in "$w"/Job-LCOL4-*.inp; do
        [ -f "$inp" ] || continue
        "$PY" - "$inp" <<'PYEOF'
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
    done
    "$PY" SpaX_CalculiX.py convert "$w" > /dev/null
    "$PY" -u SpaX_CalculiX.py solve "$w" --cpus "$CPUS" --jobs 1 \
        | sed 's/^/    /'
    SPAX_SOLVER=calculix "$PY" SpaX_PostProcess.py "$ROOT/$tag.csv" "$w" \
        "$ROOT/${tag}_drn.results.csv" \
        > "$ROOT/${tag}_drn.post.log" 2>&1 || { echo "  drn post FAILED"; }
done

echo
echo "==== 4-bridge (LCOL4) mesh convergence ===="
"$PY" - "$ROOT" <<'PYEOF'
import csv, glob, os, sys
root = sys.argv[1]

def grab(tag, state):
    p = os.path.join(root, '%s_%s.results.csv' % (tag, state))
    if not os.path.isfile(p):
        return []
    out = []
    for r in csv.DictReader(open(p)):
        try:
            out.append(float(r['E_x']) / 1e9)
        except (KeyError, ValueError):
            pass
    return out

rows = []
for d in sorted(glob.glob(os.path.join(root, 'm*'))):
    tag = os.path.basename(d)
    lm = float(list(csv.DictReader(open(os.path.join(root, tag + '.csv'))))[0]['L_mesh'])
    drn = grab(tag, 'drn')
    und = grab(tag, 'und')
    rows.append((lm, tag, drn, und))
rows.sort(reverse=True)

ref = rows[-1]
refd = sum(ref[2]) / len(ref[2]) if ref[2] else float('nan')
# undrained reference: finest mesh that actually has undrained data
und_rows = [r for r in rows if r[3]]
refu = sum(und_rows[-1][3]) / len(und_rows[-1][3]) if und_rows else float('nan')
print('%-8s %6s %12s %12s %10s %10s %10s'
      % ('L_mesh', 'seeds', 'E_drn', 'E_und', 'R', 'drn drift', 'und drift'))
for lm, tag, drn, und in rows:
    md = sum(drn) / len(drn) if drn else float('nan')
    mu = sum(und) / len(und) if und else float('nan')
    dd = (100 * (md / refd - 1)) if (drn and refd == refd) else float('nan')
    du = (100 * (mu / refu - 1)) if (und and refu == refu) else float('nan')
    print('%-8.4f %6d %12.4f %12.4f %10.4f %+9.1f%% %+9.1f%%'
          % (lm, max(len(drn), len(und)), md, mu, mu / md, dd, du))
print('\nfinest drained mesh (L_mesh=%.4f) is the drained drift reference' % ref[0])
PYEOF