#!/bin/bash -l
# Mesh convergence for 4-bridge layered cells using F-barES-FEM-T4 (c=1) for undrained
#
# Target: BRKG_n04 cells (n_bridges=4) at L_mesh = 0.0240, 0.0180, 0.0120, 0.0090, 0.0060
# Drained (K=2.2e6, nu~0.406): plain C3D4
# Undrained (K=2.2e9, nu~0.4999): F-barES-FEM-T4 c=1
#
# Requires: ccx_spax (C3D4), ccx_fbar (F-barES-FEM-T4 patched)
set -eu

cd "${SPAX_ROOT:-/home/giacomo/SpaX}"

PY=${PY:-/home/giacomo/venvs/sci/bin/python}
CCX_C3D4=${CCX_C3D4:-ccx_spax}
CCX_FBAR=${CCX_FBAR:-ccx_fbar}
ROOT=${ROOT:-out_n4_meshconv}
PARAM_CSV=${PARAM_CSV:-params/rve_bracket_n4_meshconv.csv}
CPUS=${CPUS:-8}

export SPAX_CCX="$CCX_C3D4"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SPAX_MESH_TIMEOUT=7200 SPAX_MAX_RETRIES=12
export SPAX_SEED=${SPAX_SEED:-20260828}
export CCX_ITER_TOL=${CCX_ITER_TOL:-1e-5}

mkdir -p "$ROOT/pack"

# Filter to coarsest mesh for freezing geometry
COARSEST_LM=0.0240
$PY - "$PARAM_CSV" "$ROOT/seed.csv" <<PYEOF
import csv, sys
src, dst = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(src)))
fn = list(rows[0].keys())
# Keep only coarsest mesh, one seed, drained state
for r in rows:
    if r['L_mesh'] == '0.0240' and r['run_id'].endswith('_drn_s1'):
        w = csv.DictWriter(open(dst, 'w', newline=''), fieldnames=fn)
        w.writeheader()
        w.writerow(r)
        break
PYEOF

echo "==== Freeze geometry (seed s1, L_mesh=0.0240) ===="
if [ ! -f "$ROOT/pack/BRKG_n04_drn_s1_m0p0240.npy" ]; then
    SPAX_SAVE_PACKING="$ROOT/pack" "$PY" -u SpaX_Standalone.py \
        "$ROOT/seed.csv" "$ROOT/seedgen" > "$ROOT/seed.log" 2>&1
fi
ls -la "$ROOT/pack/"

# Process each L_mesh
for lm in 0.0240 0.0180 0.0120 0.0090 0.0060; do
    tag="m${lm/./p}"
    echo "==== L_mesh = $lm ===="
    
    # Create per-mesh CSV for all states/seeds
    $PY - "$PARAM_CSV" "$ROOT/$tag.csv" "$lm" <<PYEOF
import csv, sys
src, dst, lm = sys.argv[1], sys.argv[2], sys.argv[3]
rows = list(csv.DictReader(open(src)))
fn = list(rows[0].keys())
keep = [r for r in rows if r['L_mesh'] == lm]
w = csv.DictWriter(open(dst, 'w', newline=''), fieldnames=fn)
w.writeheader()
w.writerows(keep)
PYEOF

    # Generate all decks (load frozen packing)
    echo "  Generating decks..."
    SPAX_LOAD_PACKING="$ROOT/pack" "$PY" -u SpaX_Standalone.py \
        "$ROOT/$tag.csv" "$ROOT/$tag" > "$ROOT/$tag.gen.log" 2>&1 || {
        echo "  Generation FAILED"
        tail -20 "$ROOT/$tag.gen.log"
        continue
    }
    grep -m1 "Done: BRKG" "$ROOT/$tag.gen.log" | sed 's/^ */    /'

    # Solve each run
    for run_dir in "$ROOT/$tag"/*/; do
        [ -d "$run_dir" ] || continue
        run_name=$(basename "$run_dir")
        inp="$run_dir/Job-${run_name}-utx.inp"
        [ -f "$inp" ] || continue
        
        # Determine state from run_name
        if [[ "$run_name" == *"_und_"* ]]; then
            # Undrained: use F-barES-FEM-T4 c=1
            echo "  Solving $run_name (F-bar c=1)..."
            SPAX_CCX="$CCX_FBAR" CCX_FBAR_C=1 SPAX_CCX_SIGMA_FROM_RF=1 \
                "$PY" -u SpaX_CalculiX.py solve "$run_dir" --cpus "$CPUS" --jobs 1 2>&1 | sed 's/^/    /'
        else
            # Drained: use plain C3D4
            echo "  Solving $run_name (C3D4)..."
            SPAX_CCX="$CCX_C3D4" \
                "$PY" -u SpaX_CalculiX.py solve "$run_dir" --cpus "$CPUS" --jobs 1 2>&1 | sed 's/^/    /'
        fi
    done

    # Post-process
    echo "  Post-processing..."
    "$PY" -u SpaX_PostProcess.py "$ROOT/$tag.csv" "$ROOT/$tag" "$ROOT/$tag.csv.out" \
        > "$ROOT/$tag.post.log" 2>&1
done

echo
echo "==== Mesh convergence summary ===="
"$PY" - "$ROOT" <<'PYEOF'
import csv, glob, os, re, sys
root = sys.argv[1]

def extract(lm, state):
    tag = 'm' + lm.replace('.', 'p')
    out_csv = os.path.join(root, f'{tag}.csv.out')
    if not os.path.isfile(out_csv):
        return None, None
    vals = []
    for r in csv.DictReader(open(out_csv)):
        try:
            E = float(r['E_x'])
            vals.append(E)
        except (KeyError, ValueError):
            pass
    if vals:
        return sum(vals)/len(vals), (max(vals)-min(vals))/sum(vals)*100 if sum(vals) else 0
    return None, None

print(f"{'L_mesh':>8}  {'E_drn (GPa)':>12}  {'E_und (GPa)':>12}  {'R=E_und/E_drn':>12}  {'drn spread%':>10}  {'und spread%':>10}")
for lm in ['0.0240', '0.0180', '0.0120', '0.0090', '0.0060']:
    Ed, sd = extract(lm, 'drn')
    Eu, su = extract(lm, 'und')
    if Ed and Eu:
        print(f"{float(lm):8.4f}  {Ed/1e9:12.4f}  {Eu/1e9:12.4f}  {Eu/Ed:12.4f}  {sd:10.2f}  {su:10.2f}")
    else:
        print(f"{float(lm):8.4f}  {'--':>12}  {'--':>12}  {'--':>12}  {'--':>10}  {'--':>10}")

PYEOF