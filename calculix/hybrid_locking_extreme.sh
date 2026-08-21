#!/bin/bash -l
# The severe end of the near-incompressibility question.
#
# hybrid_locking_test.sh uses the deck brine, K=2.2 GPa / G=44.3 MPa, which is
# nu = 0.4900. analysis/benchmark_micromechanics.py describes the brine phase as
# K=2.2 GPa / G=0.44 MPa instead -- a hundred times softer in shear, and
# nu = 0.49993. That is where a displacement element without the hybrid
# formulation is actually expected to lock, so the reassuring result at 0.490
# does not transfer to it by assumption.
#
# Reuses the frozen packing from hybrid_locking_test.sh, so this is the same
# geometry as every case there and the numbers are directly comparable.
set -eu
cd "$(dirname "$0")/.."

PY=${PY:-/home/giacomo/venvs/sci/bin/python}
ROOT=${ROOT:-out_lock}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SPAX_MESH_TIMEOUT=1800 SPAX_MAX_RETRIES=12 SPAX_SEED=20260821

if [ ! -f "$ROOT/pack/LOCK.npy" ]; then
    echo "run hybrid_locking_test.sh first: this reuses its frozen packing" >&2
    exit 1
fi

sed 's/,4\.43e7$/,4.43e5/' "$ROOT/brine.csv" > "$ROOT/xbrine.csv"
sed 's/,0.30,0.050,/,0.30,0.035,/' "$ROOT/xbrine.csv" > "$ROOT/xbrine_fine.csv"
tail -1 "$ROOT/xbrine.csv" | awk -F, '{print "  G_inclusion = "$NF}'

run_case () {
    local name=$1 csv=$2 order=$3
    echo "==== $name (order $order) ===="
    SPAX_LOAD_PACKING="$ROOT/pack" SPAX_MESH_ORDER="$order" \
        "$PY" -u SpaX_Standalone.py "$csv" "$ROOT/$name" > "$ROOT/$name.gen.log" 2>&1
    grep -m1 'Element types' "$ROOT/$name.gen.log" || true
    python3 SpaX_CalculiX.py convert "$ROOT/$name" > /dev/null
    python3 SpaX_CalculiX.py solve   "$ROOT/$name" --cpus "${CPUS:-4}" --jobs 2
    python3 SpaX_PostProcess.py "$csv" "$ROOT/$name" "$ROOT/$name.csv" \
        > "$ROOT/$name.post.log" 2>&1
}

run_case xbrine_o1      "$ROOT/xbrine.csv"      1
run_case xbrine_o2      "$ROOT/xbrine.csv"      2
run_case xbrine_o2_fine "$ROOT/xbrine_fine.csv" 2

echo
echo "==== summary (nu = 0.49993) ===="
python3 - "$ROOT" <<'PYEOF'
import csv, os, sys
root = sys.argv[1]


def val(name, col):
    p = os.path.join(root, name + '.csv')
    if not os.path.isfile(p):
        return float('nan')
    with open(p) as f:
        r = list(csv.DictReader(f))
    try:
        return float(r[0][col])
    except (IndexError, KeyError, ValueError, TypeError):
        return float('nan')


print('%-16s %14s %14s %10s' % ('case', 'E_eff', 'G_eff', 'nu_eff'))
for n in ('xbrine_o1', 'xbrine_o2', 'xbrine_o2_fine'):
    print('%-16s %14.6e %14.6e %10.4f' % (n, val(n, 'E_eff'), val(n, 'G_eff'),
                                          val(n, 'nu_eff')))


def gap(a, b, col):
    x, y = val(a, col), val(b, col)
    return 100.0 * (x - y) / y if y else float('nan')


print()
print('order 1 vs order 2:      E_eff %+7.3f %%   G_eff %+7.3f %%'
      % (gap('xbrine_o1', 'xbrine_o2', 'E_eff'),
         gap('xbrine_o1', 'xbrine_o2', 'G_eff')))
print('order-2 refinement:      E_eff %+7.3f %%   G_eff %+7.3f %%'
      % (gap('xbrine_o2', 'xbrine_o2_fine', 'E_eff'),
         gap('xbrine_o2', 'xbrine_o2_fine', 'G_eff')))
PYEOF
