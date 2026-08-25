#!/bin/bash -l
# Is the CalculiX-vs-Abaqus gap on rve_gas.csv bigger than packing noise?
#
# validate_gas.sh puts CalculiX 0.5-3.6% below the stored Abaqus E_x, with the
# gap growing monotonically with void content. The achieved phase fractions
# agree to 0.1-2%, so the two cells hold the same amount of each phase -- but
# not in the same ARRANGEMENT, because the campaign's generation seed is not
# recorded and the packer is random.
#
# A difference is only evidence of a solver discrepancy if it exceeds the
# spread you get from re-packing the same specification. So: the worst-case row
# (GAS_v10, ~10% voids) at several seeds, everything else held fixed. If the
# spread over seeds covers the gap, this comparison cannot resolve a solver
# difference and should not be read as showing one.
set -eu
cd "$(dirname "$0")/.."

PY=${PY:-/home/giacomo/venvs/sci/bin/python}
CCX=${CCX:-ccx_spax}
ROOT=${ROOT:-out_gasseed}
ROW=${ROW:-GAS_v10}
SEEDS=${SEEDS:-"11 22 33 44"}
CPUS=${CPUS:-8}

export SPAX_CCX="$CCX"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SPAX_MESH_TIMEOUT=7200 SPAX_MAX_RETRIES=12 SPAX_MESH_ORDER=2
export CCX_ITER_TOL=${CCX_ITER_TOL:-1e-5}

mkdir -p "$ROOT"
head -1 params/rve_gas.csv > "$ROOT/row.csv"
grep "^$ROW," params/rve_gas.csv >> "$ROOT/row.csv"
# Only the uniaxial-X case is needed; drop the second mode to halve the solves.
$PY - "$ROOT/row.csv" <<'PYEOF'
import csv, sys
p = sys.argv[1]
rows = list(csv.DictReader(open(p)))
fn = list(rows[0].keys())
rows[0]['Mode2'] = ''
w = csv.DictWriter(open(p, 'w', newline=''), fieldnames=fn)
w.writeheader(); w.writerows(rows)
PYEOF

for s in $SEEDS; do
    d="$ROOT/s$s"
    echo "==== seed $s ===="
    if [ ! -f "$d/Job-$ROW-utx.inp" ]; then
        SPAX_SEED=$s "$PY" -u SpaX_Standalone.py "$ROOT/row.csv" "$d" \
            > "$ROOT/s$s.gen.log" 2>&1 || { echo "  generation FAILED"; continue; }
    fi
    grep -m1 "Done: " "$ROOT/s$s.gen.log" | sed 's/^ */  /'
    python3 SpaX_CalculiX.py convert "$d" > /dev/null
    python3 SpaX_CalculiX.py solve "$d" --cpus "$CPUS" --jobs 1 | sed 's/^/  /'
    python3 SpaX_PostProcess.py "$ROOT/row.csv" "$d" "$ROOT/s$s.csv" > "$ROOT/s$s.post.log" 2>&1
done

echo
echo "==== packing scatter, $ROW ===="
python3 - "$ROOT" "$ROW" <<'PYEOF'
import csv, glob, os, sys, statistics
root, row = sys.argv[1], sys.argv[2]
E, P = [], []
for f in sorted(glob.glob(os.path.join(root, 's*.csv'))):
    for r in csv.DictReader(open(f)):
        try:
            E.append(float(r['E_x'])); P.append(float(r['phi_soft_total']))
        except (KeyError, ValueError):
            pass
if not E:
    print('no results'); raise SystemExit(1)
m = statistics.mean(E)
# Population s.d. (ddof=0), the convention this repository quotes everywhere.
sd = statistics.pstdev(E)
print('  packings      : %d' % len(E))
print('  E_x           : ' + ', '.join('%.4g' % v for v in E))
print('  phi_soft      : ' + ', '.join('%.4f' % v for v in P))
print('  mean          : %.4g' % m)
print('  s.d. (ddof=0) : %.4g  (%.2f%% of mean)' % (sd, 100 * sd / m))
print('  full spread   : %.2f%% of mean' % (100 * (max(E) - min(E)) / m))
ABQ = 7638133908.0        # results/results_gas.csv, GAS_v10, Abaqus
print()
print('  Abaqus E_x    : %.4g' % ABQ)
print('  CalculiX mean is %+.2f%% against it' % (100 * (m - ABQ) / ABQ))
print('  each packing  : ' + ', '.join('%+.2f%%' % (100 * (v - ABQ) / ABQ) for v in E))
PYEOF
