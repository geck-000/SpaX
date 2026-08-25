#!/bin/bash -l
# Mesh convergence on a cell whose brine slab is RESOLVED.  See make_slabconv.py
# for why the campaign cells cannot answer this (0.35-1.04 elements through the
# slab; no asymptotic regime exists there).
#
#   NS="10 15 20 30 40" KG=500 elements_ccx/tests/slabconv.sh
#
# Per mesh: plain C3D4 drained and undrained, and F-barES-FEM-T4 c=1 undrained.
# R = C1111(und)/C1111(drn), denominator always plain C3D4 -- the same
# convention report_abaqus_ratio.py uses, so the arms differ only in the
# element under test.  The _abq.inp decks are for Abaqus C3D4H on Roihu and are
# written but not solved here.
set -eu
cd "$(dirname "$0")/../.."
PY=/home/giacomo/venvs/sci/bin/python3
NS=${NS:-"10 20 30 40"}
KG=${KG:-500}
JIT=${JIT:-0.3}
BRIDGE=${BRIDGE:-one}
LOAD=${LOAD:-x}
CONF=${CONF:-sym}
ROOT=${ROOT:-out_slabconv/kg${KG}_${BRIDGE}_${LOAD}$CONF}
mkdir -p "$ROOT"

for n in $NS; do
  for st in und drn; do
    d="$ROOT/n$n/$st"; mkdir -p "$d"
    $PY elements_ccx/tests/make_slabconv.py "$d/m" "$n" "$st" "$KG" "$JIT" "$BRIDGE" "$LOAD" "$CONF" \
        > "$d/gen.log" 2>&1
  done
  # F-barES-FEM-T4 c=1 on the undrained cell only
  w="$ROOT/n$n/und_fbar1"; mkdir -p "$w"
  $PY elements_ccx/fbares.py "$ROOT/n$n/und/m_ccx.inp" "$w/m_ccx.inp" \
      --elset Sphere_Only --cycles 1 > "$w/gen.log" 2>&1

  for arm in und/m_ccx drn/m_ccx und_fbar1/m_ccx; do
    j="$ROOT/n$n/$arm"
    ( cd "$(dirname "$j")" && CCX_FBAR_C=1 SPAX_CCX_REAL=ccx_fbar \
        SPAX_CCX_MEMMAX=${MEM:-24G} OMP_NUM_THREADS=${NT:-8} \
        ccx_capped "$(basename "$j")" > solve.log 2>&1 ) || true
  done
done

$PY - "$ROOT" "$KG" "$LOAD" $NS <<'PYEOF'
import os, re, sys
root, kg, load, ns = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
EPS = 1.0e-3
AX = {'x': 0, 'y': 1, 'z': 2}[load]      # the DRIVEN component

def c1111(dat):
    """sum RF_x over the driven face / eps.  Also returns the traction check."""
    if not os.path.isfile(dat):
        return float('nan'), float('nan')
    blocks, cur, kind = {}, None, None
    for ln in open(dat):
        m = re.match(r'\s*(forces|displacements).*for set (\S+)', ln)
        if m:
            kind, cur = m.group(1), m.group(2).upper()
            blocks.setdefault((kind, cur), [])
            continue
        if cur and ln.strip() and not ln.lstrip().startswith(('*', 'S T')):
            f = ln.split()
            if len(f) >= 4:
                try:
                    blocks[(kind, cur)].append([float(x) for x in f[1:4]])
                except ValueError:
                    cur = None
    f1 = blocks.get(('forces', 'X1'), [])
    f0 = blocks.get(('forces', 'X0'), [])
    if not f1:
        return float('nan'), float('nan')
    rx = sum(r[AX] for r in f1)
    # the driven face and its opposite must carry equal and opposite load
    bal = abs(rx + sum(r[AX] for r in f0)) / abs(rx) if rx else float('nan')
    return rx / EPS, bal

print('\nK/G = %s  load along %s   R = C1111(und)/C1111(drn), '
      'denominator always plain C3D4' % (kg, load))
print('%-5s %-7s %14s %14s %14s %9s %9s'
      % ('n', 'el/slab', 'C1111 drn', 'C1111 und C3D4', 'C1111 und fbar',
         'R C3D4', 'R fbar'))
for n in ns:
    b = os.path.join(root, 'n' + n)
    d, bd = c1111(os.path.join(b, 'drn', 'm_ccx.dat'))
    u, bu = c1111(os.path.join(b, 'und', 'm_ccx.dat'))
    fb, bf = c1111(os.path.join(b, 'und_fbar1', 'm_ccx.dat'))
    print('%-5s %-7.0f %14.6e %14.6e %14.6e %9.4f %9.4f'
          % (n, 0.2 * int(n), d, u, fb, u / d if d else float('nan'),
             fb / d if d else float('nan')))
    for tag, v in (('drn', bd), ('und', bu), ('fbar', bf)):
        if v == v and v > 1e-6:
            print('      WARNING %s: sum RF_x over both faces is %.2e of the '
                  'driven reaction -- not in equilibrium' % (tag, v))
PYEOF
