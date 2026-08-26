#!/bin/bash
# One increment against ten, on a workstation, for nothing.
#
# The ramp campaign solves its cells in a single increment. Every one of them is
# nlgeom OFF with linear elastic phases, and the first-order extractor reads the
# last frame only (SpaX_PostProcess.py:1381), so the nine intermediate field
# frames the old step card wrote were never opened by anything. Collapsing them
# is not an approximation: the response is proportional to the imposed
# displacement, so one increment lands on the same answer.
#
# That is an argument about linearity, and arguments about linearity are the
# ones worth checking against a solve. The check does NOT need the production
# mesh. Linearity is a property of the equations, not of the discretisation, so
# a coarse copy of the same cell answers the same question -- and answers it
# here, in minutes, rather than costing six cells of cluster billing for a
# yes/no. A fine mesh would buy a more accurate modulus, which is not what is
# being asked; both solves see the SAME mesh, and the only thing compared is
# one against the other.
#
# Two meshes rather than one, because a single coarse cell leaves open the
# reading that the agreement is itself a coarse-mesh artefact. If the two
# increments agree at 1.2 elements across the brine layer and again at 1.6,
# nothing about the mesh is carrying the result.
#
#   ./control_local.sh [workdir]
#
# Needs Abaqus on PATH (C:\SIMULIA\Commands on this machine) and the same
# Python environment SpaX_Standalone.py normally runs in.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=${1:-${CLAUDE_JOB_DIR:-/tmp}/local_control}
mkdir -p "$WORK"
cd "$WORK"

ABQ=${ABAQUS_CMD:-/c/SIMULIA/Commands/abaqus.bat}
CPUS=${ABQ_CPUS:-4}

# The deck: the phi = 0.104 condition of rve_rampctl.csv at two coarse meshes,
# drained and undrained. Undrained matters here more than drained -- the mesh
# gate (rve_layermesh) measured the near-incompressible brine in a thin layer as
# needing roughly twice the resolution the drained response does, so if any part
# of this is going to be sensitive to something other than the increment size,
# it is that one.
python - "$HERE" <<'PY'
import os, sys
import pandas as pd
here = sys.argv[1]
d = pd.read_csv(os.path.join(here, 'params', 'rve_rampctl.csv'))
d = d[d.run_id.str.contains('_s1')]
out = []
for lm, tag in ((0.0120, 'c'), (0.0090, 'm')):
    e = d.copy()
    e['L_mesh'] = lm
    e['run_id'] = e.run_id.str.replace('RAMPC_p104', 'LCTL_%s' % tag, regex=False)
    out.append(e)
pd.concat(out).to_csv('rve_local_ctl.csv', index=False)
print('deck: rve_local_ctl.csv')
for lm in (0.0120, 0.0090):
    print('  L_mesh %.4f  %.0fk cube-equivalents  %.1f elements across the layer'
          % (lm, (0.5 / lm) ** 3 / 1e3, 0.01473 / lm))
PY

# Generated ONCE, with the flag on. The ten-increment twin is then made by
# editing the increment line of the finished deck rather than by regenerating,
# so the mesh, the packing and the periodic equations are bit-identical between
# the pair and the increment size is the only thing that differs. Regenerating
# would have repacked from a different seed and confounded the comparison with
# seed scatter -- which, at these cells, is the same size as the effect being
# looked for.
if ! ls Job-LCTL_*.inp >/dev/null 2>&1; then
  echo "generating (slow: this is the Gmsh half)..."
  OMP_NUM_THREADS=1 SPAX_LINEAR_ONE_STEP=1 \
    python "$HERE/SpaX_Standalone.py" rve_local_ctl.csv . > gen.log 2>&1
fi
echo "generated: $(ls Job-LCTL_*.inp 2>/dev/null | wc -l) decks"

for f in Job-LCTL_*.inp; do
  case "$f" in *TEN*) continue;; esac
  g="${f%.inp}_TEN.inp"
  [ -e "$g" ] && continue
  sed 's/^1\., 1\., 1e-10, 1\./0.1, 1., 1e-10, 0.1/' "$f" > "$g"
  p="${f%.inp}_periodic_pairs.csv"
  [ -e "$p" ] && cp "$p" "${g%.inp}_periodic_pairs.csv"
  grep -q '^0\.1, 1\., 1e-10, 0\.1' "$g" || echo "  WARNING: $g did not take the edit"
done

for f in Job-LCTL_*.inp; do
  j="${f%.inp}"
  [ -e "${j}.odb" ] && { echo "skip $j (odb exists)"; continue; }
  echo "solving $j ..."
  /usr/bin/time -f "  %e s wall" "$ABQ" job="$j" input="$f" cpus=$CPUS interactive \
    > "${j}.solvelog" 2>&1 || echo "  FAILED -- see ${j}.solvelog"
done

echo
echo "comparing:"
"$ABQ" python "$HERE/analysis/local_control_compare.py" "$WORK"
