"""Extract E_eff and nu_eff from a POC ODB (run under Abaqus Python).

  abaqus python poc_extract.py poc_fluidcavity.odb
  abaqus python poc_extract.py poc_void.odb

E_eff : macro axial stress / macro axial strain, where
        sigma_x = (sum of RF1 on XMAX) / (L*L)   [total face reaction / area]
        eps_x   = disp / L
nu_eff: -<E22>/<E11> volume-averaged over the meshed matrix (solid only),
        the same convention as the RVE post-processor.

Also prints the cavity pressure/volume change for the fluid-cavity case so you
can sanity-check that the brine actually pressurised under load.
"""
from odbAccess import openOdb
import sys

L = 1.0
disp = 0.01 * L

odb = openOdb(path=sys.argv[1], readOnly=True)
step = odb.steps[list(odb.steps.keys())[0]]
frame = step.frames[-1]
inst = odb.rootAssembly.instances[list(odb.rootAssembly.instances.keys())[0]]

# --- E_eff from total reaction on XMAX ---
xmax = odb.rootAssembly.nodeSets['XMAX'] if 'XMAX' in odb.rootAssembly.nodeSets.keys() \
    else inst.nodeSets['XMAX']
rf = frame.fieldOutputs['RF'].getSubset(region=xmax)
Fx = sum(v.data[0] for v in rf.values)
sigma_x = Fx / (L * L)
eps_x = disp / L
E_eff = sigma_x / eps_x

# --- nu_eff from volume-averaged strain (solid only) ---
E = frame.fieldOutputs['E']
EVOL = frame.fieldOutputs['EVOL']
vol = {v.elementLabel: v.data for v in EVOL.values}
sE11 = sE22 = sE33 = sV = 0.0
for v in E.getSubset(position=4).values if False else E.values:   # CENTROID/IP avg
    w = vol.get(v.elementLabel, 0.0)
    if w == 0.0:
        continue
    sE11 += v.data[0] * w; sE22 += v.data[1] * w; sE33 += v.data[2] * w; sV += w
E11 = sE11 / sV; E22 = sE22 / sV; E33 = sE33 / sV
nu_eff = -0.5 * (E22 + E33) / E11 if abs(E11) > 1e-30 else float('nan')

print("file      : %s" % sys.argv[1])
print("E_eff     : %.4e Pa  (%.3f GPa)" % (E_eff, E_eff / 1e9))
print("nu_eff    : %.4f   (from <E22>,<E33>/<E11>)" % nu_eff)
print("<E11>     : %.4e   <E22>: %.4e   <E33>: %.4e" % (E11, E22, E33))

# --- cavity pressure / volume, if present ---
for hr in step.historyRegions.values():
    for key, ho in hr.historyOutputs.items():
        if key in ('PCAV', 'CVOL'):
            print("%-6s    : %s" % (key, ho.data[-1][1]))
odb.close()
