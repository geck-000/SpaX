"""Stress-concentration extractor (abaqus python).

For one uniaxial-tension (utx) ODB, compute the per-element maximum principal
stress in the MATRIX, normalised by the macroscopic applied stress
S11_bar = (1/V_RVE) * sum(S11 * Vel) over all solid elements (voids contribute
zero), i.e. the local stress-concentration factor SCF = sigma_1 / S11_bar.

Reports the distribution (mean, P50/P90/P99, max) and the matrix volume fraction
above SCF thresholds. Appends one row to the output CSV.

Usage:
  abaqus python scf_extract.py <odb_path> <L> <run_id> <out_csv> [dump.npz]

If a fifth argument is given, the per-element (scf, volume) pairs are written to
that .npz instead of being discarded after the percentiles are formed. A Weibull
weakest-link measure needs the whole field, not a percentile of it,

    P_f = 1 - exp[ -(1/V0) * INT (sigma_1/sigma_0)^m dV ]

and the integral is just sum(scf^m * vol) over the matrix elements. Dumping the
arrays keeps the one expensive pass (opening the ODB, which needs an Abaqus
licence) separate from the m-sweep, which is then plain offline numpy --- see
analysis/weibull_sensitivity.py.

NOTE: matrix elements are linear C3D4 (one constant stress per element), so the
absolute peak is a mesh-limited lower bound; the P99 and the cross-case ranking
are the robust, comparable measures (identical mesh resolution across cases).
The same caveat applies to the Weibull integral, which for large m is dominated
by the few most-stressed elements and so inherits their mesh sensitivity; that
is a property of the measure, and weibull_sensitivity.py reports it as a
function of m rather than hiding it in a single number.
"""
import sys, os, csv
import numpy as np
from odbAccess import openOdb

def max_principal(s):
    s11, s22, s33, s12, s13, s23 = [float(x) for x in s]
    M = np.array([[s11, s12, s13], [s12, s22, s23], [s13, s23, s33]])
    w = np.linalg.eigvalsh(M)          # ascending
    return w[2]                         # maximum principal

def main():
    odb_path, L, run_id, out_csv = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4]
    dump_path = sys.argv[5] if len(sys.argv) > 5 else None
    odb = openOdb(odb_path, readOnly=True)
    step = max(odb.steps.values(), key=lambda s: len(s.frames))
    fr = step.frames[-1]
    S = fr.fieldOutputs['S']
    keys = fr.fieldOutputs.keys()
    EVOL = fr.fieldOutputs['EVOL'] if 'EVOL' in keys else None

    inst = None
    for n in odb.rootAssembly.instances.keys():
        if 'PBC' in n.upper() or n.upper() == 'ASSEMBLY':
            continue
        inst = odb.rootAssembly.instances[n]; break

    volmap = {}
    if EVOL is not None:
        for v in EVOL.getSubset(region=inst).values:
            volmap[v.elementLabel] = abs(float(v.data))

    # macro applied stress: volume-average S11 over the full RVE (L^3, voids -> 0)
    num = 0.0
    for v in S.getSubset(region=inst).values:
        ve = volmap.get(v.elementLabel, 1.0)
        num += float(v.data[0]) * ve
    S11_bar = num / (L ** 3)

    mat = inst.elementSets['MATRIX_ONLY']
    sp1, wts = [], []
    for v in S.getSubset(region=mat).values:
        sp1.append(max_principal(v.data))
        wts.append(volmap.get(v.elementLabel, 1.0))
    sp1 = np.array(sp1); wts = np.array(wts)
    scf = sp1 / S11_bar

    def pct(q):
        return float(np.percentile(scf, q))
    row = {
        'run_id': run_id,
        'S11_bar': '%.6e' % S11_bar,
        'n_matrix': len(scf),
        'SCF_mean': '%.3f' % float(scf.mean()),
        'SCF_p50': '%.3f' % pct(50),
        'SCF_p90': '%.3f' % pct(90),
        'SCF_p99': '%.3f' % pct(99),
        'SCF_max': '%.3f' % float(scf.max()),
        'volfrac_SCF_gt2': '%.4f' % (float(wts[scf > 2].sum()) / float(wts.sum())),
        'volfrac_SCF_gt3': '%.4f' % (float(wts[scf > 3].sum()) / float(wts.sum())),
    }
    if dump_path is not None:
        d = os.path.dirname(dump_path)
        if d and not os.path.isdir(d):
            os.makedirs(d)
        # scf is per matrix element, vol the element volume (the dV of the
        # Weibull integral); S11_bar and L carry the normalisation so the
        # offline analyzer can rebuild absolute stresses if it needs to.
        np.savez_compressed(dump_path, scf=scf, vol=wts,
                            S11_bar=S11_bar, L=L, run_id=run_id)
        print('  dumped %d matrix elements -> %s' % (len(scf), dump_path))

    odb.close()
    fields = ['run_id', 'S11_bar', 'n_matrix', 'SCF_mean', 'SCF_p50', 'SCF_p90',
              'SCF_p99', 'SCF_max', 'volfrac_SCF_gt2', 'volfrac_SCF_gt3']
    exists = os.path.isfile(out_csv)
    f = open(out_csv, 'a')
    w = csv.DictWriter(f, fieldnames=fields)
    if not exists:
        w.writeheader()
    w.writerow(row)
    f.close()
    print('%s: S11_bar=%.4e  SCF p99=%s  max=%s  (n=%d)' % (
        run_id, S11_bar, row['SCF_p99'], row['SCF_max'], len(scf)))

if __name__ == '__main__':
    main()
