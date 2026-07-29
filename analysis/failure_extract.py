"""Failure-onset extractor (abaqus python).

Extends scf_extract.py. For one uniaxial-tension (utx) ODB of a column slice,
compute per-MATRIX-element stress measures normalised by the macroscopic applied
stress  S11_bar = (1/V_RVE) * sum(S11 * Vel)  (voids contribute zero), and report
percentile distributions for two failure criteria:

  (1) Maximum-principal (tensile) stress concentration
        SCF = sigma_1 / S11_bar
      -> tensile cracking; first-failure macro stress = sigma_t / SCF_p99.

  (2) Mohr-Coulomb demand, normalised by S11_bar
        MCnorm = [ (sigma_1 - sigma_3) - (sigma_1 + sigma_3) * sin(phi) ] / S11_bar
      MC failure when  (sigma_1 - sigma_3) - (sigma_1 + sigma_3) sin(phi) >= 2 c cos(phi),
      so first-failure macro stress = 2 c cos(phi) / MCnorm_p99.
      phi (friction angle, deg) via env SPAX_MC_PHI_DEG (default 30).

Both measures are strength-independent: the cohesion c and tensile strength
sigma_t are applied OFFLINE at analysis time, so the cross-slice ranking (which
depth fails first) needs no strength assumption. P99 is the robust, mesh-stable
measure (matrix elements are linear C3D4 -> one constant stress each, so the
absolute peak SCF_max/MCnorm_max is a mesh-limited lower bound).

Usage:
  abaqus python failure_extract.py <odb_path> <L> <run_id> <out_csv>
"""
import sys, os, csv, math
import numpy as np
from odbAccess import openOdb

PHI_DEG = float(os.environ.get('SPAX_MC_PHI_DEG', '30'))
SIN_PHI = math.sin(math.radians(PHI_DEG))


def principals(s):
    s11, s22, s33, s12, s13, s23 = [float(x) for x in s]
    M = np.array([[s11, s12, s13], [s12, s22, s23], [s13, s23, s33]])
    w = np.linalg.eigvalsh(M)            # ascending: w[0]=sigma_3, w[2]=sigma_1
    return w[2], w[0]                     # (sigma_1, sigma_3)


def main():
    odb_path, L, run_id, out_csv = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4]
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
    s1, s3, wts = [], [], []
    for v in S.getSubset(region=mat).values:
        a, b = principals(v.data)
        s1.append(a); s3.append(b)
        wts.append(volmap.get(v.elementLabel, 1.0))
    s1 = np.array(s1); s3 = np.array(s3); wts = np.array(wts)

    scf = s1 / S11_bar                                        # tensile (max-principal)
    mc = ((s1 - s3) - (s1 + s3) * SIN_PHI) / S11_bar          # Mohr-Coulomb demand

    def pct(arr, q):
        return float(np.percentile(arr, q))

    row = {
        'run_id': run_id,
        'S11_bar': '%.6e' % S11_bar,
        'n_matrix': len(scf),
        'SCF_p50': '%.3f' % pct(scf, 50),
        'SCF_p90': '%.3f' % pct(scf, 90),
        'SCF_p99': '%.3f' % pct(scf, 99),
        'SCF_max': '%.3f' % float(scf.max()),
        'MC_phi_deg': '%.1f' % PHI_DEG,
        'MCnorm_p50': '%.3f' % pct(mc, 50),
        'MCnorm_p90': '%.3f' % pct(mc, 90),
        'MCnorm_p99': '%.3f' % pct(mc, 99),
        'MCnorm_max': '%.3f' % float(mc.max()),
        'volfrac_SCF_gt2': '%.4f' % (float(wts[scf > 2].sum()) / float(wts.sum())),
        'volfrac_SCF_gt3': '%.4f' % (float(wts[scf > 3].sum()) / float(wts.sum())),
    }
    odb.close()
    fields = ['run_id', 'S11_bar', 'n_matrix', 'SCF_p50', 'SCF_p90', 'SCF_p99',
              'SCF_max', 'MC_phi_deg', 'MCnorm_p50', 'MCnorm_p90', 'MCnorm_p99',
              'MCnorm_max', 'volfrac_SCF_gt2', 'volfrac_SCF_gt3']
    exists = os.path.isfile(out_csv)
    f = open(out_csv, 'a')
    w = csv.DictWriter(f, fieldnames=fields)
    if not exists:
        w.writeheader()
    w.writerow(row)
    f.close()
    print('%s: S11_bar=%.4e  SCF p99=%s max=%s  MCnorm(phi=%.0f) p99=%s  (n=%d)' % (
        run_id, S11_bar, row['SCF_p99'], row['SCF_max'], PHI_DEG, row['MCnorm_p99'], len(scf)))


if __name__ == '__main__':
    main()
