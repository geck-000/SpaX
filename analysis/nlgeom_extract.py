"""Large-deformation (nlgeom) extractor for study #8 -- run under `abaqus python`.

For each uniaxial-X ODB it reads the driven reference point RP-1's imposed
displacement U1 and reaction force RF1 at EVERY increment and forms the
unambiguous homogenized NOMINAL stress-strain path
    eps_nom(t) = U1(t) / L        (engineering strain, reference length L)
    sigma_nom(t) = RF1(t) / L^2   (1st-Piola/nominal macro stress, reference area)
The RP-1 reaction is the resultant traction on the driven face, so this is
boundary-based and free of the volume-averaging ambiguity that afflicts a
finite-strain Cauchy average.

Outputs two CSVs:
  <summary>  one row per RVE: E0 (initial tangent = linear modulus), E_sec at max
             strain, tangent at max strain, and the secant/E0 ratio -- the
             geometric-nonlinearity signature.
  <curves>   the full (run_id, eps, sigma) path, one row per increment, for plots.

Usage:  abaqus python nlgeom_extract.py <rve_csv> <workdir> <summary_csv> <curves_csv>
"""
import sys, csv, os
from odbAccess import openOdb


def _rp_region(odb, rp_name):
    sets = odb.rootAssembly.nodeSets
    for key in (rp_name, rp_name.upper(), rp_name.replace('-', '_').upper()):
        if key in sets.keys():
            return sets[key]
    return None


def _val(frame, field_key, region, dof):
    """Component `dof` of `field_key` at the single RP node in `frame`, or None."""
    try:
        f = frame.fieldOutputs[field_key].getSubset(region=region)
        if len(f.values) == 0:
            return None
        d = f.values[0].data            # tuple/array of components
        return float(d[dof])
    except Exception:
        return None


def _slope(xs, ys):
    """Least-squares slope of ys vs xs (n>=2), no intercept assumption."""
    n = len(xs)
    sx = sum(xs); sy = sum(ys)
    sxx = sum(x * x for x in xs); sxy = sum(x * y for x, y in zip(xs, ys))
    den = n * sxx - sx * sx
    return (n * sxy - sx * sy) / den if den != 0 else 0.0


def extract(odb_path, L):
    """Return (eps[], sigma[]) nominal path over all increments for a utx ODB."""
    odb = openOdb(path=odb_path, readOnly=True)
    try:
        step = odb.steps[list(odb.steps.keys())[0]]
        rp = _rp_region(odb, 'RP-1')
        eps, sig = [0.0], [0.0]
        for fi, frame in enumerate(step.frames):
            if fi == 0:
                continue
            u1 = _val(frame, 'U', rp, 0)
            rf1 = _val(frame, 'RF', rp, 0)
            if u1 is None or rf1 is None:
                continue
            eps.append(u1 / L)
            sig.append(rf1 / (L * L))
    finally:
        odb.close()
    return eps, sig


def moduli(eps, sig):
    """E0 (initial tangent over |eps|<=25% of max), secant at max, tangent at max."""
    n = len(eps)
    emax = max(abs(e) for e in eps) if n else 0.0
    if emax == 0 or n < 3:
        return 0.0, 0.0, 0.0
    lo = [(e, s) for e, s in zip(eps, sig) if 0 < abs(e) <= 0.25 * emax]
    if len(lo) >= 2:
        E0 = _slope([e for e, _ in lo], [s for _, s in lo])
    else:
        # fall back to the first non-zero increment
        e1 = next(e for e in eps if e != 0.0)
        s1 = sig[eps.index(e1)]
        E0 = s1 / e1
    E_sec = sig[-1] / eps[-1] if eps[-1] != 0 else 0.0
    E_tan = _slope(eps[-3:], sig[-3:]) if n >= 3 else E_sec
    return E0, E_sec, E_tan


def main():
    rve_csv, workdir, summary_csv, curves_csv = sys.argv[1:5]
    rows = list(csv.DictReader(open(rve_csv)))
    srows, crows = [], []
    for r in rows:
        rid = r['run_id']
        L = float(r.get('L', 0.5) or 0.5)
        nlg = r.get('nlgeom_flag', 'OFF')
        odb = os.path.join(workdir, 'Job-%s-utx.odb' % rid)
        if not os.path.exists(odb):
            print('MISSING %s' % odb); continue
        eps, sig = extract(odb, L)
        E0, E_sec, E_tan = moduli(eps, sig)
        emax = eps[-1] if eps else 0.0
        srows.append(dict(run_id=rid, nlgeom=nlg, L=L, eps_max=emax,
                          E0=E0, E_sec=E_sec, E_tan_end=E_tan,
                          sec_over_E0=(E_sec / E0 if E0 else 0.0)))
        for e, s in zip(eps, sig):
            crows.append(dict(run_id=rid, nlgeom=nlg, eps=e, sigma=s))
        print('%-13s nlgeom=%-3s eps_max=%+.4f  E0=%.4e  E_sec=%.4e  sec/E0=%.4f'
              % (rid, nlg, emax, E0, E_sec, (E_sec / E0 if E0 else 0.0)))

    with open(summary_csv, 'w') as f:
        w = csv.DictWriter(f, fieldnames=['run_id', 'nlgeom', 'L', 'eps_max',
                                          'E0', 'E_sec', 'E_tan_end', 'sec_over_E0'])
        w.writeheader(); w.writerows(srows)
    with open(curves_csv, 'w') as f:
        w = csv.DictWriter(f, fieldnames=['run_id', 'nlgeom', 'eps', 'sigma'])
        w.writeheader(); w.writerows(crows)
    print('wrote %s and %s' % (summary_csv, curves_csv))


if __name__ == '__main__':
    main()
