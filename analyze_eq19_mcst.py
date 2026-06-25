#!/usr/bin/env python3
"""MCST length-scale extraction by the Choi/Lee/Sim Eq.19 method, plus a direct
bending-vs-first-order stiffening view.

Background. Weak/soft inclusions SOFTEN the effective modulus -- a FIRST-ORDER
(classical) effect: the uniaxial RVEs already report E* < E_matrix. The MODIFIED
COUPLE STRESS (MCST) size effect is something ON TOP of that: under bending the
RVE mobilises a strain GRADIENT, and IF the microstructure carries a length scale
l, the apparent bending modulus exceeds the first-order modulus by

    E_app(h) = E_classical + 12 mu (l/h)^2            (Choi et al., Mater.&Des. 2022, Eq.19)

i.e. it rises linearly in 1/h^2. We fit this across RVE sizes (h = L):
  * SLOPE  = 12 mu l^2   -> l  (reference-free: the intercept absorbs the
    plate/beam classical modulus, so l does NOT depend on that choice);
  * INTERCEPT = the classical (h->infinity) bending modulus -> should match the
    first-order modulus (E* for a beam, E*/(1-nu^2) for a plate), a consistency
    check on whether the RVE bends as a plate or a beam.

The "MCST stiffening" the question asks for is then made explicit per RVE:
    E_bend_material = (D_rve*12/L^4) * (1 - nu^2)      # plate factor removed
    ratio = E_bend_material / E_first_order            # >1 == bending stiffer
A ratio that is ~1 (no trend) means the bending response is fully explained by the
SOFTENED first-order modulus -> no MCST stiffening. A ratio >1 that GROWS as L
shrinks is the couple-stress signature.

    python analyze_eq19_mcst.py results.csv [results2.csv ...]

Each CSV needs D_rve, L, and (for the per-RVE view) E_eff, nu_eff, G_eff.
"""
import csv, sys, math, collections
import numpy as np

R_AVG = 0.04; D_INC = 2 * R_AVG

def fnum(s):
    try: return float(s)
    except: return float('nan')

def load(fn):
    return list(csv.DictReader(open(fn, newline='')))

def eq19_fit(rows, Gstar_fallback):
    """Returns (sizes, per-size dict, fit dict). per-size dict[L] = (Eapp_mean,
    Ebend_mat_mean, Efo_mean, ratio_mean, ratio_std, n)."""
    by = collections.defaultdict(list)
    for r in rows:
        D = fnum(r.get('D_rve', 'nan')); L = round(fnum(r['L']), 3)
        if not (D > 0): continue
        nu = fnum(r.get('nu_eff', 'nan')); E = fnum(r.get('E_eff', 'nan'))
        G = fnum(r.get('G_eff', 'nan'))
        nu = nu if nu == nu else 0.348
        Eapp = D * 12.0 / L**4                 # apparent (plate) bending modulus
        Ebend_mat = Eapp * (1.0 - nu**2)       # plate factor removed -> material
        Efo = E if E > 0 else float('nan')
        ratio = Ebend_mat / Efo if Efo > 0 else float('nan')
        by[L].append((Eapp, Ebend_mat, Efo, ratio, G if G > 0 else Gstar_fallback))
    sizes = sorted(by)
    per = {}
    for L in sizes:
        a = np.array([(x[0], x[1], x[2], x[3], x[4]) for x in by[L]], dtype=float)
        per[L] = (np.nanmean(a[:,0]), np.nanmean(a[:,1]), np.nanmean(a[:,2]),
                  np.nanmean(a[:,3]), np.nanstd(a[:,3], ddof=1) if len(a) > 1 else 0.0,
                  np.nanmean(a[:,4]), len(a))
    # Eq.19 regression: Eapp vs 1/L^2
    x = np.array([1.0/L**2 for L in sizes])
    y = np.array([per[L][0] for L in sizes])
    Gm = np.nanmean([per[L][5] for L in sizes])
    A = np.vstack([np.ones_like(x), x]).T
    (b0, b1) = np.linalg.lstsq(A, y, rcond=None)[0]
    l2 = b1 / (12.0 * Gm)
    fit = dict(intercept=b0, slope=b1, Gm=Gm, l2=l2,
               l=(math.sqrt(l2) if l2 > 0 else float('nan')))
    return sizes, per, fit

def main():
    if len(sys.argv) < 2:
        print(__doc__); return 1
    # converged moduli for the consistency check (from the first-order study)
    ls = load("results_lscale.csv") if __import__('os').path.exists("results_lscale.csv") else []
    if ls:
        E=[fnum(r['E_eff']) for r in ls if fnum(r['E_eff'])>0]
        NU=[fnum(r['nu_eff']) for r in ls if fnum(r['E_eff'])>0]
        G=[fnum(r['G_eff']) for r in ls if fnum(r['E_eff'])>0]
        Estar, nustar, Gstar = np.mean(E), np.mean(NU), np.mean(G)
    else:
        Estar, nustar, Gstar = float('nan'), 0.348, float('nan')

    for fn in sys.argv[1:]:
        rows = load(fn)
        sizes, per, fit = eq19_fit(rows, Gstar)
        print("="*78)
        print("FILE: %s" % fn)
        print("="*78)
        print("%-5s %4s %4s | %-9s %-12s %-10s | %s" % (
            "L","L/d","n","E_app(GPa)","E_bend_mat","E_1storder","bend/1st (CoV%)"))
        for L in sizes:
            Eapp,Ebm,Efo,rat,rsd,Gm,n = per[L]
            cov = rsd/rat*100 if rat==rat and rat>0 else float('nan')
            print("%-5.2f %4.1f %4d | %8.3f  %8.3f     %8.3f   | %.3f (%4.1f%%)" % (
                L, L/D_INC, n, Eapp/1e9, Ebm/1e9, Efo/1e9, rat, cov))
        # Eq.19 verdict
        l2 = fit['l2']
        print("\nEq.19 fit  E_app = E0 + 12*mu*(l/L)^2   (mu=G*=%.3g GPa):" % (fit['Gm']/1e9))
        print("  slope = %+.3e   intercept E0 = %.3f GPa" % (fit['slope'], fit['intercept']/1e9))
        if not math.isnan(Estar):
            print("  intercept vs first-order:  E* = %.3f GPa   E*/(1-nu^2) = %.3f GPa  -> RVE bends as a %s"
                  % (Estar/1e9, Estar/(1-nustar**2)/1e9,
                     "PLATE" if abs(fit['intercept']-Estar/(1-nustar**2)) < abs(fit['intercept']-Estar) else "BEAM"))
        if l2 > 0:
            print("  => l = %.4f  (l/d = %.2f)   [POSITIVE slope -> MCST stiffening present]" % (
                math.sqrt(l2), math.sqrt(l2)/D_INC))
        else:
            print("  => slope <= 0  -> l imaginary: NO MCST stiffening (bending NOT stiffer than first-order)")
        # direct stiffening answer
        ratios = [per[L][3] for L in sizes]
        print("\nMCST stiffening (bending material modulus / first-order modulus):")
        print("  by size:", "  ".join("%.3f" % r for r in ratios),
              "->", "all ~1 or <1: bending is NOT stiffer than the softened first-order E"
              if max(ratios) < 1.03 else
              "some >1: possible bending stiffening (check size trend & noise)")
        print()
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
