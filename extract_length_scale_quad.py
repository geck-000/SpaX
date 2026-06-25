#!/usr/bin/env python3
"""Extract the MCST length scale l from the QUADRATIC (C3D10H) bending sweep,
calibrated by the quadratic homogeneous-cube baseline.

Rationale (see memory bending-length-scale-artifact):
  * Linear C3D4H bending shear-locks (+12..28% spurious stiffness, worst at
    small L). Quadratic C3D10H removes almost all of it.
  * A homogeneous cube has true l=0, so its quadratic D_ratio = f_quad(L) is the
    residual cube-vs-thin-plate continuum shape factor (small, ~1.0..1.08).
  * For each porous RVE the couple-stress excess uses the CALIBRATED reference
    D_classical = f_quad(L) * E*_plate * L^4/12, so:
        l^2 = (D_rve - D_classical) / (G* * L^2)
    A REAL length scale is SIZE-INDEPENDENT (l const across L).  l ~ L means a
    residual artifact;  l ~ 0 means no measurable length scale.
  * Equivalently S-1 = D_rve/D_classical - 1 should scale as (l/L)^2; the slope
    of (S-1) vs (1/L^2) gives l^2 * (E*_plate/(12 G*)).
First-order E*,nu*,G* come from the LINEAR study (uniaxial is constant-strain,
exact for C3D4H, so linear is fine and already converged).

    python extract_length_scale_quad.py [porous_q_results.csv]
"""
import csv, sys, os, math
import numpy as np

R_AVG = 0.04; D_INC = 2*R_AVG

def fnum(s):
    try: return float(s)
    except: return float('nan')

def load(fn):
    return list(csv.DictReader(open(fn, newline=''))) if os.path.exists(fn) else []

def main():
    por_fn = sys.argv[1] if len(sys.argv) > 1 else "results_porous_q.csv"
    por = load(por_fn)
    if not por:
        raise SystemExit("no porous quadratic results at %s (run solve+post first)" % por_fn)

    # converged first-order moduli (LINEAR study, size-independent)
    ls = load("results_lscale.csv")
    Es=[fnum(r['E_eff']) for r in ls if fnum(r['E_eff'])>0]
    NU=[fnum(r['nu_eff']) for r in ls if fnum(r['E_eff'])>0]
    G =[fnum(r['G_eff']) for r in ls if fnum(r['E_eff'])>0]
    Estar, nustar, Gstar = np.mean(Es), np.mean(NU), np.mean(G)
    Eplate = Estar/(1-nustar**2)
    print("First-order (linear, converged): E*=%.3g GPa  nu*=%.3f  G*=%.3g GPa"%(
        Estar/1e9, nustar, Gstar/1e9))

    # quadratic homogeneous baseline f_quad(L). Default to the xz homogeneous
    # runs; override with one or more CSVs as extra args (e.g. an xy baseline
    # for an xy-bending study — equal to xz for an isotropic cube, but pass it
    # explicitly to verify rather than assume).
    base_files = sys.argv[2:] if len(sys.argv) > 2 else [
        "results_homog_q_small.csv", "results_homog_q.csv"]
    fq = {}
    base_rows = []
    for bf in base_files:
        base_rows += load(bf)
    for r in base_rows:
        D=fnum(r['D_rve']); E=fnum(r['E_eff']); nu=fnum(r['nu_eff']); L=fnum(r['L'])
        if D>0 and E>0:
            fq[round(L,3)] = D/(E/(1-nu**2)*L**4/12.0)
    if not fq:
        raise SystemExit("no quadratic homogeneous baseline (results_homog_q*.csv)")
    print("Quadratic homogeneous baseline f_quad(L):")
    for L in sorted(fq): print("   L=%.2f (L/d=%.1f)  f=%.4f"%(L, L/D_INC, fq[L]))
    print()

    # per-RVE l with the calibrated reference. Use each RVE's OWN first-order
    # moduli when present (the -utx/-ss13 decks were solved) so the reference
    # tracks that realization's actual stiffness/VoF — this removes the VoF
    # confound at small L (cap rejection lowers achievable VoF differently per
    # size). Fall back to the converged ensemble moduli if a row lacks them.
    n_perrve = 0
    bysize = {}
    for r in por:
        D=fnum(r.get('D_rve','nan')); L=round(fnum(r['L']),3)
        if not (D>0): continue
        f = fq.get(L)
        if f is None:  # nearest baseline size (f is ~scale-invariant for quad)
            f = fq[min(fq, key=lambda x: abs(x-L))]
        Er=fnum(r.get('E_eff','nan')); nur=fnum(r.get('nu_eff','nan')); Gr=fnum(r.get('G_eff','nan'))
        if Er>0 and Gr>0 and nur==nur:
            Ep, Gv = Er/(1-nur**2), Gr; n_perrve += 1
        else:
            Ep, Gv = Eplate, Gstar
        Dcl = f*Ep*L**4/12.0
        l2 = (D - Dcl)/(Gv*L**2)
        l = math.sqrt(l2) if l2>0 else 0.0
        S = D/Dcl
        bysize.setdefault(L, []).append((l, S, D))
    print("Per-RVE moduli used for %d/%d RVEs (rest fall back to converged E*,G*).\n"
          % (n_perrve, sum(len(v) for v in bysize.values())))
    sizes = sorted(bysize)

    print("%-6s %5s %4s | %-22s | %-20s | %s"%(
        "L","L/d","n","D_rve mean (CoV%)","l (mean std CoV%)","S-1 mean"))
    Ldd=[]; Sm1=[]; invL2=[]
    for L in sizes:
        a=bysize[L]; ls_=np.array([x[0] for x in a]); Ss=np.array([x[1] for x in a]); Ds=np.array([x[2] for x in a])
        lm, lsd = ls_.mean(), (ls_.std(ddof=1) if len(ls_)>1 else 0.0)
        print("%-6.2f %5.1f %4d | %.3e (%4.1f%%)   | %6.4f %6.4f %5.1f | %+.4f"%(
            L, L/D_INC, len(a), Ds.mean(), Ds.std(ddof=1)/Ds.mean()*100 if len(Ds)>1 else 0,
            lm, lsd, lsd/lm*100 if lm>0 else float('nan'), Ss.mean()-1))
        Ldd.append(L); Sm1.append(Ss.mean()-1); invL2.append(1.0/L**2)

    # verdict: is l size-independent (real) or tracking L (artifact)?
    print("\nl/L by size:", "  ".join("%.3f"%(np.mean([x[0] for x in bysize[L]])/L) for L in sizes))
    lmeans=np.array([np.mean([x[0] for x in bysize[L]]) for L in sizes])
    if len(sizes)>=3:
        # regress (S-1) vs 1/L^2 through origin: slope = l^2 * Eplate/(12 G*)
        x=np.array(invL2); y=np.array(Sm1)
        slope = np.dot(x,y)/np.dot(x,x)
        l2_fit = slope*12*Gstar/Eplate
        l_fit = math.sqrt(l2_fit) if l2_fit>0 else 0.0
        # also check flatness of l(L)
        spread = (lmeans.max()-lmeans.min())/lmeans.mean() if lmeans.mean()>0 else float('nan')
        print("\nRegression (S-1) ~ slope/L^2  (couple-stress form):")
        print("   slope=%.3e  ->  l_fit = %.4f  (l_fit/d = %.2f)"%(slope, l_fit, l_fit/D_INC))
        print("   l(L) relative spread across sizes: %.0f%%"%(spread*100))
        if l_fit < 0.3*R_AVG and abs(np.mean(Sm1)) < 0.03:
            print("   => l ~ 0: NO measurable microstructural length scale (S-1 within noise).")
        elif spread < 0.25:
            print("   => l SIZE-INDEPENDENT (spread <25%%): real length scale l ~ %.4f."%lmeans.mean())
        else:
            print("   => l still size-dependent: inspect (S-1) vs 1/L^2 linearity / add seeds.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
