#!/usr/bin/env python3
"""Homogeneous-cube calibration of the RVE bending size effect.

The raw bending study found D_rve/D_classical ~ 1.11, FLAT in box size L, so the
extracted MCST length scale l just tracked L (l ~ 0.167 L) -- the signature of a
NON-microstructural artifact, not a couple-stress size effect.

A homogeneous cube (no inclusions) has, by construction, NO intrinsic length
scale: its true l = 0. Run through the IDENTICAL pipeline (same Lesicar bending
BC, same D_rve extraction, same coarse C3D4H bulk mesh), whatever D_ratio it
produces is pure artifact -- the cube-vs-thin-plate continuum shape correction
PLUS any linear-tet bending stiffening of the coarse bulk. This script:

  1. SIZE sweep  (results_homog.csv, L_mesh=0.033 matching the porous runs):
     D_ratio_homog(L). If ~constant -> scale-invariant artifact (expected).
  2. MESH sweep  (results_homog_meshconv.csv, fixed L, L_mesh refined):
     D_ratio_homog -> 1 as mesh refines  => the offset is linear-tet LOCKING
     (a fixable mesh artifact);  plateaus > 1  => genuine continuum shape effect.
  3. CORRECTION: use D_ratio_homog as the bending reference for the POROUS RVEs
     instead of thin-plate D_classical, and recompute l. The microstructural
     couple-stress signal is the EXCESS  (D_ratio_porous - D_ratio_homog).

    python analyze_homog_calib.py
"""
import csv, math, os
import numpy as np

R_AVG = 0.04; D_INC = 2*R_AVG

def fnum(s):
    try: return float(s)
    except: return float('nan')

def load(fn):
    return list(csv.DictReader(open(fn, newline=''))) if os.path.exists(fn) else []

def dratio(r):
    """D_rve/D_classical for a fully-solved row (needs E_eff, nu_eff, D_rve)."""
    D = fnum(r.get('D_rve','nan')); E = fnum(r.get('E_eff','nan'))
    nu = fnum(r.get('nu_eff','nan')); L = fnum(r['L'])
    if not (D > 0 and E > 0): return None
    Dcl = E/(1-nu**2) * L**4/12.0
    return D/Dcl

def main():
    homog = load("results_homog.csv")
    mesh  = load("results_homog_meshconv.csv")

    # ---- 1. size sweep: D_ratio_homog(L) at L_mesh=0.033 ----
    print("="*70)
    print("1. HOMOGENEOUS cube, size sweep (L_mesh=0.033, true l=0 by construction)")
    print("="*70)
    print("%-6s %5s | %-12s %-12s | %s"%("L","L/d","E_eff(tens)","D_ratio","E_bend/E_tens"))
    homog_by_L = {}
    for r in sorted(homog, key=lambda x: fnum(x['L'])):
        dr = dratio(r)
        if dr is None:
            print("%-6.2f %5.1f | UNSOLVED (E_eff or D_rve missing)"%(fnum(r['L']), fnum(r['L'])/D_INC))
            continue
        L = round(fnum(r['L']),3); homog_by_L[L] = dr
        print("%-6.2f %5.1f | %.3e   %.4f       | %.4f"%(
            L, L/D_INC, fnum(r['E_eff']), dr, dr))
    if homog_by_L:
        vals = list(homog_by_L.values())
        print("  -> D_ratio_homog mean %.4f, spread %.4f (max-min). %s"%(
            np.mean(vals), max(vals)-min(vals),
            "scale-INVARIANT artifact" if (max(vals)-min(vals))<0.02 else "size-dependent (!)"))

    # ---- 2. mesh sweep: does D_ratio_homog -> 1 as we refine? ----
    print("\n"+"="*70)
    print("2. HOMOGENEOUS cube, MESH convergence (fixed L=0.64)")
    print("="*70)
    print("%-10s %-9s | %-10s %s"%("run_id","L_mesh","D_ratio","verdict-trend"))
    mc = []
    for r in sorted(mesh, key=lambda x: -fnum(x['L_mesh'])):
        dr = dratio(r)
        if dr is None:
            print("%-10s %-9s | UNSOLVED"%(r['run_id'], r['L_mesh'])); continue
        mc.append((fnum(r['L_mesh']), dr))
        print("%-10s %-9s | %.4f"%(r['run_id'], r['L_mesh'], dr))
    if len(mc) >= 2:
        mc.sort(reverse=True)  # coarse -> fine
        trend = mc[-1][1] - mc[0][1]
        print("  -> coarsest %.4f  finest %.4f  (delta %+.4f)"%(mc[0][1], mc[-1][1], trend))
        if mc[-1][1] < 1.02:
            print("  => converges toward 1.0: the offset is LINEAR-TET LOCKING (mesh artifact, fixable).")
        elif abs(trend) < 0.01:
            print("  => plateaus above 1.0: GENUINE cube-vs-thin-plate continuum shape effect.")
        else:
            print("  => still moving; refine further to separate locking from shape.")

    # ---- 3. correct the POROUS l using the homogeneous baseline ----
    print("\n"+"="*70)
    print("3. POROUS length scale, CORRECTED by homogeneous baseline")
    print("="*70)
    lscale = load("results_lscale.csv"); bend = load("results_bending.csv")
    # converged porous first-order moduli (size-independent)
    Es=[fnum(r['E_eff']) for r in lscale if fnum(r['E_eff'])>0]
    NU=[fnum(r['nu_eff']) for r in lscale if fnum(r['E_eff'])>0]
    G =[fnum(r['G_eff']) for r in lscale if fnum(r['E_eff'])>0]
    Estar, nustar, Gstar = np.mean(Es), np.mean(NU), np.mean(G)
    Eplate = Estar/(1-nustar**2)
    if not homog_by_L:
        print("  (no homogeneous results yet -- run the solve/post first)"); return 0
    # use a single representative baseline: mean over sizes if scale-invariant
    f_base = np.mean(list(homog_by_L.values()))
    print("  porous moduli: E*=%.3g GPa  nu*=%.3f  G*=%.3g GPa"%(Estar/1e9,nustar,Gstar/1e9))
    print("  homogeneous artifact factor f = %.4f (thin-plate over-stiffness baseline)\n"%f_base)

    # gather porous D_rve per size
    pby = {}
    for r in lscale+bend:
        D=fnum(r.get('D_rve','nan')); L=fnum(r['L'])
        if D>0: pby.setdefault(round(L,3),[]).append(D)

    def lcalc(D, L, fref):
        Dcl = fref * Eplate * L**4/12.0     # artifact-corrected classical reference
        l2 = (D - Dcl)/(Gstar*L**2)
        return (math.sqrt(l2) if l2>0 else 0.0), D/(Eplate*L**4/12.0)

    print("%-6s %5s %4s | %-18s | %-22s"%("L","L/d","n","RAW l (l/L)","CORRECTED l (l/L)"))
    for L in sorted(pby):
        Ds = np.array(pby[L])
        # raw: reference = thin-plate (f=1); corrected: reference = f_base
        raw = np.array([lcalc(D,L,1.0)[0] for D in Ds])
        cor = np.array([lcalc(D,L,f_base)[0] for D in Ds])
        drp = np.mean([lcalc(D,L,1.0)[1] for D in Ds])
        print("%-6.2f %5.1f %4d | %6.4f (%.3f)    | %6.4f (%.3f)   [D_ratio_por=%.3f vs f=%.3f]"%(
            L, L/D_INC, len(Ds), raw.mean(), raw.mean()/L, cor.mean(), cor.mean()/L, drp, f_base))
    print("\n  If CORRECTED l ~ 0 and flat -> no measurable microstructural length scale.")
    print("  If CORRECTED l plateaus at a size-independent value -> that is the real l.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
