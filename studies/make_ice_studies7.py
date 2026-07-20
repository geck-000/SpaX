"""Seventh battery of sea-ice studies (2026-07-06).

  #8 nlgeom -> rve_nlgeom_lin.csv / rve_nlgeom_ten.csv / rve_nlgeom_cmp.csv
     (one geometrically-nonlinear / large-deformation study)

Probes whether GEOMETRIC nonlinearity (finite-strain kinematics, nlgeom=ON)
changes the homogenized response of the sea-ice RVE relative to the small-strain
linear result -- and whether the porous warm base shows a tension/compression
asymmetry (pores elongate in tension, close in compression). The matrix stays
linear-elastic; only the kinematics are made finite.

METHOD. The standard extractor fits the INITIAL slope of a volume-averaged Cauchy
stress, which by construction hides large-strain curvature. Instead we read the
driven reference point's history (U, RF) directly and form the unambiguous
homogenized NOMINAL stress-strain path
    sigma_nom(t) = RF_RP(t) / L^2 ,   eps_nom(t) = U_RP(t) / L
over the whole 0 -> +-4% ramp (see analyze_nlgeom.py). The geometric-nonlinearity
signature is the deviation of the secant modulus at max strain from the initial
tangent (= the linear modulus), read off this reaction-based curve.

DESIGN (paired, geometry-controlled). Three representative column slices
(z25 stiff/low-porosity, z65 transition, z95 warm channelled base) are written to
three CSVs in the SAME row order, so generating all three with the SAME SPAX_SEED
gives byte-identical morphology per slice. The files differ ONLY in nlgeom_flag
and Disp:
    lin : nlgeom OFF, Disp = +0.020  (4% tension, straight reference line E_lin)
    ten : nlgeom ON,  Disp = +0.020  (4% tension,  finite-strain path)
    cmp : nlgeom ON,  Disp = -0.020  (4% compression, finite-strain path)
Single-axis (utx only; Mode2 blanked) to keep it cheap and focused.
"""
from make_ice_studies import row, write, phi_brine, E_matrix, temperature, ZS

S_CSHAPE = [7.0, 5.5, 4.8, 4.5, 4.3, 4.3, 4.5, 5.0, 6.0, 8.0]
SLICES = (0.25, 0.65, 0.95)          # stiff top / transition / warm channelled base
DISP = 0.020                          # |applied strain| = Disp/L = 0.02/0.50 = 4%

def _column_row(run_id, z):
    """Build one FY C-shape column slice (same physics as study_brineK), utx only."""
    S = S_CSHAPE[ZS.index(z)]
    T = temperature(z, T_top=-20.0)
    if T > -0.5:
        T = -0.5
    phi_b = phi_brine(S, T)
    gas    = 0.012
    spher  = 0.86 - 0.26 * z
    growth = 0.40 + 0.32 * z
    r_avg  = 0.030 + 0.016 * z
    ch = 0.40 if phi_b > 0.05 else 0.0
    r = row(run_id, E_matrix(T), phi_b, gas, spher, growth,
            channels_frac=ch, r_avg=r_avg)
    r['Mode2'] = ''                   # single-axis: suppress the utz deck
    r['Disp2'] = ''
    return r

# ---------------------------------------------------------------------------
# Study 8: LARGE-DEFORMATION (nlgeom) homogenization on 3 representative slices.
def study_nlgeom():
    cases = (
        ('rve_nlgeom_lin.csv', 'OFF', +DISP),   # linear reference (straight line)
        ('rve_nlgeom_ten.csv', 'ON',  +DISP),   # finite-strain tension
        ('rve_nlgeom_cmp.csv', 'ON',  -DISP),   # finite-strain compression
    )
    for fname, nl, disp in cases:
        rows = []
        for z in SLICES:
            tag = {'OFF': 'LIN', 'ON': 'TEN' if disp > 0 else 'CMP'}[nl]
            r = _column_row(f'NLG{tag}_z{int(z*100):02d}', z)
            r['nlgeom_flag'] = nl
            r['Disp'] = f'{disp:+.3f}'
            rows.append(r)
        write(fname, rows)

if __name__ == '__main__':
    study_nlgeom()
