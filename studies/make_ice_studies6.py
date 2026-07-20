"""Sixth battery of sea-ice studies (2026-07-06).

  #5 brineK -> rve_brineKconst.csv + rve_brineKtemp.csv
     (temperature-dependent brine bulk modulus K(T) down the FY C-shape column)

Closes the material model: the trapped brine has so far been modelled with a
FIXED bulk modulus K = 2.2 GPa at every depth. Physically the brine gets colder
and more concentrated with depth toward the surface (equilibrium brine salinity
rises as T drops), and its bulk modulus grows with salt content (K = rho c^2).
This study layers a physically-varying K(T) onto the column and measures how much
E(z) and the anisotropy E_z/E_x move -- expected small, since an earlier sweep
found E only ~3-7% sensitive to brine K over a 100x range, and here the physical
K span is only ~1.3x.

DESIGN (single-mesh, geometry-controlled). Mesh generation is NON-deterministic
(inclusions are randomly placed), so two separate generation runs do NOT share
geometry -- a matched-seed "paired CSV" scheme would leave the K(T) signal buried
under run-to-run mesh scatter. Instead this file only DEFINES the column: it emits
rve_brineKconst.csv (fixed K=2.2 GPa, used to generate ONE mesh per slice, utx+utz)
and rve_brineKtemp.csv (identical rows but K_inclusion = K(T), used only to carry
the per-slice K values). build_brineK_decks.py then stamps the K(T) twin onto each
Kconst mesh by rewriting ONLY the brine *Elastic card, so the paired decks share
byte-identical geometry, mesh, PBC and matrix and differ solely in brine K. The
per-slice E difference is then the pure K(T) sensitivity. First-order
(utx + utz -> E_x, E_z), linear C3D4H-hybrid tets.
"""
from make_ice_studies import row, write, phi_brine, E_matrix, temperature, ZS

# FY C-shape bulk salinity (Cox & Weeks 1983), same profile as the column studies.
S_CSHAPE = [7.0, 5.5, 4.8, 4.5, 4.3, 4.3, 4.5, 5.0, 6.0, 8.0]

def brine_salinity(T):
    """Equilibrium brine salinity (ppt) from freezing-point depression
    (Cox & Weeks 1983 / Assur polynomial), valid -22.9 < T < -2 degC. Clamped at
    the warm base (T -> -2, seawater ~35 ppt) and cold surface (T -> -22.9)."""
    Tc = min(max(T, -22.9), -2.0)
    return 1.725 - 18.756 * Tc - 0.3964 * Tc * Tc

def K_brine(T):
    """Brine bulk modulus (Pa), rising with concentration as the brine cools.
    Linear in brine salinity: pure-water-like 2.2 GPa at seawater salinity
    (35 ppt), up to ~2.8 GPa for the ~220 ppt cold concentrated surface brine.
    Slope 3.2e6 Pa/ppt (K = rho c^2 grows with dissolved salt)."""
    return 2.2e9 + 3.2e6 * (brine_salinity(T) - 35.0)

# ---------------------------------------------------------------------------
# Study 5: TEMPERATURE-DEPENDENT BRINE MODULUS K(T) down the FY C-shape column.
def study_brineK():
    T_TOP = -20.0
    rows_const, rows_temp = [], []
    for z, S in zip(ZS, S_CSHAPE):
        T = temperature(z, T_top=T_TOP)
        if T > -0.5:                       # Frankenstein-Garner validity floor
            T = -0.5
        phi_b = phi_brine(S, T)
        gas    = 0.012
        spher  = 0.86 - 0.26 * z
        growth = 0.40 + 0.32 * z
        r_avg  = 0.030 + 0.016 * z
        E_mat  = E_matrix(T)
        ch = 0.40 if phi_b > 0.05 else 0.0
        # constant-K reference (BASE bakes K_inclusion = 2.2 GPa)
        rc = row(f'BKC_z{int(z*100):02d}', E_mat, phi_b, gas, spher, growth,
                 channels_frac=ch, r_avg=r_avg)
        # K(T) twin: identical geometry inputs, only the brine bulk modulus differs
        rt = dict(rc)
        rt['run_id'] = f'BKT_z{int(z*100):02d}'
        rt['K_inclusion'] = f'{K_brine(T):.6e}'
        rows_const.append(rc)
        rows_temp.append(rt)
    write('rve_brineKconst.csv', rows_const)
    write('rve_brineKtemp.csv', rows_temp)

if __name__ == '__main__':
    study_brineK()
    # quick echo of the K(T) profile applied
    print('\n z/H    T(C)   S_brine(ppt)   K_brine(GPa)')
    for z in ZS:
        T = temperature(z, T_top=-20.0)
        if T > -0.5:
            T = -0.5
        print(f' {z:4.2f}  {T:6.1f}   {brine_salinity(T):8.1f}      '
              f'{K_brine(T)/1e9:6.3f}')
