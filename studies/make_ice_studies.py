"""Generate the four sea-ice column study CSVs from physical parameters.

All four reuse the same physics as parametric_sea_ice_column.csv:
  brine volume   phi_b = S * (-49.185/T + 0.532) / 1000   (Frankenstein & Garner 1967)
  matrix modulus E(T) = 9.36e9 + 0.0129e9 * (-2 - T)      (ice stiffens ~13 MPa/degC)
  brine = K/G soft solid (K=2.2 GPa, G=0.44 MPa, nu~0.4999 near-fluid)
  gas    = true void (VoF_void_sphere)
  total meshed sphere VoF_sphere = VoF_void_sphere + VoF_incl_sphere
  total brine phi_b = VoF_incl_sphere + channel_vof_target (channels carry part of
                      brine once percolated, phi_b > 0.05 -- rule of fives, Golden 1998)

To keep the sweeps cheap each RVE solves only utx + utz (full_tensor=No):
that yields E_x (=E_eff) and E_z, hence the vertical/horizontal anisotropy E_z/E_x,
which is the headline of every column study, at 1/3 the cost of the full 6x6 tensor.

Studies:
  1 morphology  -> rve_morphology.csv    (fixed phi_b, sweep sphericity x growth-conc)
  2 percolation -> rve_percolation.csv   (sweep phi_b 0.03-0.10, channels off/on)
  3 monotonic   -> rve_monotonic.csv     (column but monotonic salinity, beam-eff matrix)
  4 seasonal    -> rve_seasonal.csv      (column at several T_top, winter->spring)
"""
import csv

# ---- shared constants -------------------------------------------------------
G_BRINE = 440029.33528897085   # Pa, soft-solid shear (nu~0.4999 with K=2.2 GPa)
K_BRINE = 2.2e9
NU_MAT  = 0.33
L, LMESH = 0.50, 0.033
BEAM_FACTOR = 0.49             # matrix -> vibrating-beam-effective modulus (Marchenko)

def phi_brine(S, T):
    """Frankenstein & Garner (1967) brine volume fraction; T in degC (<0), S in ppt."""
    return S * (-49.185 / T + 0.532) / 1000.0

def E_matrix(T, scale=1.0):
    return (9.36e9 + 0.0129e9 * (-2.0 - T)) * scale

BASE = {
    'L': f'{L:.2f}', 'L_mesh': f'{LMESH:.3f}', 'Is_Porous': 'Composite',
    'nu_matrix': f'{NU_MAT}', 'r_avg': '0.035', 'r_std': '0.010',
    'Mode': 'Uniaxial Tension X', 'Disp': '0.005',
    'Mode2': 'Uniaxial Tension Z', 'Disp2': '0.005',
    'E_sphere_inclusion': '2.2e9', 'nu_sphere_inclusion': '0.48',
    'sphericity_std': '0.10', 'min_distance': '0.002', 'max_iterations': '200000',
    'nlgeom_flag': 'OFF', 'PBC_Method': 'Gmsh', 'Kappa': '0',
    'Bending_Plane': 'xz', 'Bending_PBC_Type': 'Lesicar',
    'r_channel_avg': '0.022', 'r_channel_std': '0.005',
    'Growth_Direction': 'Z', 'Inclusion_Type': 'Liquid',
    'K_inclusion': f'{K_BRINE}', 'G_inclusion': f'{G_BRINE}',
    'full_tensor': 'No',
}
COLS = ['run_id','L','L_mesh','Is_Porous','E_matrix','nu_matrix','VoF_sphere','r_avg',
        'r_std','Mode','Disp','Mode2','Disp2','VoF_void_sphere','VoF_incl_sphere',
        'E_sphere_inclusion','nu_sphere_inclusion','sphericity_avg','sphericity_std',
        'min_distance','max_iterations','nlgeom_flag','PBC_Method','Kappa',
        'Bending_Plane','Bending_PBC_Type','generate_channels','channel_vof_target',
        'r_channel_avg','r_channel_std','Growth_Direction','Growth_Concentration',
        'Inclusion_Type','K_inclusion','G_inclusion','full_tensor']

def row(run_id, E_mat, phi_b, gas, spher, growth, channels_frac=0.0, r_avg=0.035):
    """Build a CSV row dict. channels_frac = fraction of brine carried by the
    percolating vertical channel network (rest goes to meshed pockets)."""
    r = dict(BASE)
    chan_vof = round(channels_frac * phi_b, 4)
    incl     = round(phi_b - chan_vof, 4)
    r.update({
        'run_id': run_id,
        'E_matrix': f'{E_mat:.4e}',
        'VoF_sphere': f'{round(gas + incl, 4)}',
        'VoF_void_sphere': f'{round(gas, 4)}',
        'VoF_incl_sphere': f'{incl}',
        'sphericity_avg': f'{spher:.2f}',
        'Growth_Concentration': f'{growth:.2f}',
        'generate_channels': 'Yes' if channels_frac > 0 else 'No',
        'channel_vof_target': f'{chan_vof}' if channels_frac > 0 else '0',
        'r_avg': f'{r_avg:.3f}',
    })
    return r

def write(path, rows):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print(f'wrote {path}  ({len(rows)} RVEs, {2*len(rows)} solves)')

# ===========================================================================
# depth grid shared by column-type studies (slice mid-points)
ZS = [0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95]

def temperature(z, T_top, T_bot=-1.8):
    return T_top + (T_bot - T_top) * z

# ---------------------------------------------------------------------------
# Study 1: MORPHOLOGY-FROM-CT  -- fixed phi_b, sweep sphericity x growth-conc
# Isolate how pocket SHAPE + ORIENTATION (the CT-measured quantities) drive
# vertical anisotropy E_z/E_x, with brine fraction and matrix held constant.
def study_morphology():
    phi_b = 0.05          # at the percolation threshold, isolated pockets
    gas   = 0.012
    E_mat = E_matrix(-10.0)
    rows = []
    for spher in (0.60, 0.70, 0.80, 0.90):
        for growth in (0.30, 0.60, 0.90):
            rid = f'MORF_s{int(spher*100)}g{int(growth*100)}'
            rows.append(row(rid, E_mat, phi_b, gas, spher, growth))
    write('rve_morphology.csv', rows)

# ---------------------------------------------------------------------------
# Study 2: PERCOLATION TRANSITION -- sweep phi_b 0.03..0.10, channels off & on.
# Locate the mechanical-anisotropy onset at the rule-of-fives (phi_b ~ 0.05).
# 'off' = brine all in meshed pockets; 'on' = 35% of brine in vertical channels.
def study_percolation():
    gas   = 0.012
    E_mat = E_matrix(-8.0)
    spher, growth = 0.70, 0.60
    levels = [0.03,0.04,0.05,0.06,0.07,0.08,0.10]
    rows = []
    for p in levels:                                   # pockets only (no channels)
        rows.append(row(f'PERC_p{int(p*1000):03d}off', E_mat, p, gas, spher, growth))
    for p in [x for x in levels if x >= 0.05]:         # percolated: channels on
        rows.append(row(f'PERC_p{int(p*1000):03d}on', E_mat, p, gas, spher, growth,
                        channels_frac=0.35))
    write('rve_percolation.csv', rows)

# ---------------------------------------------------------------------------
# Study 3: MONOTONIC SALINITY  -- column shape match to Marchenko.
# Salinity increases monotonically with depth (no C-shape interior minimum) so
# phi_b(z) and hence E(z) decrease monotonically top->base like Marchenko's fit.
# Matrix uses the vibrating-beam-effective modulus (x0.49) for absolute match.
def study_monotonic():
    rows = []
    for z in ZS:
        T = temperature(z, T_top=-20.0)
        S = 4.0 + 4.0 * z                 # monotonic 4.2 -> 7.8 ppt (skeletal base)
        phi_b = phi_brine(S, T)
        gas = 0.012
        spher  = 0.86 - 0.26 * z          # 0.85 -> 0.62
        growth = 0.40 + 0.32 * z          # 0.40 -> 0.70
        r_avg  = 0.030 + 0.016 * z        # 0.030 -> 0.045
        E_mat = E_matrix(T, scale=BEAM_FACTOR)
        ch = 0.40 if phi_b > 0.05 else 0.0
        rid = f'MONO_z{int(z*100):02d}'
        rows.append(row(rid, E_mat, phi_b, gas, spher, growth,
                        channels_frac=ch, r_avg=r_avg))
    write('rve_monotonic.csv', rows)

# ---------------------------------------------------------------------------
# Study 4: SEASONAL WARMING  -- the full column at several surface temperatures.
# As T_top rises winter->spring, phi_b crosses the percolation threshold higher
# up the column, so the warm soft (channelled) base layer thickens.  Uses the
# real FY C-shape salinity profile.  3 scenarios x 10 slices.
def study_seasonal():
    S_cshape = [7.0,5.5,4.8,4.5,4.3,4.3,4.5,5.0,6.0,8.0]   # Cox & Weeks 1983
    rows = []
    for tag, T_top in (('w20', -20.0), ('w12', -12.0), ('w06', -6.0)):
        for z, S in zip(ZS, S_cshape):
            T = temperature(z, T_top=T_top)
            if T > -0.5:        # Frankenstein-Garner validity floor
                T = -0.5
            phi_b = phi_brine(S, T)
            gas = 0.012 + 0.008 * (1 - z)       # slightly gassier near surface
            spher  = 0.86 - 0.26 * z
            growth = 0.40 + 0.32 * z
            r_avg  = 0.030 + 0.016 * z
            E_mat = E_matrix(T)
            ch = 0.40 if phi_b > 0.05 else 0.0
            rid = f'SEAS_{tag}_z{int(z*100):02d}'
            rows.append(row(rid, E_mat, phi_b, gas, spher, growth,
                            channels_frac=ch, r_avg=r_avg))
    write('rve_seasonal.csv', rows)

if __name__ == '__main__':
    study_morphology()
    study_percolation()
    study_monotonic()
    study_seasonal()
