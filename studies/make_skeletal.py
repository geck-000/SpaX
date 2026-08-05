"""Skeletal basal layer: resolve the bottom few percent of the sheet.

Motivation
----------
The ten-slice column averages the lowest 5% of the sheet into one RVE at
z/H=0.95 with phi_b=0.15, giving E_base=4.47 GPa against E_top=8.74 GPa, i.e.
a bottom-to-top modulus ratio alpha = 0.51. Four-point bending tests on
floating Baltic beams (Kujala et al. 1990, IAHR, Table 2) instead report
alpha = 0.12-0.19 and a neutral axis at z0/H = 0.37-0.39, against 0.45 for
alpha=0.51. For a linearly graded plate,

    z0/H = (1+2a)/(3(1+a))          a = E_bot/E_top

which reproduces their measured z0/H from their measured alpha, and inverting
it for the measured neutral axis gives a = 0.163. Their flexural-to-tensile
correction, E_top/E_flex = 3(1+a)/(a^2+4a+1), is 2.1 at a=0.163 but only 1.32
at our a=0.51 -- which is why our column cannot presently explain the 2.3x
offset against vibrating-beam moduli as a gradient artefact.

The missing physics is the skeletal layer: the lowest few centimetres are an
unconsolidated mush of ice platelets separated by brine lamellae, whose brine
volume approaches unity at the ice-water interface and which carries almost no
load. Frankenstein-Garner evaluated at a slice-mean temperature cannot produce
it, and a single 5%-thick slice cannot resolve it.

Two decks
---------
rve_skeletal.csv          E(phi_b) into the high-porosity regime, so that the
                          basal modulus is measured rather than assumed. Two
                          morphologies at matched phi_b -- channel-dominated
                          (brine in vertical lamellae, the physical skeletal
                          form) and pocket-dominated (control) -- to test
                          whether morphology matters once phi_b is large.

rve_skeletal_laminae.csv  the four sub-laminae of the bottom 5% that feed the
                          CLT re-assembly (analysis/skeletal_clt.py). Solves
                          utx + ss12, which is what classical lamination theory
                          needs (in-plane E and G); E_z is not used by CLT and
                          in-plane isotropy E_y=E_x holds to 0.2% above the
                          base and 1% at it.

Meshing ceiling
---------------
Random sequential adsorption of near-spherical pockets jams well below the
phi_b -> 1 limit of the true interface, so the sweep is capped at phi_b=0.50
and the laminae at 0.465. The residual, genuinely unconsolidated skin is not
meshable as a periodic cell and is handled in the CLT assembly as a
zero-stiffness lamina rather than pretended to be resolved here.
"""
import csv, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_ice_studies import phi_brine, E_matrix, BASE, COLS, row

PARAMS = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, 'params')

# --- basal conditions, matching the ICE_z95 slice of the production column ----
T_BASE   = -1.9          # degC at the base of the sheet
GAS_BASE = 0.020         # gas void fraction carried down from ICE_z95
PHI_95   = 0.15          # brine fraction of the existing lowest slice
PHI_MAX  = 0.50          # meshing ceiling (see module docstring)

# channel-dominated = brine in vertical lamellae between platelets (skeletal);
# pocket-dominated  = the ICE_z95 partition, kept as a morphology control.
CH_SKELETAL = 0.85
CH_POCKET   = 0.40


def write(name, rows):
    path = os.path.join(PARAMS, name)
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print('wrote params/%s  (%d RVEs)' % (name, len(rows)))


def basal_row(rid, phi_b, channels_frac, seed_tag=''):
    """A basal RVE at fixed T, gas and cell size, varying only phi_b and how the
    brine is partitioned between lamellae and pockets."""
    # sphericity falls and Z-elongation rises with phi_b: at high brine fraction
    # the residual ice is platelet-like rather than a matrix with round holes.
    t      = max(0.0, min(1.0, (phi_b - PHI_95) / (PHI_MAX - PHI_95)))
    spher  = 0.62 - 0.22 * t
    growth = 0.70 + 0.25 * t
    r = row(rid + seed_tag, E_matrix(T_BASE), phi_b, GAS_BASE, spher, growth,
            channels_frac=channels_frac, r_avg=0.045)
    # Widen the lamellae with phi_b. A vertical cylinder of radius rc in an
    # L-cell carries pi*rc^2/L^2 of volume each, so holding rc fixed at 0.022
    # would need ~70 channels to reach channel_vof=0.43 in an L=0.50 cell --
    # 70 near-touching cylinders is the sliver regime that crashes Gmsh. Fatter,
    # fewer lamellae reach the same volume with meshable ligaments, and is also
    # the physical trend: skeletal brine layers coarsen toward the interface.
    r['r_channel_avg'] = '%.3f' % (0.022 + 0.023 * t)
    r['r_channel_std'] = '%.3f' % (0.005 + 0.004 * t)
    # the mesh must resolve the thinning ligaments between lamellae
    r['L_mesh'] = '0.028'
    r['Mode'], r['Disp']   = 'Uniaxial Tension X', '0.005'
    r['Mode2'], r['Disp2'] = 'Uniaxial Tension Z', '0.005'
    return r


def study_phi_sweep():
    """E(phi_b) from the existing base up to the meshing ceiling, 3 packings each."""
    rows = []
    for phi in (0.15, 0.22, 0.29, 0.36, 0.43, 0.50):
        for s in (1, 2, 3):
            rows.append(basal_row('SKEL_c%03d' % round(phi * 1000),
                                  phi, CH_SKELETAL, '_s%d' % s))
    for phi in (0.22, 0.36, 0.50):                      # morphology control
        for s in (1, 2, 3):
            rows.append(basal_row('SKEL_p%03d' % round(phi * 1000),
                                  phi, CH_POCKET, '_s%d' % s))
    write('rve_skeletal.csv', rows)


def study_laminae():
    """The four sub-laminae of the bottom 5%, for the CLT re-assembly.

    Solid fraction is taken to rise linearly from the interface across the
    skeletal thickness, the standard first approximation for a mush, so
    phi_b(z) = PHI_95 + (PHI_MAX - PHI_95) * (z - 0.95)/0.05 truncated at the
    meshing ceiling. Each lamina solves utx + ss12 for the in-plane pair CLT
    needs.
    """
    rows = []
    for z in (0.955, 0.970, 0.985, 0.995):
        t     = (z - 0.95) / 0.05
        phi_b = round(PHI_95 + (PHI_MAX - PHI_95) * t, 3)
        for s in (1, 2, 3):
            r = basal_row('SKLM_z%03d' % round(z * 1000), phi_b,
                          CH_SKELETAL, '_s%d' % s)
            r['Mode'],  r['Disp']  = 'Uniaxial Tension X', '0.005'
            r['Mode2'], r['Disp2'] = 'Simple Shear S12',   '0.005'
            rows.append(r)
    write('rve_skeletal_laminae.csv', rows)


def study_steep_column():
    """A steeply monotonic column carried down to a resolved skeletal base.

    analysis/skeletal_clt.py shows that inserting a skeletal lamina under the
    existing column moves the neutral plane only from z0/H=0.466 to 0.452,
    nowhere near the measured 0.37-0.39: a thin compliant layer at the plate
    surface is too thin to shift the centroid of E(z). What does control the
    neutral plane is the SHAPE of the profile. A linear E(z) reaches the
    measured band at alpha~0.16, while our C-shape-salinity profile -- flat
    through the cold interior, softening only in the bottom fifth -- stays at
    z0/H=0.438 even with its endpoints forced to the same alpha.

    So the discriminating experiment is not a softer base but a steeper
    interior. This deck imposes a steeply monotonic salinity, which removes the
    interior plateau, and appends the four skeletal laminae so the base is
    resolved as well. If the neutral plane then falls into the measured band,
    the C-shape salinity is the source of the disagreement; if it does not, the
    disagreement is with Kujala's linear-E(z) inversion rather than with our
    microstructure.
    """
    rows = []
    ZS = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95]
    T_TOP = -15.0
    for z in ZS:
        S = 3.0 + 9.0 * z                       # steeply monotonic, no C-shape
        T = T_TOP + (T_BASE - T_TOP) * z
        phi_b = phi_brine(S, T)
        gas = 0.012 + 0.008 * (1 - z)
        spher = 0.86 - 0.26 * z
        growth = 0.40 + 0.32 * z
        ch = 0.40 if phi_b > 0.05 else 0.0
        for s in (1, 2, 3):
            r = row('STEEP_z%02d' % round(z * 100), E_matrix(T), phi_b, gas,
                    spher, growth, channels_frac=ch, r_avg=0.030 + 0.016 * z)
            r['run_id'] += '_s%d' % s
            r['Mode'], r['Disp'] = 'Uniaxial Tension X', '0.005'
            r['Mode2'], r['Disp2'] = 'Simple Shear S12', '0.005'
            rows.append(r)
    # The same four skeletal laminae, so the base is resolved too. The ramp has
    # to start from THIS column's lowest slice, not from the C-shape column's
    # PHI_95, or the first lamina would be less porous than the slice above it.
    phi_start = phi_brine(3.0 + 9.0 * 0.95, T_TOP + (T_BASE - T_TOP) * 0.95)
    for z in (0.955, 0.970, 0.985, 0.995):
        t = (z - 0.95) / 0.05
        phi_b = round(phi_start + (PHI_MAX - phi_start) * t, 3)
        for s in (1, 2, 3):
            r = basal_row('STEEP_z%03d' % round(z * 1000), phi_b,
                          CH_SKELETAL, '_s%d' % s)
            r['Mode2'], r['Disp2'] = 'Simple Shear S12', '0.005'
            rows.append(r)
    write('rve_steep_column.csv', rows)


if __name__ == '__main__':
    study_phi_sweep()
    study_laminae()
    study_steep_column()
    a = 0.163
    print('\ntarget: alpha=%.3f -> z0/H=%.3f, E_top/E_flex=%.2f'
          % (a, (1 + 2 * a) / (3 * (1 + a)), 3 * (1 + a) / (a * a + 4 * a + 1)))
    print('need E_base ~ %.2f GPa against E_top = 8.74 GPa' % (a * 8.739))
