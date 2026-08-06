"""Decks for case study 3: the Novik Bay cantilever beams of Gogolaze et al. (2026).

Why these runs
--------------
The case study drives their layered-composite beam construction with an
RVE-computed E(z) instead of an empirical E(v_b) regression.  Scoping against
the existing database (analysis/case_study_3_scoping.py) left two gaps, and
only the second of them is what we first thought.

1. COVERAGE is adequate.  Indexed on the TRUE total soft fraction,
   phi_tot = VoF_sphere + channel_vof_target -- note VoF_sphere counts only the
   meshed pockets and excludes the channel network -- the database already
   reaches phi_tot = 0.52, well past the phi_tot = 0.24 the base of these beams
   needs.

2. MORPHOLOGY, not coverage, is the binding uncertainty.  At phi_tot = 0.24 the
   existing cells split into two populations that differ by half again in
   modulus at identical soft fraction: channel-dominated (SKEL_c220,
   channels_frac 0.85) gives E/E_matrix = 0.47, pocket-dominated (SKEL_p220,
   channels_frac 0.40) gives 0.30.  Above percolation the knockdown is therefore
   not a function of phi_tot alone, and pooling the two is what produced the
   +-0.089 scatter in the scoping curve.  Since the beam is percolated at every
   depth (v_b >= 64 per-mille throughout, ~12x the Golden threshold), this
   ambiguity propagates into every slice, not just the base.

Two decks follow.

rve_gogo_column.csv    The beam itself: twelve slices whose brine fraction is
                       taken from their measured profile (their eq. 14), five
                       packings each, at the baseline channel partition of the
                       production column so that the result is commensurable
                       with Section sec:macro.  This is what E(z) is read from.

rve_gogo_chanfrac.csv  The sensitivity that decides the error bar: the brine
                       partition swept from pocket-dominated to skeletal at the
                       base slice, five packings each.  Cheap, and it converts
                       the morphology ambiguity from an unknown into a quoted
                       band.

Mesh resolution
---------------
These decks hold ELEMENTS-PER-INCLUSION constant rather than L_mesh, which is
the criterion docs/RVE_STUDY_README.md states ("mesh resolution fixed ... so
elements-per-inclusion is constant") but which only bites once r_avg varies with
depth, as it does here.  The target is 3.21 elements across an inclusion
diameter, i.e. the ratio of the skeletal campaign (L_mesh = 0.028 at
r_avg = 0.045), so

    L_mesh(z) = 2 * r_avg(z) / 3.21 .

Pocket radius follows the production column's growth law, r_avg = 0.030 + 0.016 z,
which happens to reach 0.0453 at the base -- so at the base these cells recover
the skeletal mesh almost exactly, and higher up they carry production-column
pocket sizes at the same relative resolution.  Cost is ~1.8x a flat
L_mesh = 0.028 deck, almost all of it in the cold slices where the pockets are
small and the mesh therefore fine.

Channel lamellae still widen with phi_b, as make_skeletal.py established: the
sliver regime that defeats Gmsh is many near-touching thin cylinders, not a few
fat ones.  Narrow matrix ligaments are handled by the mesher's own local
refinement balls (SpaX_Standalone._collect_gap_balls), not by the global size.

Usage:  python studies/make_cantilever.py
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_ice_studies import COLS, E_matrix, row

PARAMS = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, 'params')

# --- the measured beam ------------------------------------------------------
H_BEAM = 0.32                 # m, Beams 3 and 4
VB_POLY = (0.29315, -5.124, 85.977)   # their eq. (14): z in cm -> v_b in per-mille
T_TOP, T_BOT = -7.5, -1.8     # degC, their Fig. 6 (Beam 3)
GAS = 0.020                   # gas void fraction, as carried by the basal cells

N_SLICES = 12
SEEDS = (1, 2, 3, 4, 5)

# Baseline brine partition. 0.40 is the production column's value above the
# percolation threshold; the sweep below brackets it.
CH_BASE = 0.40
CH_SWEEP = (0.20, 0.40, 0.60, 0.85)

PHI_REF = 0.15                # reference brine fraction for the morphology ramp
PHI_MAX = 0.50                # meshing ceiling established by make_skeletal.py

# Mesh criterion: elements across an inclusion diameter, held fixed at every
# depth.  3.21 is the skeletal campaign's ratio, L_mesh = 0.028 at r_avg = 0.045.
EL_PER_DIAMETER = 2 * 0.045 / 0.028

# Production-column pocket growth: r_avg = R_AVG_0 + R_AVG_GROWTH * z/H.
R_AVG_0, R_AVG_GROWTH = 0.030, 0.016


def pocket_radius(z_over_h):
    # rounded here, because row() writes r_avg at 3 dp and L_mesh must be
    # derived from the value that actually reaches the generator
    return round(R_AVG_0 + R_AVG_GROWTH * z_over_h, 3)


def mesh_size(r_avg):
    """L_mesh that puts EL_PER_DIAMETER elements across a pocket diameter."""
    return 2.0 * r_avg / EL_PER_DIAMETER


def brine_fraction(z_over_h):
    """Their eq. (14), evaluated at a depth fraction, returned as a fraction."""
    a, b, c = VB_POLY
    zc = z_over_h * H_BEAM * 100.0
    return (a * zc ** 2 + b * zc + c) / 1000.0


def temperature(z_over_h):
    return T_TOP + (T_BOT - T_TOP) * z_over_h


def beam_row(rid, phi_b, T, channels_frac, z_over_h, seed_tag=''):
    """One slice of the Novik Bay beam.

    Morphology follows the basal trend of make_skeletal.py: as brine fraction
    rises the residual ice becomes platelet-like, so sphericity falls and
    Z-elongation rises.  Pocket size follows the production column, and the
    element size follows the pocket, so elements-per-inclusion is constant.
    """
    t = max(0.0, min(1.0, (phi_b - PHI_REF) / (PHI_MAX - PHI_REF)))
    spher = 0.62 - 0.22 * t
    growth = 0.70 + 0.25 * t
    r_avg = pocket_radius(z_over_h)
    r = row(rid + seed_tag, E_matrix(T), phi_b, GAS, spher, growth,
            channels_frac=channels_frac, r_avg=r_avg)
    # Fatter, fewer lamellae as brine rises -- the sliver regime that crashes
    # Gmsh is many near-touching thin cylinders, not a few wide ones.
    r['r_channel_avg'] = '%.3f' % (0.022 + 0.023 * t)
    r['r_channel_std'] = '%.3f' % (0.005 + 0.004 * t)
    r['L_mesh'] = '%.4f' % mesh_size(r_avg)
    r['Mode'], r['Disp'] = 'Uniaxial Tension X', '0.005'
    r['Mode2'], r['Disp2'] = 'Uniaxial Tension Z', '0.005'
    return r


def write(name, rows):
    path = os.path.join(PARAMS, name)
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print('wrote params/%s  (%d RVEs, %d solves)' % (name, len(rows), 2 * len(rows)))


def study_column():
    """Twelve slices of the measured beam, five packings each."""
    rows = []
    for i in range(N_SLICES):
        z = (i + 0.5) / N_SLICES
        phi_b = round(brine_fraction(z), 4)
        T = temperature(z)
        for s in SEEDS:
            rows.append(beam_row('GOGO_z%03d' % round(z * 1000), phi_b, T,
                                 CH_BASE, z, '_s%d' % s))
    write('rve_gogo_column.csv', rows)
    return rows


def study_chanfrac():
    """Brine partition swept at the base slice: the dominant uncertainty."""
    z = (N_SLICES - 0.5) / N_SLICES
    phi_b = round(brine_fraction(z), 4)
    T = temperature(z)
    rows = []
    for ch in CH_SWEEP:
        for s in SEEDS:
            rows.append(beam_row('GOCH_f%03d' % round(ch * 100), phi_b, T,
                                 ch, z, '_s%d' % s))
    write('rve_gogo_chanfrac.csv', rows)
    return rows


def summary(col_rows):
    from make_ice_studies import L
    print('\nbeam profile (12 slices, h = %.2f m), %.2f elements per inclusion '
          'diameter at every depth:' % (H_BEAM, EL_PER_DIAMETER))
    print('  %-12s %6s %8s %8s %8s %8s %8s %9s'
          % ('slice', 'z/H', 'T degC', 'v_b', 'phi_tot', 'r_avg', 'L_mesh', '~cells'))
    seen = set()
    cells = 0.0
    for r in col_rows:
        rid = r['run_id'].rsplit('_s', 1)[0]
        if rid in seen:
            continue
        seen.add(rid)
        z = int(rid.split('_z')[1]) / 1000.0
        phi_tot = float(r['VoF_sphere']) + float(r['channel_vof_target'])
        lm = float(r['L_mesh'])
        n = (float(L) / lm) ** 3
        cells += n * len(SEEDS)
        print('  %-12s %6.3f %8.2f %8.4f %8.4f %8.4f %8.4f %9.0f'
              % (rid, z, temperature(z), brine_fraction(z), phi_tot,
                 float(r['r_avg']), lm, n))
    print('  (~cells is (L/L_mesh)^3, a shape-independent proxy, not an element count)')
    return cells


if __name__ == '__main__':
    col = study_column()
    ch = study_chanfrac()
    summary(col)
    print('\ntotal: %d RVEs, %d solves' % (len(col) + len(ch),
                                           2 * (len(col) + len(ch))))
