#!/usr/bin/env python3
r"""Build the two layered-morphology decks the column campaigns are missing.

Every depth-column campaign in this study meshes isolated pockets and channels
at every depth, including the slices below the percolation threshold where
Section 4.4 shows that description fails. The slab primitive exists but has only
ever been exercised in the bracket studies, which are single-morphology cells at
one fixed composition. These two decks close that gap.

LAYERCOL -- a layered cell at each brine fraction a percolated column slice
    actually reaches, phi = 0.06 to 0.15. The bridge fraction is set from
    Assur's geometry, b = 1 - sqrt(phi/phi_0), so the cells test whether a
    column built that way reproduces the closure rather than assuming it. Both
    load directions are solved, which gives the ANISOTROPY of the layered
    morphology -- a quantity we do not currently have for any column, since the
    E_z/E_x = 1.13 reported at the base is a pocket-and-channel value and the
    confinement study suggests the layered figure is far larger.

LAYERB -- b swept independently of phi at two fixed brine fractions. The closure
    uses a single exponent n at every phi; if n drifts with phi that assumption
    is wrong, and bracket_bridge cannot say, having swept b at one phi only.

Both decks solve drained and undrained, because Section 4.4.3 measured that to
be worth up to 19x in this morphology and nothing in a column campaign has
carried the distinction before.

    python3 make_layer_decks.py <outdir>
"""
import csv
import os
import sys

PHI_0 = 0.20
PHI_C = 0.05

# Template taken from rve_bracket_bridge.csv so the two are directly
# comparable: same cell, same mesh density, same matrix, same pocket content.
# L_mesh is halved from the bracket decks' 0.024. Those cells put well under one
# element across a brine layer -- t = slab_vof*L/(n_slabs*(1-b)) runs 0.014 to
# 0.022 here against a 0.024 element -- which was tolerable for measuring a
# trend in b but is not good enough for the absolute moduli a column needs.
# 0.012 buys one to two elements across the layer at roughly eight times the
# element count. See the thickness table printed by this script.
BASE = {
    'L': '0.50', 'L_mesh': '0.0120', 'Is_Porous': 'Composite',
    'E_matrix': '9.37e9', 'nu_matrix': '0.33',
    'VoF_sphere': '0.0325', 'r_avg': '0.030', 'r_std': '0.008',
    'Mode': 'Uniaxial Tension X', 'Disp': '0.005',
    'Mode2': 'Uniaxial Tension Z', 'Disp2': '0.005',
    'VoF_void_sphere': '0.0100', 'VoF_incl_sphere': '0.0225',
    'E_sphere_inclusion': '2.2e9', 'nu_sphere_inclusion': '0.48',
    'sphericity_avg': '0.85', 'sphericity_std': '0.10',
    'min_distance': '0.002', 'max_iterations': '200000',
    'nlgeom_flag': 'OFF', 'PBC_Method': 'Gmsh', 'Kappa': '0',
    'Bending_Plane': 'xz', 'Bending_PBC_Type': 'Lesicar',
    'generate_channels': 'No', 'channel_vof_target': '0.0',
    'r_channel_avg': '0.025', 'r_channel_std': '0.006',
    'Growth_Direction': 'Z', 'Growth_Concentration': '0.40',
    'Inclusion_Type': 'Liquid',
    'G_inclusion': '440029.33528897085', 'full_tensor': 'No',
    'n_slabs': '4', 'n_bridges': '2', 'slab_axis': 'x',
    'bridge_correlation': '0.0000',
}
COLS = ['run_id', 'L', 'L_mesh', 'Is_Porous', 'E_matrix', 'nu_matrix',
        'VoF_sphere', 'r_avg', 'r_std', 'Mode', 'Disp', 'Mode2', 'Disp2',
        'VoF_void_sphere', 'VoF_incl_sphere', 'E_sphere_inclusion',
        'nu_sphere_inclusion', 'sphericity_avg', 'sphericity_std',
        'min_distance', 'max_iterations', 'nlgeom_flag', 'PBC_Method',
        'Kappa', 'Bending_Plane', 'Bending_PBC_Type', 'generate_channels',
        'channel_vof_target', 'r_channel_avg', 'r_channel_std',
        'Growth_Direction', 'Growth_Concentration', 'Inclusion_Type',
        'K_inclusion', 'G_inclusion', 'full_tensor', 'n_slabs', 'slab_vof',
        'bridge_fraction', 'n_bridges', 'slab_axis', 'bridge_correlation']

# drained releases the brine bulk modulus by three decades; undrained seals it.
K = {'drn': '2.2e+06', 'und': '2.2e9'}


def assur_b(phi, phi_0=PHI_0):
    return max(0.0, 1.0 - (min(phi, phi_0) / phi_0) ** 0.5)


def row(rid, phi_slab, b, state):
    r = dict(BASE)
    r['run_id'] = rid
    r['slab_vof'] = '%.4f' % phi_slab
    r['bridge_fraction'] = '%.4f' % b
    r['K_inclusion'] = K[state]
    return [r[c] for c in COLS]


def thickness(phi_slab, b, L=0.5, n_slabs=4):
    """Brine-layer thickness in model units, from phi_layer = t(1-b)/a_0."""
    return phi_slab * L / (n_slabs * max(1.0 - b, 1e-6))


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    os.makedirs(out, exist_ok=True)

    lm = float(BASE['L_mesh'])
    print('layer thickness vs element size (L_mesh = %.4f)' % lm)
    print('  %-7s %-7s %-9s %s' % ('phi', 'b', 't', 'elements across'))
    worst = 1e9
    for q in (0.06, 0.08, 0.10, 0.12, 0.15):
        t = thickness(q, assur_b(q))
        worst = min(worst, t / lm)
        print('  %-7.2f %-7.3f %-9.4f %.1f' % (q, assur_b(q), t, t / lm))
    if worst < 1.0:
        print('  WARNING: thinnest layer is under one element across')
    print()

    # ---- LAYERCOL: the brine fractions percolated column slices reach ------
    # phi is the slab target; the mesh realises somewhat more once the pockets
    # are counted, and every analysis reads phi_inclusion back rather than
    # trusting this number.
    p = os.path.join(out, 'rve_layercol.csv')
    n = 0
    with open(p, 'w', newline='', encoding='utf8') as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        for phi in (0.06, 0.08, 0.10, 0.12, 0.15):
            b = assur_b(phi)
            for state in ('drn', 'und'):
                for s in (1, 2, 3):
                    w.writerow(row('LCOL_p%03d_%s_s%d' % (round(phi * 1000),
                                                          state, s),
                                   phi, b, state))
                    n += 1
    print('wrote %s  (%d cells)' % (p, n))
    print('   phi -> b :', ', '.join('%.2f->%.3f' % (q, assur_b(q))
                                     for q in (0.06, 0.08, 0.10, 0.12, 0.15)))

    # ---- LAYERB: b independent of phi, to test a phi-independent exponent --
    p = os.path.join(out, 'rve_layerb.csv')
    n = 0
    with open(p, 'w', newline='', encoding='utf8') as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        for phi in (0.08, 0.15):
            for b in (0.10, 0.20, 0.30, 0.45):
                for state in ('drn', 'und'):
                    for s in (1, 2, 3):
                        w.writerow(row('LB_p%03d_b%03d_%s_s%d'
                                       % (round(phi * 1000), round(b * 100),
                                          state, s), phi, b, state))
                        n += 1
    print('wrote %s  (%d cells)' % (p, n))
    print('   Assur b at phi=0.08 is %.3f, at phi=0.15 is %.3f -- both inside'
          ' the swept range, so the sweep brackets the assumed value'
          % (assur_b(0.08), assur_b(0.15)))

    # ---- LAYERMESH: is a layer this thin actually resolved? ----------------
    # One condition at three element sizes. Only one to two elements span a
    # brine layer at the production setting, so the moduli above are worth
    # nothing until this shows them flat under refinement. The bracket decks
    # never ran this check, having put under one element across the layer.
    p = os.path.join(out, 'rve_layermesh.csv')
    n = 0
    phi, b = 0.10, assur_b(0.10)
    with open(p, 'w', newline='', encoding='utf8') as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        for lm in ('0.0240', '0.0120', '0.0080', '0.0060'):
            for state in ('drn', 'und'):
                for s in (1, 2):
                    r = dict(BASE)
                    r['L_mesh'] = lm
                    r['run_id'] = 'LMESH_m%s_%s_s%d' % (lm.replace('.', 'p'),
                                                        state, s)
                    r['slab_vof'] = '%.4f' % phi
                    r['bridge_fraction'] = '%.4f' % b
                    r['K_inclusion'] = K[state]
                    w.writerow([r[c] for c in COLS])
                    n += 1
    print('wrote %s  (%d cells)' % (p, n))
    print('   phi=%.2f, b=%.3f, t=%.4f -> %.1f/%.1f/%.1f/%.1f elements across'
          % (phi, b, thickness(phi, b),
             thickness(phi, b) / 0.024, thickness(phi, b) / 0.012,
             thickness(phi, b) / 0.008, thickness(phi, b) / 0.006))


if __name__ == '__main__':
    main()
