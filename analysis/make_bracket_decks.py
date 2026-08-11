# -*- coding: utf-8 -*-
"""Production decks for the drained/undrained bracket.

The bracket is reported as a relation E(phi) for each morphology rather than
slice by slice down a column, for two reasons. The relation is what gets
compared against the empirical laws, which are themselves written in phi; and a
layer carrying only a few percent of the cell volume is a few percent of L
thick, so resolving it at the cold end would need elements an order of
magnitude finer than the column decks use. Layers are the WARM-ice morphology
in any case, so the layered sweep starts at phi = 0.08 and the pocket sweep
covers the whole range.

Each cell is written twice, at K = 2.2 GPa (brine sealed in, the undrained
limit the cells have always computed) and at K = 2.2 MPa (the fill carries
nothing, drained). Everything else is held identical between the pair, so the
difference is drainage and nothing else.

Three realisations per point, because every scatter quoted in these papers is a
population s.d. over realisations and a bracket whose ends are single cells
would not support one.

    python3 analysis/make_bracket_decks.py          # writes into params/
"""
import os
import sys

HDR = ('run_id,L,L_mesh,Is_Porous,E_matrix,nu_matrix,VoF_sphere,r_avg,r_std,'
       'Mode,Disp,Mode2,Disp2,VoF_void_sphere,VoF_incl_sphere,'
       'E_sphere_inclusion,nu_sphere_inclusion,sphericity_avg,sphericity_std,'
       'min_distance,max_iterations,nlgeom_flag,PBC_Method,Kappa,Bending_Plane,'
       'Bending_PBC_Type,generate_channels,channel_vof_target,r_channel_avg,'
       'r_channel_std,Growth_Direction,Growth_Concentration,Inclusion_Type,'
       'K_inclusion,G_inclusion,full_tensor,'
       'n_slabs,slab_vof,bridge_fraction,n_bridges,slab_axis')

K_UND, K_DRN = 2.2e9, 2.2e6
G_BRINE = 440029.33528897085
VOID = 0.010
E_ICE = '9.37e9'
LAYER_SHARE = 0.85          # of the brine that sits in layers, warm ice
BRIDGE = 0.030
PHI_LAYER = (0.08, 0.10, 0.13, 0.16, 0.20, 0.23)
PHI_POCKET = (0.02, 0.04, 0.06, 0.08, 0.10, 0.13, 0.16, 0.20, 0.23)
SEEDS = (1, 2, 3)


def row(run_id, phi, K, layered, L_mesh):
    if layered:
        pockets, slab_vof, n_slabs = phi * (1 - LAYER_SHARE), \
            phi * LAYER_SHARE, 1
    else:
        pockets, slab_vof, n_slabs = phi, 0.0, 0
    return ','.join([
        run_id, '0.50', '%.4f' % L_mesh, 'Composite', E_ICE, '0.33',
        '%.4f' % (pockets + VOID), '0.030', '0.008',
        'Uniaxial Tension X', '0.005', 'Uniaxial Tension Z', '0.005',
        '%.4f' % VOID, '%.4f' % pockets, '2.2e9', '0.48', '0.85', '0.10',
        '0.002', '200000', 'OFF', 'Gmsh', '0', 'xz', 'Lesicar',
        'No', '0.0', '0.025', '0.006', 'Z', '0.40', 'Liquid',
        '%.6g' % K, '%.11f' % G_BRINE, 'No',
        str(n_slabs), '%.4f' % slab_vof, '%.4f' % BRIDGE, '2', 'x'])


def build(layered):
    rows = []
    phis = PHI_LAYER if layered else PHI_POCKET
    pre = 'BRKL' if layered else 'BRKP'
    # 0.024 resolves a layer down to phi = 0.08; pockets need no such margin
    # and stay on the column decks' own size.
    L_mesh = 0.024 if layered else 0.033
    for phi in phis:
        for tag, K in (('und', K_UND), ('drn', K_DRN)):
            for s in SEEDS:
                rows.append(row('%s_p%03d_%s_s%d'
                                % (pre, int(round(1000 * phi)), tag, s),
                                phi, K, layered, L_mesh))
    return rows


def build_nlayers():
    """Layer count at FIXED porosity -- separating two variables that were
    confounded in the exploratory runs.

    Series theory says the number of layers cannot matter at fixed total layer
    thickness: the compliances add and only the sum enters. Three dimensions
    need not agree, because a thinner layer is more strongly confined by the
    ice on either side, and confinement is exactly what makes an undrained
    layer stiff. The exploratory cells hint that it does matter -- two layers at
    phi = 0.179 came out at 5.23 GPa against one layer at phi = 0.150 at 2.87,
    stiffer on more brine, which cannot be a porosity effect.

    If the count matters undrained and not drained, that is confinement and it
    belongs in the model as the lamellar spacing, a real measurable length.
    If it matters at both ends it is a discretisation artefact and the bracket
    numbers need re-basing. Either way it has to be resolved before the layered
    moduli can be quoted, so phi, b, L and the mesh are all held fixed and only
    the count moves.
    """
    rows = []
    phi = 0.15
    for n in (1, 2, 3, 4):
        for tag, K in (('und', K_UND), ('drn', K_DRN)):
            for s in SEEDS:
                r = row('BRKN_n%d_%s_s%d' % (n, tag, s), phi, K, True, 0.024)
                p = r.split(',')
                p[36] = str(n)                                  # n_slabs
                p[37] = '%.4f' % (phi * LAYER_SHARE)             # unchanged total
                rows.append(','.join(p))
    return rows


def build_spacing():
    """Fixed plate spacing, varying cell size: the homogenisation check.

    The layer pitch is L/n_slabs, so the earlier count sweep at fixed L varied
    the plate SPACING from 0.50 to 0.125 model units. The modulus moved because
    the material changed, not because the cell failed to converge: each layer
    contributes a constriction compliance where load funnels through its
    bridges, independent of that layer's thickness, so n constrictions per cell
    means compliance grows with n. Series theory misses it by summing only
    thickness terms.

    Spacing is therefore a physical parameter of the model, and the cell edge
    maps to a few millimetres, which puts the tested range squarely on the
    0.5-1 mm plate spacing of real sea ice. What has NOT been shown is that the
    modulus depends on the spacing alone rather than on the cell size as well.
    Holding a0 fixed and moving L over a factor 2.5, with n = L/a0 following,
    tests exactly that. A flat result means the cell is a genuine RVE for this
    morphology and a0 is its length parameter; a drifting one means it is not
    homogenising and no spacing can be quoted.
    """
    rows = []
    phi, a0 = 0.15, 0.125
    for L in (0.250, 0.375, 0.500, 0.625):
        n = int(round(L / a0))
        for tag, K in (('und', K_UND), ('drn', K_DRN)):
            for s in SEEDS:
                r = row('BRKS_L%03d_%s_s%d' % (int(round(1000 * L)), tag, s),
                        phi, K, True, 0.024).split(',')
                r[1] = '%.3f' % L
                r[36] = str(n)
                rows.append(','.join(r))
    return rows


def build_bridge():
    """Drained bridge-fraction sweep: the exponent that decides the mechanism.

    The earlier b-sweep was run undrained and showed no trend, because sealed
    brine resists at its bulk modulus and the bridges never control anything.
    In the drained limit they control everything, and E(b) has never been
    measured there.

    The exponent is the point. Weeks & Assur write E = 9.5 (1 - sqrt(v))^4,
    which is a load-bearing area raised to the FOURTH power: at phi = 0.227
    that area is 0.52, a large fraction, and the softness comes from the
    exponent rather than from a near-severed plane. Our series-and-constriction
    picture reaches the same modulus from b ~ 0.03 with a roughly linear
    dependence. Those are different mechanisms agreeing by arithmetic, and
    fitting E ~ b^n separates them: n near 4 means the cells reproduce the
    empirical law and b can then be taken from Assur geometry with nothing left
    free; n near 1 means the agreement was coincidence and the layered match
    should not be claimed as a prediction.

    Spacing is held at the physical value (four layers in a 0.5 cell) and phi
    at the column's mid-range, so b is the only thing moving.
    """
    rows = []
    phi = 0.15
    for b in (0.02, 0.04, 0.08, 0.15, 0.28, 0.50):
        for tag, K in (('drn', K_DRN), ('und', K_UND)):
            for s in SEEDS:
                r = row('BRKB_b%03d_%s_s%d' % (int(round(1000 * b)), tag, s),
                        phi, K, True, 0.024).split(',')
                r[36] = '4'                       # n_slabs, physical spacing
                r[38] = '%.4f' % b                # bridge_fraction
                rows.append(','.join(r))
    return rows


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(root, 'params')

    rows = build_bridge()
    p = os.path.join(out, 'rve_bracket_bridge.csv')
    with open(p, 'w') as f:
        f.write(HDR + '\n')
        for r in rows:
            f.write(r + '\n')
    print('wrote %s : %d cells (%d solves)' % (p, len(rows), 2 * len(rows)))

    rows = build_spacing()
    p = os.path.join(out, 'rve_bracket_spacing.csv')
    with open(p, 'w') as f:
        f.write(HDR + '\n')
        for r in rows:
            f.write(r + '\n')
    print('wrote %s : %d cells (%d solves)' % (p, len(rows), 2 * len(rows)))

    rows = build_nlayers()
    p = os.path.join(out, 'rve_bracket_nlayers.csv')
    with open(p, 'w') as f:
        f.write(HDR + '\n')
        for r in rows:
            f.write(r + '\n')
    print('wrote %s : %d cells (%d solves)' % (p, len(rows), 2 * len(rows)))

    for layered, name in ((True, 'rve_bracket_layer.csv'),
                          (False, 'rve_bracket_pocket.csv')):
        rows = build(layered)
        p = os.path.join(out, name)
        with open(p, 'w') as f:
            f.write(HDR + '\n')
            for r in rows:
                f.write(r + '\n')
        print('wrote %s : %d cells (%d solves)' % (p, len(rows), 2 * len(rows)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
