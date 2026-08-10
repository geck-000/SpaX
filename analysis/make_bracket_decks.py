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


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(root, 'params')

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
