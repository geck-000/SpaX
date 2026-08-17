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


ELEM_ACROSS = 2.5      # elements across a brine layer, set by the mesh gate
LM_MIN, LM_MAX = 0.005, 0.012


def mesh_for(phi_slab, b):
    """Element size giving ELEM_ACROSS elements across this cell's layer.

    The gate (rve_layermesh) measured convergence against the finest mesh at
    phi = 0.10, b = 0.293, where the layer spans t = 0.0177:

        elements across   drained err   undrained err
              3.0             0.0%          0.0%
              2.2             0.3%          0.3%
              1.5             1.0%         20.7%
              0.7             8.7%         35.0%

    So the undrained response needs roughly twice the resolution the drained
    one does -- near-incompressible brine in a thin layer -- and a single fixed
    element size cannot serve a deck whose layer thickness varies by a factor
    of two. Sizing per cell puts every one at the same resolution instead.
    """
    t = thickness(phi_slab, b)
    return min(max(t / ELEM_ACROSS, LM_MIN), LM_MAX)


def row(rid, phi_slab, b, state):
    r = dict(BASE)
    r['run_id'] = rid
    r['slab_vof'] = '%.4f' % phi_slab
    r['bridge_fraction'] = '%.4f' % b
    r['K_inclusion'] = K[state]
    r['L_mesh'] = '%.4f' % mesh_for(phi_slab, b)
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

    # ---- LAYERSKEL: the skeletal layer, with the SAME primitive ------------
    # The closure stops at phi_0 because Assur's b(phi) reaches zero there, and
    # the paper has been saying the skeletal layer "needs its own description".
    # That conflates two things. The slab primitive is a set of brine layers
    # between ice platelets held apart by bridges, which is exactly what the
    # skeletal layer IS; nothing about it breaks above phi_0. What breaks is
    # the RELATION tying b to phi, since b is an independent parameter of the
    # generator and only the closure ties it to Assur.
    #
    # So these cells decouple them: high brine fraction with bridges imposed
    # directly. If they carry load, the floor above phi_0 is premature and what
    # is missing is a measured b(phi) there rather than a model. Both load
    # directions, because the interesting quantity is the anisotropy: along the
    # platelets the ice is continuous and across them it is not, so E_z/E_x
    # should diverge where the pocket cells have it saturating near 1.7.
    p = os.path.join(out, 'rve_layerskel.csv')
    n = 0
    with open(p, 'w', newline='', encoding='utf8') as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        for phi in (0.22, 0.28, 0.35):
            for b in (0.03, 0.07, 0.13):
                for state in ('drn', 'und'):
                    for s in (1, 2):
                        w.writerow(row('LSK_p%03d_b%03d_%s_s%d'
                                       % (round(phi * 1000), round(b * 100),
                                          state, s), phi, b, state))
                        n += 1
    print('wrote %s  (%d cells)' % (p, n))
    print('   all above phi_0 = %.2f, where Assur gives b = 0' % PHI_0)
    for q in (0.22, 0.28, 0.35):
        t = thickness(q, 0.03)
        print('     phi %.2f: layer t = %.4f (%.1f elements), thinnest bridge '
              'r = %.4f (%.1f elements)'
              % (q, t, t / float(BASE['L_mesh']), (0.03 * 0.25 / 2 / 3.14159)
                 ** 0.5, ((0.03 * 0.25 / 2 / 3.14159) ** 0.5)
                 / float(BASE['L_mesh'])))

    # ---- WBLLAYER: the localisation study, in the right morphology ---------
    # rve_weibull.csv compares six morphologies for stress concentration, and
    # its two warmest cases (WBL_CHAN at phi ~ 0.09, WBL_BASE at ~ 0.21) sit at
    # or above the in-plane threshold, so both are meshed as pockets and
    # channels where the layered description belongs. These are the layered
    # counterparts at the same brine fractions, so the SCF comparison can be
    # made morphology against morphology rather than assumed.
    #
    # The expectation worth testing: in a layered cell the entire transverse
    # load path is the ice bridges, so the stress should concentrate there and
    # the tail of the distribution should be far longer than any pocket case.
    # Section 4.6.1 finds the ranking of morphologies to be m-independent, and
    # this is the case most likely to disturb that.
    p = os.path.join(out, 'rve_weibull_layer.csv')
    n = 0
    with open(p, 'w', newline='', encoding='utf8') as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        for tag, phi in (('chan', 0.09), ('base', 0.21)):
            b = assur_b(phi) if phi < PHI_0 else 0.05
            for state in ('drn', 'und'):
                for s in range(1, 6):
                    w.writerow(row('WBLL_%s_%s_s%d' % (tag, state, s),
                                   phi, b, state))
                    n += 1
    print('wrote %s  (%d cells)' % (p, n))
    print('   phi 0.09 -> b %.3f (Assur);  phi 0.21 -> b 0.05 (imposed, past '
          'phi_0)' % assur_b(0.09))

    # ---- TORLAYER: shear of a layered cell, swept in size ------------------
    # rve_torsion.csv sweeps cell size at phi_soft ~ 0.142, again as pockets
    # and channels. The layered counterpart has to hold the lamellar spacing
    # AND the bridge density fixed while the cell grows, which is the trap
    # Section 4.1 documents: holding counts instead changes the microstructure.
    # Cell edges are therefore chosen so that n_slabs comes out integer at a
    # fixed spacing a_0 = 0.125, and the bridge count is scaled with L^2 from a
    # density of 32 per unit area, which is high enough that integer rounding
    # costs a few percent rather than a factor.
    p = os.path.join(out, 'rve_torsion_layer.csv')
    n = 0
    phi, b = 0.14, assur_b(0.14)
    with open(p, 'w', newline='', encoding='utf8') as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        for L in (0.24, 0.36, 0.48, 0.60, 0.72):
            n_sl = int(round(L / 0.12))
            n_br = max(1, int(round(32.0 * L * L)))
            for s in (1, 2, 3, 4, 5, 6):
                r = dict(BASE)
                r['run_id'] = 'TORL_L%03d_s%d' % (round(L * 1000), s)
                r['Kappa'] = '0.11'
                r['Bending_Plane'] = 'torsion'
                r['L'] = '%.3f' % L
                r['n_slabs'] = str(n_sl)
                r['n_bridges'] = str(n_br)
                r['slab_vof'] = '%.4f' % phi
                r['bridge_fraction'] = '%.4f' % b
                r['K_inclusion'] = K['drn']
                r['full_tensor'] = 'No'       # the -tor load case carries it
                t = phi * L / (n_sl * max(1.0 - b, 1e-6))
                r['L_mesh'] = '%.4f' % min(max(t / ELEM_ACROSS, LM_MIN), LM_MAX)
                w.writerow([r[c] for c in COLS])
                n += 1
    print('wrote %s  (%d cells)' % (p, n))
    print('   a_0 held at 0.120, bridge density at 32/unit area, full tensor on')
    for L in (0.24, 0.36, 0.48, 0.60, 0.72):
        n_sl = int(round(L / 0.12)); n_br = max(1, int(round(32.0 * L * L)))
        print('     L=%.3f: %d slabs (a_0=%.4f), %d bridges (density %.1f)'
              % (L, n_sl, L / n_sl, n_br, n_br / (L * L)))

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

    # ---- ERGLAYER: the size sweep of Section 4.4.1, on layered cells --------
    # rve_eringen.csv bends pocket-and-channel cells at six sizes and regresses
    # the apparent modulus on 1/L^2. Its conclusion cannot be carried to the
    # layered regime, because the bound it produces is expressed in inclusion
    # diameters and a layered material's length is the lamellar spacing. This
    # deck repeats the measurement where that length is the one in play.
    #
    # a_0 = 0.12 is held exactly (L is an integer multiple at every size), so
    # the layer thickness -- and hence the element count across it -- does not
    # drift with cell size. Bridges are scaled with L^2 from a density of 32
    # per unit area. Drained, since the base is drained at every depth where
    # the layered description applies.
    A0_SWEEP = 0.12
    SWEEP_L = (0.24, 0.36, 0.48, 0.60, 0.72)
    phi, b = 0.14, assur_b(0.14)
    t_sweep = phi * A0_SWEEP / max(1.0 - b, 1e-6)
    lm_sweep = 0.0167          # nodes/face <= 5300 out to L = 0.60

    p = os.path.join(out, 'rve_eringen_layer.csv')
    n = 0
    with open(p, 'w', newline='', encoding='utf8') as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        for L in SWEEP_L:
            n_sl = int(round(L / A0_SWEEP))
            n_br = max(1, int(round(32.0 * L * L)))
            for s_ in (1, 2, 3, 4, 5, 6):
                r = dict(BASE)
                r['run_id'] = 'ERGL_L%03d_s%d' % (round(L * 1000), s_)
                r['Mode2'] = ''
                r['Disp2'] = ''
                r['L'] = '%.3f' % L
                r['L_mesh'] = '%.4f' % lm_sweep
                r['Kappa'] = '0.11'          # bending, as the pocket sweep
                r['n_slabs'] = str(n_sl)
                r['n_bridges'] = str(n_br)
                r['slab_vof'] = '%.4f' % phi
                r['bridge_fraction'] = '%.4f' % b
                r['K_inclusion'] = K['drn']
                w.writerow([r[c] for c in COLS])
                n += 1
    print('wrote %s  (%d cells)' % (p, n))
    print('   a_0 held at %.3f exactly; t=%.4f, L_mesh=%.4f -> %.1f elements '
          'across at every size' % (A0_SWEEP, t_sweep, lm_sweep,
                                    t_sweep / lm_sweep))
    for L in SWEEP_L:
        n_sl = int(round(L / A0_SWEEP)); n_br = max(1, int(round(32.0 * L * L)))
        print('     L=%.2f: %d slabs (a_0=%.4f), %d bridges (density %.1f), '
              '%d elem/edge' % (L, n_sl, L / n_sl, n_br, n_br / (L * L),
                                round(L / lm_sweep)))

    # The matched control. A phi=0 cell has no microstructure of either kind,
    # so this measures the extraction bias alone -- the artefact of imposing a
    # plate-like curvature on a cubic cell, which the pocket sweep found to be
    # as large as the trend it would otherwise have been read as.
    p = os.path.join(out, 'rve_eringen_layer_homog.csv')
    n = 0
    with open(p, 'w', newline='', encoding='utf8') as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        for L in SWEEP_L:
            r = dict(BASE)
            r['run_id'] = 'ERGLH_L%03d' % round(L * 1000)
            r['Mode2'] = ''
            r['Disp2'] = ''
            r['L'] = '%.3f' % L
            r['L_mesh'] = '%.4f' % lm_sweep
            r['Kappa'] = '0.11'
            r['VoF_sphere'] = '0.0'
            r['VoF_void_sphere'] = '0.0'
            r['VoF_incl_sphere'] = '0.0'
            r['n_slabs'] = '0'
            r['n_bridges'] = '0'
            r['slab_vof'] = '0.0000'
            r['bridge_fraction'] = '0.0000'
            r['K_inclusion'] = K['drn']
            w.writerow([r[c] for c in COLS])
            n += 1
    print('wrote %s  (%d cells)' % (p, n))

    # ---- NLGLAYER: the base slice of Section 4.4.2, on its own geometry -----
    # The pocket sweep reloaded three slices; the deepest sits at phi = 0.150,
    # above both thresholds, and was reloaded as a channelled cell. This is the
    # same composition with the geometry the rest of this work assigns to that
    # depth. Tension and compression to 2%, plus a linear reference on the same
    # meshes so that a loss of convergence can be read as kinematic rather than
    # as a property of the discretisation.
    p = os.path.join(out, 'rve_nlgeom_layer.csv')
    n = 0
    phi, b = 0.150, assur_b(0.150)
    for_lm = mesh_for(phi, b)
    with open(p, 'w', newline='', encoding='utf8') as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        for tag, disp, flag in (('TEN', '+0.020', 'ON'),
                                ('CMP', '-0.020', 'ON'),
                                ('LIN', '+0.020', 'OFF')):
            for s_ in (1, 2, 3):
                r = dict(BASE)
                r['run_id'] = 'NLGL_%s_s%d' % (tag, s_)
                r['Mode2'] = ''
                r['Disp2'] = ''
                r['Disp'] = disp
                r['nlgeom_flag'] = flag
                r['L_mesh'] = '%.4f' % for_lm
                r['slab_vof'] = '%.4f' % phi
                r['bridge_fraction'] = '%.4f' % b
                r['K_inclusion'] = K['drn']
                w.writerow([r[c] for c in COLS])
                n += 1
    print('wrote %s  (%d cells)' % (p, n))
    print('   phi=%.3f, b=%.3f (the layer plane is %.0f%% ice), t=%.4f, '
          'L_mesh=%.4f -> %.1f elements across'
          % (phi, b, 100 * b, thickness(phi, b), for_lm,
             thickness(phi, b) / for_lm))

    # ---- LCOL p060: the one layercol condition that would not mesh ---------
    # At phi = 0.06 the layer is at its thinnest (t = 0.0137) and the bridges at
    # their widest (b = 0.45), and a pocket sphere straddling a layer plane
    # leaves a sliver facet Gmsh reports as an overlapping boundary and cannot
    # repair. Keeping the pockets clear of the planes is enough: min_distance
    # goes to 0.005, still a hundredth of the cell edge, with fresh seeds so the
    # packer starts from a different configuration. This condition is the low
    # end of the layered branch, so leaving it out would calibrate the takeover
    # from phi >= 0.08 alone.
    p = os.path.join(out, 'rve_layercol_p060.csv')
    n = 0
    phi, b = 0.06, assur_b(0.06)
    with open(p, 'w', newline='', encoding='utf8') as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        for state in ('drn', 'und'):
            for s_ in (4, 5, 6):
                r = dict(BASE)
                r['run_id'] = 'LCOL_p060_%s_s%d' % (state, s_)
                r['slab_vof'] = '%.4f' % phi
                r['bridge_fraction'] = '%.4f' % b
                r['K_inclusion'] = K[state]
                r['L_mesh'] = '%.4f' % mesh_for(phi, b)
                r['min_distance'] = '0.005'
                w.writerow([r[c] for c in COLS])
                n += 1
    print('wrote %s  (%d cells)' % (p, n))
    print('   phi=%.2f, b=%.4f, t=%.4f, min_distance 0.002 -> 0.005'
          % (phi, b, thickness(phi, b)))

    # ---- TORLAYER control: the matched phi=0 cells for the torsion sweep ----
    # Section 4.4.1's own finding is that a size sweep without a matched control
    # measures the extraction bias as much as the material, and the bias there
    # was as large as the trend. A torsion sweep is no different, so the layered
    # torsion cells get phi=0 twins at the same five edges and the same element
    # size. Without these the layered sweep is uninterpretable.
    p = os.path.join(out, 'rve_torsion_layer_homog.csv')
    n = 0
    phi, b = 0.14, assur_b(0.14)
    with open(p, 'w', newline='', encoding='utf8') as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        for L in (0.24, 0.36, 0.48, 0.60, 0.72):
            n_sl = int(round(L / 0.12))
            t = phi * L / (n_sl * max(1.0 - b, 1e-6))
            r = dict(BASE)
            r['run_id'] = 'TORLH_L%03d' % round(L * 1000)
            r['Kappa'] = '0.11'
            r['Bending_Plane'] = 'torsion'
            r['L'] = '%.3f' % L
            r['L_mesh'] = '%.4f' % min(max(t / ELEM_ACROSS, LM_MIN), LM_MAX)
            r['VoF_sphere'] = '0.0'
            r['VoF_void_sphere'] = '0.0'
            r['VoF_incl_sphere'] = '0.0'
            r['n_slabs'] = '0'
            r['n_bridges'] = '0'
            r['slab_vof'] = '0.0000'
            r['bridge_fraction'] = '0.0000'
            r['K_inclusion'] = K['drn']
            r['full_tensor'] = 'No'
            w.writerow([r[c] for c in COLS])
            n += 1
    print('wrote %s  (%d cells)' % (p, n))


if __name__ == '__main__':
    main()
