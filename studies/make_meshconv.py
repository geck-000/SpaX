"""Does the meshed brine fraction converge on the placed one as the mesh refines?

The volume of material the mesh assigns to the brine phase exceeds the volume
the packer places, by about 1.5x at the cold surface and 2.4x at the warm base.
The packer is accurate -- it reports ~105% of target, and the gas voids, which
are absent from the mesh rather than assigned to it, hit their targets exactly.
The suspicion is therefore discretisation: at the resolution used throughout,
basal inclusions are only two to three elements across their diameter, and an
element straddling a phase boundary is assigned whole.

If that is right, refining the mesh at fixed geometry must drive the assigned
fraction down toward the placed one. If it is not, the assigned fraction will
sit still and the explanation is wrong.

Two slices bracket the effect -- z05, where it is mildest, and z95, where it is
worst -- at three mesh sizes each, two packings apiece. Nothing needs solving:
the assigned volume is a property of the deck, so generation alone answers it.

    cd params && python3 ../studies/make_meshconv.py   # -> rve_meshconv.csv
"""
import csv
import os

from make_ice_studies import BASE, COLS, write

RES = os.path.join(os.path.dirname(__file__), '..', 'results')
SRC = 'results_marchenko.csv'
SLICES = (5, 95)
MESHES = (0.033, 0.0231, 0.0165)      # baseline, 0.7x, 0.5x
N_SEED = 2

ECHOED = ['L', 'Is_Porous', 'E_matrix', 'nu_matrix', 'VoF_sphere', 'r_avg',
          'VoF_void_sphere', 'VoF_incl_sphere', 'E_sphere_inclusion',
          'sphericity_avg', 'PBC_Method', 'Bending_PBC_Type',
          'Growth_Direction', 'generate_channels', 'channel_vof_target']


def depth_of(run_id):
    tail = run_id.rsplit('_z', 1)[-1]
    return int(round(float(''.join(c for c in tail if c.isdigit()))))


def study_meshconv():
    path = os.path.join(RES, SRC)
    rows = []
    for x in csv.DictReader(open(path)):
        if not x.get('E_eff') or float(x['E_eff']) <= 0:
            continue
        z = depth_of(x['run_id'])
        if z not in SLICES:
            continue
        base = dict(BASE)
        for k in ECHOED:
            if k in x and x[k] != '':
                base[k] = x[k]
        base['Growth_Concentration'] = '%.2f' % (0.40 + 0.32 * (z / 100.0))
        for lm in MESHES:
            for s in range(1, N_SEED + 1):
                r = dict(base)
                r['L_mesh'] = '%.4f' % lm
                r['run_id'] = 'MCONV_z%02d_m%03d_s%d' % (z, round(lm * 10000), s)
                rows.append(r)
    write('rve_meshconv.csv',
          [{c: row.get(c, BASE.get(c, '')) for c in COLS} for row in rows])
    print('wrote rve_meshconv.csv: %d slices x %d meshes x %d packings = %d RVEs'
          % (len(SLICES), len(MESHES), N_SEED, len(rows)))
    for lm in MESHES:
        print('   L_mesh %.4f  -> ~%.1f elements across a 0.083 basal inclusion'
              % (lm, 0.083 / lm))


if __name__ == '__main__':
    study_meshconv()
