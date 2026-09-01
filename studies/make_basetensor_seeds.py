"""Full-6x6-tensor replicates of the warm base slice.

The paper reports the base stiffness tensor -- and with it the in-plane
isotropy statement E_x = E_y -- from a single solved cell (CTEN_z95). That one
cell splits E_x = 4.85 against E_y = 5.02 GPa, a 3.3% in-plane difference. It
is about 1.7 standard deviations of the in-plane replicate scatter, so it is
consistent with the small number of channels one cell can hold rather than a
real in-plane texture, but with n = 1 the question cannot be settled.

This deck replicates the base slice N_SEED times with full_tensor=Yes, so each
packing is solved in all six load cases and the in-plane split becomes an
ensemble statement: if E_y/E_x scatters about unity across packings the
isotropy holds and the single-cell split was a realisation effect; if it sits
systematically above unity the base is genuinely orthotropic in plane.

The slice physics is reconstructed exactly as in make_colseeds.py -- from the
echoed parameters of results_column.csv -- so these replicates are directly
comparable both to the existing five-packing colseeds ensemble (E_x mean
4.849 GPa) and to the single full-tensor slice CTEN_z95 (E_x 4.853 GPa).

Not to be confused with the earlier `rve_basetensor.csv` deck, which sweeps four
*different* depths (BTEN_z65..z95) with one packing each to show where the
anisotropy switches on. This one holds the depth fixed at the base and varies
the packing instead, which is the question that deck cannot answer.

    cd params && python3 ../studies/make_basetensor_seeds.py   # -> rve_basetensor_seeds.csv

Then, from a directory holding SpaX_Standalone.py (same seed as colseeds, so
row index -> packing is reproducible):

    SPAX_SEED=20260723 python3 SpaX_Standalone.py params/rve_basetensor_seeds.csv out_basetensor_seeds/

full_tensor=Yes: six load cases per RVE (utx,uty,utz,ss12,ss13,ss23), so
N_SEED x 6 = 30 solves.
"""
import csv
import os

from make_ice_studies import BASE, COLS, write

N_SEED = 5
# Which depth slice is replicated. The base, z/H = 0.95, is the one the paper
# reports; SPAX_BT_ZSLICE=85 replicates the slice above it, which carries the
# other half of the anisotropy claim (Table 4 shows the scatter rising in the
# bottom *two* slices, not one) and is asked for by referee comment M9.
Z_SLICE = int(os.environ.get('SPAX_BT_ZSLICE', 95))
COLUMN_CSV = os.path.join(os.path.dirname(__file__), '..', 'results',
                          'results_column.csv')

# Deck fields echoed verbatim in results_column.csv (same list as make_colseeds).
ECHOED = ['L', 'L_mesh', 'Is_Porous', 'E_matrix', 'nu_matrix', 'VoF_sphere',
          'r_avg', 'VoF_void_sphere', 'VoF_incl_sphere', 'E_sphere_inclusion',
          'sphericity_avg', 'PBC_Method', 'Bending_PBC_Type', 'Growth_Direction',
          'generate_channels', 'channel_vof_target']


def study_basetensor(L=None, prefix='BTEN', out=None, full_tensor=True):
    """Replicate the base slice N_SEED times, each solved for the full 6x6.

    With `L` given, the cell edge is overridden while every inclusion size,
    volume fraction and the mesh size are left untouched -- so the cell holds
    proportionally more inclusions at the same physical resolution, which is
    how the size sweep (rve_sizechan.csv) grows its cells. This is the knob
    that matters for in-plane isotropy at the base: the scatter is set by how
    few channels fit, so only a bigger cell reduces it.
    """
    src = None
    for x in csv.DictReader(open(COLUMN_CSV)):
        if int(round(float(x['run_id'].split('z')[-1]))) == Z_SLICE:
            src = x
            break
    if src is None:
        raise SystemExit(f'no z{Z_SLICE} slice in {COLUMN_CSV}')

    r = dict(BASE)
    r.update({k: src[k] for k in ECHOED})
    r['Growth_Concentration'] = f'{0.40 + 0.32 * (Z_SLICE / 100.0):.2f}'
    # 6 load cases (utx,uty,utz,ss12,ss13,ss23) or just utx+utz
    r['full_tensor'] = 'Yes' if full_tensor else 'No'
    if L is not None:
        r['L'] = f'{float(L):.2f}'               # L_mesh deliberately unchanged

    rows = []
    for s in range(1, N_SEED + 1):
        rr = dict(r)
        rr['run_id'] = f'{prefix}_z{Z_SLICE:02d}_s{s}'
        rows.append(rr)

    default = ('rve_basetensor_seeds.csv' if Z_SLICE == 95 else
               'rve_basetensor%d_seeds.csv' % Z_SLICE)
    write(out or default,
          [{c: row.get(c, BASE.get(c, '')) for c in COLS} for row in rows])


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:                        # e.g. `... 0.80 BT80 [first]`
        L = float(sys.argv[1])
        pref = sys.argv[2] if len(sys.argv) > 2 else 'BT%02d' % round(L * 100)
        ft = not (len(sys.argv) > 3 and sys.argv[3].lower().startswith('first'))
        study_basetensor(L=L, prefix=pref, full_tensor=ft,
                         out='rve_basetensor_%s.csv' % pref.lower())
    else:
        study_basetensor()
