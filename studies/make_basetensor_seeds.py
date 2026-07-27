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
Z_SLICE = 95                                     # warm base, z/H = 0.95
COLUMN_CSV = os.path.join(os.path.dirname(__file__), '..', 'results',
                          'results_column.csv')

# Deck fields echoed verbatim in results_column.csv (same list as make_colseeds).
ECHOED = ['L', 'L_mesh', 'Is_Porous', 'E_matrix', 'nu_matrix', 'VoF_sphere',
          'r_avg', 'VoF_void_sphere', 'VoF_incl_sphere', 'E_sphere_inclusion',
          'sphericity_avg', 'PBC_Method', 'Bending_PBC_Type', 'Growth_Direction',
          'generate_channels', 'channel_vof_target']


def study_basetensor():
    """Replicate the base slice N_SEED times, each solved for the full 6x6."""
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
    r['full_tensor'] = 'Yes'                     # 6 load cases

    rows = []
    for s in range(1, N_SEED + 1):
        rr = dict(r)
        rr['run_id'] = f'BTEN_z{Z_SLICE:02d}_s{s}'
        rows.append(rr)

    write('rve_basetensor_seeds.csv', [{c: row.get(c, BASE.get(c, '')) for c in COLS}
                                 for row in rows])


if __name__ == '__main__':
    study_basetensor()
