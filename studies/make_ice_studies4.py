"""Brine-modulus sensitivity: base meshes only.
Two microstructures (isolated pockets vs channelled base). The brine modulus
sweep is then applied by patch_brine.py, which copies these meshed decks and
rewrites ONLY the inclusion *Elastic card -- so geometry is byte-identical across
the sweep and the response difference is purely the soft-phase modulus.
"""
from make_ice_studies import row, write, E_matrix

def study_brine_base():
    E_mat = E_matrix(-8.0); gas = 0.012
    rows = [
        row('BRINE_iso',  E_mat, 0.05, gas, 0.75, 0.55),                    # isolated pockets
        row('BRINE_chan', E_mat, 0.08, gas, 0.65, 0.65, channels_frac=0.40),# channelled base
    ]
    write('rve_brine_base.csv', rows)

if __name__ == '__main__':
    study_brine_base()
