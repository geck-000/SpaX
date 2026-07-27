"""Replicate (seeded) column deck for the depth-profile scatter envelopes.

Reviewer request (Jani #12 / Reviewer 2): the depth profiles show a single packing
per slice; add mean/scatter envelopes. This builds the seeded deck by reconstructing
the *exact* per-slice microstructure of the paper's reference column
(results_column.csv -> Table 2 / Fig. column) -- its echoed VoF, gas, sphericity,
channel and matrix parameters -- and replicating each slice N_SEED times. Building
from results_column.csv (rather than re-deriving the physics) guarantees the
replicate means are directly comparable to the single-packing profile, so the
per-slice mean and +/-1 s.d. can re-centre the profile consistently.

The original ICE_z column deck is not in studies/ (it predates the study
generators), but results_column.csv preserves every deck parameter except the
per-slice Growth_Concentration (pocket-orientation strength), which barely affects
the effective-modulus magnitude and is taken from the columnar-growth formula
grw = 0.40 + 0.32 z used by the sibling column studies.

    cd params && python3 ../studies/make_colseeds.py      # -> rve_colseeds.csv

Then, from a directory holding SpaX_Standalone.py (fixed seed = reproducible):

    SPAX_SEED=20260723 python3 SpaX_Standalone.py params/rve_colseeds.csv out_colseeds/

First-order (full_tensor=No): two load cases per RVE (Uniaxial Tension X and Z),
so 10 slices x N_SEED replicates x 2 solves.
"""
import csv, os
from make_ice_studies import BASE, COLS, write

N_SEED = 5
COLUMN_CSV = os.path.join(os.path.dirname(__file__), '..', 'results',
                          'results_column.csv')

# Deck fields echoed verbatim in results_column.csv (results-col -> deck-col).
ECHOED = ['L', 'L_mesh', 'Is_Porous', 'E_matrix', 'nu_matrix', 'VoF_sphere',
          'r_avg', 'VoF_void_sphere', 'VoF_incl_sphere', 'E_sphere_inclusion',
          'sphericity_avg', 'PBC_Method', 'Bending_PBC_Type', 'Growth_Direction',
          'generate_channels', 'channel_vof_target']


def study_colseeds():
    """Reconstruct the reference column slice-by-slice from results_column.csv and
    replicate each slice N_SEED times (distinct run_id -> distinct row index ->
    independent packing under a fixed SPAX_SEED)."""
    rows = []
    for x in csv.DictReader(open(COLUMN_CSV)):
        z = int(round(float(x['run_id'].split('z')[-1])))   # ICE_z95 -> 95
        r = dict(BASE)                                       # constant defaults
        r.update({k: x[k] for k in ECHOED})                 # exact per-slice physics
        r['Growth_Concentration'] = f'{0.40 + 0.32 * (z / 100.0):.2f}'
        for s in range(1, N_SEED + 1):
            rr = dict(r)
            rr['run_id'] = f'CSEED_z{z:02d}_s{s}'
            rows.append(rr)
    # write() expects rows keyed by COLS
    write('rve_colseeds.csv', [{c: row.get(c, BASE.get(c, '')) for c in COLS}
                               for row in rows])


if __name__ == '__main__':
    study_colseeds()
