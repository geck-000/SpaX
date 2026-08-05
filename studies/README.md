# studies/ — RVE parameter-deck generators

Scripts that build the parameter CSVs (one row = one RVE) consumed by
`SpaX_Standalone.py`, and the deck "stamping" helpers used for paired studies
that must share a single mesh.

| File | Role |
|------|------|
| `make_ice_studies.py` | Base module: shared helpers (`row`, `write`, `phi_brine`, `E_matrix`, `temperature`, `ZS`, `COLS`, `G_BRINE`) and the first batch of studies. All other `make_*` scripts import from it. |
| `make_ice_studies2.py` … `make_ice_studies7.py` | Later study generators (percolation, morphology, orientation, brine `K(T)` [#5], nlgeom [#8], …), each writing one or more `rve_*.csv` decks into `../params/`. |
| `make_seaice_2nd.py` | Second-order (quadratic, channelled) RVE deck generator. |
| `make_colseeds.py` | Seeded first-year C-shape column (`rve_colseeds.csv`): the 10 depth slices replicated 5× each for the statistical scatter envelopes on the depth profiles. Independent packings via `SPAX_SEED` per-row reseeding. |
| `build_brineK_decks.py` | Stamps a temperature-dependent brine `K(T)` twin deck onto a single shared base mesh per slice (isolates the `K(T)` effect from mesh noise). |
| `build_nlgeom_decks.sh` | Stamps the linear / tension / compression nlgeom decks from one shared base mesh. |
| `make_colseeds_extra.py` | Additional packings per column slice, extending the replicate ensemble from five to ten. |
| `make_fieldseeds.py` | Replicate decks for the three field-comparison columns. |
| `make_basesweep.py` | Re-verification of the RVE-size convergence check at the base slice. |
| `make_basetensor_seeds.py` | Full 6×6-tensor replicates of the warm base slice. |
| `make_eringen.py` | Nonlocal length-scale study: extends the bending sweep across cell sizes. |
| `make_skeletal.py` | Skeletal basal layer — resolves the bottom few percent of the sheet as sub-laminae. |
| `make_weibull.py` | Weibull / weakest-link sensitivity of the stress-localisation measure. |
| `patch_brine.py` | Utility to patch brine material cards in existing decks. |

**Running:** the `make_*` scripts import each other, so run them from *this*
folder (`cd studies`). They write parameter CSVs; point `SpaX_Standalone.py`
(in the repo root) at the resulting `../params/*.csv`.
