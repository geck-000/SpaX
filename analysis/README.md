# analysis/ — result analyzers & field extractors

Post-solve analysis: turn the `../results/*.csv` produced by
`SpaX_PostProcess.py` into figures and derived quantities. These are plain
`python3` (matplotlib / numpy) except where noted — no Abaqus needed.

| File | Role |
|------|------|
| `analyze_studies.py` | General study analyzer / summary tables across the campaign CSVs. |
| `analyze_brineK.py` | Brine `K(T)` study [#5] figure (`study_brineK.png`). |
| `analyze_nlgeom.py` | Large-deformation study [#8] figure (`study_nlgeom.png`). |
| `analyze_tilt.py` | Channel-tilt study [#6] figure (`study_tilt.png`): `E_z/E_x` vs tilt. |
| `make_rev_figs.py` | The manuscript depth-profile figures (`study_scfdepth`, `ice_column_profiles`, `study_coltensor`) as PDF + PNG; `--out DIR` to redirect. |
| `figstyle.py` | Shared figure conventions (Okabe-Ito palette, enlarged fonts, `z/H` on the vertical axis) imported by `make_rev_figs.py` and the `analyze_*` study figures so every depth profile matches Fig.2. |
| `aggregate_coltensor.py` | Aggregates per-slice 6×6 elasticity tensors into a depth profile. |
| `macro_plate.py` | Laminated-plate (classical lamination theory) whole-sheet macro assembly. |
| `nlgeom_extract.py` | Extracts nominal stress–strain from the nlgeom reference-point reaction. |
| `scf_extract.py` | Stress-concentration-factor percentiles from matrix principal stress. |
| `failure_extract.py` | Failure-onset extractor (`abaqus python`). Extends `scf_extract.py`: Mohr–Coulomb and max-principal criteria over a slice's `utx` ODB. |
| `analyze_failure.py` | Failure-onset figure (`study_failure.png`) from `results_failure.csv`. |
| `plot_ice_column.py` | Depth profiles and cross-plots of the graded column. |
| `plot_tensor.py` | Per-slice 6×6 elasticity-tensor figure. |
| `plot_marchenko.py`, `plot_marchenko_match.py` | Comparison against the Marchenko (2024) vibrating-beam field data — raw, and after the matrix/salinity matching. |
| `verify_column.py`, `verify_sizeeffect.py` | Self-checking harnesses. Re-derive every quoted number — for the depth-graded column and for the bending size effect respectively — and fail loudly if one drifts. Run from `../results/`. |
| `build_ensemble_column.py` | Builds the production depth column as the ensemble mean over independent packings. |
| `recentre_column.py` | Builds the re-centred depth column used by the laminated-plate macro model. |
| `skeletal_clt.py` | Neutral plane and flexural modulus of a graded sheet with a skeletal base. |
| `gradient_correction.py` | How much of the RVE-versus-field modulus offset is a through-thickness gradient artefact? |
| `plot_field_decomp.py` | Separates level from shape in the field comparison, so the profile ranking can be judged independently of the overall offset. |
| `plot_kujala.py` | Figure: the computed depth profile against the four-point bending decomposition. |
| `plot_sizeeffect.py` | Figure: the bending size effect against a matched no-microstructure baseline. |
| `fit_nonlocal.py` | Fits the bending size sweep to both nonclassical families — couple-stress and integral nonlocal — bounding each length scale. |
| `plot_scf_field.py` | Figure: where the stress concentrates inside the cell, not just how much. |
| `weibull_sensitivity.py` | Weakest-link sensitivity of the stress-localisation measure to the Weibull modulus. |
| `aggregate_basetensor_seeds.py` | Aggregates the base-slice full-tensor replicates into an ensemble statement. |
| `compare_basetensor_sizes.py` | Compares the warm-base full-tensor ensembles across cell size. |
| `check_channel_isotropy.py` | Tests whether the channel generator is isotropic in the RVE plane. |
| `prose_audit.py` | Locates the hardest-to-read prose in a LaTeX manuscript. Editorial aid, not part of the analysis chain. |

**Running:** most analyzers read result CSVs by bare filename (e.g.
`results_tilt00.csv`). Run them with those files on the path — either copy the
relevant `../results/*.csv` alongside, or invoke from a directory where they are
reachable. Figures are written next to the CSVs by default.
Use a Python that has matplotlib and pandas available (a conda environment, for
instance); a bare system `python3` may lack them.
