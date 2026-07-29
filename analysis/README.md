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

**Running:** most analyzers read result CSVs by bare filename (e.g.
`results_tilt00.csv`). Run them with those files on the path — either copy the
relevant `../results/*.csv` alongside, or invoke from a directory where they are
reachable. Figures are written next to the CSVs by default.
Use a Python that has matplotlib and pandas available (a conda environment, for
instance); a bare system `python3` may lack them.
