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
| `aggregate_coltensor.py` | Aggregates per-slice 6×6 elasticity tensors into a depth profile. |
| `macro_plate.py` | Laminated-plate (classical lamination theory) whole-sheet macro assembly. |
| `nlgeom_extract.py` | Extracts nominal stress–strain from the nlgeom reference-point reaction. |
| `scf_extract.py` | Stress-concentration-factor percentiles from matrix principal stress. |

**Running:** most analyzers read result CSVs by bare filename (e.g.
`results_tilt00.csv`). Run them with those files on the path — either copy the
relevant `../results/*.csv` alongside, or invoke from a directory where they are
reachable. Figures are written next to the CSVs by default.
Use the conda/matplotlib Python (`~/miniconda3/bin/python`); the system
`python3` may lack matplotlib.
