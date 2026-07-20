# results/ — homogenisation outputs & figures

Effective properties produced by `SpaX_PostProcess.py` from the solved ODBs, the
nlgeom stress–strain curves, and the analysis figures rendered by the scripts in
`../analysis/`.

**Effective-property tables** (`results_*.csv`) — one row per RVE, with effective
moduli, anisotropy ratios, and (where computed) tensor / SCF quantities:

| File | Study |
|------|-------|
| `results_column.csv` | Depth-graded column. |
| `results_coltensor.csv` | Full 6×6 tensor down the column. |
| `results_percolation` / `sizechan` / `salfamily` … | Corresponding `../params/` sweeps. |
| `results_brineK{const,temp}.csv` | Brine `K(T)` study [#5]. |
| `results_nlgeom_{lin,ten,cmp}.csv` | Large-deformation study [#8]. |
| `results_tilt{00,15,30}.csv` | Channel-tilt study [#6]. |
| `results_scf.csv`, `results_failure.csv` | Stress-concentration / failure-onset. |
| `results_macro_plate.csv` | Laminated-plate macro assembly. |

**Curves** (`curves_nlgeom_*.csv`) — per-frame nominal stress–strain for the
large-deformation study.

**Figures** (`study_*.png`) — the analysis plots (`study_brineK`, `study_nlgeom`,
`study_tilt`, `study_coltensor`, `study_scf`, `study_failure`,
`study_macro_plate`), regenerable from the CSVs via `../analysis/`.
