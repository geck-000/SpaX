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
| `results_column_ensemble.csv` | The production depth column: the ensemble mean over independent packings, and the file the manuscript numbers are drawn from. |
| `results_colseeds.csv`, `results_colseeds_extra.csv`, `results_colseeds_all.csv` | The column replicate ensemble — first five packings, the additional five, and the merged ten. |
| `results_fieldseeds.csv` | Replicates of the three field-comparison columns. |
| `results_eringen.csv`, `results_eringen_homog.csv` | The bending size sweep and its matched inclusion-free (`φ=0`) control. The control is what makes the sweep interpretable — see `../docs/USER_DOCS.md` §6. |
| `results_basesweep*.csv`, `results_steep_column.csv` | Base-slice size convergence, and the steeply monotonic salinity column. |
| `results_skeletal.csv`, `results_skeletal_laminae.csv` | Skeletal basal layer, whole-slice and resolved into sub-laminae. |
| `results_weibull.csv` | Weakest-link sensitivity of the localisation measure. |

**Curves** (`curves_nlgeom_*.csv`) — per-frame nominal stress–strain for the
large-deformation study.

**Figures** (`study_*.png`) — the analysis plots (`study_brineK`, `study_nlgeom`,
`study_tilt`, `study_coltensor`, `study_scf`, `study_failure`,
`study_macro_plate`), regenerable from the CSVs via `../analysis/`. The
manuscript figures carry their own names — `ice_column_profiles`,
`kujala_comparison`, `field_decomposition`, `stress_field_3d` — and are written
as PDF + PNG by `../analysis/make_rev_figs.py` and the plotting scripts beside
it.

## Legacy campaign tables

The tables below predate this folder — they are the first sea-ice campaign, run
when everything still lived in the repository root.
`results_bending`, `results_brine`, `results_bten`, `results_chan_ft`,
`results_channel`, `results_channels_q`, `results_fymy`, `results_gas`,
`results_homog*`, `results_lscale`, `results_marchenko`, `results_mono{,2}`,
`results_morph`, `results_old`, `results_orient`, `results_perc`,
`results_porous_q3{,_interim}`, `results_seas`, `results_seeds`.

> **Scatter convention.** The current manuscript declares and uses the *population* standard deviation
> (`ddof=0`).

## Stress-field dumps held outside the repository

`scf_extract.py` keeps the full per-element stress array rather than discarding
it once the percentiles are formed, because `weibull_sensitivity.py`,
`weibull_mnorm.py` and `plot_scf_field.py` all re-read it. Those dumps are
`.npz` and they are large.

Thirty of them are tracked here under `weibull_dumps/` and `scf_fields/`,
totalling about 100 MB. The twenty **layered** ones, `WBLL_*.npz`, are not: at
1.1 GB they would have to live in the history permanently. They are archived
outside the repository at

    ~/SpaX/spax_scratch_archive/weibull_dumps/

having been pulled off `/scratch` on 2026-08-19 and verified against the remote
byte sizes and by loading each array. That copy is now the only one — the
scratch originals remain for the moment but scratch is a shared allocation with
no backup, so anything needed from `WBLL_*` should be taken from the archive
rather than assumed to be on the cluster.
