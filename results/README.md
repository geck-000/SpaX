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

## Legacy campaign tables

The tables below predate this folder — they are the first sea-ice campaign, run
when everything still lived in the repository root, and they are the only
surviving data for the studies written up in `claude.tex` (the older campaign
report on Overleaf). They were untracked on one machine until 2026-07-29:
`results_bending`, `results_brine`, `results_bten`, `results_chan_ft`,
`results_channel`, `results_channels_q`, `results_fymy`, `results_gas`,
`results_homog*`, `results_lscale`, `results_marchenko`, `results_mono{,2}`,
`results_morph`, `results_old`, `results_orient`, `results_perc`,
`results_porous_q3{,_interim}`, `results_seas`, `results_seeds`. Their input
decks are the matching `../params/rve_*.csv` restored in the same commit.

> **Scatter convention.** These legacy tables were summarised in `claude.tex`
> using the *sample* standard deviation (`ddof=1`). The current manuscript
> (`main_rev.tex` §4.1) declares and uses the *population* standard deviation
> (`ddof=0`), which with five packings is ~12% narrower. `results_seeds.csv` is
> the provenance for claude.tex's "0.3%, 0.5%, 1.9%" three-configuration
> replicate study — under the population convention those are 0.2%, 0.5%, 1.7%.
