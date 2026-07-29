# Handover — revision campaign state

Written 2026-07-27. Everything below is the state at SpaX `main` = `9fc3487`
and Overleaf `main` = `3e7f8bd`. Read this first if you are picking the work up
on another machine.

---

## 1. Where things live

| What | Where |
|---|---|
| Toolkit | `SpaX_Standalone.py`, `SpaX_GmshPeriodic.py`, `SpaX_PostProcess.py` (top level) |
| Deck builders | `studies/` — one script per campaign, writes a parameter CSV |
| Parameter CSVs | `params/rve_*.csv` |
| HPC submit/post | `hpc/submit_*.sh`, `hpc/postprocess_*.sh` |
| Result CSVs | `results/results_*.csv` |
| Stiffness matrices | `post_coltensor/`, `post_basetensor_seeds/`, `post_bt80/` |
| Analysis + figures | `analysis/` (plain `python3`, no Abaqus needed) |
| Renders | `viz/render_rve.py`, `viz/odb_to_vtk.py` |

The two manuscripts are **not** in this repo. They live in the Overleaf project
(`main_rev.tex`, `software.tex`, `revision_notes.tex`, `claude.tex`), cloned via
the git bridge. `claude.tex` is the older campaign report and is the only record
of six studies whose CSVs predate `results/` — do not delete it.

`out_*/` directories hold generated Abaqus decks. They are **derived and
untracked**: ~2.8 GB in total, regenerable from `params/`. Exclude them from any
zip unless you specifically need them.

---

## 2. What state the papers are in

Both compile clean, no undefined references, every float cited.

- `main_rev.tex` — 21 pp. Revision of the modelling paper, responding to three
  reviewers. All reviewer rows are actioned; `revision_notes.tex` (9 pp) is the
  response document and is the record of traceability.
- `software.tex` — 12 pp, SoftwareX. Every number re-verified against
  `results/` on 2026-07-27.

**No published number needed revision as a result of this campaign.** Two
scares were investigated and both resolved as artefacts — see §4.

---

## 3. To resume on a new machine

### Environment

```bash
pip install numpy gmsh                 # generation + meshing
pip install matplotlib pandas          # analysis figures
pip install pyvista                    # optional, 3-D renders
```

Abaqus is needed only to *solve* decks and to read ODBs. Generation and the
whole analysis stage run without it. On this project the solving is done on CSC
Roihu, not locally.

### Regenerating the manuscript figures

```bash
cd results && PYTHONPATH=../analysis python3 ../analysis/make_rev_figs.py
```

Writes `study_scfdepth`, `ice_column_profiles`, `study_coltensor` as PDF+PNG.
Copy the PNGs into the Overleaf `figures/` directory.

### Running a campaign end to end

```bash
cd params && python3 ../studies/make_<campaign>.py          # -> rve_<campaign>.csv
cd .. && SPAX_SEED=20260723 python3 SpaX_Standalone.py params/rve_<campaign>.csv out_<campaign>/
rsync -az out_<campaign>/ <cluster>:/scratch/project_XXXXXX/test_rve/ --include='Job-*' --exclude='*'
ssh <cluster> "cd /scratch/project_XXXXXX/test_rve && bash submit_<campaign>.sh"
# ... then pull results_*.csv or post_*/ back
```

`SPAX_SEED` fixes the packing: the seed is resolved per *row index*, so the same
seed with a different number of rows gives different packings. Two campaigns are
only comparable packing-for-packing if their CSVs have the same shape.

---

## 4. Findings from this campaign

### Confirmed, no change needed

- **The warm-base modulus is box-size converged.** Four cell sizes at base
  microstructure (`L=0.50/0.65/0.80/1.00`, 5 packings each): `E_x` = 4.79, 4.77,
  4.70, 4.66 GPa — −2.6% across a doubling of the cell edge. The published
  sweep (`results_sizechan.csv`) only ever ran at a third of the base
  soft-phase fraction, so this closed a real gap in the argument.
- **In-plane isotropy holds at the base.** A `L=0.50` cell holds only 3–5
  channels, too few for the two in-plane directions to be equivalent within one
  realisation; it shows an apparent ~1% split. At `L=0.80` (10–11 channels)
  `E_y/E_x = 0.998 ± 0.010` with the packings straddling unity.
- **The channel generator is unbiased in plane** —
  `analysis/check_channel_isotropy.py`, three seed streams, all CIs spanning
  zero.

### Corrected

- **The warm base was re-centred** from a single packing (4.47 GPa, a ~6σ low
  outlier) onto the five-packing mean, 4.85 GPa. This propagated to Table 2,
  §4.2, the study-overview table, the Conclusions, the summary schematic — and
  to the laminated-plate model, which had been left on the outlier
  (`analysis/recentre_column.py`; `B/sqrt(AD)` 0.12 → 0.11).
- **§4.4's base stiffness matrix** was still the pre-recentring run and
  contradicted its own shear ratio. Replaced with the current full-tensor values.

### Two scares that were artefacts

- A "systematic" in-plane split at `L=0.50` — a small-sample coincidence, see
  above.
- A "39% non-convergence" of the base modulus at `L=0.80` — **a bug in this
  repo**, now fixed. See §5.

---

## 5. Traps — read before running anything

- **`extract_elasticity_tensor` used to assume the applied strain.** It
  defaulted to `0.01`, while the decks prescribe a fixed *displacement*
  (`Disp=0.005`), so the true strain is `Disp/L`. At `L=0.50` those coincide
  exactly, so it was silently right for every campaign until the cell size
  changed, then silently rescaled the whole tensor. **Fixed** — it now reads the
  deck and raises rather than guessing. Ratios were always immune, since a
  common factor cancels.
- **`csc_solve_array.sh` skips a deck if its `.odb` file exists**, testing for
  the file rather than a completed solve. A job killed at walltime leaves a
  truncated ODB, so a naive resubmit silently skips exactly the jobs that
  failed. Identify failures by the absence of an `Abaqus exit:` line in the
  log, and delete those ODBs before resubmitting.
- **`module load abaqus` is broken in batch** on the CSC v2026_03 stack. Source
  `~/abaqus_env.sh` instead. Three scripts still carry the broken form —
  `postprocess_firstorder.sh`, `postprocess_nlgeom.sh`,
  `postprocess_coltensor.sh` — and three are correct:
  `postprocess_basesweep.sh`, `_basetensor.sh`, `_bt80.sh`. Copy from a
  correct one when writing a new campaign's post-processor.
- **Generation forks one mesher per core** unless `SPAX_GEN_WORKERS` says
  otherwise. At ~7×10⁵ elements each worker takes ~3 GB, so on a 16 GB laptop
  set `SPAX_GEN_WORKERS=2`, and `=1` above ~10⁶ elements. Killing the parent
  leaves **orphaned `SpaX_GmshPeriodic` children** still holding the memory —
  kill those explicitly, and note their command line says `GmshPeriodic`, not
  `Standalone`, so an obvious `pkill` pattern misses them.
- **Periodic meshing retries are normal at high inclusion counts.** `L=1.00`
  needed 3 re-packs of a 6-attempt budget. A traceback in the log is usually the
  retry mechanism working, not a failure. Check `GENERATION COMPLETE: N RVEs`
  and the per-deck node counts before assuming the decks are bad.
- **Scatter convention**: everything quoted is the *population* standard
  deviation (`np.std`, `ddof=0`), now stated explicitly in §4.1.
- **CSC scratch is shared** — delete decks and ODBs after pulling results.

---

## 6. Open items

Nothing is blocking. What remains is editorial, and is listed in full in
Section G of `revision_notes.tex`:

- author/affiliation/funding fields in both papers
- a citable Zenodo/Code Ocean release for the software paper's C3 metadata row
- whether to merge Chapters 2–3, and whether to merge the two depth-profile
  figures — both deliberately left as author decisions
- the three-configuration replicate study behind "0.3%, 0.5%, 1.9%" in §4.1 has
  no surviving CSV (it predates `results/`). It is left as published; the
  ten-slice campaign gives the same base quantity as `1.041 ± 0.010` and is
  reproducible, if you would rather quote that.
