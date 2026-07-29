# SpaX — user guide

Operating notes for running SpaX beyond a single deck: how the pipeline fits
together, how to keep a study reproducible, and what goes wrong in practice.

The repository `README.md` is the reference — parameter columns, environment
variables, output columns. This document is the part that reference material
does not tell you. `RUNBOOK.md` covers running on a cluster.

---

## 1. Installation

```bash
pip install numpy gmsh          # generation + periodic meshing
pip install matplotlib pandas   # analysis and figures
pip install pyvista             # optional, 3-D renders
```

Abaqus is needed only to **solve** decks and to read ODBs. Deck generation and
the entire analysis stage run without it, so a machine with no licence can still
build a campaign and reproduce every figure from the result tables.

Two Python contexts are in play, and mixing them is the most common first
stumble:

| Run with | For |
|----------|-----|
| `python3` | generation, all analysers, the `analyze` subcommands |
| `abaqus python` | anything that opens an ODB — the extraction pass, `principal`, `elasticity` |

---

## 2. Where things live

| What | Where |
|---|---|
| The toolkit | `SpaX_Standalone.py`, `SpaX_GmshPeriodic.py`, `SpaX_PostProcess.py` |
| Deck builders | `studies/` — one script per campaign, writes a parameter CSV |
| Parameter CSVs | `params/rve_*.csv` |
| Cluster submit/post | `hpc/submit_*.sh`, `hpc/postprocess_*.sh` |
| Result tables | `results/results_*.csv` |
| Stiffness matrices | `tensors/`, `post_coltensor/`, `post_basetensor_seeds/`, `post_bt80/` |
| Analysers and figures | `analysis/` (plain `python3`) |
| Renderers | `viz/render_rve.py`, `viz/odb_to_vtk.py` |

`out_*/` directories hold generated Abaqus decks. They are derived and
untracked — several GB, and regenerable from `params/`. Exclude them from any
archive unless you specifically need them.

---

## 3. A campaign end to end

```bash
cd studies && python3 make_<campaign>.py           # -> ../params/rve_<campaign>.csv
cd .. && SPAX_SEED=<seed> python3 SpaX_Standalone.py \
         params/rve_<campaign>.csv out_<campaign>/
# ... solve the decks with Abaqus (see RUNBOOK.md for the cluster path) ...
abaqus python SpaX_PostProcess.py params/rve_<campaign>.csv \
              out_<campaign>/ results/results_<campaign>.csv
cd results && python3 ../analysis/<analyzer>.py
```

Analysers read their inputs by bare filename, so run them from `results/`.

---

## 4. Reproducibility

### Seeding is per row index

`SPAX_SEED` fixes the packing, but the seed is resolved **per row index**. The
same seed with a different number of rows gives different packings. Two
campaigns are only comparable packing-for-packing if their CSVs have the same
shape.

### Mesh generation is not deterministic across runs

Inclusions are placed randomly, so two runs of the same deck do not produce the
same mesh. Any study that compares a *material* change must therefore control
the geometry, or the mesh-discretisation scatter will bury the signal.

The way to do that is to generate one mesh and stamp the variants onto it —
rewrite only the material card, leaving geometry, mesh and boundary conditions
byte-identical. `studies/build_brineK_decks.py` and `studies/build_nlgeom_decks.sh`
are the worked examples. An early separate-mesh attempt at the same comparison
reported scatter of random sign an order of magnitude larger than the real
effect; that scatter was the discretisation floor, not the physics.

Alternatively, freeze the packing and re-mesh it:

```bash
SPAX_SAVE_PACKING=pack/ python3 SpaX_Standalone.py one_row.csv m_033/   # freeze
SPAX_LOAD_PACKING=pack/ python3 SpaX_Standalone.py row_020.csv m_020/   # reuse
```

### Scatter convention

Every scatter quoted in this project is the **population** standard deviation
(`ddof=0`) over the packings of an ensemble. Note that pandas `.std()` and
`statistics.stdev` both default to the *sample* convention, which at five
packings is ~12% wider. Standard errors and Welch statistics correctly use
`ddof=1` — see `analysis/check_channel_isotropy.py` and
`analysis/compare_basetensor_sizes.py`.

---

## 5. Traps

Read this section before running anything long.

### Generation forks one mesher per core

Unless `SPAX_GEN_WORKERS` says otherwise. At ~7×10⁵ elements each worker holds
about 3 GB, so on a 16 GB machine set `SPAX_GEN_WORKERS=2`, and `=1` above ~10⁶
elements.

Killing the parent process leaves **orphaned mesher children** still holding
that memory. Kill them explicitly — and note their command line reads
`SpaX_GmshPeriodic`, not `SpaX_Standalone`, so the obvious `pkill` pattern
misses them.

### Periodic meshing retries are normal

At high inclusion counts the packer may need several attempts within its
`SPAX_MAX_RETRIES` budget; a large cell can use half of it. **A traceback in the
generation log is usually the retry mechanism working, not a failure.** Confirm
with the `GENERATION COMPLETE: N RVEs` line and the per-deck node counts before
concluding the decks are bad.

### Coarsening `L_mesh` is not a reliable speed-up

For porous RVEs the mesh-in-gap refinement triggers on a threshold proportional
to `lc_fine`, so a coarser base mesh flags *more* gaps and can add elements
back. Leave `SPAX_FORCE_HYBRID` at its default instead, for a clean solve saving
from the non-hybrid matrix.

### Applied strain comes from the deck, not a default

The decks prescribe a fixed *displacement*, so the engineering strain is
`Disp/L` — not a constant. `extract_elasticity_tensor` reads this from the deck
and raises rather than guessing. If you write your own extractor, do the same:
assuming a fixed strain is silently correct only at the one cell size where the
two coincide, and silently rescales the whole tensor everywhere else. Ratios are
immune, since a common factor cancels.

### A solve array skips on ODB *existence*

Not on a completed solve. A job killed at walltime leaves a truncated ODB, so a
naive resubmit skips exactly the jobs that failed. Identify failures by the
absence of an `Abaqus exit:` line in the log and delete those ODBs first.

### `module load abaqus` can fail in batch

On some cluster software stacks the module system is not initialised for
non-interactive shells, and the Abaqus modulefile cannot be parsed — every solve
dies immediately. Sourcing a snapshot of a working interactive environment
avoids it. In `hpc/`, `postprocess_basesweep.sh`, `postprocess_basetensor.sh`
and `postprocess_bt80.sh` use the robust form; the other five still carry the
plain `module load`. Copy from a robust one when writing a new post-processor.

### Clean up scratch

Delete decks and ODBs once the result CSV is pulled. Cluster scratch is shared
and usually purged on a timer.

---

## 6. Choosing load cases

For **isotropic effective moduli**, one uniaxial plus one shear is enough:
`Mode=Uniaxial Tension X` with `Mode2=Simple Shear S13`, adding `Kappa>0` if you
also want the bending response.

For **anisotropy** — the transverse isotropy of a vertical channel network, say —
set `full_tensor=Yes` to solve all six load cases and obtain
`E_x,E_y,E_z,G_xy,G_xz,G_yz` and the full 6×6 matrix.

A cell must be large enough to hold enough features for the directions being
compared to be equivalent. A cell holding only three to five channels is too
small for the two in-plane directions to be equivalent within one realisation,
and will show an apparent in-plane split that is a sampling artefact rather than
a material property. Check convergence by growing the cell and by adding
packings, not by trusting a single realisation.
