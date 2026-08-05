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

## 2. What SpaX models

A periodic representative volume element of a **two-phase microstructure**: a
continuous matrix holding a population of discrete inclusions. The inclusions
are packed by random sequential addition, meshed periodically so that opposite
faces correspond node-for-node, and loaded through reference-point-driven
periodic constraints. The macroscopic stress is recovered as a volume average,
with the Hill–Mandel condition satisfied identically by the periodic
constraints.

That covers a broad class of materials:

| Class | How it is set up |
|---|---|
| **Particle-reinforced composites** | Stiff inclusions in a compliant matrix — the classic case. Set `E_sphere_inclusion` above `E_matrix`. |
| **Porous media and foams** | Voids, via `VoF_void_sphere`. True voids are left unmeshed rather than filled with a near-zero stiffness. |
| **Compliant-inclusion composites** | A second phase softer than the matrix, via its own bulk and shear moduli. |
| **Anisotropic microstructures** | Elongated inclusions aligned by `Growth_Direction` and sharpened by `Growth_Concentration`, and optionally a percolating network of parallel channels. |

Voids and a soft second phase can coexist in one cell; they are distinct
phases, not two names for the same thing.

### Modelling a compliant phase

Give a compliant inclusion phase **its own bulk and shear moduli** rather than
treating it as a cavity. A nearly incompressible fluid-filled inclusion is well
represented by a soft solid with a realistic `K` and a very small `G`, giving
`ν → 0.5` from below. A true cavity constraint is both physically wrong for a
phase that carries pressure and numerically fragile in a periodic cell, where
the constraint interacts badly with the face-tying equations.

The distinction matters most when the inclusion phase is far more compressible
than intuition suggests: it is the *ratio* of the phase moduli to the matrix,
not the label "fluid", that sets the effective response.

### What controls the effective anisotropy

Two mechanisms are easy to confuse. Inclusion **shape and orientation** — the
aspect ratio of aligned ellipsoids — produce an anisotropy that saturates
quickly. Inclusion **connectivity** — whether the second phase percolates —
produces a much stronger one. Where both are present, connectivity dominates,
so a study sweeping shape at fixed volume fraction will underestimate the
anisotropy of a microstructure whose second phase is actually connected.

### Units are free

Classical effective moduli are scale-invariant, so the cell edge and the
inclusion radii are in arbitrary model units. What matters is their ratio, and
the elements-per-inclusion the mesh size implies. Map the cell onto a physical
size only when reporting.

---

## 3. Where things live

| What | Where |
|---|---|
| The toolkit | `SpaX_Standalone.py`, `SpaX_GmshPeriodic.py`, `SpaX_PostProcess.py` |
| Deck builders | `studies/` — one script per campaign, writes a parameter CSV |
| Parameter CSVs | `params/rve_*.csv` |
| Cluster submit/post | `hpc/submit_*.sh`, `hpc/postprocess_*.sh` |
| Result tables | `results/results_*.csv` |
| Stiffness matrices | `tensors/` (earlier campaigns at the top level; the manuscript ensembles in `tensors/column/`, `tensors/basetensor_seeds/`, `tensors/bt80/`) |
| Analysers and figures | `analysis/` (plain `python3`) |
| Renderers | `viz/render_rve.py`, `viz/odb_to_vtk.py`, `viz/render_stress_field.py` |

`out_*/` directories hold generated Abaqus decks. They are derived and
untracked — several GB, and regenerable from `params/`. Exclude them from any
archive unless you specifically need them.

---

## 4. A campaign end to end

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

## 5. Reproducibility

### Seeding is per row index

`SPAX_SEED` fixes the packing, but the seed is resolved **per row index**. The
same seed with a different number of rows gives different packings. Two
campaigns are only comparable packing-for-packing if their CSVs have the same
shape.

### Split generation must seed from the *global* row index

When generation is spread across parallel array tasks, each task receives a
slice of the CSV. If the per-task seed is derived from the row's index *within
its own slice*, every task starts from the same index — and a campaign whose
rows differ only in replicate number comes back with near-identical packings
that look like an implausibly well-converged ensemble.

The symptom is a median packing coefficient of variation one to two orders of
magnitude below what serial generation of the same deck gives. Derive the seed
from the task's **global** row index instead; `hpc/generate_array.sh` has the
worked form. If an ensemble's scatter looks too good, check this before
believing it.

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

### Report ensembles, not realisations

A single packing can sit several standard deviations off the ensemble mean, and
a production number taken from one realisation will move when the campaign is
repeated. Quote the ensemble mean with its scatter, and treat a difference
smaller than the replicate scatter as no difference at all.

---

## 6. Traps

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

### A bending size effect is not by itself a length scale

Solving the cell in bending at several sizes and regressing the apparent modulus
on `1/L²` looks like a direct measurement of an intrinsic length. It is not,
because the measurement carries a systematic error of the same shape: imposing
plate-like bending kinematics on a *cubic* cell is itself size-dependent.

A geometrically identical **inclusion-free** cell — no microstructure, and
therefore no possible intrinsic length — shows a size dependence of the same
order as the microstructured one. Any bending size-effect measurement must
therefore be referred to such a `φ = 0` control and the control divided out.

Two consequences. **Do not solve larger bending cells hoping the length scale
will converge** — the cost grows as `(L/d)³` and the trend being chased is not
microstructural. And do not read the *sign* of the raw slope as a verdict: a
positive slope is diagnostic of couple-stress stiffening, but a negative one is
ambiguous between nonlocal softening, first-order dilution and the extraction
bias itself.

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
and usually purged on a timer. `RUNBOOK.md` §5 has the procedure and, more
importantly, the order the steps have to happen in.

---

## 7. Choosing load cases

For **isotropic effective moduli**, one uniaxial plus one shear is enough:
`Mode=Uniaxial Tension X` with `Mode2=Simple Shear S13`, adding `Kappa>0` if you
also want the bending response.

For **anisotropy** — the transverse isotropy of an aligned inclusion population
or a parallel channel network, say — set `full_tensor=Yes` to solve all six load
cases and obtain `E_x,E_y,E_z,G_xy,G_xz,G_yz` and the full 6×6 matrix.

A cell must be large enough to hold enough features for the directions being
compared to be equivalent. A cell holding only three to five channels is too
small for the two in-plane directions to be equivalent within one realisation,
and will show an apparent in-plane split that is a sampling artefact rather than
a material property. Check convergence by growing the cell and by adding
packings, not by trusting a single realisation.
