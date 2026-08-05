# SpaX — user guide

SpaX builds periodic representative volume elements (RVEs) of a two-phase
microstructure, meshes them so that opposite faces correspond node-for-node,
solves them under periodic boundary conditions, and returns the effective
(homogenised) stiffness. It was written for sea ice — an ice matrix holding
brine pockets, gas voids and, near the warm base, percolating vertical brine
channels — but nothing in the core is specific to ice.

This file is the *what and why*. For the *how to run it*, see
[`RUNBOOK.md`](RUNBOOK.md).

---

## 1. What the toolkit does

A single run takes one row of a parameter CSV — one RVE — and carries it
through four stages:

1. **Pack.** Inclusions are placed by random sequential addition into a cube of
   edge `L`, honouring a minimum separation so that the mesher can resolve the
   gap between neighbours. Packing is periodic: an inclusion crossing a face
   re-enters through the opposite one.
2. **Mesh.** Gmsh meshes the cube under a strict periodicity constraint, so
   every boundary node has an exact partner on the opposite face. This is what
   makes the periodic boundary conditions exact rather than approximate.
3. **Solve.** Abaqus/Standard applies a macroscopic strain through
   reference-point-driven constraints tying each pair of opposite faces, one
   load case per independent macro-strain component.
4. **Homogenise.** The macroscopic stress is recovered as the volume average
   over the cell, and the effective stiffness follows. Energetic consistency
   between the scales is enforced by the Hill–Mandel condition, which the
   periodic constraints satisfy identically.

Only stage 3 needs Abaqus. Packing, meshing and the whole analysis stage are
plain Python.

### Two phases, not a cavity

Brine is modelled as a **soft solid** with its own bulk and shear moduli
(`K = 2.2 GPa`, `G = 0.44 MPa`, i.e. `ν ≈ 0.4999`), not as a fluid cavity and
not as a void. Brine at `K = 2.2 GPa` is compressible on the scale that matters
here, and a genuinely incompressible cavity is both wrong and numerically
fragile in a periodic cell. Gas inclusions *are* true voids and are left
unmeshed — they are a distinct phase from the brine.

---

## 2. Repository layout

| Path | Contents |
|---|---|
| `SpaX_Standalone.py` | Generation: packing, meshing, deck writing. |
| `SpaX_GmshPeriodic.py` | The periodic mesher. Forked per RVE during generation. |
| `SpaX_PostProcess.py` | ODB extraction and homogenisation. |
| `studies/` | Deck builders — one script per campaign, each writing a parameter CSV. |
| `params/` | The generated parameter CSVs (`rve_*.csv`), one row per RVE. |
| `hpc/` | Slurm batch scripts. Site-specific archival examples, not a dependency. |
| `results/` | Homogenisation output (`results_*.csv`), curves and figures. |
| `analysis/` | Analyzers and field extractors. Plain Python, no Abaqus licence. |
| `viz/` | RVE and field rendering. |
| `tensors/`, `post_coltensor/`, `post_basetensor_seeds/`, `post_bt80/` | Per-RVE 6×6 elasticity tensors, one CSV each. |
| `docs/` | This guide and the runbook. |

`out_*/` directories hold generated Abaqus decks. They are derived and
untracked — regenerable from `params/` — and should be excluded from any
archive.

---

## 3. The parameter CSV

One row per RVE. The columns fall into groups; the root `README.md` documents
every column individually, and the summary here is only to orient you.

- **Identity and geometry** — `run_id`, box edge `L`, target mesh size
  `L_mesh`. Elements-per-inclusion is what actually matters for accuracy, so
  hold `L_mesh` fixed when sweeping `L`.
- **Matrix material** — `E_matrix`, `nu_matrix`.
- **Inclusion population** — volume fractions for the soft phase and the voids,
  mean and spread of the radius, and `sphericity_avg`, which controls how far
  the ellipsoids depart from spheres.
- **Anisotropic growth** — `Growth_Direction` and `Growth_Concentration`,
  which align and sharpen the inclusion orientation distribution. Columnar ice
  grows with its long axis vertical, so this is `Z` throughout the ice work.
- **Channels** — an optional percolating network of vertical tubes, switched on
  by `generate_channels` with its own volume fraction target. Channels may be
  tilted off-axis.
- **Load cases** — which macro-strain components to solve, and whether to
  extract the full 6×6 tensor.

---

## 4. What the toolkit has been used to study

### 4.1 The depth-graded sea-ice column

The flagship application. A first-year ice sheet is sliced into ten horizontal
layers, each becoming one RVE that carries its own local temperature. That
temperature sets the brine volume fraction, which in turn sets the pocket
morphology, so a single depth coordinate drives the whole microstructure.

```
   air / snow    T ≈ -20 °C   ── cold, low brine, isolated near-spherical pockets
   ┌────────────────────────┐  z05   φb ≈ 0.022   sphericity 0.85   no channels
   │  z05                   │
   │  z15 … z45             │  interior, low salinity (C-shape minimum)
   │  z55, z65              │
   │  z75                   │  ── percolation threshold, φb ≈ 0.05 ──
   │  z85   channels on     │  z85   φb ≈ 0.068   sphericity 0.68   tubes + channels
   │  z95   channels on     │  z95   φb ≈ 0.150   sphericity 0.62   skeletal layer
   └────────────────────────┘
   ocean         T ≈ -1.8 °C  ── warm, brine-rich, connected vertical channels
```

The physical inputs and where they come from:

- **Temperature.** A linear conductive winter profile,
  `T(z) = T_top + (T_bot − T_top)·(z/H)`, from a cold surface to a base pinned
  at the seawater freezing point, ≈ −1.8 °C.
- **Salinity.** The classic first-year **C-shape**: elevated at top and bottom,
  with an interior minimum near 4.3 ppt. Bulk salinity, not brine salinity.
- **Brine volume fraction.** Frankenstein & Garner,
  `φb = S·(−49.185/T + 0.532)/1000`, valid for −0.5 ≥ T ≥ −22.9 °C. The sharp
  rise toward the warm base is the physical driver of the soft skeletal layer.
- **Gas.** Drained and entrapped air at 1–2 % by volume, slightly higher near
  the surface where bubbles are retained. Modelled as true voids.
- **Percolation.** Brine becomes connected above `φb ≈ 0.05` — the "rule of
  fives". The profile crosses it between the eighth and ninth slice, which is
  why channels switch on there. It is a microstructural transition, not a
  tuning knob.
- **Morphology.** Cold ice holds near-spherical pockets; toward the warm base
  they merge into vertical tubes and lamellae. The sphericity trend and the
  vertical orientation both come from X-ray micro-CT.
- **Matrix.** Ice stiffens as it cools, roughly 13 MPa/°C, so `E_matrix` is a
  function of the slice temperature rather than a constant.

**Units are free.** Classical effective moduli are scale-invariant, so the RVE
edge is in model units. Map it to a physical cell of a few millimetres — sub-mm
pockets, ~1 mm channels — when reporting.

The output is the depth profile `E(z)`, `G(z)` and the anisotropy
`E_z(z)/E_xy(z)`, which then feeds a laminated-plate model of the whole sheet.

### 4.2 Statistical representativeness and cell size

How large must the cell be, relative to the inclusions, for the homogenised
properties to be representative? The study holds everything fixed except the
box edge `L` and sweeps it across several sizes with independent packings at
each, so that the mean and the packing-to-packing scatter can be separated.

Two thresholds decide the answer together: the mean must have stopped moving,
and the scatter must have fallen below tolerance. The first-order moduli
satisfy both at a modest `L/d`.

**Report scatter as the population standard deviation (`ddof=0`).** This is a
project-wide convention. The `pandas` and `statistics` defaults are `ddof=1`
and mixing them has been a source of drift.

### 4.3 The bending size effect — a resolved question

An earlier line of work solved the cell in bending across a range of sizes,
looking for the intrinsic length scale `l` that a couple-stress continuum would
imply. The apparent `l` kept climbing with cell size and never plateaued, which
was initially read as a convergence problem to be solved with bigger cells.

**It was not.** Solving a geometrically identical *inclusion-free* cell — same
sizes, same kinematics, no microstructure at all, and therefore no possible
intrinsic length — produces a size dependence of the same order. The apparent
length scale is an artefact of imposing plate-like bending kinematics on a
cubic cell, not a property of the microstructure. Once the control is divided
out, the residual trend is null.

Two consequences for anyone using this toolkit:

- Do not run larger bending cells hoping the length scale will converge. It
  will not, and the cost grows as `(L/d)³`.
- Any bending size-effect measurement must be referred to a matched `φ = 0`
  control. A sign test on its own is not enough: a *positive* slope is
  diagnostic of couple-stress stiffening, but a *negative* one is ambiguous
  between nonlocal softening, first-order dilution and the extraction bias
  itself.

### 4.4 Other directions the toolkit supports

- Refining the percolation transition to locate the onset of mechanical
  anisotropy, which pairs naturally with fluid permeability.
- Sweeping pocket shape and orientation at fixed brine fraction, to isolate how
  the directly tomography-measured quantities control the anisotropy.
- First-year versus multi-year ice: C-shape against a desalinated profile, with
  raised gas content for drained summer ice.
- Seasonal sweeps, tracing how the soft base thickens as the sheet warms.

---

## 5. Output

`results_*.csv` carries one row per RVE with the directional moduli, Poisson
ratios, the achieved volume fractions (which differ slightly from the targets,
since packing is stochastic), and the anisotropy ratio. When the full tensor is
requested, the 6×6 stiffness is written per RVE into the `post_*/` directories.

A note on the achieved fractions: always analyse against what was actually
packed, not what was requested. The difference is small but systematic, and it
is the reason composition columns are carried through post-processing.

---

## 6. Conventions

- **Scatter is the population standard deviation**, `ddof=0`, everywhere.
- **Seeding is resolved per row index.** The same seed with a different number
  of rows gives different packings, so two campaigns are comparable
  packing-for-packing only if their CSVs have the same shape.
- **Ensemble means, not single realisations.** A single packing can sit several
  standard deviations off the mean; production numbers are ensemble means with
  the scatter quoted alongside.
