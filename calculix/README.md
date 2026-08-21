# Solving the SpaX decks with CalculiX

`SpaX_CalculiX.py` (top level) translates the generator's Abaqus `.inp` into a
CalculiX deck, runs `ccx`, and reads the results back. The post-processor then
produces the same `results.csv` it produces from an ODB — with **plain
`python3` and no Abaqus licence anywhere in the chain**.

```bash
python3 SpaX_Standalone.py  params.csv out/      # unchanged
python3 SpaX_CalculiX.py    all        out/      # translate + solve
python3 SpaX_PostProcess.py params.csv out/ results.csv
```

Converted decks are written beside the originals as `Job-<id>-<mode>-ccx.inp`,
so an Abaqus deck and its CalculiX twin coexist in one directory. The
post-processor prefers an `.odb` when both are present, and records which
solver produced each row in a new `solver` column. `SPAX_SOLVER=calculix`
reverses that preference, which is how to read both out of one directory.

## What this folder is

Verification runs, not part of the toolkit.

| Script | Question it answers |
|---|---|
| `validate_ccx.sh` | Does CalculiX reproduce the stored Abaqus results? |
| `compare_ccx.py` | Row-by-row relative difference between two results tables. |
| `hybrid_locking_test.sh` | What does CalculiX's lack of hybrid elements cost? |
| `hybrid_locking_extreme.sh` | …and at the softest brine the repo describes? |

## Validation: the homogeneous cube

`params/rve_homog_qxy.csv` is the right reference because it has no inclusions:
`E_eff`, `nu_eff` and `G_eff` are then closed-form and mesh-independent, so a
disagreement there is a defect rather than discretisation.

Against the stored `results/results_homog_qxy.csv` (Abaqus, `SPAX_MESH_ORDER=2`):

| Quantity | Agreement |
|---|---|
| `E_eff`, `nu_eff`, `G_eff` | **3e-6 %** — eight significant figures |
| `D_classical` | 3e-6 % (derived from the above) |
| `D_rve`, `E_bending` | CalculiX 1.0–1.4 % higher |
| `l` | 13–49 % higher |

The first-order agreement is conclusive: the periodic-BC translation, the
volume averaging and the modulus fit are all correct.

**The bending gap is not attributable to the solver.** Re-meshing the same cell
at `L_mesh` 0.028 / 0.033 / 0.040 moves `D_rve` by 1.27 % *within CalculiX
alone*, while `E_eff`/`nu_eff`/`G_eff` do not move at all. The gap against a
table built on a mesh this tree cannot reproduce is therefore the same size as
the mesh noise, and the two cannot be separated from the stored table.
Separating them needs one Abaqus run on a deck generated here — the `.inp` the
converter reads is the same file Abaqus would solve, so that comparison is
exact when someone has a licence.

`l` amplifies the gap because it comes from the small difference
`D_rve - D_classical`; for the homogeneous cube the true length scale is zero
and the whole excess (0.5–3.3 % in Abaqus, 1.9–4.3 % in CalculiX) is numerical.
`analyze lengthscale` divides it out through the homogeneous baseline
`f_quad(L)`, so this calibrates away **provided the baseline and the porous
cells come from the same solver**. Do not mix solvers inside one length-scale
analysis.

## Validation: a real campaign deck

`params/rve_gas.csv` is the campaign's dominant configuration — L=0.50 at
L_mesh=0.033, the size shared by 859 of the deck rows — and its Abaqus results
are stored in `results/results_gas.csv`. `validate_gas.sh` reruns it through
CalculiX on the shipped default solver. Six cells, 46k–187k elements, up to
~1.16M equations.

Unlike the homogeneous cube this cell is randomly packed, and the campaign's
generation seed is not recorded in the tree, so a fresh run cannot reproduce
its exact microstructure. The comparison therefore bundles packing, meshing and
solver, and the achieved phase fractions are what say how much of each:

| run | total soft fraction, Abaqus / CalculiX | E_x [GPa] | Δ | E_z [GPa] | Δ | ν_x |
|---|---|---|---|---|---|---|
| GAS_v00 | 0.0202 / 0.0198 | 9.157 / 9.109 | −0.53 % | 9.179 / 9.157 | −0.24 % | 0.3303 / 0.3308 |
| GAS_v02 | 0.0397 / 0.0395 | 8.840 / 8.734 | −1.20 % | 8.890 / 8.845 | −0.51 % | 0.3303 / 0.3310 |
| GAS_v04 | 0.0593 / 0.0592 | 8.533 / 8.378 | −1.82 % | 8.626 / 8.530 | −1.12 % | 0.3304 / 0.3310 |
| GAS_v06 | 0.0793 / 0.0789 | 8.220 / 8.055 | −2.00 % | 8.353 / 8.222 | −1.57 % | 0.3303 / 0.3309 |
| GAS_v08 | 0.0990 / 0.0985 | 7.919 / 7.705 | −2.71 % | 8.065 / 7.950 | −1.42 % | 0.3304 / 0.3309 |
| GAS_v10 | 0.1183 / 0.1185 | 7.638 / 7.366 | −3.57 % | 7.805 / 7.632 | −2.22 % | 0.3304 / 0.3310 |

The cells hold the same amount of each phase — total soft fraction agrees to
0.1–2 % — but not in the same arrangement. Every solve converged
(`equilibrium_gap` 1.9e-7 to 1.7e-6).

ν agrees to 0.2 % and is flat across the series, while E is 0.5–3.6 % low with
the gap **growing monotonically with void content**. A solver discrepancy would
not be porosity-dependent; that shape belongs to geometry or discretisation.

Two things are worth stating about the reference before reading anything into
that. The current `results_gas.csv` is the **2026-08-07 re-run**, made after the
2026-08-06 mesher-geometry corrections (`98226e6`, `4aeb689`, `2f689f8`), so it
is a modern Abaqus reference and not the legacy table — the numbers it replaced
were 1–2 % different and carried no achieved-phase columns at all. But **the
element count behind it is not recorded anywhere in the tree.** (`results_scf.csv`
lists ~20 700 matrix elements for a run named GAS_v00, but that file dates from
2026-07-20, before those same fixes, so it describes a different mesh and cannot
be used to characterise this one.)

### It was quadratic against linear

The table above compares the wrong things. It was generated with
`SPAX_MESH_ORDER=2`, and the reference was solved with **linear** tets.

The campaign's own submitter records this. Every first-order campaign in the
table in `hpc/spax_submit.sh` carries mesh order 1 — `weibull`,
`weibull_layer`, `nlgeom_layer`, `layercol`, all `-utx` — and only the bending
and torsion campaigns carry 2; `hpc/generate_array.sh` documents the variable
as "`SPAX_MESH_ORDER` (2 for bending)". `rve_gas.csv` is a first-order deck
(`Kappa=0`, `Mode` utx, `Mode2` utz), so it was solved on C3D4.

And the size was already measured here: `hybrid_locking_test.sh` found order 1
reads **+4.03 %** stiffer than order 2 on one frozen geometry. That is the
direction and magnitude of the "gap", and its growth with void content is what
linear-tet over-stiffening does when there is more geometry to resolve.

Rerunning the same deck at order 1 (`validate_gas_order1.sh`):

| run | E_x: Abaqus | ccx order 2 | ccx order 1 |
|---|---|---|---|
| GAS_v00 | 9.157e9 | −0.53 % | **+0.01 %** |
| GAS_v02 | 8.840e9 | −1.20 % | **−0.06 %** |
| GAS_v04 | 8.533e9 | −1.82 % | **−0.24 %** |
| GAS_v06 | 8.220e9 | −2.00 % | **−0.11 %** |
| GAS_v08 | 7.919e9 | −2.71 % | **+0.05 %** |
| GAS_v10 | 7.638e9 | −3.57 % | **−0.54 %** |

E_z behaves the same way: −0.24…−2.22 % at order 2, +0.44…−0.06 % at order 1.
ν agrees to 0.04 %.

The deficit collapses to **±0.54 %**, and — the part that matters — the
monotonic trend with void content is gone. What remains scatters around zero.

`packing_scatter.sh` measures what that residual should be: GAS_v10 re-packed
at four seeds, everything else fixed, gives a population s.d. of **0.29 %** and
a full spread of **0.77 %**. The residual sits inside its own packing noise.

The same run also shows the order effect is systematic rather than scatter: at
order 2 the four packings land at −3.33 %, −3.90 %, −3.15 % and −3.61 % against
Abaqus — a 3.5 % offset with 0.77 % of spread around it. Re-packing moves the
answer by a third of a percent; changing the element order moves it by three
and a half.

(That scatter was measured at order 2 while the matched comparison is at order
1. The spread is a property of the geometry rather than the element, so it
should carry over, but it was not measured at order 1.)

**Nothing is left over for the solver.** Which is the same answer the
homogeneous cube gave at 6–8 significant figures, now confirmed on a real
microstructure at campaign scale.

The lesson is about the comparison, not the code: match the element order
before reading anything into a modulus difference. A 4 % discrepancy that grows
with porosity looked like physics and was a deck setting.

### The bug this comparison found

`phi_inclusion` initially came back at exactly **2×** the campaign's value while
the porosity printed beside it was right. ccx emits the EVOL block with every
element listed twice under one header; anything reading volumes through the
label→volume dict was immune — the volume-averaged stress, `V_solid`,
`porosity` and hence `E_eff` were all correct — so only the raw sum behind
`phi_inclusion` doubled. Total meshed volume now sums to L³ exactly on a
void-free cell.

The defect was in the reader, not the generator, and nothing about the moduli
would have revealed it. That is the argument for comparing achieved geometry
alongside results rather than results alone.

## Hybrid elements

CalculiX has none, and says so rather than guessing:

```
*ERROR reading *ELEMENT:
C3D4H    is an unknown element type
```

There is no user-defined route to one. `*USER ELEMENT` exists in `ccx` but is
the substructure/superelement interface, not a way to code a mixed
displacement/pressure formulation from the deck; adding a genuine hybrid tet
means changing ccx's Fortran element routines and rebuilding. So the converter
strips the `H`, and the phase the generator marked hybrid — the brine, at
`nu >= SPAX_HYBRID_NU` — runs on the plain displacement element that is
supposed to volumetrically lock as `nu -> 0.5`.

The CalculiX community thread on this
([discourse 1509](https://calculix.discourse.group/t/incompressible-material-models/1509))
confirms the limitation from the manual — *"Perfectly incompressible materials
require hybrid finite elements. CalculiX does not provide such elements"* — and
offers two remedies. **Neither transfers to this problem**, and it is worth
saying why, because both look applicable at first glance:

* **Reduced integration (`C3D20R`), the manual's advice for isochoric
  behaviour.** CalculiX has reduced integration only for hexahedra — `C3D8R`
  and `C3D20R`. There is no reduced-integration tetrahedron (`C3D4`, `C3D10`,
  `C3D15` only), and the periodic mesher produces tets for these geometries.
  Not available.

* **Accept `nu ≈ 0.475`** (ccx's own default when a hyperelastic `D1` is zero).
  That advice is for `*HYPERELASTIC`, where incompressibility is a modelling
  choice. Here the brine is `*ELASTIC` and its `nu` is not a knob — it encodes
  a measured bulk modulus. Because `G` is only 0.44 MPa, `K` is violently
  sensitive to `nu` near 0.5:

  | ν | K | vs the specified 2.2 GPa |
  |---|---|---|
  | 0.49990 | 2.2 GPa | as specified |
  | 0.499 | 220 MPa | 0.10× |
  | 0.490 | 21.9 MPa | 0.010× |
  | 0.475 | 8.65 MPa | **0.0039×** |

  Capping at 0.475 cuts the brine's bulk modulus by **254×**, turning the
  undrained cell into something close to its drained twin. Since K = 2.2 GPa
  against the ice's 9.25 GPa is exactly what makes the confined brine a load
  path, that does not work around the numerics — it deletes the physics under
  test.

Which leaves the two routes below, and the cheap alternative of solving those
cells at order 2.

**Measured, it does not.** One frozen packing (`SPAX_SAVE_PACKING`), the same
geometry throughout, inclusion fraction 0.295:

| Case | `E_eff` order 1 → 2 | `G_eff` order 1 → 2 |
|---|---|---|
| brine, `nu` = 0.4900 | +3.62 % | +4.59 % |
| brine, `nu` = 0.49993 | +3.78 % | +4.82 % |
| **compressible twin, `nu` = 0.300** | **+4.03 %** | **+4.51 %** |

The twin has the *same* inclusion Young's modulus (`E = 9KG/(3K+G) = 1.320e8`)
and differs only in compressibility. Its order-1 → order-2 change is the same
as the brine's, so that change is ordinary linear-tet over-stiffness, not
locking. Had the missing hybrid formulation been biting, the near-incompressible
rows would have moved much further than the twin.

Refining the order-2 brine mesh (`L_mesh` 0.050 → 0.035) moves `E_eff` by
0.22 % at `nu` = 0.490 and 0.23 % at `nu` = 0.49993 — converged, and equally so
at both.

The reason is physical, and worth stating because it means the result is
structural rather than luck. The brine's ν = 0.49 comes from a large K/G ratio,
not from being volumetrically stiff:

| | K | G | E | ν |
|---|---|---|---|---|
| ice | 9.25 GPa | 3.55 GPa | 9.43 GPa | 0.330 |
| brine | 2.20 GPa | 0.044 GPa | 0.132 GPa | 0.490 |

The brine is **four times more compressible in bulk than the ice around it**,
and 80× softer in shear. Volumetric locking bites when a phase has to represent
isochoric deformation *while carrying load* — an incompressible phase acting as
a constraint on its surroundings. This one does the opposite: it is a soft,
nearly-void inclusion that the matrix flows around. It neither carries the load
nor constrains the matrix volumetrically, so the element's ability to represent
isochoric motion inside it barely enters the effective stiffness. The
incompressibility *is* being represented — `nu_eff` reads 0.336 with the brine
against 0.307 with the compressible twin — it simply is not what governs.

That also says exactly when the conclusion would stop holding: a soft phase
that **percolates** and carries load, or geometric **confinement**. Both are
present in the layered decks — see below.

### The layered cells: where it does show up

`layered_incompressible.sh` tests the case the spherical result does not cover.
The layered decks carry drained/undrained pairs:

| | K | G | ν | Abaqus element |
|---|---|---|---|---|
| drained | 2.2 MPa | 0.44 MPa | 0.406 | C3D4 (below the 0.45 threshold) |
| **undrained** | **2.2 GPa** | 0.44 MPa | **0.49993** | **C3D4H — hybrid** |

The undrained cell has everything the spherical one lacked: the brine is a
cell-spanning slab, confined between ice plates, with the ice bridges as
constrictions — and K = 2.2 GPa is only 4× below the ice, so it genuinely
carries load. Pressure transmission through that confined brine is the
mechanism the layered closure rests on. The campaigns also run it on linear
tets (`5196ff1`, "Run the layered sweeps with linear elements to fit the
face-constraint limit").

One geometry per element order; the drained twin made by rewriting the single
inclusion elastic card on the **identical mesh**, so drainage is the only
difference:

| order 1 vs order 2 | E_x (across layers) | E_z (in-plane) |
|---|---|---|
| drained, ν = 0.406 | **+9.08 %** | +0.75 % |
| undrained, ν = 0.49993 | **+11.19 %** | +1.06 % |
| **attributable to incompressibility** | **+2.11 points** | +0.32 points |
| *spherical brine vs its twin, for contrast* | *−0.41 points (none)* | — |

**Read the control before the headline.** The raw +11.19 % looks like locking
and mostly is not: the drained cell, which Abaqus also meshes without hybrid
elements, loses +9.08 % on the same geometry. That is ordinary linear-tet
stiffness amplified by thin slabs and narrow bridges, and **Abaqus carries it
too** — it is not a CalculiX-versus-Abaqus difference at all.

What *is* CalculiX-specific is the excess: **+2.1 points** in E_x, the part
Abaqus recovers with C3D4H and CalculiX cannot. Small next to the geometry
term, but real, and qualitatively unlike the spherical case where the
incompressibility penalty was zero (−0.41 points, i.e. the brine tracked its
compressible twin exactly). It is also directional — 2.1 points across the
layers against 0.3 in-plane — which is the signature of the confined phase
having to deform near-isochorically under pressure.

So the expected CalculiX-vs-Abaqus discrepancy on an undrained layered cell at
order 1 is of order **2 % in the across-layer modulus**, not 11 %.

**Exposure: 204 undrained cells across 14 decks** (`rve_bracket_*`,
`rve_layerb`, `rve_layercol*`, `rve_layermesh`, `rve_layerskel`,
`rve_weibull_layer`). The drained-only decks (`rve_eringen_layer`,
`rve_nlgeom_layer`, `rve_torsion_layer`) sit at ν = 0.406, below the hybrid
threshold, so Abaqus used no hybrid there either and CalculiX matches it.

### Do not rely on the numbers above: that cell was under-resolved

The table was measured on a cell whose brine slab is barely one element thick.
The generator refines to `lc_fine = 0.4 x L_mesh` near an inclusion, so:

| | slab thickness | `lc_fine` | elements through the slab |
|---|---|---|---|
| the test above, L_mesh=0.020 | 0.0075 | 0.0080 | **0.9** |
| the test above, L_mesh=0.011 | 0.0075 | 0.0044 | **1.7** |
| campaign `rve_layerb` | 0.0100 | 0.0020 | **5.0** |
| campaign `rve_layercol_p060` | 0.0075 | 0.0022 | **3.4** |

With one element across it, the confined brine layer — the whole feature under
test — is not represented at all. Both the +9.08 % geometric term and the
+2.11-point incompressibility excess were measured on a cell that cannot carry
the mechanism they are about, and the +7.40 % order-2 drift under refinement is
at least as likely to be the unresolved slab as locking. **These numbers do not
transfer to the campaign's cells.**

`layered_incompressible.sh` now takes `NSLABS` and `SLABVOF` so the local
resolution can be held at the campaign's while the cell is shrunk to something
this machine can solve — `L=0.15 NSLABS=2 LMESH=0.0045` gives 4.2 elements
through the same 0.0075 slab. That is the measurement to trust.

Two further caveats that stand either way. This infers the hybrid benefit from
an order-1-vs-order-2 comparison rather than measuring C3D4 against C3D4H
directly — that needs one Abaqus run on an undrained layered deck, and the deck
the converter reads is the same file Abaqus would solve. And it is one
morphology: `n_bridges=2`, `bridge_fraction=0.29`. Narrower bridges mean tighter
constrictions, so the penalty should be re-checked at the low-`bridge_fraction`
end of `rve_bracket_bridge.csv` before any number is relied on.

**One limitation to keep in view.** This validated `E_eff`, `G_eff` and
`nu_eff` — homogenised quantities, which average over the inclusion interiors.
Volumetric locking distorts the *local* stress field well before it moves a
homogenised modulus. The per-slice principal-stress and SCF extraction is not
ported to CalculiX (see *Not ported*), and if it ever is, this measurement does
not carry over to it: the fields inside a near-incompressible phase would need
checking on their own terms.

**So: use `SPAX_MESH_ORDER=2` for CalculiX and the missing hybrid element costs
nothing measurable.** Order 1 is a separate 4 % error that Abaqus has too. This
was measured at one inclusion fraction on one packing; a cell where the soft
phase percolates and carries load could behave differently, and the twin
comparison above is the cheap way to re-check it.

### If it ever does need building

Two routes, and they are not the same size of job.

**A true mixed displacement/pressure element** — an extra pressure unknown per
element entering the global system — is the invasive one. `lakon`, the element
type string, is consumed in 38 source files; a new global DOF also touches the
matrix structure (`mastruct.c`), the solvers and the results recovery. This is
not a weekend.

**B-bar (mean dilatation) is the tractable one, and it is the standard fix.**
`e_c3d.f` assembles the isotropic stiffness directly from the Lamé constants,
with the volumetric and deviatoric parts already written as separate terms:

```fortran
s(ii1,jj1) = s(ii1,jj1) + (al*w(1,1) +
&              um*(2.d0*w(1,1)+w(2,2)+w(3,3)))*weight
```

Locking lives entirely in the `al` (λ) term, which blows up as ν → 0.5. B-bar
replaces the discrete divergence operator in that term by its element average —
accumulate it through the integration loop, then add the volumetric
contribution once afterwards instead of at every point. No extra DOF, no change
to the matrix structure, no solver changes; roughly 50–100 lines in `e_c3d.f`
plus registering a type name, and a matching change in `resultsmech.f` so the
recovered volumetric stress uses the same operator.

**What B-bar would and would not buy here.** For the spherical decks, nothing:
the 4 % order-1 error is not volumetric — the compressible twin shows the same
4 % — so B-bar would not let those cells drop from quadratic to linear
elements, and there is no incompressibility penalty to recover.

For the **undrained layered** decks it would buy the +2.1 points measured
above, which is the part Abaqus gets from C3D4H. That is the one place in this
repository where the missing hybrid element costs something measurable. It is
worth weighing against the alternative, which costs nothing to try: solve those
cells at `SPAX_MESH_ORDER=2`, where the penalty is by construction absent.
The campaigns moved to linear elements for a reason (`5196ff1`, the
face-constraint limit), so that alternative may not be free in practice — but
it should be priced before writing Fortran.

## Two ccx incompatibilities worth knowing

**Node sets in `*EQUATION`.** CalculiX resolves set names everywhere else but
wants a node *number* inside `*EQUATION`. Untranslated it stops with one
`*ERROR reading *EQUATION` per periodic-BC equation — loud, at least.

**The 20-character field limit.** CalculiX copies each comma-separated field
into a fixed 20-character buffer and **truncates** anything longer. Python's
`repr` overruns it routinely: a coordinate of `3.8163916471489756E-17` is 22
characters and truncates to `3.8163916471489756E-`, which fails to parse.

That failure is the lucky case. Truncation only errors when it lands
mid-exponent — a field that truncates to a still-valid number is accepted
silently at the wrong value. An `*Equation` coefficient of
`3.469446951953614e-18` (21 characters) truncates to `3.469446951953614e-1`, a
periodicity constraint off by seventeen orders of magnitude in a deck that
solves without complaint. The converter therefore re-renders **every** numeric
field at 12 significant digits, not only the ones observed to break.

## Scale, and why the default is the iterative solver

A stock `ccx` is usually linked against SPOOLES only — check with
`ldd $(which ccx)`. SPOOLES is a *direct* solver, and at production size that
is not a trade-off but a wall:

| equations | SPOOLES | ITERATIVE CHOLESKY |
|---|---|---|
| 167 306 | 23.8 s, 1276 MB | 20.6 s, **370 MB**, 338 iters |
| 488 390 | 149.8 s, 4402 MB | 64.8 s, **1072 MB**, 306 iters |
| ~1.4 M | still factorising at 19 GB after ~25 min, stopped | — |

Two things matter in that table. The iterative memory is **linear** in the
model (2.90× for 2.92× the equations) while SPOOLES' is super-linear and its
time worse still (6.3× for 2.92×). And the **iteration count did not grow** —
338 → 306 — so the incomplete-Cholesky preconditioner is holding up as the cell
grows, which is what makes the method usable at all here.

So `SPAX_CCX_SOLVER` now defaults to `ITERATIVE CHOLESKY`. Set it to `SPOOLES`
for a small-cell reference, or `DEFAULT` to write no `SOLVER` parameter and let
ccx choose.

### The tolerance, which is not optional

Out of the box the iterative solvers return the **wrong answer, quietly**. The
convergence test in `pcgsolver.c` stops once the residual max-norm falls to
0.5 % of the mean load (`c1 = 0.005`, a local constant no input deck reaches),
which leaves `E_eff` 0.15 % below the direct answer on the same deck.

`patches/0001-iterative-tolerance.patch` exposes it as `CCX_ITER_TOL`:

| `CCX_ITER_TOL` | iters | wall | peak | `E_eff` error | equilibrium gap |
|---|---|---|---|---|---|
| 0.005 (stock) | 70 | 14.0 s | 370 MB | **−0.1477 %** | 1.8e-3 |
| 1e-3 | 94 | 14.6 s | 370 MB | +0.0133 % | 3.3e-6 |
| 1e-5 | 338 | 20.9 s | 370 MB | −0.0000 % | 1.8e-7 |
| SPOOLES | — | 23.2 s | 1275 MB | reference | 7.4e-8 |

`1e-5` is the default `SpaX_CalculiX.solve()` exports. Memory is flat in the
tolerance — tightening costs iterations, not storage.

Beware the decoy: the `eps` argument threaded down from `preiter.c` looks like
the tolerance and is documented as "required accuracy", but is never read — PCG
and CG overwrite it with the final residual and return it. Patching `eps`
changes nothing.

An unpatched ccx silently ignores `CCX_ITER_TOL`, so `solve()` checks the log
for the tolerance the patched binary echoes and warns plainly when it is
absent.

### The check that survives into production

`equilibrium_gap`, recorded per RVE in `results.csv`, is the relative
disagreement between the volume-averaged stress and the reference-point
reaction — two independent measurements of the same macroscopic stress, linked
by nothing but the model being right and the system being solved. It tracked
the true error across the whole sweep above (1.8e-3 when the solve was bad,
1.8e-7 when good) and it needs **no reference solve**, which is the point: on a
cell too large to solve directly even once, it is the only convergence evidence
available. Treat anything above ~1e-4 as a solve not to be believed.

## Not ported

`extract_principals` — the per-slice principal-stress and SCF field CSVs behind
`results_scf.csv` and the stress-field figures. `SpaX_CalculiX.extract_averages`
supplies the volume-averaged stress and strain tensors that the 6x6 elasticity
route needs, but not the slice fields.
