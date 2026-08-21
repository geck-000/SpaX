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

The reason is physical: the brine is a *soft* inclusion, ~70x more compliant
than the ice. Locking hurts when an incompressible phase has to carry the
deformation; here the matrix governs. The incompressibility is still being
represented — `nu_eff` comes out 0.336 with the brine against 0.307 with the
compressible twin — it simply is not sensitive to the element technology.

**So: use `SPAX_MESH_ORDER=2` for CalculiX and the missing hybrid element costs
nothing measurable.** Order 1 is a separate 4 % error that Abaqus has too. This
was measured at one inclusion fraction on one packing; a cell where the soft
phase percolates and carries load could behave differently, and the twin
comparison above is the cheap way to re-check it.

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

## Scale

A stock `ccx` is often linked against SPOOLES only — check with
`ldd $(which ccx)` for PARDISO or PaStiX. SPOOLES is a direct solver whose
memory grows much faster than the model, and the production cells in this
campaign run to millions of elements. `SPAX_CCX_SOLVER=ITERATIVE CHOLESKY`
switches to a built-in iterative solver that needs far less memory, but this is
exactly the sort of ill-conditioned system (70x phase contrast, a
near-incompressible phase) that iterative solvers struggle on. **The largest
cell verified here is 35k elements.** Anything approaching campaign scale needs
a ccx built against PaStiX or PARDISO, and a check that it converges.

## Not ported

`extract_principals` — the per-slice principal-stress and SCF field CSVs behind
`results_scf.csv` and the stress-field figures. `SpaX_CalculiX.extract_averages`
supplies the volume-averaged stress and strain tensors that the 6x6 elasticity
route needs, but not the slice fields.
