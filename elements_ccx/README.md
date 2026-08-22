# New CalculiX elements

Fortran sources added to `ccx` to give it the mixed displacement/pressure
element it does not ship — the one Abaqus calls `C3D4H`. Copy into
`<ccx_2.23>/src`, add them to `SCCXF` in `Makefile.inc`, apply the dispatcher
wiring in `../patches_ccx/`, and rebuild.

| File | Role |
|---|---|
| `u4mat.f` | the element operator: builds and bubble-condenses the mixed blocks, in closed form |
| `e_c3d_u4.f` | turns that operator into the element stiffness (`mafillsm.f` calls it) |
| `resultsmech_u4.f` | turns the *same* operator into stresses and internal forces (`resultsmech.f` calls it) |

Both consumers go through `u4mat`, deliberately. Splitting them is what leaves
`../patches_ccx/0002-bbar-mean-dilatation.patch` reporting an `equilibrium_gap`
of 1.0: patch the stiffness, recover the forces with the unpatched operator,
and the reaction cross-check — the only reference-free convergence evidence
this repository has — stops meaning anything.

## Using it

```
*USER ELEMENT,TYPE=U4,NODES=4,INTEGRATIONPOINTS=1,MAXDOF=4
*ELEMENT,TYPE=U4,ELSET=Brine
...
*STATIC,SOLVER=SPOOLES
```

`MAXDOF=4` is what makes it work at all: `allocation.f:1110` raises `mi(2)` to
it, which gives `mastruct.c` room for a fourth nodal DOF. DOFs are node-major —
1,2,3 displacement, 4 pressure.

**`SOLVER=SPOOLES` is not optional in this build.** The pressure block enters
negative, so the assembled system is symmetric *indefinite*, and
incomplete-Cholesky PCG — this repository's production solver, and the only one
whose memory fits a large cell — requires positive definiteness.

And SPOOLES is worse on a saddle point than its SPD numbers suggest, because
pivoting drives fill-in: the campaign cell at 218k equations passed **5.9 GB
without finishing**, against 1.28 GB for a 167k-equation SPD system in
`../calculix/README.md`. This is the binding constraint on using U4 at
production scale, not the element.

## Getting an efficient solver

Three routes, in increasing order of effort.

**1. Compile in PARDISO. ccx already has the interface and it is already set up
for exactly this matrix.** `src/pardiso.c` sets `mtype=-2` — MKL PARDISO's code
for *real symmetric indefinite* — and `src/pastix.c` exists beside it. Both are
guarded by `-DPARDISO` / `-DPASTIX`, and this build sets neither:

```
$ grep -oE '\-D[A-Z]+' Makefile | sort -u
-DARCH -DARPACK -DMATRIXSTORAGE -DNETWORKOUT -DSPOOLES
```

So no source needs writing — install MKL (it is not on this machine;
`ldconfig -p | grep mkl` is empty), add `-DPARDISO` and the MKL libraries to
the Makefile, rebuild, and `*STATIC,SOLVER=PARDISO` handles U4 with proper
indefinite ordering and threading. **This is the cheapest path by a wide
margin** and should be tried before anything below.

**2. Eliminate the pressure and recover a positive-definite system.** Because
`A_lb = 0` (see `u4mat.f`), the assembled system is exactly

```
[ A    Bᵀ ]        C' = C + B_b A_bb⁻¹ B_bᵀ
[ B   −C' ]
```

`C` is SPD for any finite `K` and the added term is PSD, so `C'` is invertible
and the Schur complement

```
A_eff = A + Bᵀ C'⁻¹ B
```

is symmetric **positive definite** — incomplete-Cholesky PCG applies, memory
returns to the ~370 MB regime, and there are no extra global DOFs. Lumping `C'`
to its row sums makes `C'⁻¹` diagonal and the assembly direct.

The catch is sparsity. `B` couples a pressure node to the displacement nodes of
its own element, so `Bᵀ C'⁻¹ B` couples displacement nodes that share a
*pressure node* — the node patch, wider than the element stencil. `mastruct.c`
has to build that larger pattern. That is the same structural change
`../calculix/README.md` identifies for nodal-averaged B-bar, and the two turn
out to be the same piece of work.

**3. Keep the mixed form and precondition it properly** (block-diagonal or
Uzawa). Most code, least reuse of what ccx has. Not recommended while route 1
is uncompiled.

## Closed-form integration, and why A_lb vanishes

`u4mat.f` uses no quadrature. Every integral is a monomial in the barycentric
coordinates over a straight tet, where
`∫ L₁^a L₂^b L₃^c L₄^d dV = 6V·a!b!c!d!/(a+b+c+d+3)!`, so with `gᵢ = ∇Lᵢ`
constant and `Σᵢ gᵢ = 0`:

| integral | value |
|---|---|
| `∫ Lᵢ dV` | `V/4` |
| `∫ Lᵢ Lⱼ dV` | `V/10` (i=j), `V/20` (i≠j) |
| `∫ Lᵢ ∂b/∂x_d dV` | `−(256V/840)·gᵢ[d]` |
| `∫ ∂b/∂x_c ∂b/∂x_d dV` | `(256²V/15120)·Σᵢ gᵢ[c]gᵢ[d]` |

**The linear–bubble block `A_lb` is exactly zero.** The bubble vanishes on every
face, so `∫∇b dV = ∮ b n dS = 0`, and every linear–bubble deviatoric term
carries a factor of it. Checked against the earlier 15-point implementation:
`max|A_lb| / max|A_ll| = 7e-16`, and the two agree on `A_ll` to 13 digits.

That collapses the condensation to one line — the bubble's *only* effect is a
stabilisation added to the pressure block:

```
A_ll' = A_ll ,   B_l' = B_l ,   C' = C + B_b A_bb⁻¹ B_bᵀ
```

which is the standard characterisation of MINI as P1/P1 plus a parameter-free
stabilisation, and it is what makes route 2 above possible.

Two practical consequences. The element declares **`INTEGRATIONPOINTS=1`**,
which matters more than it sounds: ccx sizes `sti`, `eme`, `xstiff` and `stx`
at `mi(1)` for *every* element in the model, so a 15-point brine element taxed
the whole mesh — **1.74 GB against 0.12 GB** on a 345k-element cell. And the
single output point carries the *exact* element volume average, because the
bubble contributes nothing to the mean strain; that removes a real trap, since
ccx's `.dat` reader collapses integration points by arithmetic mean, which
equals the volume average only for equal quadrature weights — true for C3D4 and
C3D10, false for the 15-point tet rule.

## Why MINI (P1⊕bubble/P1) and not the P1/P0 that Abaqus documents

`C3D4H` is documented as a linear tet with *constant* pressure, and the obvious
implementation is to make that pressure element-internal and statically
condense it, which needs no new global DOF and no solver change. That is the
recipe in the note this work was checked against, and **it produces an element
bit-for-bit identical to plain `C3D4`.**

The reason is short. A linear tet has one integration point, so the strain is
constant over the element and an element-constant pressure represents the
divergence *exactly*. Condensing it is then algebraically the mean-dilatation
B-bar operator, and B-bar on a one-point tet is the identity. This is not a
prediction: `../patches_ccx/0002-bbar-mean-dilatation.patch` implements exactly
that operator and measures `E_x` on an undrained layered cell as bit-identical
with it on and off.

So a faithful P1/P0 clone cannot close any gap. Making a linear tet
locking-free requires enriching it:

* **P1/P1 unenriched** violates the inf-sup condition — checkerboard pressure.
* **MINI (P1⊕bubble/P1)** restores inf-sup parameter-free. The bubble is
  element-internal and condenses locally; the pressure stays nodal, continuous
  and global. Implemented here.
* **P1/P1 with Brezzi–Pitkäranta stabilisation** would also work and is less
  code, but introduces a mesh-dependent tuning constant — unwanted when the
  whole purpose is validating against Abaqus.

The cost of the choice is exactly the thing the condensation route was avoiding:
a global pressure DOF, hence an indefinite system, hence SPOOLES.

## What is verified, and what is not

**Verified: the patch test passes exactly, at brine incompressibility.** See
`tests/`. A uniform strain state prescribed on every node of one tet has a
constant-stress exact solution, and U4 reproduces it to roundoff at all 15
integration points:

| | σxx | σyy = σzz | shear |
|---|---|---|---|
| E = 9.43e9, ν = 0.33, unit tet — analytic | 1.397192e+07 | 6.881690e+06 | 0 |
| U4 | **1.397192E+07** | **6.881690E+06** | ≤ 2e-10 |
| K = 2.2 GPa, G = 0.44 MPa, ν = 0.4999, distorted tet — analytic | 2.200587e+06 | 2.199707e+06 | 0 |
| U4 | **2.200587E+06** | **2.199707E+06** | ≤ 7e-12 |

The formulation is the standard mixed one, `σ = 2μ dev(ε) + p I` with
`div u − p/K = 0`, assembled as `[[A, Bᵀ],[B, −C]]`. λ never appears, which is
the point — it is the term that blows up as ν → ½.

**Not verified:**

* **Inf-sup stability.** A patch test is passed by plenty of unstable elements;
  it establishes consistency, not stability. The bubble is there to supply
  stability, and showing it needs a genuinely constrained problem.
* **Anything at RVE scale.** No layered cell has been run with U4, so nothing
  here yet says whether it closes the 8–12 % gap against Abaqus `C3D4H` that
  `../calculix/README.md` measures.
* **Periodic pressure coupling** — see below. On the layered decks this has to
  be solved before a U4 run means anything.

One bug worth recording because the obvious test misses it: `shape4tet` returns
at `iflag=2` *before* applying the inverse Jacobian, so `shp(1:3,i)` are
derivatives with respect to the local coordinates. Only `iflag=3` gives global
ones. On a unit axis-aligned tet the two coincide, so the first patch test
passed with the bug in place; the distorted tet in `tests/u4brine.inp` is what
exposes it.

## Known gap for the layered decks

The periodic `*EQUATION` constraints the generator writes cover DOFs 1–3 only.
In the spherical decks the brine is interior and that is fine. In the **layered**
decks the brine slab spans the cell and reaches the periodic faces, so its
pressure field needs periodic coupling too. Running U4 on those cells without
adding pressure equations across the face pairs would leave the pressure
unconstrained where it matters most.

# STATUS: U4 fails the Abaqus validation. Do not use it yet.

Verified: consistency (patch test exact at ν = 0.4999 on a distorted tet) and
absence of *shear*-regime locking (cantilever, C3D4 collapses 36×, U4 holds).

**Not correct on a layered RVE.** Against the campaign cell `LMESH_m0p0240`,
using the drained twin that is C3D4 in both codes:

| | R = E_x(und)/E_x(drn) | vs Abaqus |
|---|---|---|
| Abaqus C3D4H / C3D4 | 2.4701 (seed spread 0.61 %) | — |
| CalculiX C3D4 / C3D4 | 2.5094 | +1.59 % |
| **CalculiX U4 / C3D4** | **1.2751** | **−48.4 %** |

`equilibrium_gap` was 1.9e-07, so this is a clean solve of a wrong model.

**The diagnosis, and one wrong turn on the way to it.** The stabilisation is
439× the compressibility term in a single matrix entry, which looks like the
answer and is not: it annihilates uniform pressure *exactly* (row-sum ratio
0.99999999999999534), which is precisely what MINI should do. The element is
structurally right.

What is wrong is magnitude. The bubble penalises pressure *gradients* with a
coefficient ~h²/μ, derived for Stokes where μ is O(1). Relative to the physical
compressibility ~h³/K the ratio is

```
S/C  ~  0.05 K/mu        -- independent of h
```

and the brine has K/μ = 5000. The consequence, measured on the small cell:

| | E_x (across layers) |
|---|---|
| C3D4 undrained, K = 2.2 GPa | 6.349e9 |
| C3D10 undrained, K = 2.2 GPa | 4.738e9 |
| **U4 undrained, K = 2.2 GPa** | **3.600e9** |
| C3D4 **drained**, K = 2.2 MPa | 3.446e9 |

**U4's undrained answer is within 4.5 % of a genuinely drained cell.** A factor
of a thousand in bulk modulus contributes essentially nothing.

**Why the earlier tests missed it, which is the lesson worth keeping.** Both
passing tests are structurally blind to this:

* the **patch test** prescribes a uniform field, so the pressure is uniform,
  `B_bᵀp = 0`, the bubble amplitude is identically zero, and the stabilisation
  is not exercised at all;
* the **cantilever** is deviatoric-dominated, so the bulk modulus barely enters.

Consistency and shear behaviour were verified; the bulk response under a
*varying* pressure field was not, and that is the only regime the confined
brine occupies. A test that passes is not evidence for behaviour it cannot see.

**The remedy under test.** `CCX_U4_STAB=CAPPED` scales the stabilisation by
`θ = tr(C)/(tr(C)+tr(S))`, so it can never exceed the physical compressibility
it perturbs: `θ → 1` when the material is compressible (recovering textbook
MINI exactly), and `S_eff → C` when it is not. Parameter-free, and it leans on
the pressure mass term being nonzero — true here, because the brine is
ν = 0.49993 rather than exactly ½. Default remains unscaled MINI.

This is a deliberate departure from the textbook element and is **not yet
validated**. It must be measured against Abaqus C3D4H before use, exactly as
the unscaled version was.

# CAPPED also fails, for the opposite reason. Use C3D10.

`CCX_U4_STAB=CAPPED` was meant to stop the stabilisation swamping the physics.
It does — and it breaks stability instead.

| element | L_mesh | R | vs converged (1.9897) |
|---|---|---|---|
| Abaqus C3D4H | 0.0240 | 2.4701 | +24.1 % |
| Abaqus C3D4H | 0.0120 | 2.3786 | +19.5 % |
| CalculiX C3D4 | 0.0240 | 2.5094 | +26.1 % |
| CalculiX C3D4 | 0.0120 | 2.5852 | +29.9 % |
| **CalculiX C3D10** | **0.0120** | **2.1249** | **+6.8 %** |
| CalculiX U4 MINI | 0.0240 | 1.2751 | −35.9 % |
| CalculiX U4 CAPPED | 0.0240 | 2.2627 | +13.7 % |
| CalculiX U4 CAPPED | 0.0120 | 1.7020 | −14.5 % |

CAPPED passes *through* the right answer between the two meshes and keeps
softening. The coarse-mesh +13.7 % was a crossing, not agreement.

## The pressure field says why — and the first metric for it was wrong

The obvious checkerboard measure, RMS edge-to-edge jump over RMS field
variation, gives 0.991 for MINI and 1.281 for CAPPED, which looks conclusive
and is not. Controls on the same mesh: `u_y`, a physical displacement
fluctuation, scores **1.316** — higher than either pressure. On an unstructured
tet mesh that ratio is dominated by edge length, not oscillation.

The measure that does discriminate is each node's deviation from the mean of
its own neighbours, relative to the field variation:

| | u_x | u_y | u_z | **pressure** |
|---|---|---|---|---|
| MINI | 0.367 | 0.794 | 0.780 | **0.754** |
| CAPPED | 0.614 | 0.825 | 0.761 | **1.034** |

MINI's pressure (0.754) sits *inside* the range of the physical fluctuation
fields — no checkerboard, exactly as inf-sup stability promises. CAPPED's
(1.034) exceeds every displacement control on the same mesh and is 37 % above
MINI's. **CAPPED is under-stabilised**, which is why it is too soft and why it
gets worse under refinement.

So the two failures are complementary and neither is fixable by tuning:

* **MINI** — pressure stable, but the Stokes-scaled stabilisation destroys the
  brine's bulk stiffness (R = 1.28, within 4.5 % of a drained cell).
* **CAPPED** — bulk stiffness restored, stability lost.

## What to use instead

**Plain C3D10, and nothing else needs building.** At `L_mesh` = 0.0120 it is
+6.8 % from the converged answer, against Abaqus's own C3D4H at +19.5 % and
CalculiX C3D4 at +29.9 % on the same mesh. No new element, no MKL, no mixed
formulation — order 2 simply does not lock enough to matter here.

The cost is 8× the equations (601k → 4.9M on this cell) and it is not free, but
it is correct, and it is the recommendation until someone builds a properly
inf-sup stable element scaled for K/μ ~ 5000. `pressure_check.py` and the
neighbour-deviation measure are the acceptance test any such element must pass,
alongside the refinement sequence.

# Next: nodal-averaged B-bar, and why it escapes the trap

## Why U4 cannot be rescued

`S = B_b A_bb⁻¹ B_bᵀ` is *simultaneously* the stabilisation and the spurious
compliance. `CAPPED` scaled S down; `STIFFB` stiffened `A_bb`, which also
shrinks S. Two different-looking fixes, one identical lever:

| variant | u_x | u_y | u_z | pressure | R @0.0240 |
|---|---|---|---|---|---|
| MINI | 0.367 | 0.794 | 0.780 | **0.754** OK | 1.2751 |
| CAPPED | 0.614 | 0.825 | 0.761 | **1.034** above controls | 2.2627 |
| STIFFB | 0.653 | 0.834 | 0.781 | **1.061** above controls | 2.3241 |

Any pressure stabilisation for a P1/P1-type space scales as h²/μ on dimensional
grounds, while the physical compressibility scales as h²/K. At K/μ = 5000 the
stabilisation exceeds the physics by ~K/μ **regardless of mesh**. Shrink it to
protect the bulk stiffness and you fall below inf-sup. No window exists at this
contrast — a property of the linear tet, not of this implementation, and why
Abaqus's own C3D4H is only +19.5 % at `L_mesh` = 0.0120.

## Why nodal-averaged B-bar is different

Stability comes from **patch averaging**, a geometric construction, not a
1/μ-scaled penalty, so the scaling argument does not apply. It carries **no
pressure DOF**: a pure displacement method, so the matrix stays symmetric
**positive definite** and incomplete-Cholesky PCG works. PARDISO becomes an
option rather than a requirement.

```
V_a = sum over elements at node a of V_e/4            (nodal volume)
θ_a = (1/V_a) sum over those elements of (V_e/4) div(u)|_e
K   = K_dev (element-wise, as now) + K_vol (built from θ)
```

## Cost, measured on `LMESH_m0p0120` (1 211 410 elements, 214 539 nodes)

| | equations | nnz |
|---|---|---|
| C3D4 | 601k | 26.5 M |
| **nodal B-bar** | **601k** | **126.5 M** (4.8× C3D4) |
| C3D10 | 4.9M | 343 M |

2.7× fewer nonzeros than C3D10 and 8× fewer equations. Elements per node 22.6;
1-ring 14.7 mean / **27 max**; 2-ring 70.1 mean / 132 max.

## Implementation route — no `mastruct.c` change needed

The objection to nodal B-bar is that `K_vol` couples a whole patch, widening
the stencil past the element graph. Not needed here: `*USER ELEMENT` accepts up
to **255 nodes** (`userelements.f:83`) and `mastruct.c:137` reads the node count
from the element label, so a patch *is* expressible as an element. Our worst
patch is 27 nodes.

Two pure-displacement user elements, both `MAXDOF=3`:

| type | nodes | role |
|---|---|---|
| `U5` | 4 | linear tet, **deviatoric only** — C3D4 minus its volumetric term |
| `U6` | the 1-ring of one node (≤27 here) | that node's volumetric patch stiffness |

One `U6` per mesh node, one `U5` per tet. Every piece of plumbing already
exists and is exercised: label-driven `nope`, `mafillsm` dispatch,
`resultsmech_u` recovery, and the `printoutelem.f` volume fix. A generator
computes the 1-rings and writes the `U6` connectivity, as `u4ify.py` writes the
periodic pressure equations.

Acceptance tests unchanged — they caught every failure so far: patch test
exact, refinement tracking Abaqus toward R = 1.9897, and the *volumetric strain*
field checked for oscillation against the displacement controls on the same
mesh, the way `pressure_check.py` checks pressure.
