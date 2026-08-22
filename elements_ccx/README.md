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

# U5 + nodal-averaged B-bar

The route that replaces U4. Two pieces, and only one of them is an element:

| file | role |
|---|---|
| `e_c3d_u5.f`, `resultsmech_u5.f` | `U5` — linear tet carrying the **deviatoric** stiffness only |
| `nodalbbar.py` | builds the volumetric half out of existing ccx features |

```
K = K_dev   per tet (U5)
  + K_vol   per node: theta_a tied to the surrounding displacements by an
            *EQUATION, given energy by a grounded SPRING1
```

The bulk modulus never enters an element-local operator, so it cannot lock.
Stability comes from averaging over the node patch — geometry, with no 1/μ
scaling anywhere — which is exactly what U4 could not have: there the
stabilisation and the spurious compliance were the same term.

And there is **no pressure unknown**, so the matrix stays symmetric positive
definite: `SOLVER=ITERATIVE CHOLESKY` works, PARDISO optional rather than
mandatory.

## Why no patch element was needed

`U6` was going to be a patch element — one per node, spanning its 1-ring.
`*USER SECTION` can attach constants to user elements, but only the *same*
constants to a whole set, so per-patch data would have meant ~36 000
single-element sections. Carrying θ as one extra DOF on a dummy node instead
needs no new element at all, and ccx's MPC cascade handles the widened stencil
that `mastruct.c` would otherwise have to build.

Two details that make it practical:

* **Every spring is identical.** With `t_a = √(K·V_a)·θ_a` the energy is exactly
  `½ t_a²`, so one `*SPRING` of unit stiffness covers every patch and `V_a`
  lives in the `*EQUATION` coefficients. Otherwise each node needs its own
  element set.
* **The θ node's unused DOFs must be constrained.** It carries three DOFs and
  only DOF 1 is used; left free the matrix is singular and SPOOLES *stops
  silently* after "Factoring the system of equations" — no error, no
  `Job finished`, an empty `.dat`.

## Stress output uses a different operator from the stiffness, deliberately

`resultsmech_u5` reports `σ = 2μ dev(ε) + K div(u)` from the **element's** own
divergence, while the stiffness contains no volumetric term at all. That is
exact for the homogenisation because θ_a is the V-weighted mean of the
surrounding element divergences and each element feeds four nodes:

```
sum_a V_a theta_a  =  sum_e V_e div(u)|_e
```

so the volume-*integrated* volumetric stress is identical computed either way.
Internal forces stay deviatoric, matching the stiffness; the springs supply the
rest of the reaction, so `equilibrium_gap` stays a real check rather than the
~1.0 it reads under the B-bar patch.

## Verification

| test | U5 + nodal B-bar | for contrast |
|---|---|---|
| patch test, ν = 0.33, unit tet | **exact** 1.397192E+07 | — |
| patch test, ν = 0.4999, distorted | **exact** 2.200587E+06 | — |
| locking sweep, vs EB at ν = 0.49999 | **0.5632**, settling | C3D4 0.0155 |
| volumetric strain smoothness | **0.456**, controls 0.412/0.620/0.608 | U4 pressure: MINI 0.754, CAPPED 1.034, STIFFB 1.061 |

θ is smoother than two of the three displacement controls, which is what a
patch average should be. This is the measure that failed both U4 repairs.

**Still open:** the layered RVE against Abaqus at both mesh levels. CAPPED
cleared the coarse point and then shot past the converged answer under
refinement, so a single good R proves nothing — the 0.0120 point is the test.
Also untested: the phase interface. Only the brine is split into U5 + patches;
the ice stays ordinary C3D4, and patches at interface nodes cover brine
elements only, so no 1000× modulus contrast is smeared. That is a design
assertion, not yet a measurement, and a bad interface treatment would show up
as an offset in R rather than in any of the four tests above.

## STATUS: the formulation is right, the ccx delivery is not

The RVE match looked perfect and is not to be believed. On `LMESH_m0p0240`,
`R = 2.4677` against Abaqus's 2.4701 — 0.10 %, inside the 0.61 % seed spread —
but `equilibrium_gap = 6.93e-01` and `E_z` came out 3.6 % *stiffer* than C3D4,
which replacing volumetric stiffness with an average can never do.

A confined-compression block, where the answer is closed form
(`M = K + 4G/3`), isolates it:

| mesh | free displacement DOFs | reaction / exact |
|---|---|---|
| 1×1×1 (6 tets, 8 nodes) | none | **1.0000** |
| 2×2×2 (48 tets, 27 nodes) | few | 1.5485 |
| 4×4×4 (384 tets, 125 nodes) | many | 3.0515 |

The error appears **only when displacement DOFs are free** and grows with the
interior-node count.

**Every ingredient is exact in isolation:**

| checked | result |
|---|---|
| U5 deviatoric stiffness, confined block | reaction ratio **1.0000** |
| nodal volumes | Σ V_a = mesh volume exactly |
| divergence operator vs an analytic field | θ_a = 1e-4 = div(u) at every node |
| `SPRING1` convention | k = 1, u = 1 → RF = 1.0 |
| **the whole operator, assembled in Python** | **reaction ratio 1.0000** |

That last row is the one that settles it. The same `K_dev + Σ_a K V_a b⊗b`
assembled directly reproduces the closed form exactly; routed through ccx's
`*EQUATION` + `SPRING1` it does not. **Nodal B-bar is sound here — the delivery
mechanism is at fault.**

It is also provably a bug rather than a limitation: nodal averaging satisfies
`Σ_a V_a θ_a² ≤ Σ_e V_e θ_e²` by Cauchy–Schwarz, so this method can only ever be
*softer* than element-wise volumetric stiffness. Coming out 3× stiffer is
impossible for a correct assembly.

**Prime suspect**: the θ DOF carries the spring *and* is the dependent DOF of
its own `*EQUATION`, so its stiffness entry has both indices MPC-dependent.
`mafillsm.f` has a distinct branch for that case and it is the one path none of
the passing tests exercises.

**Routes out**, in order of preference:

1. Give the volumetric term to a real element so no MPC is involved — the
   original `U6` patch element. `*USER SECTION` blocks the obvious per-patch
   data route, but ccx's substructure path (`matrix2userelem.f`,
   `writesubmatrix.f`) reads an externally assembled stiffness matrix as a user
   element, and the Python assembly above already produces exactly that matrix.
2. Keep the MPC but move the spring off the dependent DOF.
3. Fix the both-DOFs-dependent branch in `mafillsm.f`, if that is genuinely
   where it is.

## What the verification ladder did and did not catch

Patch test, locking sweep and volumetric-strain smoothness (0.456, inside
controls) **all passed** on the broken assembly. The patch test cannot see it
because prescribing every displacement means the springs never carry load; the
other two are insensitive to an overall volumetric scale error.

Only **reaction against a closed-form modulus on a confined block** found it,
and that test should run before any RVE, not after. It is now
`tests/oedometer_m2.inp`.

### Root cause: overlapping MPCs, not the method

Every candidate was eliminated by test, not by argument:

| hypothesis | test | verdict |
|---|---|---|
| U5 stiffness wrong | confined block, deviatoric only | exact, 1.0000 |
| nodal volumes wrong | Σ V_a vs mesh volume | exact |
| divergence operator wrong | θ_a vs analytic field | exact, 1e-4 everywhere |
| emitted equations wrong | θ reconstructed from the deck | exact, 1e-4 everywhere |
| `SPRING1` convention wrong | k=1, u=1 | RF = 1.0 |
| MPC + spring condensation broken | 2 nodes, θ = 3u, k = 1, F = 1 | u = 1/9 **exact** |
| equations over-constrain u | equations with springs removed | deviatoric exactly, 1.0000 |
| **the operator itself** | **same operator assembled in Python** | **exact, 1.0000** |

The discriminator:

| patches | Python | ccx |
|---|---|---|
| 1 | 7.069748e+01 | **7.069748e+01** |
| 27 (all) | 2.200587e+05 | 3.407531e+05 |

**One patch is exact; the full set is 1.55× too stiff.** ccx handles a single
`*EQUATION` + `SPRING1` pair perfectly and degrades once many equations share
the same displacement DOFs — which is intrinsic here, because adjacent node
patches overlap by construction. Each mesh DOF appears as an independent term
in up to ~24 equations.

This cannot be fixed in the generator. The volumetric term has to reach the
matrix without going through ccx's MPC machinery.

### The fix: make the patch a real element after all

Give `K_vol^a = K V_a b⊗b` to a `U6` element spanning the patch's 1-ring —
the original design, which the `*USER SECTION` per-set-constants limit pushed
me away from. Two ways to supply the per-patch data:

1. **Compute it inside the element.** Make node 1 of the `U6` connectivity the
   patch centre and have `e_c3d_u6.f` build a node→element map once (cached,
   built in a serial phase to stay safe under the OpenMP assembly), then form
   `b` from the surrounding tets exactly as `nodalbbar.py` does. No deck
   plumbing at all.
2. **Feed the matrices in.** ccx's substructure path (`matrix2userelem.f`,
   `writesubmatrix.f`) reads an externally assembled stiffness matrix as a user
   element, and the Python assembly above already produces exactly these
   matrices — but one file per patch makes this impractical at 36 000 patches.

Route 1 is the one to build. `*USER ELEMENT` already permits 255 nodes
(worst patch here is 27), `mastruct.c` reads the count from the label, and U5
plus the dispatcher and volume fixes are all in place — the remaining work is
one element routine.

## Per-material patches: the design, and an unresolved deck-format failure

Applying U5+U6 to the brine alone leaves interface brine nodes with one-sided
patches, and in a slab 2-3 elements thick most brine nodes *are* interface
nodes -- so the averaging has little to average over. Measured: R goes
2.3843 -> 2.5602 under refinement, i.e. back to plain C3D4, while the *same
element* gains 8.9x -> 50.9x over C3D4 on a homogeneous cantilever as the mesh
refines. The element is fine; treating one phase is not.

**The ice needs patches too.** Locking follows isochoric *deformation*, not the
material's own ν: ice beside a near-incompressible slab is forced to deform
almost isochorically, so C3D4 ice locks despite ν = 0.33.

**But patches must not span the interface.** The stress recovery relies on
`Σ_a V_a θ_a = Σ_e V_e div(u)|_e`, which holds for ONE K. Across a 1000×
modulus jump it fails, and the reported stress stops matching the transmitted
force. So: one patch per (node, material). `u6patch` filters contributing tets
by the patch element's own material; the generator emits one `*USER ELEMENT`
type per (material, ring size).

### The deck-format failure: found and fixed

**Cause: ccx's built-in element dispatch claims type names containing
digits.** `elements.f:361` has `elseif(label(4:4).eq.'4') then nope=4`, `:358`
has `'10' -> nope=10`, `:320` has `'20' -> nope=20`. A user element named
`U614` therefore gets `nope=4` from the digit rule *before* the `*USER ELEMENT`
lookup runs, so ccx reads 4 nodes instead of 18 and treats the continuation
line as a fresh element card.

Only 3 of ~36 000 continuation-needing patches reported an error, because
`"element N is already defined"` fires only when the misread line's first field
collides with a real element id. The rest were corrupted silently.

The fix is to name types from **letters only** (`U6AA`…`U6ZZ`, 676 available
against 58 needed), which avoids every digit-keyed rule.

Ruled out on the way, all by test rather than argument: a missing or extra
trailing comma on continued cards; more than 16 fields overflowing
`textpart(16)`; interleaving `*USER ELEMENT` with `*ELEMENT` across a
keyword-chain boundary; a cap on the number of user element types (`nuel_` is
dynamic); the deck being internally malformed (a validator found 0 mismatches
over all 93 904 elements); and continuation being broken in general (a minimal
18-node deck whose continuation line deliberately starts with a colliding
element id parses correctly).

A second real bug found alongside it: patch element ids started at 1e8, so
ccx's `ne = max(ne, id)` sized every element array for 100 million elements.
Ids now continue from the deck's real maximum.

### The equilibrium gap: NOT the constraints -- the internal forces vanish

Corrects an earlier conclusion in this file. With **no constraints of any kind**
-- the real RVE mesh, uniaxial displacement on the two x-faces, lateral faces
traction-free, so the average-stress theorem gives `<σxx>·A = F` exactly:

| | `<σxx>·A` | reaction | ratio |
|---|---|---|---|
| C3D4 | 5.8196e3 | 5.8179e3 | **1.000291** |
| U5+U6 | 5.1345e3 | **−1.2985e2** | **−39.5** |

The U5+U6 reaction is essentially zero while its volume-averaged stress is
sane. Internal forces are not reaching the constrained nodes. This reproduces
without periodic constraints, so everything ruled out earlier about MPCs,
three-term equations and 3D chaining was ruled out correctly but was never the
issue -- the failure needs only the real mesh and free DOFs.

What is now known:

| condition | result |
|---|---|
| structured block, free DOFs, any constraint form tested | **exact** |
| real RVE mesh, every node prescribed | stress **exact**, but no free DOFs to exercise the stiffness |
| real RVE mesh, free DOFs, no constraints | **reaction ≈ 0** |

So the discriminator is the *mesh*, not the boundary conditions. The structured
block has patches built from a regular 6-tets-per-hex subdivision; the RVE mesh
is unstructured, with rings from 4 to 32 nodes and a wide spread of element
sizes. A reaction near zero with a plausible stress means `fn` is being
accumulated for only a small part of the model.

### ROOT CAUSE: patches larger than 20 nodes overflow ccx's element matrix

`mafillsm.f:69` and `e_c3d_u.f:68` declare **`s(60,60)`** — 60 DOFs, i.e. 20
nodes at 3 DOF/node. U6 patch rings reach **32 nodes on the small cell and 37
on the campaign cell**, so `e_c3d_u6` writes `s(ii,jj)` with indices up to 96
into a 60×60 array. Fortran does not bounds-check: the writes land in other
columns or past the array, and `mafillsm` reads the same out-of-bounds entries
when assembling.

It accounts for every observation:

| observation | explanation |
|---|---|
| structured block exact under every BC | its rings are ≤15 nodes = 45 DOFs, inside the limit |
| real RVE mesh fails with free DOFs | rings to 32 nodes = 96 DOFs, past the limit |
| stress output exact under prescribed strain | stress comes from `u6patch`, never from `s` |
| stiffness and recovery return **bit-identical** `va`, `kva`, `b` | both call `u6patch`, which is correct — the corruption is downstream, in `s` |
| reaction ≈ 0 with a plausible `<σxx>` | the assembled stiffness is corrupted, so the solution is too soft, while U5 still reports `K_e·div(u)` at full stiffness |

Two things this clears, both of which looked like suspects and were not:

* **U5-only giving RF ≈ 1.1e-11 is physically correct**, not a missing
  assembly: deviatoric-only means `K = 0`, so `E = 9KG/(3K+G) = 0` and a
  laterally-free bar has genuinely zero axial stiffness.
* **`va = 2.8e-10` on a 4-node ring is correct too.** Its centre node sits in 7
  tets but only one of that material, so per-material patching is working as
  designed.

**The fix** is to enlarge the element matrix along the whole user-element path —
`s`, `sm` and `ff` in `mafillsm.f`, `e_c3d_u.f` and the U5/U6 routines — to at
least `3 × max_ring`. Capping the ring instead is not an option: dropping nodes
from a patch breaks `Σ_a V_a θ_a = Σ_e V_e div(u)|_e`, and the rank-one patch
stiffness cannot be split across two elements.

A cheap interim check, before touching ccx: regenerate with only patches of ≤20
nodes given to U6 and the rest left as plain C3D4. That is not the right
method, but if the equilibrium gap collapses it confirms the diagnosis end to
end.

### Two earlier claims to correct

* The shared-patch (both-phase) run gave `equilibrium_gap` 2.5e-3 at 0.0240 and
  2.6e-1 at 0.0120, and I attributed that to the K-identity breaking across the
  interface. **Those decks were also malformed** by the same continuation issue,
  so the attribution is unproven -- the algebra stands, the demonstration does
  not.
* `R = 2.1568` (+8.4 % from converged) from that run was reported as a large
  improvement before its equilibrium gap was checked. It is not trustworthy.
