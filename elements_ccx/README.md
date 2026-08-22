# New CalculiX elements

Fortran sources added to `ccx` to give it the mixed displacement/pressure
element it does not ship — the one Abaqus calls `C3D4H`. Copy into
`<ccx_2.23>/src`, add them to `SCCXF` in `Makefile.inc`, apply the dispatcher
wiring in `../patches_ccx/`, and rebuild.

| File | Role |
|---|---|
| `u4mat.f` | the element operator: builds and bubble-condenses the mixed blocks |
| `e_c3d_u4.f` | turns that operator into the element stiffness (`mafillsm.f` calls it) |
| `resultsmech_u4.f` | turns the *same* operator into stresses and internal forces (`resultsmech.f` calls it) |

Both consumers go through `u4mat`, deliberately. Splitting them is what leaves
`../patches_ccx/0002-bbar-mean-dilatation.patch` reporting an `equilibrium_gap`
of 1.0: patch the stiffness, recover the forces with the unpatched operator,
and the reaction cross-check — the only reference-free convergence evidence
this repository has — stops meaning anything.

## Using it

```
*USER ELEMENT,TYPE=U4,NODES=4,INTEGRATIONPOINTS=15,MAXDOF=4
*ELEMENT,TYPE=U4,ELSET=Brine
...
*STATIC,SOLVER=SPOOLES
```

`MAXDOF=4` is what makes it work at all: `allocation.f:1110` raises `mi(2)` to
it, which gives `mastruct.c` room for a fourth nodal DOF. DOFs are node-major —
1,2,3 displacement, 4 pressure.

**`SOLVER=SPOOLES` is not optional.** The pressure block enters negative, so
the assembled system is symmetric *indefinite*. Incomplete-Cholesky PCG — this
repository's production solver, and the only one whose memory fits a large cell
— requires positive definiteness and cannot be used. That caps the reachable
model size at whatever SPOOLES can factorise, which `../calculix/README.md`
measures at roughly 500k equations in 4.4 GB.

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
