# F-barES-FEM-T4: formulation, small-strain reduction, and what it commits us to

Source: Y. Onishi, R. Iida, K. Amaya, *F-barES-FEM-T4 for viscoelastic problems*,
Int. J. Comput. Methods **15**(7) 1845003 (2018).  Equation numbers below are the
paper's.  `~/Downloads/onishi2019.pdf`; text dump kept alongside the run logs.

The method it builds on is Onishi, Iida & Amaya (2017), which is where the
recommended `c` per Poisson ratio actually lives.  We do not have that paper.
What the 2018 paper states about the envelope is quoted verbatim in section 5.

---

## 1. Notation

| symbol | meaning |
|---|---|
| `ᵉ(·)` | element quantity, constant over a T4 |
| `ⁿ(·)` | node quantity |
| `ᵸ(·)` | edge quantity, constant over an edge smoothing domain |
| `ᵉV`   | element volume (reference) |
| `ᵸV`   | edge smoothing volume, `Σ_{e∈ᵸE} ᵉV/6` |
| `ⁿV`   | node smoothing volume, `Σ_{e∈ⁿE} ᵉV/4` |
| `ᵸE`   | elements attached to edge `h`; `ⁿE` elements attached to node `n`; `ᵉN` the 4 nodes of `e` |

The `/6` and `/4` are the edge and node counts of a T4, so each element
distributes its whole volume across its 6 edges and its 4 nodes.  Every
smoothing operator below therefore has **row sums equal to 1**, which is what
makes the small-strain reduction in section 3 exact rather than approximate.

## 2. The finite-strain chain

**(a) Edge-smoothed deformation gradient — ES-FEM, applied once.** Eq. (1):

```
ᵸF̃ = Σ_{e∈ᵸE} ᵸw_e ᵉF ,        ᵸw_e = (ᵉV/6) / ᵸV
```

equivalently through smoothed shape-function gradients, eq. (3):
`ᵸF̃_ij = ᵸÑ_{P,j} x_{P:i}`.

**(b) Isovolumetric split.** Eqs. (4)-(5):

```
ᵸJ̃ = det(ᵸF̃) ,      ᵸF̃^iso = (ᵸJ̃)^(-1/3) ᵸF̃
```

**(c) c-time cyclic smoothing of J — NS-FEM, applied c times.** Eqs. (6)-(7):

```
ⁿJ̄ = (1/ⁿV) Σ_{e∈ⁿE} ᵉJ̃ (ᵉV/4)          (6)   element -> node
ᵉJ̄ = (1/4) Σ_{n∈ᵉN} ⁿJ̄                   (7)   node -> element
```

repeated `c` times, with `ᵉJ̄` fed back in as `ᵉJ̃` on the second and later
passes.  Then once, to the edges, eq. (8):

```
ᵸJ̄ = (1/ᵸV) Σ_{e∈ᵸE} ᵉJ̄ (ᵉV/6)           (8)   element -> edge
```

**(d) Volumetric gradient and recombination.** Eqs. (9), (11):

```
ᵸF^vol = (ᵸJ̄)^(1/3) I
ᵸF̄     = ᵸF^vol · ᵸF̃^iso = (ᵸJ̄ / ᵸJ̃)^(1/3) ᵸF̃
```

`det(ᵸF̄) = ᵸJ̄` by construction: the edge keeps ES-FEM's shape change and takes
its volume change from the cyclically smoothed field.

**(e) Stress.** Eqs. (12)-(13), Hencky strain `H̄` of `ᵸF̄`:

```
T^hyd = K tr(H̄) I ,     T^dev = 2G₀ (H̄^dev - Σᵢ gᵢ H^vᵢ)
```

For us the Prony terms are absent, so `T^dev = 2G H̄^dev`.  Two identities
matter: `tr(H̄) = ln det(ᵸF̄) = ln ᵸJ̄`, and `H̄^dev = H̃^dev` because `F̄` and `F̃`
differ only by a spherical factor.  **So the deviatoric stress never sees the
cyclic smoothing, and the hydrostatic stress never sees anything else.**

**(f) Internal force — and the one line that defines F-bar.** Eq. (17), with
the paper's own note:

> "Note that the stretching tensor in this equation, `D̃`, is **not** the
> deformation rate of `ᵸF̄` in Eq. (11) but that of `ᵸF̃` in Eq. (1) due to the
> adoption of the F-bar method."

Stress from the **modified** gradient, virtual work paired with the
**unmodified** one.  This is Petrov-Galerkin and the tangent, eq. (19), is
**non-symmetric**.  That is not an incidental detail — see section 4.

## 3. Small-strain reduction (what we actually implement)

Let `ε_e = sym ∇u|_e` and `θ_e = tr ε_e = div u|_e`, both constant per T4.
Collect the three smoothing operators as sparse matrices:

```
Q : elem -> node ,  Q[n,e] = (ᵉV/4)/ⁿV
P : node -> elem ,  P[e,n] = 1/4
E : elem -> edge ,  E[h,e] = (ᵉV/6)/ᵸV
A = P Q             one cycle, elem -> elem
```

All three have unit row sums, so on `ᵉJ ≈ 1 + θ_e` they act affinely and the
constant passes through untouched: the whole chain is **linear in θ**.

```
ᵸθ̃ = (E θ)_h                 edge strain trace        (ES-FEM, once)
ᵸθ̄ = (E A^c θ)_h             cyclically smoothed      (eqs. 6-8)
ᵸε̃ = (E ε)_h
ᵸε̄ = dev(ᵸε̃) + (1/3) ᵸθ̄ I    the F-bar strain, eq. (11)
```

Stress at the edge: `σ_h = 2G dev(ᵸε̃) + K ᵸθ̄ I`.  Virtual work paired with
`δᵸε̃` (section 2f).  With `B̃_h = Σ_e ᵸw_e B_e` the edge strain-displacement
matrix and `D_div` the element divergence operator:

```
K = Σ_h ᵸV B̃_hᵀ D_dev B̃_h            deviatoric: symmetric ES-FEM
  + K_bulk (E D_div)ᵀ W (E A^c D_div)   volumetric: NON-SYMMETRIC
```
with `W = diag(ᵸV)`.  Test space `E D_div`, trial space `E A^c D_div`.

At `c = 0` the two spaces coincide and the scheme collapses to plain selective
ES-FEM-T4, symmetric.  That is the only `c` for which symmetry is legitimate.

## 4. The Petrov-Galerkin form, and what it does not change

The obvious-looking Galerkin assembly `B̄ᵀ D B̄` is not what section 2f
specifies, and the difference is measurable.  What it is *not* is a change in
the constraint count: `K_vol = G_testᵀ W G_trial` has rank `min(rank G_test,
rank G_trial)` and the same right null space `ker(G_trial)` in either form, so
both admit exactly the same set of volumetric-free displacement modes.  The
counts, measured directly by `verify_fbar.py` V4 on a brine sphere with 1104
elements, 311 nodes and 1630 edges:

| c | rank `E` (test) | rank `E A^c` (trial) | constraints vs 3·n_node = 933 |
|---|---|---|---|
| 0 | 1011 | 1011 | ratio 0.92 — over-constrained, still locking territory |
| 1 | 1011 |  308 | ratio 3.03 — the Q1P0 / MINI count |
| 2 | 1011 |  308 | ratio 3.03 |

So the cyclic smoothing does exactly what it is designed to do: it relieves the
edge-volumetric constraint from over-constrained at `c = 0` to the textbook
count at `c ≥ 1`.  An ideal constraint *count* is not inf-sup stability,
though, which is the open question in section 6.

What the Petrov-Galerkin form changes is the distribution of internal force and
energy across those modes, not which modes exist.  It is required by eq. (17)
regardless of how large the difference turns out to be.

**Consequence for a CalculiX implementation.** ccx calls PARDISO with
`mtype = -2`, symmetric indefinite, and stores one triangle.  A correct
F-barES-FEM-T4 cannot be assembled into that path.  It needs `mtype = 11` and
a full-matrix storage scheme — a solver-level change, not a `*USER ELEMENT`.
This is a second, independent blocker on top of the `lakon(8:8)` 255-node
connectivity cap, which the `E A^c` stencil (~2c+1 element rings) breaches at
`c ≥ 2` and is already tight at `c = 1`.

## 5. The envelope the paper claims for itself

Verbatim, section 3.1:

> "According to our previous research for hyperelastic materials [Onishi et al.
> (2017)], the recommended c of F-barES-FEM-T4 for this material (**Poisson's
> ratio ν = 0.49 at most**) is 1 or 2."

and the validation material is `E = 1 MPa, ν = 0.3` instantaneous relaxing to
`ν ≈ 0.49` long-term.  **ν = 0.49 is K/G ≈ 50.**  Our operating point is
K/G = 500 (ν = 0.499); brine as specified is K/G = 5000 (ν = 0.4999).  So we
are asking the method for one to two orders of magnitude beyond anything the
authors published, with `c` stated to be ν-dependent and no formula given for
it.  That is a reason to measure our own operating point rather than to assume
either outcome; section 7 does.  It is not, on the evidence there, a wall.

## 6. Verification

`elements_ccx/tests/verify_fbar.py`, all passing:

* **V1** unit row sums on `Q`, `P`, `E`, `R` and on the composed chain `S`, to
  1e-15, for `c = 0..3` and both readings of eq. (6).  A constant `J` field
  survives the chain untouched, which is what makes section 3 exact.
* **V2** patch test: homogeneous block, uniform strain, `C1111 = K + 4G/3`
  recovered to 4e-16 relative for every `c` and both readings.
* **V3** `S = E` identically at `c = 0`, so the scheme collapses to selective
  ES-FEM-T4 exactly, as it must.
* **V4** the constraint-count table in section 4.

`elements_ccx/tests/verify_fbar_nl.py` checks the finite-strain element of
section 2 against closed-form answers, all passing:

* **N1** `f(0) = 0`; **N2** `f(rigid translation) = 0`.
* **N3a** unjittered box, free surfaces, uniform stretch: interior nodal
  residual vanishes to 1e-13 relative, and the resultant on the `x = 1` face
  matches the analytic Hencky Cauchy stress `(K + 4G/3) ln λ` to 1e-13, for
  `c = 0..3` at both `λ = 1.001` and `λ = 1.2`.  A 20% stretch, so this
  exercises eqs. (1)-(18) well away from the linear regime.
* **N3b** jittered periodic cell, homogeneous: the affine field is recovered
  with fluctuation 2e-12.
* **N4** the finite-strain element converges to its own small-strain reduction
  at O(eps): 6e-6 relative in C1111 at `eps = 1e-5`, for `c = 0, 1, 2`.

The reading ambiguity in eq. (6) — the input to the first cycle is written
`ᵉJ̃`, with the tilde the paper otherwise reserves for smoothed quantities, so
it is either the raw element Jacobian `det(ᵉF)` or an element restriction of
the already edge-smoothed `ᵸJ̃` — is implemented both ways and switchable with
`SPAX_FBAR_JIN=elem|edge`.  It changes the numbers by a few per cent and
changes no conclusion below.  `elem` is the default.

## 7. What it does, measured

**A retraction first.**  An earlier version of this section reported that the
method comes apart above K/G = 50 — fluctuation growing monotonically with `c`
and 29% of C1111 lost by `c = 3` — and concluded that it could not reach our
operating point.  That was wrong, and the cause was not the element.  At
`jitter = 0.3` the prototype's structured mesh **tangles**: one inverted tet at
`n = 6`, fourteen at `n = 16`.  `grads()` silently flipped their node order,
which restored positive volumes and hid it, but the tets still overlapped their
neighbours, so the mesh was not a partition of the cell.  On that mesh nothing
passes its own patch test — C3D4 included, at `fluc` 1.4e-3 with a residual of
0.21 against the exact affine field, at interior nodes.  `mesh_box` now shrinks
the jitter amplitude until no tet is inverted and the worst is at least
`SPAX_MESH_QMIN` (0.15) of the mean, the same guard `make_block.py` has carried
since the acceptance suite hit this.  Every number below is remeasured.

Brine sphere `r = 0.30` in ice, periodic cell, requested jitter 0.3 (shrunk to
0.64 of that at `n = 8` by the quality guard), `n = 8`.  `fluc` is max
displacement fluctuation per unit applied strain; `p-sch` is the scheme's own
pressure field, jump across brine-brine faces over mean `|p|`.

| K/G | scheme | C1111 | fluc | p-sch |
|---|---|---|---|---|
| 500  | c3d4   | 2.6926e+08 | 0.44 | 0.193 |
| 500  | ns_vol | 2.5067e+08 | 0.45 | 0.009 |
| 500  | fbar_0 | 2.5788e+08 | 0.54 | 0.049 |
| 500  | fbar_1 | 2.4946e+08 | 0.51 | 0.005 |
| 500  | fbar_2 | 2.4841e+08 | 0.44 | 0.002 |
| 500  | fbar_3 | 2.4814e+08 | 0.41 | 0.001 |
| 5000 | c3d4   | 2.5110e+09 | 0.38 | 0.139 |
| 5000 | ns_vol | 2.4141e+09 | 0.92 | 0.004 |
| 5000 | fbar_0 | 2.4527e+09 | 0.67 | 0.035 |
| 5000 | fbar_1 | 2.4102e+09 | 1.07 | 0.002 |
| 5000 | fbar_2 | 2.4046e+09 | 0.61 | 0.001 |
| 5000 | fbar_3 | 2.4032e+09 | 0.47 | 0.000 |

Both of the paper's claims now hold at our operating point and beyond it:

1. **Pressure.**  `p-sch` falls monotonically with `c` at every K/G, by two
   orders of magnitude from C3D4 by `c = 3`.  This is Fig. 6 of the paper.
2. **Displacement.**  `fluc` is not degraded by `c`; at K/G = 5000 it is
   non-monotone and `c = 3` (0.47) is *better* than `c = 0` (0.67) and better
   than NS-FEM's 0.92.  C1111 moves only −2.0% from `c = 0` to `c = 3` at
   K/G = 5000, against −29% on the tangled mesh.  This is the paper's "regardless
   of the number of cyclic smoothings".

At K/G = 5000, `fbar_3` is the best arm in the prototype: it is 4.3% softer
than C3D4 (which locks, so softer is the right direction), it carries the
smallest fluctuation of any smoothed scheme, and its pressure field is
effectively oscillation-free.  It beats the incumbent `ns_vol` on all three.

The paper's stated envelope (section 5, `ν = 0.49 at most`) is a statement
about where the authors validated `c`, not a wall we have found.

## 8. Where this leaves a CalculiX implementation

The method works at K/G = 500 and at 5000, so the case for building it is now
open rather than closed.  Two obstacles are real and neither is about accuracy:

* **The tangent is non-symmetric.**  ccx calls PARDISO with `mtype = -2`
  (symmetric indefinite, one triangle stored).  A correct F-barES-FEM-T4 needs
  `mtype = 11` and full-matrix storage — a solver-path change, not something a
  `*USER ELEMENT` can do on its own.
* **The stencil is wide.**  `E A^c` spans roughly `2c+1` element rings, which
  breaches the `lakon(8:8)` 255-node connectivity cap at `c ≥ 2` and is already
  tight at `c = 1`.  The measurements above want `c = 2` or `3`.

So the honest position is: the physics is validated in the prototype and the
blockers are structural to ccx's assembly and solver paths.  Sizing that work
is the next decision, and it should be taken against the incumbent — the
unstabilised U5+U6 nodal B-bar at K/G = 500, whose known defect is a 2-3 point
over-softening that grows with refinement.
