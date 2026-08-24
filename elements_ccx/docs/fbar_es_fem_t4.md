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

## 4. Why the symmetric form fails, and by how much

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
energy across those modes, not which modes exist.  Measured on the brine-sphere
cell (n=8, jitter 0.3, K/G 500), displacement fluctuation per unit applied
strain:

| c | `B̄ᵀDB̄` (wrong) | `B̃ᵀDB̄` (paper) |
|---|---|---|
| 0 | 0.51 | 0.51 (identical by construction) |
| 1 | 2.29 | 1.77 |
| 2 | 3.19 | 2.38 |

The correction is real and in the right direction, and it is required by
eq. (17) regardless of its size.  It is not on its own sufficient to reverse
the trend with `c`.  See section 6.

**Consequence for a CalculiX implementation.** ccx calls PARDISO with
`mtype = -2`, symmetric indefinite, and stores one triangle.  A correct
F-barES-FEM-T4 cannot be assembled into that path.  It needs `mtype = 11` and
a full-matrix storage scheme — a solver-level change, not a `*USER ELEMENT`.
This is a second, independent blocker on top of the `lakon(8:8)` 255-node
connectivity cap, which the `E A^c` stencil (~2c+1 element rings) breaches at
`c ≥ 2` and is already tight at `c = 1`.

## 5. The envelope the paper actually claims

Verbatim, section 3.1:

> "According to our previous research for hyperelastic materials [Onishi et al.
> (2017)], the recommended c of F-barES-FEM-T4 for this material (**Poisson's
> ratio ν = 0.49 at most**) is 1 or 2."

and the validation material is `E = 1 MPa, ν = 0.3` instantaneous relaxing to
`ν ≈ 0.49` long-term.  **ν = 0.49 is K/G ≈ 50.**  Our operating point is
K/G = 500 (ν = 0.499); brine as specified is K/G = 5000 (ν = 0.4999).  Every
published result for this method sits one to two orders of magnitude below
where we need it, and `c` is stated to be ν-dependent with no formula given.

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

The reading ambiguity in eq. (6) — the input to the first cycle is written
`ᵉJ̃`, with the tilde the paper otherwise reserves for smoothed quantities, so
it is either the raw element Jacobian `det(ᵉF)` or an element restriction of
the already edge-smoothed `ᵸJ̃` — is implemented both ways and switchable with
`SPAX_FBAR_JIN=elem|edge`.  It changes the numbers by a few per cent and
changes no conclusion below.  `elem` is the default.

## 7. What it does, measured

Brine sphere `r = 0.30` in ice, periodic cell, jitter 0.3, `n = 8`.  `fluc` is
max displacement fluctuation per unit applied strain; `p-sch` is the scheme's
own pressure field, jump across brine-brine faces over mean `|p|`.

| K/G | ν | `fluc` c=1 / 2 / 3 | `p-sch` c=0 → 3 | C1111 loss c=0 → 3 |
|---|---|---|---|---|
| 10   | 0.458  | 0.41 0.41 0.41 | 0.081 → 0.018 | −7% |
| 25   | 0.481  | 0.41 0.41 0.41 | 0.076 → 0.010 | −7% |
| 50   | 0.4901 | 0.42 0.41 0.42 | 0.071 → 0.006 | −6% |
| 100  | 0.495  | 0.55 0.68 0.75 | 0.066 → 0.004 | −6% |
| 250  | 0.498  | 1.05 1.38 1.61 | 0.060 → 0.003 | −7% |
| 500  | 0.499  | 1.77 2.38 2.86 | 0.057 → 0.002 | −9% |
| 1000 | 0.4995 | 3.03 4.08 5.11 | 0.054 → 0.002 | −12% |
| 5000 | 0.4999 | 10.67 13.26 16.05 | 0.043 → 0.004 | −29% |

Two things are true at once, and both are reproduced faithfully:

1. **The pressure claim holds everywhere.**  `p-sch` falls monotonically with
   `c` at every K/G, which is the paper's Fig. 6.  Nothing about the
   checkerboard suppression fails at high K/G.
2. **The displacement claim holds only inside the paper's envelope.**  Up to
   K/G = 50 — precisely the `ν = 0.49 at most` of section 5 — `fluc` is flat in
   `c`, which is the paper's "regardless of the number of cyclic smoothings".
   From K/G = 100 it starts to drift, and by K/G = 5000 it is an order of
   magnitude out with 29% of the macroscopic stiffness gone.

A 29% loss of C1111 in a cell that is ~12% brine by volume is not locking
relief; no redistribution of compliance within the soft phase can produce it.
The oscillation has not been removed, it has been **moved**: out of the
pressure, which is now averaged over a `2c+1`-ring stencil and therefore looks
smooth by construction, and into the displacement, where `ker(G_trial)` — 625
of 933 modes at `c ≥ 1`, section 4 — does not penalise it.  A smooth pressure
field is a necessary condition for stability, not a sufficient one, and this is
what the difference looks like.

## 8. Conclusion for the campaign

Implemented to the paper and verified against it, F-barES-FEM-T4 reproduces
every published claim inside `ν ≤ 0.49` (K/G ≤ 50) and fails outside it.  Our
operating point is K/G = 500 and the physics target is 5000.  Dropping to
K/G = 50 to enter the envelope costs a full floor (+1.73%/+1.90% on R, measured
earlier in the campaign), which is the thing we refused for exactly this
reason.

Two further blockers stand behind that one, both from section 4, and neither is
worth paying down until the envelope problem is solved:

* the tangent is non-symmetric, and ccx calls PARDISO with `mtype = -2`
  (symmetric indefinite, one triangle stored) — a solver-level change, not a
  `*USER ELEMENT`;
* the `E A^c` stencil spans ~`2c+1` element rings, which breaches the
  `lakon(8:8)` 255-node connectivity cap at `c ≥ 2` and is already tight at
  `c = 1`.

Recommendation: do not implement F-barES-FEM-T4 in CalculiX.  The incumbent
unstabilised U5+U6 nodal B-bar at K/G = 500 remains the better candidate; its
known defect is a 2-3 point over-softening that grows with refinement, which is
a smaller and better-characterised error than what is above.
