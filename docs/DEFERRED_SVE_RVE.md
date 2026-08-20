# Deferred: putting a number on how far our cells are from an RVE

**Status: NOT STARTED. Deliberately queued behind the closure work.**
Raised 2026-08-20. Do not begin until the bridge-factor question is settled and
written up — it answers a different question and would compete for the same
cells and the same attention.

Two papers, both Ostoja-Starzewski:

- **[A]** *Material spatial randomness: from statistical to representative
  volume element*, Probab. Eng. Mech. **21** (2006) 112–132.
- **[B]** *A probabilistic measure of SVE-to-RVE convergence*, Probab. Eng.
  Mech. **85** (2026) 103979.

**[A] is the better entry point for us**, and the reason is statistical cost,
set out below. [B] is the refinement.

## Where our cells actually sit

[A] Eq. (1.7) defines the mesoscale parameter `delta = L/d`, L the cell edge and
d the microscale, and notes that the usual recipes put an RVE at delta of 10–100
while warning that it "strongly depends on the type of problem studied."

For the layered cells the microscale that matters is the **lamellar spacing**,
because the bridges live in those planes — not the pocket diameter:

| | d | delta |
|---|---|---|
| layered cells, L = 0.50 | a0 = 0.125 | **4.0** |
| same cells, by pocket diameter | 0.060 | 8.3 |
| Eringen sweep, 2 / 3 / 4 lamellae | 0.12 | 2.0 / 3.0 / 4.0 |

**delta = 4 is the honest number**, and it is well under the heuristic. That is
not a footnote — it is the explanation for the size effect we measured. The
Eringen sweep, normalised by its matched `_homog` control, runs 0.3581, 0.2729,
0.2454 at two, three and four lamellae: still moving 11% on the last doubling
and still falling. A cell at delta = 4 is not expected to be converged, and ours
is not.

We currently defend L = 0.5 by arguing the residual is **common-mode** — LCOL,
LAYERB, RAMP and SUBC all carry it, so it cancels from any comparison across
them. That is sound for ratios, silent about absolute level, and it is an
argument rather than a measurement. See `analysis/layered_cell_size.py`.

## Why [A] is cheap and [B] is not

**[A] gives bounds.** Its Eq. (2.21) is

    <S^t_delta>^-1  <=  C^eff  <=  <C^d_delta>

the effective stiffness lies between the harmonic average of the Neumann
(uniform traction) moduli and the arithmetic average of the Dirichlet (uniform
displacement) moduli, over the same ensemble. Eq. (4.12) extends this to a
hierarchy that tightens monotonically as the window grows. So one **brackets**
C^eff at whatever size one can afford, without ever reaching the RVE and without
extrapolating a trend — which is precisely what our 5%-above-the-asymptote claim
currently rests on.

Statistically this is modest: a bracket needs an ensemble mean, so N of 5–10 is
already useful. **We have that today.**

**[B] gives a distance.** It models the Dirichlet and inverse-Neumann tensors as
random elements of the SPD cone under `T = log C`, fits a Gaussian to each, and
takes the Hellinger distance, with the RVE criterion `H < eps`. Richer — it uses
the whole distribution, not the mean — but it needs a covariance. [B] is **2D
throughout**; 3D extension is its own future direction #2. In 3D `log C` carries
21 independent components, so its covariance has 231 free parameters and wants
roughly 200 realisations per size. The largest ensemble in this repository is 50
rows.

Reductions, with the middle one honest here rather than convenient — our cells
are transversely isotropic by construction, lamellae normal to x, growth along z:

| representation | dim | Sigma params | realisations |
|---|---|---|---|
| full 3D | 21 | 231 | ~200+ |
| transversely isotropic (5 constants) | 5 | 15 | ~50 |
| scalar, `E_x` alone ([B] §3.3, closed form) | 1 | 1 | ~25–30 |

## The one blocker, common to both

**SpaX has no Dirichlet or Neumann boundary condition.** The generator is
periodic-only: Gmsh periodic pairs plus Lesicar Eq. 14 constraints for the
bending cases. Both methods are built on the gap between uniform-displacement
and uniform-traction responses, and periodic cannot stand in for either — [A]
Eq. (2.14) lists `pp` as an admissible Hill-condition BC, so our apparent moduli
are legitimate, but periodic sits between the bounds and converges faster than
either, so it **cannot produce the gap**.

Adding `dd` (Eq. 2.11) and `tt` (Eq. 2.12) is real generator work: uniform
displacement is straightforward, uniform traction needs a consistent surface
load and rigid-body suppression. Both then need the full six load cases, since
the bounds are on tensors.

There is one route that avoids it entirely: **[B]'s appendix protocol** compares
*sizes*, `H(p_L, p_2L) < eps`, and never specifies a boundary condition, so our
existing periodic cells qualify. It certifies convergence rather than bracketing
the answer, which is weaker, but it needs no generator work at all.

## Recommended order, when the time comes

1. **`delta = 4` in the paper, now.** It costs nothing, it is already computed
   above, and it reframes the size effect from an embarrassment into an expected
   consequence of the geometry. Worth a sentence in the limitations regardless
   of whether anything below gets done.
2. **[A] bounds at one size.** Add `dd` and `tt`, solve the existing L = 0.5
   layered conditions under both with `full_tensor='Yes'`, take Eq. (2.21).
   Turns "roughly 5% above the asymptote" into an interval containing C^eff.
   Modest ensembles suffice.
3. **[A] bounds at two sizes**, to show the hierarchy tightening — L = 0.25 and
   0.5. Sub-RVE cells are cheap: 1/8 and 1/64 the elements.
4. **[B]**, scalar or transversely isotropic, only if a referee wants a
   probabilistic criterion rather than a bracket.

## What none of this settles

Nothing about the bridge factor. `g = b^n` is unsupported for reasons unrelated
to cell size — see `results/control/bt_grid.txt`, where n varies with b by
0.15–0.17 at **fixed** t, across three rows. An RVE criterion would certify the
size of the cells that produced that result; it would not change the result.
