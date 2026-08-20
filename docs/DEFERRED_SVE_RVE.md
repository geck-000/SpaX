# Deferred: a probabilistic SVE-to-RVE criterion for the layered cells

**Status: NOT STARTED. Deliberately queued behind the closure work.**
Raised 2026-08-20. Do not begin this until the bridge-factor question is
settled and written up — it answers a different question and would compete for
the same cells and the same attention.

## Why it is worth doing at all

Our layered cells are not size-converged and we know it. The Eringen size sweep
(`results_eringen_layer.csv`, with its matched `_homog` control) gives, after
normalising by the control:

| L | n_slabs | normalised | vs L=0.48 |
|---|---|---|---|
| 0.24 | 2 | 0.3581 | +45.9% |
| 0.36 | 3 | 0.2729 | +11.2% |
| 0.48 | 4 | 0.2454 | — |

The control moves only 2.2% across that range, so the trend is material rather
than extraction bias. The normalised ratios are still falling at the largest
size, which puts the L = 0.5 we use everywhere roughly 5% above the asymptote.

We currently defend L = 0.5 by arguing the residual is **common-mode**: LCOL,
LAYERB, RAMP and SUBC all carry it, so it cancels out of any comparison across
them. That argument is sound for ratios and silent about absolute level — and
it is an argument, not a measurement. See `analysis/layered_cell_size.py`.

## The method

Ostoja-Starzewski, *A probabilistic measure of SVE-to-RVE convergence*,
Probab. Eng. Mech. **85** (2026) 103979.

Solve two canonical BVPs on a finite cell: Dirichlet (uniform displacement)
gives a stiffness `C^d`, Neumann (uniform traction) gives a compliance `S^n`.
Dirichlet is always the stiffer; the two coincide only in the RVE limit. Invert
the Neumann response so both are stiffnesses, treat each as a random element of
the SPD cone, map `T = log C` into a flat space where ordinary probability
applies, fit a Gaussian over realisations, and take the **Hellinger distance**
between the two laws. The RVE size is the smallest `L` with `H < ε`.

What that buys over our current size studies is a criterion that is
coordinate-invariant and that uses the whole distribution rather than the mean
— which matters here, because the quantity we would be certifying is a tensor
with `E_z/E_x` around 2.3, not a scalar.

## Two routes, and only one of them is cheap

**Route B — the appendix protocol. Feasible with the code as it stands.**
It compares *sizes* rather than boundary conditions: find the smallest `L` with
`H(p_L, p_2L) < eps_H`. The protocol never specifies a boundary condition, so
our existing periodic cells qualify. This needs no generator work at all.

**Route A — the paper's core, its Eq. (5).** Needs Dirichlet and Neumann
responses. **SpaX has neither.** The generator is periodic-only: Gmsh periodic
pairs, plus Lesicar Eq. 14 constraints for the bending cases. Periodic cannot
stand in for either bound — it sits between them and converges faster, so it
cannot produce the gap the method measures. Adding KUBC and SUBC is a real
piece of generator work and roughly doubles the load cases on top of the six a
full tensor already needs.

## The binding constraint is dimensionality, not compute

The paper is **two-dimensional** throughout: `C` is 2x2 SPD, `log C` lives in
R^3, and the Gaussian needs a 3x3 covariance. Extension to 3D is listed as the
author's own future direction #2 — it is not done in the paper.

In 3D, `C` is 6x6 in Voigt form, so `log C` carries 21 independent components
and its covariance has 231 free parameters. At the usual `N >= 10 x dim` that
is about 210 realisations per size. The largest ensemble anywhere in this
repository is 50 rows.

Two reductions make it tractable, and the middle one is physically honest here
rather than a convenience:

| representation | dim | Sigma params | realisations needed |
|---|---|---|---|
| full 3D | 21 | 231 | ~200+ |
| transversely isotropic (5 constants) | 5 | 15 | ~50 |
| scalar, `E_x` alone (paper §3.3) | 1 | 1 | ~25-30 |

Our cells *are* transversely isotropic by construction — lamellae normal to x,
growth along z — so the five-constant form is the natural one. The paper also
gives the scalar case in closed form, in modified Bessel functions.

## What we would need

1. `full_tensor='Yes'` as standard rather than exceptional. The capability
   already exists (`rve_coltensor.csv`, `rve_failure_rep.csv`,
   `rve_basetensor*.csv`); it costs six load cases instead of two.
2. Ensembles of ~30-50 realisations at two or more sizes. Cheaper than it
   sounds: an SVE study *wants* sub-RVE cells, and L = 0.125 and L = 0.25 carry
   1/64 and 1/8 the elements of the cells we have been solving.
3. A correlation length to define `delta = L / l_c`. We have one: the lamellar
   spacing `a0 = 0.75 mm`.
4. For Route A only: KUBC and SUBC in the generator.

## Recommended first cut, when the time comes

**Route B, scalar, on `E_x`** — roughly 30 periodic cells each at L = 0.25 and
L = 0.5, no generator changes, runnable on a workstation. That alone would
convert "we believe the size bias cancels" into a number. Escalate to the
five-constant transversely isotropic form only if the scalar answer sits near
the tolerance, and to Route A only if a referee asks for the Dirichlet-Neumann
criterion specifically.

## What this does NOT settle

Nothing about the bridge factor. `g = b^n` is unsupported for reasons that have
nothing to do with cell size — see `results/control/bt_grid.txt`, where n varies
with b by 0.15-0.17 at *fixed* t. An SVE-to-RVE criterion would certify the size
of the cells that produced that result; it would not change the result.
