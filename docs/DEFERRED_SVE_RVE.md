# Deferred: putting a number on how far our cells are from an RVE

**Status: NOT STARTED. Queued behind the closure work.** Raised 2026-08-20,
revised the same day after reading both papers properly.

- **[A]** Ostoja-Starzewski, *Material spatial randomness: from statistical to
  representative volume element*, Probab. Eng. Mech. **21** (2006) 112–132.
- **[B]** Ostoja-Starzewski, *A probabilistic measure of SVE-to-RVE
  convergence*, Probab. Eng. Mech. **85** (2026) 103979.

## Correction to the first version of this note

The first draft recommended computing the Dirichlet/Neumann bracket of [A]
Eq. (2.21). **That is not viable at our contrast and the paper says so.** Our
phase contrast is

    ice 9.37 GPa / drained brine 1.24 GPa e-3  =  7572          (undrained: 7099)

and [A] §4.4.2 states that bringing the gap between `<S^n_d>^-1` and `<C^e_d>`
down to about 30% needs delta = 10 at contrast 1e2 and **delta = 50 at contrast
1e4**. We sit essentially at 1e4. delta = 50 is fifty lamellae, roughly 1950x
the elements we solve now, and it buys a *30% bracket* — wider than the effects
we are arguing about. The two-sided bounds route is dead on arrival here; do not
spend generator work on `dd`.

## What actually governs: the bridges, not the lamellae

[A] Eq. (1.7) defines `delta = L/d`. The choice of d is not obvious here and it
changes the conclusion.

| candidate microscale | d | delta |
|---|---|---|
| lamellar spacing a0 | 0.125 | 4.0 |
| pocket diameter | 0.060 | 8.3 |
| **bridge diameter, b=0.30, n=2** | **0.219** | **2.3** |
| bridge diameter, b=0.20, n=2 | 0.178 | 2.8 |
| bridge diameter, b=0.10, n=2 | 0.126 | 4.0 |

The lamellar spacing is the wrong answer, and [A] §1.1 says why: an RVE is
*exactly* defined for the unit cell of a periodic microstructure. Our lamellae
are deterministic and evenly spaced by construction, so they need no
homogenising. The randomness is the pocket packing and, above all, the **bridge
placement** — and the bridges are the load-bearing feature, since the transverse
path crosses the planes through them.

At `n_bridges = 2` a bridge is **13–22% of the cell edge**, and there are two per
plane over four planes. `delta_bridge = 2.3–4.0` is the smallest number in the
problem and the one that governs. That is far below any RVE criterion, and it is
the honest explanation for the size effect we measured: the Eringen sweep,
normalised by its matched control, runs 0.358 / 0.273 / 0.245 at two, three and
four lamellae and is still falling.

It also predicts something testable that we are *already* running: raising
`n_bridges` from 2 to 4 at fixed b raises `delta_bridge` from 2.3 to 3.2. The
`NBR_p095_n4` cells in the gap deck were built to test whether the sharp feature
at phi ~ 0.095 is a two-bridge artefact — if bridge count matters there, it is
the same finite-size story.

## What [A] gives us that is good news

**Fig. 7(b), contrast 1000.** Six boundary conditions against delta: `dd`,
`tt`, `pp`, `dp`, `dt`, `tp`. `dd` starts near 37 at delta = 4 and is still ~5 at
delta = 48 — hopeless. **`pp`, `tt`, `dt` and `tp` sit together at the bottom,
essentially flat from delta = 4 onward.**

So at high contrast the displacement-controlled condition is the outlier, and
periodic is close to converged where `dd` is nowhere near. SpaX being
periodic-only stops being a limitation and becomes the reason our delta ~ 3 cells
are usable at all. This is an argument **for** what we already do, and it belongs
in the paper.

## The tension, unresolved

Fig. 7(b) says `pp` is flat from delta = 4. Our own sweep says 46% at two
lamellae, 11% at three, still falling at four. Two candidate explanations, and I
cannot presently choose between them:

1. **Morphology.** Fig. 7 is disk-matrix — isolated soft inclusions. Ours is a
   *percolating* soft phase spanning whole planes, with load crossing through
   bridges. [A] §4.4.2 notes convergence is "relatively much slower (!)" for soft
   inclusions, and a spanning soft phase is the extreme of that.
2. **Measure.** Our sweep is `E_bending` from the Eringen decks. Uniaxial may
   converge faster, and uniaxial is what the closure uses.

Resolving (2) is cheap and is the first thing to do.

## Plan, in order

**Phase 0 — free, do it whether or not the rest happens.**
Put `delta_bridge = 2.3–4.0` and contrast 7572 in the limitations, and cite
Fig. 7(b) to justify periodic BCs at modest delta. This converts the size effect
from an awkward number into an expected consequence of two bridges per plane.

**Phase 1 — 9 local cells, no generator work.**
Repeat the Eringen size sweep under **uniaxial** loading: L = 0.24, 0.36, 0.48
at 2, 3, 4 lamellae with a0 held at 0.12, plus matched phi=0 controls. Decides
whether the 11% residual is a bending artefact or real for the quantity the
closure actually uses. Cells are small (delta 2–4), all local, zero billing.

**Phase 2 — extend delta, if Phase 1 shows the effect survives.**
The honest way to raise `delta_bridge` is more bridges per plane, not a bigger
cell: `n_bridges = 4, 6, 8` at fixed b and fixed L. That holds the areal ice
fraction — the physics the closure depends on — while subdividing it, and it
costs nothing in element count. A bigger cell at fixed mesh goes as L^3 and is
not affordable locally past L = 0.6.

**Phase 3 — `tt` only, if a rigorous bound is wanted.**
Add the uniform-traction condition ([A] Eq. 2.12) and *skip* `dd`. At our
contrast `tt` converges fast and gives a rigorous lower bound on C^eff, which
can be set beside the periodic value. One BC of generator work rather than two,
and the one that is actually informative here. Note the risk: with a percolating
soft phase the traction problem is close to the "holes" limit and may be poorly
conditioned.

**Phase 4 — [B], only on request.**
The Hellinger criterion, scalar (`E_x`, closed form, ~25–30 realisations) or
transversely isotropic (5 constants, ~50). [B] is 2D throughout and its own
future direction #2 is the 3D extension; full 3D wants ~200 realisations for a
21-component covariance. [B]'s appendix protocol — `H(p_L, p_2L) < eps` — needs
no new boundary condition and would work on Phase 1/2 output directly.

## What none of this settles

Nothing about the bridge factor. `g = b^n` is unsupported for reasons unrelated
to cell size — see `results/control/bt_grid.txt`, where n varies with b by
0.15–0.17 at **fixed** t across three rows. An RVE criterion would certify the
size of the cells that produced that result; it would not change it.

Though there is one connection worth watching: if `n_bridges` turns out to
matter (Phase 2, or the `NBR` cells already running), then the b-dependence of n
and the finite-size story are the same phenomenon seen twice, and both would need
restating together.
