# Figures for ice_rve.tex

Ten figures, mapped to sections. Status is one of **port** (exists in main_fix,
carries over), **have** (already generated in this campaign), **build** (script
needed), **blocked** (waiting on a running campaign).

The guiding rule: each figure should carry one claim, and the claim should be
one that could have come out otherwise. Figures 6 and 7 are the two where a
prediction was tested and could have failed, so they get the most space.

---

## §2 Microstructure and cell generation

### Fig 1 — The sheet, the column and the cell  ·  **port**
Schematic in tikz: floating sheet, the depth slices, one cell extracted, the
coordinate convention. Unchanged from main_fix; the preamble in ice_rve.tex is
already identical so it compiles without edits.

### Fig 2 — The three morphologies, as meshed  ·  **have (c), build (a,b)**
Three columns, each showing the brine phase read back from a solved `.inp`
rather than from the geometry that was requested, so what is shown is what was
solved.

| panel | morphology | source |
|---|---|---|
| (a) | pockets | `plot_slab_mesh.py` on a `BRKP_*` deck |
| (b) | channels | same, on a channelled column deck |
| (c) | layers with bridges | `results/slab_mesh.png` — **have** |

Panel (c) already exists and its layer-plane slice, where the bridges appear as
clean gaps, is the single most useful image in the paper: it makes the whole
constriction argument visible before any of the numbers arrive.

---

## §3 Homogenisation and verification

### Fig 3 — The layered cell is a representative volume  ·  **build**
Two panels, and the second is what makes it worth a figure rather than a table.

- (a) `E_x` against cell edge at **fixed microstructure** (spacing 0.125,
  density 64), drained and undrained. Flat: CV 1.9% and 0.8%.
- (b) the same sweep with the bridge **count** held instead, giving `L^-1.14`
  and 39%. Both drawn on the same axes.

Data: `results/results_bracket_density.csv`, `results/results_bracket_spacing.csv`.
Claim: the cell homogenises, and the earlier apparent failure was a
microstructure changing under the test.

---

## §4 The pocket column and the limit of its range

### Fig 4 — The pocket column against the field inversions  ·  **port + update**
main_fix Fig 18 rebuilt on re-run data. Left: `E(z)` against Kujala's inferred
profile and Marchenko's Kerr–Palmer fit. Right: the three shape metrics,
`alpha`, `z0/H`, `Et/Eflex`, computed and measured side by side.
Blocked on the column campaigns finishing.

### Fig 5 — What does not explain the gap  ·  **build**
One panel, five bars: aspect ratio (<3%), channels (null), pocket drainage
(1.04x), the porosity route (needs phi = 0.52 against 0.23 available), and the
inversion artefact (5.7%). Each against the factor of four to seven that has to
be explained, drawn as a band across the figure.
Claim: the failure is not parametric. This is the figure that earns the rest of
the paper.

---

## §5 Three mechanisms

### Fig 6 — Confinement  ·  **have (partial), build**
- (a) `E_x` against the fill's bulk modulus over three decades, at fixed
  geometry. Collapses 6.9x.
- (b) `E_z` on the same cells, flat to 0.6%. This is the panel that identifies
  the mechanism rather than merely reporting it: a layer parallel to the load
  is unconfined and cannot care about the fill's bulk modulus, so the
  invariance is the signature.
- (c) the same release applied to a pocket cell: 1.04x.

Partly in `results/layered_bracket.png` panel (a); needs rebuilding as a
standalone figure with the `E_z` control given equal weight.

### Fig 7 — Constriction  ·  **build**
`E_x` against bridge count at **fixed total bridge area**, log–log, drained and
undrained, with the `N^0.5` spreading-compliance prediction drawn as a line
rather than fitted.

Measured `N^0.458` drained, `N^0.017` undrained.

The caption must state what the alternatives would have looked like: flat would
have meant area-limited, falling would have meant bending-limited. It came out
rising, which is neither, and matches a prediction made before the campaign ran.

Data: `results/results_bracket_nbridges.csv`.

### Fig 8 — Lamellar spacing  ·  **build**
`E_x` against lamellar spacing `a0`, drained and undrained, with Pringle's
measured 200–500 um marked as a band and the solved cells marked as points.
Shows that the physical spacing lies beyond the soft end of what was solved, so
the extrapolation direction is stated rather than hidden.

Data: `results/results_bracket_nlayers.csv`, re-read as a spacing sweep.

---

## §6 A closure for the depth profile

### Fig 9 — The closure and its uncertainty  ·  **have**
`results/ez_closure.png`, four panels: `E(phi)` with the exponent band against
Weeks & Assur, Marchenko and the pocket cells; `E(z)` for a stated column
against Kujala; the two calibration datasets; and the sensitivity apportioned
between ingredients.

Panel (d) is the one to keep: the **assumed** `phi_0` is worth x3.5, the
**measured** spacing x2.5, the **calibrated** exponent only x1.3. The largest
uncertainty is the parameter nobody has measured, which inverts the usual
worry about calibrated models.

### Fig 10 — The three comparison cases  ·  **have**
`results/match_ez.png`: Marchenko on his own inverted porosity, Gogolaze's
required exponent against what theory allows, and the porosity Kujala's beams
imply against the synthetic column.

---

## Blocked on running campaigns

- Fig 4 needs the column campaigns from the full re-run (at campaign 9 of 36).
- Fig 7's exponent claim is corroborated by `bracket_bridge`, running now,
  which measures the b-dependence at fixed N. If it disagrees with
  `sqrt(b nu)`, Fig 7 stands but §6's calibration has to be revisited.

## Not planned, and why

No figure for the strength correlations. They sit below the isotropic bounds
and so contain damage this elastic closure does not model, and putting them on
the same axes invites the reader to treat them as a target.
