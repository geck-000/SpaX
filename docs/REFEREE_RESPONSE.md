# Referee report — Mechanics of Advanced Materials and Structures

**Manuscript:** Sea ice as a depth-graded particle-weakened composite: periodic
homogenisation from micro-CT-constrained RVEs

**Recommendation:** Major revision

The manuscript is careful, unusually well replicated, and refreshingly explicit
about its own limits (Section 5.1 is a model of its kind). The depth-resolved
stiffness, the separation of amount from arrangement, and the matched-control
treatment of the bending size effect are all solid contributions. My concerns are
concentrated in one area: the manuscript is written as a geophysics paper and
submitted to a mechanics-of-composites journal, and the verification and
benchmarking a composites readership expects are largely absent.

---

## Major points

**M1. No comparison with classical micromechanics.**
The central result is a knockdown law, Eq. (6), for a two-phase particulate
composite with spherical-to-ellipsoidal inclusions. Nowhere is it compared with
Mori–Tanaka, the Hashin–Shtrikman bounds, or the dilute Eshelby estimate. For
this journal that is the first question a reader will ask. Either the computed
law agrees with the standard estimates — in which case the authors have a
verification result they are not claiming — or it does not, in which case the
discrepancy is itself the finding. As it stands the reader cannot tell which.

**M2. The homogenisation implementation is not verified.**
The Hill–Mandel condition is invoked five times but never checked numerically. No
patch test is reported. There is no mesh-convergence study: cell-size convergence
(Section 4.1) is a different thing and does not bound discretisation error. Given
that the paper's headline claims rest on effects of order 0.4% in $E_z/E_x$, the
reader needs to know the discretisation error is smaller than that.

**M3. The two per-phase coefficients are not physically defensible as stated.**
Eq. (6) gives $1.68\,\phi_\mathrm{brine}$ and $1.64\,\phi_\mathrm{gas}$, and the
text states the two phases are "indistinguishable per unit volume", quoting
standard errors of 0.01 that make the 0.04 difference look resolved at 3 sigma.
But a void is strictly more compliant than a pocket filled with a phase of finite
bulk modulus, so the gas coefficient must exceed the brine one. Mori–Tanaka puts
them at 2.00 and 1.84 respectively. The fit orders them the other way. Whatever
is causing that — different volume-fraction estimators for meshed and non-meshed
phases is the obvious candidate — it means the per-phase split is at or beyond
the resolution of the data, and the claim of indistinguishability is being made
on the wrong side of the physics.

**M4. The second-order scheme is implemented but not described.**
[WITHDRAWN on inspection of the source. My original point asserted that no
second-order boundary value problem was solved. That is wrong: the bending
studies impose a quadratic macroscopic field driven by a curvature reference
point, with Lesicar Eq. 14 integral constraints on the face-mean fluctuations,
and Abaqus steps named accordingly. The real defect is that Section 3.2
described only the first-order scheme, so a reader could not tell. The revision
adds Eq. (5) and the accompanying description, which also strengthens the
size-effect null: it comes from a scheme that carries a macroscopic curvature
degree of freedom and could have returned a length scale.]

---

## Minor points

**m5.** Twenty-two figures is a great deal. Several are schematics that a
specialist reader does not need (workflow, homogenisation schematic, bending
method).

**m6.** The Weibull modulus is swept from 1 to 50 but never calibrated against a
measured value for sea ice. State that the sweep is a sensitivity analysis, not a
prediction.

**m7.** The finite-strain study uses one packing per slice. The instability
strains (2.8%, 2.0%, 0.24%) carry no replication, yet an order-of-magnitude
ordering is claimed from them. Either replicate or attach the caveat.

**m8.** Section 4.4.2 compares against four beams from a single 1990 campaign.
That is a thin observational base for a claim about the domain of validity of the
whole method.
