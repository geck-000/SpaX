# Referee report, second round — MAMS

The first-round points are adequately addressed. The new verification section is
exactly what was needed, and the retraction of the two-phase coefficient claim is
handled with more candour than most authors manage. The description of the
second-order bending scheme (Eq. 5) closes what was, on reflection, my own
misreading rather than an omission of substance — though the manuscript did
invite it by describing only the first-order scheme.

Four further points, one of which I consider serious.

---

## R2-1 (serious). The percolation-threshold agreement is circular.

Section 2.2 states that the transition from isolated pockets to a connected
network "is **imposed** at the rule-of-fives threshold $\phi_b \approx 0.05$".
The Conclusions then state that "vertical connectivity switches on at the 'rule
of fives' threshold of Golden et al. (1998), **which the step in our anisotropy
recovers without having been fitted to it**".

These cannot both stand. The anisotropy steps where the channels are switched
on, and the channels are switched on at $\phi_b=0.05$ by construction. The
agreement is an identity, not a prediction, and presenting it as independent
confirmation of Golden's threshold is the kind of claim that damages an otherwise
careful paper. Either remove it, or — better — say what the study genuinely
shows, which is that *given* a network switched on at the observed threshold, the
mechanical consequence is a step in anisotropy of a specific size.

## R2-2. The micromechanics benchmark tests the wrong quantity.

The new Section 4 benchmarks the isotropic knockdown against spherical
Mori–Tanaka. But the paper's headline claim is an *anisotropy* produced by
*aligned, non-spherical* pockets, and spherical estimates are blind to it by
construction — they return $E_z/E_x = 1$ identically. The verification therefore
does not touch the result that most needs it, particularly now that the authors'
own discretisation study bounds the anisotropy error at 0.007, which is larger
than the 0.0043 they claim at the cold surface. Eshelby's solution for a spheroid
is closed form and the dilute estimate for aligned spheroids is a short
calculation; it should be done.

## R2-3. "Sphericity" is used throughout but never defined, and is not sphericity.

The term appears from Section 2.2 onward as the primary shape descriptor, with
values 0.85 down to 0.58, and is nowhere defined. Inspection suggests it is the
ratio of semi-axes, not the geometric sphericity $\pi^{1/3}(6V)^{2/3}/A$. These
differ substantially — a semi-axis ratio of 0.6 corresponds to a geometric
sphericity near 0.94, and conversely a geometric sphericity of 0.6 implies an
aspect ratio near 9. Readers comparing against tomography literature, which
reports geometric sphericity, will be misled by a large factor. Define it.

## R2-4. Two smaller matters of internal consistency.

(a) The effective Poisson ratio is never reported, yet the lamination model uses
$Q = E/(1-\nu^2)$. Which $\nu$? If the matrix value 0.33 is used while the
computed cells give 0.319 to 0.347 varying with depth, that should be stated and
its effect on $B/\sqrt{AD}$ bounded.

(b) The "41 negative eigenvalues" is, I take it, the solver's factorisation
diagnostic rather than an eigenvalue extraction. As written a reader may take it
for a buckling analysis. Say which it is.
