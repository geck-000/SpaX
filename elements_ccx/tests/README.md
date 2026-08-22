# U4 patch tests

Both prescribe a uniform strain state on **every** node of a single tet, so the
exact answer is a constant stress and any consistent element must reproduce it
to roundoff. Run with `ccx_u4 u4test` / `ccx_u4 u4brine` and read the `.dat`.

| Deck | Geometry | Material | Point |
|---|---|---|---|
| `u4test.inp` | unit tet on the axes | E = 9.43e9, ν = 0.33 | compressible baseline |
| `u4brine.inp` | distorted, off-axis | K = 2.2 GPa, G = 0.44 MPa, ν = 0.4999 | the real case, and it exercises the global-derivative path |

Results, against closed form:

| | σxx | σyy = σzz | shear |
|---|---|---|---|
| `u4test` analytic | 1.397192e+07 | 6.881690e+06 | 0 |
| `u4test` U4 | **1.397192E+07** | **6.881690E+06** | ≤ 2e-10 |
| `u4brine` analytic | 2.200587e+06 | 2.199707e+06 | 0 |
| `u4brine` U4 | **2.200587E+06** | **2.199707E+06** | ≤ 7e-12 |

Exact at all 15 integration points in both.

`u4brine` is the one that matters twice over. It is at the brine's own Poisson
ratio, where a displacement element locks; and it is deliberately not
axis-aligned, because on a unit axis-aligned tet `dN/dxi` equals `dN/dx` and the
test cannot see whether the shape-function derivatives were converted to global
coordinates. They were not, at first — `shape4tet` returns at `iflag=2` before
applying the inverse Jacobian, and only `iflag=3` gives global derivatives.
A unit-tet-only test would have shipped that bug.

**What these do not test:** the inf-sup behaviour. A patch test is satisfied by
plenty of unstable elements — it says the element is consistent, not that it is
stable. Stability is what MINI's bubble is for, and demonstrating it needs a
real constrained problem, not one tet.
