# SPAX — TODO, pending runs, and possible studies

Living backlog for the sea-ice RVE homogenization work. Committed to the repo so
it is reachable from any machine (`git pull` on the
`fix/homogenisation-strain-and-parallel-gen` branch). Last updated 2026-07-01.

Generators live in `make_ice_studies*.py` / `patch_brine.py` / `make_seaice_2nd.py`;
offline analysis (pandas only, no Abaqus) is `analyze_studies.py`. Effective moduli
are scale-invariant; matrix conventions: ~9.4 GPa high-frequency for trend studies,
~4.6 GPa "vibrating-beam-effective" (×0.49) for field comparison.

---

## A. IMMEDIATE / PENDING

### A1. Second-order (bending) sea-ice run — ✅ DONE (2026-07-01, all 3 sizes on CSC Roihu)
**Why:** `rve_seaice_2nd.csv` is a quadratic-element (C3D10H), channelled-base
bending size-sweep at L = 0.24 / 0.32 / 0.40 (L/d = 3/4/5, d=0.08), 4 seeds each,
to confirm the MCST length-scale **null** for the realistic sea-ice morphology via
the reference-free Eq.19 slope method. Locally the **L=0.24 and L=0.32 sets solved
fine (8 RVEs)**, but the **L=0.40 quadratic bending jobs are too slow on this
laptop** (each large quadratic `-ben` solve runs many minutes) and were scrapped.
Run those on CSC where the walltime/cores are available.

**What to run on CSC (the heavy part = the 4 `L=0.40` RVEs, all 3 modes each):**
- Decks are already generated in `out_si2nd/` — the 12 `Job-SI2_L400_s*-{utx,ss13,ben}.inp`
  exist. Either copy those `.inp` to CSC and solve directly, **or** regenerate there:
  ```bash
  export SPAX_SEED=0 SPAX_MESH_ORDER=2 SPAX_MAX_RETRIES=12 \
         SPAX_SLIVER_START=1 SPAX_SLIVER_MULT_Q=1.0 SPAX_OPT_PASSES=2
  python3 SpaX_Standalone.py rve_seaice_2nd.csv out_si2nd/   # quadratic, channelled
  ```
- Solve with the standard CSC pipeline (see `RVE_STUDY_README.md`): `01_generate.sh`
  → `02_solve.sh` (array over the `-ben`/`-utx`/`-ss13` jobs) → `03_postprocess.sh`.
  Quadratic bending is the hungry case — give `02_solve.sh` more `--mem-per-cpu`
  (≥8 G) and walltime.
- Post-process and analyse:
  ```bash
  abaqus python SpaX_PostProcess.py rve_seaice_2nd.csv out_si2nd/ results_si2nd.csv
  python3 SpaX_PostProcess.py analyze eq19 results_si2nd.csv     # slope E_app vs 1/L^2
  ```
**Expected result:** slope ≈ 0 → **no measurable length scale** (consistent with the
prior exhaustive null for spheres AND channels — see the `bending-length-scale-artifact`
memory). If confirmed at 3 sizes, the paper's second-order section reports a clean
negative result: classical first-order homogenization suffices; no couple-stress/MCST
model needed. **Local partial result (L=0.24,0.32 only) is in `results_si2nd.csv`** —
2-point slope already points to ~0, CSC adds the third size + scatter to nail it.
RESULT (2026-06-30, L=0.24 & 0.32, 4 seeds each): Eq.19 slope = **-2.0e7 (negative)**
-> l imaginary, **NO MCST length scale**; bend/first-order ratio 0.92/0.96 (<1, softer
not stiffer); intercept = plate modulus (RVE bends as a plate). Confirms the prior
null for the channelled sea-ice morphology. CSC just needs L=0.40 as the 3rd point.
FINAL (2026-07-01, L=0.40 added on CSC Roihu, 4 seeds; `results_si2nd.csv` now 12 rows):
3-size Eq.19 fit slope = **-4.36e7 (negative)** -> l imaginary, **NO MCST length scale**
confirmed at 3 sizes. bend/first-order ratio climbs **0.924 -> 0.960 -> 0.999** across
L=0.24 -> 0.32 -> 0.40 (the first-order softening vanishes as the RVE grows; no size-
independent excess), intercept E0=6.21 GPa = plate modulus. **Closes the paper's
second-order section: classical first-order homogenization suffices, no MCST needed.**

### A2. Paper write-up (Overleaf)  [in progress 2026-07-01]
- `claude.tex` (full study-campaign report, all studies + figures) is already pushed
  to Overleaf (`\input{claude}` in `main.tex`).
- The **second-order null** subsection can now be finalized: the 3-size Eq.19 fit is
  complete (A1), slope = -4.36e7, ratio 0.924/0.960/0.999. Being added to Overleaf.
- The **failure-onset** subsection (B1 below, figure `study_failure.png`, data
  `results_failure.csv`) is likewise ready and being added to Overleaf.
- (Still TODO) the **stress-concentration** subsection (data `results_scf.csv`,
  figure `study_scf.png`) if not already merged into the failure-onset write-up.

---

## B. POSSIBLE FUTURE STUDIES (not yet run)

Ranked roughly by value / novelty for the paper.

1. **Strength / failure-onset mapping.** ✅ **DONE (2026-07-01, CSC Roihu).** Swept a
   Mohr–Coulomb (φ=30°) and max-principal criterion down all 10 column slices
   (`ICE_z05..z95`, uniaxial, `failure_extract.py` → `results_failure.csv`;
   `analyze_failure.py` → `study_failure.png`). **First failure = the porous base
   (z95, 95% depth) for both criteria.** SCF P99 is flat ~1.8–1.9 through the top 70%
   then climbs 2.09→2.27→**4.48** at z75/z85/z95 (MCnorm P99 up to 3.31). First-failure
   macro stress σ_fail = σ_t/SCF_p99 bottoms at **0.223 MPa (tensile) / 0.314 MPa (MC)**
   at the base (σ_t=1 MPa, c=0.6 MPa; the depth *ranking* is strength-independent).
   Confirms the pre-registered expectation: the warm, channelled base cracks first.

2. **Full transverse-isotropy tensor down the WHOLE column.** ✅ **DONE
   (2026-07-05, CSC Roihu).** `study_coltensor()` (make_ice_studies2.py) generated
   `rve_coltensor.csv` — all 10 slices CTEN_z05..z95, `full_tensor=Yes`, same FY
   C-shape physics as the base-4 basetensor (a clean superset). 60 linear C3D4H
   solves (10 slices × 6 load cases) → `extract_elasticity_tensor` per slice →
   `aggregate_coltensor.py` → `results_coltensor.csv` + `study_coltensor.png`.
   **Result:** the column is elastically isotropic through the top ~70%
   (E_z/E_p and G_ax/G_xy within ±0.3% of 1 for z/H≤0.65), then anisotropy emerges
   only in the bottom ~20%: E_z/E_p = 1.003 → 1.022 → **1.031** and shear
   G_ax/G_xy = 1.007 → 1.010 → **1.036** at z/H=0.75/0.85/0.95, while all three
   Young's moduli soften together 8.74 → 4.85 GPa. Confirms, over the whole depth,
   that the base is the only anisotropic zone (both Young's and shear). Cost ~1.7 BU.
   NOTE (fix): CSC batch scripts now `export CSC_ENV_INIT_NON_INTERACTIVE=yes` +
   `source /etc/profile.d/zz-csc-env.sh` before `module load abaqus/2026` — else a
   non-interactive submit leaves the Tcl env-modules active, which cannot parse the
   .lua Abaqus modulefile ("Magic cookie missing") and every solve dies exit 127.

3. **RVE-size convergence for the channelled morphology** (first-order E, not bending).
   The generic size study (`RVE_STUDY_README.md`) used spheres; repeat for the
   percolated channel network to confirm E_z/E_x is box-size-converged. ~5 sizes×5 seeds.

4. **Salinity-profile family.** We did C-shape, monotonic, steeper-monotonic. Add a
   measured Arctic station profile and a linear profile; map profile shape → E(z) and
   anisotropy. Cheap (10 slices each, 2-mode).

5. **Temperature-dependent brine modulus.** ✅ **DONE (2026-07-06, CSC Roihu).**
   `make_ice_studies6.study_brineK()` laid a physically-varying brine bulk modulus
   K(T) onto the FY C-shape column: equilibrium brine salinity (Cox&Weeks/Assur)
   rises 50→215 ppt as T drops toward the surface, so K = 2.25→2.78 GPa (×1.236),
   slope 3.2e6 Pa/ppt. **Single-mesh, geometry-controlled design** (mesh generation
   is NON-deterministic — inclusions are randomly placed — so a matched-seed "paired
   CSV" scheme would bury the K(T) signal under run-to-run mesh scatter): generate
   ONE mesh per slice (`rve_brineKconst.csv`, utx+utz), then `build_brineK_decks.py`
   stamps the K(T) twin onto each mesh by rewriting ONLY the brine `*Elastic` card
   (E,ν from 9KG/(3K+G), (3K−2G)/(2(3K+G))). Paired decks differ by exactly one line
   → byte-identical geometry, mesh, PBC and matrix; the per-slice E diff is the PURE
   K(T) sensitivity. 40 linear C3D4H-hybrid solves (20 RVEs × utx+utz) →
   `results_brineK{const,temp}.csv` → `analyze_brineK.py` → `study_brineK.png`.
   **Result: negligible and clean.** With identical geometry the shift is a UNIFORM
   positive stiffening — all 10 slices, both E_x and E_z, **+0.05% to +0.11%** (max
   0.109%), Δ(E_z/E_x) ≤ 0.0002. (An earlier separate-mesh attempt reported ≤0.22%
   with random signs — that scatter was the mesh-discretization floor, not K(T).)
   The brine is near-incompressible (K/G ~5000–6300), so the ×1.24 K span only nudges
   its Poisson ratio 0.49990→0.49992 — the hybrid pressure DOF resolves it but the
   macro effect is ~0.1%. **Closes the material model: the sea-ice column's effective
   moduli are robust to the brine's thermal-compaction state.** Cost ~2.5 BU (2 runs).

6. **Anisotropic / inclined channels.** ✅ **FEATURE SOLVED (2026-07-06).** Wavy
   mean-vertical channels now mesh with STRICT periodicity (skipped: 0) at every
   tilt/seed tested (0/15/20/30°, 4 seeds each = all pass). Enable via env
   `SPAX_CHANNEL_TILT_DEG` (>0); default 0 keeps the straight cylinder → zero
   regression (verified: straight path unchanged, "0 manual closed").
   THREE root-cause fixes over the earlier shelved attempts:
   (a) **`occ.addPipe` is broken in this OCC build** — it makes a malformed off-centre
       blob even for a straight vertical sweep; that (not the periodic seam, as
       previously thought) was the real "degenerate geometry". Replaced by a boolean
       UNION of short `addCylinder` segments along the centreline (`_add_wavy_channel`).
   (b) **CAD-level z-periodicity** — segment joints are sampled on an L-periodic
       z-grid (spacing L/M, half-shifted) so the segments crossing z=0 and z=L are
       exact translates and the two face cuts match.
   (c) **Vertical-at-faces waveform** — off(z)=amp*(1-cos(2πz/L)) has off=0 AND
       slope=0 at both z-faces, so the cut there is a plain CIRCLE (not a tilted
       ellipse). `setPeriodic(1)` aligns circle cuts fine (as for straight channels),
       dissolving the closed-ellipse-seam problem entirely — no manual closed-curve
       stamping needed. Channel still leans up to atan(2π·amp/L) in the interior
       (bulges laterally ≤ 2·amp), amp = L·tan(tilt)/(2π).
   Azimuth of the lean is per-channel random by default (deterministic hash of the
   channel centre) so inclination has NO net in-plane direction; set
   `SPAX_CHANNEL_TILT_AZIMUTH_DEG` to force a single direction. Packing margin is
   inflated by 2·amp so wavy channels stay off the x/y faces (no XY periodic copies
   needed). NOTE (design): a straight NET tilt is impossible on a cubic periodic cell
   without sheared periodicity; the wavy (mean-vertical) form tests whether channel
   inclination *dilutes* the vertical anisotropy, not whether the axis *rotates*.
   **STUDY ✅ DONE (2026-07-06, CSC Roihu):** 0/15/30° × 4 seeds, utx+utz → E_z/E_x
   (rve_tilt{00,15,30}.csv, submit_tilt.sh, analyze_tilt.py → study_tilt.png).
   **Result — inclination DILUTES the vertical anisotropy (but does not remove it):**
   E_z/E_x = **1.041 ± 0.006 (straight) → 1.029 ± 0.004 (15°) → 1.031 ± 0.009 (30°)**,
   i.e. the anisotropy excess drops to ~70–77% of the vertical-channel value and
   saturates by ~15° (the lateral excursion breaks pure verticality; beyond 15° the
   channels still span z so the ratio flattens). ~0.9 BU. NOTE: the solves first hit
   a transient CSC outage (mid-day move to software stack v2026_03: `/appl` I/O
   errors reading Abaqus libs → exit 127); resolved when CSC recovered. Durable
   fix for the v2026_03 Lmod-in-batch breakage: scripts `source ~/abaqus_env.sh`
   (snapshot of the working interactive PATH/MODULEPATH/LD_LIBRARY_PATH) — see the
   `csc-batch-abaqus-env` memory.

7. **Two-scale / laminated-plate macro model.** Feed the depth-resolved E(z), E_z/E_x
   into a graded laminated-plate model of a whole ice sheet (macro deliverable the
   intro promises). Pure post-processing, no new RVE solves.

8. **Nonlinear / large-deformation (nlgeom).** ✅ **DONE (2026-07-06, CSC Roihu).**
   One geometrically-nonlinear study on 3 representative slices (z25 low-porosity,
   z65 transition, z95 channelled base), matrix kept linear-elastic so only the
   kinematics are finite. Mesh generation is NON-deterministic, so the linear /
   tension / compression cases per slice are stamped onto ONE shared mesh
   (`build_nlgeom_decks.sh`: edit only the *Step nlgeom flag + RP-1 disp) → identical
   geometry. Homogenised response read as the reaction-based nominal path
   σ_nom=RF_RP/L², ε_nom=U_RP/L over the whole ramp (`nlgeom_extract.py`, since the
   standard extractor fits only the initial slope and would hide curvature).
   **Result:** linear ref is exactly straight (sec/E0=1.0000, validates extraction);
   **compression stiffens uniformly +1.5–1.8%** (all slices reach 2%); **tension
   softens −1.7%** (z25) but the **percolated base is geometrically UNSTABLE in
   tension** — the solve loses convergence at 0.7% (z65) and 0.1% (z95) strain as the
   thin ice ligaments across the channels neck/buckle. So the warm base is the
   softest, most anisotropic, most stress-concentrating AND first-to-destabilise-in-
   tension zone. Geometric corrections ≲2% at 2% strain confirm small-strain
   homogenization is quantitatively adequate for the intact sheet.
   `make_ice_studies7.py` (column def), `build_nlgeom_decks.sh` (single-mesh stamp),
   `nlgeom_extract.py`, `analyze_nlgeom.py` → `study_nlgeom.png`. Cost ~1.9 BU
   (incl. a first 4% pass that over-drove the porous tension). Remaining #8 extension
   (viscoelastic/creep brine) still open + heavy.

---

## C. QUICK REFERENCE

**Run a study end-to-end (local):**
```bash
python make_ice_studies.py            # (or 2/3/4) -> rve_*.csv
SPAX_SEED=0 SPAX_MESH_ORDER=1 python SpaX_Standalone.py rve_X.csv out_X/   # decks
# solve Job-*.inp with Abaqus (throttled runner pattern in $CLAUDE_JOB_DIR/tmp/*.ps1)
abaqus python SpaX_PostProcess.py rve_X.csv out_X/ results_X.csv           # extract
python analyze_studies.py             # all figures (pandas only, no Abaqus)
```

**Mesh-order / second-order env flags (quadratic, channelled-safe):**
`SPAX_MESH_ORDER=2 SPAX_MAX_RETRIES=12 SPAX_SLIVER_START=1 SPAX_SLIVER_MULT_Q=1.0 SPAX_OPT_PASSES=2`

**Throttled solve:** 4 jobs × cpus=2 keeps all 12 cores busy without overloading
(machine crashed once at 4 parallel *generators*; keep generation sequential).
**Always delete ODBs after post-processing** (`find out_X -type f ! -name '*.inp' -delete`)
— each campaign is ~5–10 GB of ODBs; the `results_*.csv` hold everything for offline work.

**Done so far (2026-06-29 … 07-05):** 13 microstructure studies + SCF + **full 3-size
2nd-order** (L=0.24/0.32/0.40, MCST null) + **failure-onset mapping** (10 slices,
base cracks first) + **whole-column full 6×6 tensor** (10 slices, anisotropy only in
the bottom ~20%), ~490 solves; the second-order L=0.40, failure-onset, and
whole-column-tensor sets ran on **CSC Roihu** (`abaqus/2026`); earlier batches
committed (`76f257f`, `e3566b5`, `9f35833`); `claude.tex`+figures on Overleaf.
Headline: anisotropy comes ONLY from percolated vertical channels (warm bottom
~15–20%, confirmed down the full column); Marchenko matched to 8% RMS; results
robust to brine G; **no MCST length scale (confirmed at 3 sizes)**; **first failure
at the porous base**.
