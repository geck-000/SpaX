# SPAX — Periodic RVE Homogenization Toolkit

Generate, solve, and homogenize periodic representative-volume-element (RVE)
microstructures (e.g. sea ice: ice matrix + brine inclusions + voids), including
first-order effective moduli, anisotropy, and second-order (couple-stress)
bending response.

## Start here

This repository is two things at once, kept deliberately apart. Pick your path:

**"I want the tool."** Take the three modules at the top level and read
§1–§6 below. Nothing in the subfolders is needed — they are a worked example,
not part of the toolkit. Start at §5, *Worked examples*.

**"I want to reproduce the sea-ice study."** The subfolders are a complete
campaign: decks in `params/`, the outputs they produced in `results/`, and the
analyzers that turn those into figures in `analysis/`. See *Reproducing the
published figures* below.

**"I want to understand a specific number in the paper."** Every folder has its
own `README` naming what each file is for; `docs/` holds the campaign notes.

```
SpaX_Standalone.py  SpaX_GmshPeriodic.py  SpaX_PostProcess.py   <- the toolkit
────────────────────────────────────────────────────────────────────────────
studies/ params/ hpc/ results/ analysis/ viz/ tensors/          <- the campaign
```

The toolkit is **three Python files**:

| File | Run with | Role |
|------|----------|------|
| `SpaX_Standalone.py` | `python3` | packing + periodic meshing + Abaqus `.inp` generation |
| `SpaX_GmshPeriodic.py` | (called by Standalone) | the periodic Gmsh mesher |
| `SpaX_PostProcess.py` | `abaqus python` (extraction) / `python3` (analysis) | ODB extraction, effective properties, anisotropy, study analyzers |

Everything else in the repo supports the sea-ice study campaign and is grouped
into folders, each with its own `README`:

| Folder | Contents |
|--------|----------|
| `studies/` | RVE parameter-deck generators (`make_ice_studies*.py`, deck-stamping helpers). |
| `params/` | The generated input parameter decks (`rve_*.csv`) — one row per RVE. |
| `hpc/` | Slurm batch scripts (`submit_*.sh`, `postprocess_*.sh`) — worked examples from the cluster this work ran on, kept for reference. |
| `results/` | Homogenisation output tables (`results_*.csv`), curves, and figures (`study_*.png`). |
| `analysis/` | Analyzers and field extractors that turn results into figures/quantities. |
| `viz/` | RVE visualization (`render_rve.py`, `odb_to_vtk.py`). |
| `tensors/` | Per-slice 6×6 elasticity tensors, one CSV per RVE. Earlier campaigns at the top level; each manuscript ensemble in its own subdirectory (`column/`, `basetensor_seeds/`, `bt80/`). |
| `docs/` | User guide (`USER_DOCS.md`) and the cluster runbook (`RUNBOOK.md`). |

## Reproducing the published figures

All analyzers are plain `python3` (matplotlib + pandas, no Abaqus licence) and
read the CSVs already in `results/`. Run them from `results/`, which is where
they look for their inputs:

```bash
cd results
PYTHONPATH=../analysis python3 ../analysis/make_rev_figs.py   # 3 figures, PDF+PNG
python3 ../analysis/analyze_brineK.py
python3 ../analysis/analyze_nlgeom.py
python3 ../analysis/macro_plate.py
```

| Figure | Produced by |
|--------|-------------|
| `ice_column_profiles`, `study_coltensor`, `study_scfdepth` | `analysis/make_rev_figs.py` |
| `study_brineK` | `analysis/analyze_brineK.py` |
| `study_nlgeom` | `analysis/analyze_nlgeom.py` |
| `study_macro_plate` | `analysis/macro_plate.py` |
| `rve_mesh`, `rve_meshcut`, `rve_straight_channels`, `rve_tilted_channels` | `viz/render_rve.py` (needs `pyvista` and a generated deck) |

Regenerating the *inputs* rather than the figures means re-solving the decks in
`params/` with Abaqus — see §1. The `hpc/` scripts are the Slurm jobs that did
so; they are archival, site-specific examples of how the campaigns were
submitted, not a dependency. See `docs/RUNBOOK.md` for the generic procedure.

> **Scatter convention.** Every quoted scatter is the population standard
> deviation (`ddof=0`) over the packings of an ensemble. Note that pandas
> `.std()` and `statistics.stdev` default to the *sample* convention, ~12%
> wider at five packings. Standard errors and Welch statistics correctly use
> `ddof=1` — see `check_channel_isotropy.py` and `compare_basetensor_sizes.py`.

---

## 1. Workflow at a glance

```
 parameters.csv ──(generate)──> Job-<id>-<mode>.inp ──(Abaqus solve)──> Job-<id>-<mode>.odb
                                                                              │
                                                          (post-process) ─────┘──> results.csv
                                                                              │
                                                                  (analyze) ──┴──> length scale / anisotropy / ...
```

1. **Generate** decks from a parameter CSV (one row = one RVE):
   ```bash
   python3 SpaX_Standalone.py parameters.csv  output_dir/
   ```
   Produces `Job-<run_id>-<mode>.inp` for each requested load case
   (`utx`, `uty`, `utz`, `ss12`, `ss13`, `ss23`, `ben`).

2. **Solve** each deck with Abaqus (your own SLURM array or locally):
   ```bash
   abaqus job=Job-<run_id>-utx cpus=4 interactive
   ```

3. **Post-process** the ODBs into `results.csv`:
   ```bash
   abaqus python SpaX_PostProcess.py parameters.csv output_dir/ results.csv
   # one RVE only (SLURM array task i):  ... results_i.csv  i
   # union per-task partials:            python3 SpaX_PostProcess.py --merge parts_dir/ results.csv
   ```

4. **Analyze** (plain `python3`, reads CSVs only — no Abaqus):
   ```bash
   python3 SpaX_PostProcess.py analyze eq19        results.csv
   python3 SpaX_PostProcess.py analyze lengthscale porous.csv [homog_baseline.csv ...]
   python3 SpaX_PostProcess.py analyze homog-calib
   python3 SpaX_PostProcess.py analyze hybrid      firstorder.csv bending.csv
   python3 SpaX_PostProcess.py analyze rve-study   results.csv
   ```

5. **Tensors** (Abaqus python): principal quantities / full 6×6 elasticity tensor
   ```bash
   abaqus python SpaX_PostProcess.py principal  <odb> <out.csv> [L] [run_id]
   abaqus python SpaX_PostProcess.py elasticity <odb_dir> <out.csv> <L> <run_id>
   ```

---

## 2. CSV parameters (one row per RVE)

The parameter CSV has a header row; each subsequent row defines one RVE. Columns
may appear in any order (read by name). The authoritative set:

### 2.1 Identity & geometry
| Column | Type | Example | Meaning |
|--------|------|---------|---------|
| `run_id` | str | `RVE_L48_s01` | unique name; all output files are `Job-<run_id>-<mode>.*` |
| `L` | float | `0.48` | RVE cube edge length (model units) |
| `L_mesh` | float | `0.033` | target mesh size. Bulk grows to ≈`3·L_mesh`; surface near inclusions ≈`0.4·L_mesh` (see `SPAX_LC_FINE_MULT`) |

### 2.2 Matrix material
| Column | Type | Example | Meaning |
|--------|------|---------|---------|
| `Is_Porous` | str | `Composite` | `Composite` (matrix + meshed inclusions), `Porous` (matrix + voids), or `Hybrid` (both) |
| `E_matrix` | float | `9.43e9` | matrix Young's modulus |
| `nu_matrix` | float | `0.33` | matrix Poisson ratio. Drives matrix element hybrid choice via `SPAX_HYBRID_NU` |

### 2.3 Inclusion population (packing)
| Column | Type | Example | Meaning |
|--------|------|---------|---------|
| `VoF_sphere` | float | `0.20` | target total inclusion volume fraction. `0.0` ⇒ homogeneous matrix cube |
| `r_avg` | float | `0.04` | mean inclusion radius (`d = 2·r_avg`) |
| `r_std` | float | `0.01` | radius standard deviation |
| `sphericity_avg` | float | `0.75` | mean sphericity (1 = sphere; <1 = ellipsoid) |
| `sphericity_std` | float | `0.1` | sphericity spread |
| `min_distance` | float | `0.002` | minimum surface separation between inclusions |
| `max_iterations` | int | `200000` | RSA placement attempt cap |
| `VoF_void_sphere` | float | `0.10` | fraction (of the box) assigned to **voids** (not meshed) |
| `VoF_incl_sphere` | float | `0.10` | fraction assigned to **soft inclusions** (meshed) |
| `E_sphere_inclusion` | float | `2.2e9` | inclusion Young's modulus (Solid inclusions) |
| `nu_sphere_inclusion` | float | `0.48` | inclusion Poisson ratio (Solid inclusions) |
| `Inclusion_Type` | str | `Liquid` | `Solid` (use E/nu above) or `Liquid` (use K/G below; brine) |
| `K_inclusion` | float | `2.2e9` | inclusion bulk modulus (Liquid only) |
| `G_inclusion` | float | `4.43e7` | inclusion shear modulus (Liquid only). E,ν derived from K,G; ν≈0.49 ⇒ hybrid element |

> **Void vs. soft inclusion.** `VoF_void_sphere` + `VoF_incl_sphere` partitions the
> inclusion population into non-meshed voids and meshed soft phase. For a pure
> soft-phase RVE set void=0; for a pure-porous RVE use `Is_Porous=Porous`.

### 2.4 Anisotropic growth (optional)
| Column | Type | Example | Meaning |
|--------|------|---------|---------|
| `Growth_Direction` | str | `Z` | preferred elongation axis (`X`/`Y`/`Z`/`Random`) during densification |
| `Growth_Concentration`| float | `0.5` | orientation concentration toward `Growth_Direction` (0 = isotropic) |

### 2.5 Brine channels (optional, vertical Z network)
| Column | Type | Example | Meaning |
|--------|------|---------|---------|
| `generate_channels` | str | `No` | `Yes` to add vertical (Z) channels (anisotropic, percolating) |
| `channel_vof_target`| float | `0.0` | target channel volume fraction (when channels on) |
| `r_channel_avg` | float | `0.02` | mean channel radius |
| `r_channel_std` | float | `0.005`| channel radius spread |

### 2.6 Load cases & boundary conditions
| Column | Type | Example | Meaning |
|--------|------|---------|---------|
| `Mode` | str | `Uniaxial Tension X` | primary first-order load case (see list below) |
| `Disp` | float | `0.0048` | applied displacement for `Mode` (engineering strain = `Disp/L`) |
| `Mode2` | str | `Simple Shear S13` | optional secondary load case (usually a shear for G) |
| `Disp2` | float | `0.0048` | applied displacement for `Mode2` |
| `full_tensor` | str | `No` | `Yes` ⇒ generate **all six** uniaxial+shear modes (needed for anisotropy) |
| `PBC_Method` | str | `Gmsh` | periodic-BC construction method |
| `nlgeom_flag` | str | `OFF` | Abaqus geometric nonlinearity (`OFF` for small-strain homogenization) |
| `Kappa` | float | `0.11` | imposed bending curvature. `>0` ⇒ also generate the bending (`ben`) deck |
| `Bending_Plane` | str | `xz` | bending plane: `xz` / `yz` (gradient through Z) or `xy` (gradient through Y) |
| `Bending_PBC_Type` | str | `Lesicar` | bending PBC formulation |

**`Mode` / `Mode2` values** (mapped to a deck suffix):

| String | Suffix | Drives |
|--------|--------|--------|
| `Uniaxial Tension X` | `utx` | E_x, ν_x |
| `Uniaxial Tension Y` | `uty` | E_y, ν_y |
| `Uniaxial Tension Z` | `utz` | E_z, ν_z |
| `Simple Shear S12` | `ss12` | G_xy |
| `Simple Shear S13` | `ss13` | G_xz |
| `Simple Shear S23` | `ss23` | G_yz |

> **For isotropic effective moduli** use `Mode=Uniaxial Tension X` + `Mode2=Simple
> Shear S13` (+ `Kappa>0` for bending). **For anisotropy** (e.g. transverse
> isotropy of vertical channels) set `full_tensor=Yes` to solve all six and obtain
> `E_x,E_y,E_z,G_xy,G_xz,G_yz`.

---

## 3. Environment variables

All are optional and **opt-in**; defaults reproduce the legacy (linear,
isotropic-report) behaviour. Set them before `SpaX_Standalone.py` (generation)
or `SpaX_PostProcess.py` (post) as appropriate.

### 3.1 Element order & quadratic (couple-stress) bending
| Variable | Default | Effect |
|----------|---------|--------|
| `SPAX_MESH_ORDER` | `1` | `2` ⇒ quadratic `C3D10(H)` (locking-free in bending). `1` ⇒ linear `C3D4(H)` |
| `SPAX_MIN_SICN` | `0.01` | minimum scaled-inverse-condition for order-2; below ⇒ re-pack (validity gate) |
| `SPAX_SLIVER_MULT_Q`| `1.0` | order-2 boundary-cap rejection floor (× `lc_fine`), applied from attempt 0 |
| `SPAX_OPT_PASSES` | `3` | Netgen quality passes before order promotion |
| `SPAX_HYBRID_NU` | `0.45` | ν threshold: phases with ν ≥ this get **hybrid** (H) elements, else non-hybrid |
| `SPAX_FORCE_HYBRID` | `0` | `1` ⇒ legacy all-hybrid (matrix + inclusion both H) |
| `SPAX_MESH_ALGO2D` | (Gmsh) | optional Gmsh 2-D mesher id (e.g. `6` = Frontal-Delaunay) |

### 3.2 Packing density & mesh-safe gaps
| Variable | Default | Effect |
|----------|---------|--------|
| `SPAX_LC_FINE_MULT` | `0.4` | fine surface size = this × `L_mesh` |
| `SPAX_FLOOR_MULT` | `0.75/LC_FINE_MULT` | min inclusion radius floor = this × `LC_FINE_MULT` × `L_mesh` |
| `SPAX_GAP_MULT` | `1.0` | extra min-distance widening (× `lc_fine`) for inter-inclusion gaps |
| `SPAX_GAP_REFINE` | `1` | `0` ⇒ disable mesh-in-gap local refinement (legacy widen-only) |
| `SPAX_GAP_RESOLVE` | `0.5` | gap-refinement resolution factor |
| `SPAX_SLIVER_MULT` | `0.5` | max sliver-rejection strength on retries (× `lc_fine`) |
| `SPAX_SLIVER_START` | `2` | failed attempts before sliver-rejection escalation begins |
| `SPAX_OFFAXIS` | `1` | `0` ⇒ disable off-axis ellipsoid sliver repair |
| `SPAX_OFFAXIS_FRAC` | `0.6` | off-axis true-gap floor (× `lc_fine`) |
| `SPAX_OFFAXIS_CHANNEL_FRAC` | `1.0` | channel↔inclusion gap floor (× `lc_fine`) |
| `SPAX_CHANNEL_SEP` | `1.0` | extra channel-side min-distance widening (× `lc_fine`) |
| `SPAX_ZGROW_BIAS` | `0.0` | anisotropic densification strength along `Growth_Direction` (w>0 lowers sphericity) |

### 3.3 Channel inclination (wavy channels)
| Variable | Default | Effect |
|----------|---------|--------|
| `SPAX_CHANNEL_TILT_DEG` | `0` | maximum lean off the Z axis, in degrees. `0` keeps the straight vertical cylinder. A channel is built as a *mean-vertical* wave of amplitude `amp = L·tan(tilt)/2π`, returning to the same XY position **and** tangent on both z-faces, so strict periodicity is preserved. It bulges laterally by up to `2·amp`, and the packing margin is inflated by `amp` to keep channels off the x/y faces |
| `SPAX_CHANNEL_TILT_AZIMUTH_DEG` | (per-channel) | direction of the lean. Unset ⇒ each channel gets its **own** azimuth from a deterministic hash of its centre, so the inclination has no net in-plane direction — this is what isolates *dilution* of the vertical anisotropy from an artificial x/y anisotropy. Set it to force every channel to lean the same way |
| `SPAX_TILT_DEBUG` | — | dump per-channel tilt geometry while meshing |

> A straight *net* tilt is impossible on a cubic periodic cell without sheared
> periodicity. The wavy form tests whether inclination **dilutes** the vertical
> anisotropy, not whether the anisotropy axis rotates.

### 3.4 Reproducibility, decoupling & runtime
| Variable | Default | Effect |
|----------|---------|--------|
| `SPAX_SEED` | (OS entropy) | reproducible packing (deterministic per row index) |
| `SPAX_SAVE_PACKING` | — | dir; freeze each packed sphere array to `<dir>/<run_id>.npy` |
| `SPAX_LOAD_PACKING` | — | dir; reuse a frozen packing (skip packer) — re-mesh the SAME geometry at any `L_mesh` |
| `SPAX_RESUME` | (off) | `1`/`yes`/`true`/`on` ⇒ skip rows whose decks already exist |
| `SPAX_GEN_WORKERS` | cores | parallel row workers (generation) |
| `SPAX_MAX_RETRIES` | `6` | mesh attempts per RVE before giving up |
| `SPAX_MESH_TIMEOUT` | `900` | per-RVE mesh wall-clock cap (seconds) |
| `SPAX_DEBUG_SURF` | — | debug: map a surface tag to its nearest inclusion |

> **Thread pinning** (recommended for parallel generation, to avoid
> oversubscription): `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`.

### 3.5 Analysis-side
| Variable | Default | Effect |
|----------|---------|--------|
| `SPAX_MC_PHI_DEG` | `30` | Mohr–Coulomb friction angle (degrees) used by `analysis/failure_extract.py` |

---

## 4. Output columns (`results.csv`)

Per RVE, depending on which load cases were solved:

| Column | From | Meaning |
|--------|------|---------|
| `E_eff`, `nu_eff`, `G_eff` | primary `Mode`/`Mode2` | scalar effective moduli (backward-compatible) |
| `E_x`,`E_y`,`E_z`,`nu_x`,`nu_y`,`nu_z` | `utx`/`uty`/`utz` | directional Young's moduli (full-tensor runs) |
| `G_xy`,`G_xz`,`G_yz` | `ss12`/`ss13`/`ss23` | directional shear moduli |
| `E_anisotropy` | — | max(E)/min(E) (1.0 = isotropic) |
| `E_z_over_xy` | — | transverse-isotropy index E_z / ½(E_x+E_y) |
| `D_rve` | `ben` | effective bending rigidity |
| `E_bending`, `porosity` | `ben` | bending modulus, void fraction |
| `D_classical`,`D_ratio`,`l`,`l_squared`,`E_bending_material` | derived | classical-vs-RVE bending, MCST length scale |
| `N_membrane`,`B_coupling` | `ben` | membrane resultant / membrane-bending coupling |

---

## 5. Worked examples

**Isotropic effective moduli + bending (linear), one RVE:**
```csv
run_id,L,L_mesh,Is_Porous,E_matrix,nu_matrix,VoF_sphere,r_avg,r_std,Mode,Disp,Mode2,Disp2,VoF_void_sphere,VoF_incl_sphere,E_sphere_inclusion,nu_sphere_inclusion,sphericity_avg,sphericity_std,min_distance,max_iterations,nlgeom_flag,PBC_Method,Kappa,Bending_Plane,Bending_PBC_Type,generate_channels,channel_vof_target,r_channel_avg,r_channel_std,Growth_Direction,Growth_Concentration,Inclusion_Type,K_inclusion,G_inclusion
RVE_a,0.48,0.033,Composite,9.43e9,0.33,0.20,0.04,0.01,Uniaxial Tension X,0.0048,Simple Shear S13,0.0048,0.10,0.10,2.2e9,0.48,0.75,0.1,0.002,200000,OFF,Gmsh,0.11,xz,Lesicar,No,0,0.02,0.005,Z,0.5,Liquid,2.2e9,4.43e7
```
```bash
python3 SpaX_Standalone.py params.csv out/
```

**Quadratic (couple-stress) bending, locking-free:**
```bash
SPAX_MESH_ORDER=2 python3 SpaX_Standalone.py params.csv out/
```

**Full anisotropy (all six modes) + the 6×6 elasticity tensor** — set
`full_tensor=Yes` in the CSV, then:
```bash
SPAX_MESH_ORDER=2 python3 SpaX_Standalone.py params.csv out/
# ... solve all six Job-RVE_a-{utx,uty,utz,ss12,ss13,ss23}.odb ...
abaqus python SpaX_PostProcess.py params.csv out/ results.csv      # E_x,E_y,E_z,...,E_anisotropy
abaqus python SpaX_PostProcess.py elasticity out/ Cij.csv 0.48 RVE_a
```

**Mesh-convergence on a FIXED geometry** (freeze once, re-mesh at several sizes):
```bash
SPAX_SAVE_PACKING=pack/ SPAX_MESH_ORDER=2 python3 SpaX_Standalone.py one_row.csv m_033/   # freeze
for lm in 0.045 0.026 0.020; do
  sed "s/,0.033,/,$lm,/" one_row.csv > row_$lm.csv
  SPAX_LOAD_PACKING=pack/ SPAX_MESH_ORDER=2 python3 SpaX_Standalone.py row_$lm.csv m_$lm/
done
```

---

## 6. Notes & gotchas

- **`VoF_sphere=0`** ⇒ a homogeneous matrix cube (no inclusions) — useful as a
  couple-stress *calibration* baseline (true length scale = 0).
- **Quadratic meshing** of porous RVEs may re-pack a few times (the `SPAX_MIN_SICN`
  validity gate rejects inverted sliver elements); this is normal and converges.
- **Coarsening `L_mesh` is not a reliable speed-up** for porous RVEs: the
  mesh-in-gap refinement triggers on a threshold ∝ `lc_fine`, so a coarse base
  mesh flags *more* gaps and can add elements back. Use `SPAX_FORCE_HYBRID=0`
  (default, non-hybrid matrix) for the clean ~13% solve saving instead.
- **Run analyzers with plain `python3`** (CSV-only, no Abaqus license); run the
  ODB extraction (`<csv> <dir>`, `principal`, `elasticity`) with `abaqus python`.
