# Sea-ice column study — depth-resolved RVE homogenization

`parametric_sea_ice_column.csv` is a stack of 10 RVEs, one per horizontal slice
through a ~0.5–1.0 m first-year (FY) sea-ice sheet in late-winter conditions.
Each slice carries its **local temperature**, which sets the **brine volume
fraction** and the **brine-inclusion morphology** observed in X-ray micro-CT.
Solving the stack gives the effective stiffness *as a function of depth*
`E(z), G(z)` and, because `full_tensor=Yes`, the vertical-vs-horizontal
anisotropy `E_z(z)/E_xy(z)` — the property that actually changes down the column.

```
   air / snow   T ≈ -20 °C   ── cold, low brine, isolated near-spherical pockets
   ┌───────────────────────┐  z05  φb≈0.022   spher 0.85   no channels
   │  ICE_z05              │
   │  ICE_z15 ... z45     │  interior, low salinity (C-shape minimum)
   │  ICE_z55             │
   │  ICE_z65             │
   │  ICE_z75              │  ─ rule-of-fives percolation (φb ≈ 0.05) ─
   │  ICE_z85  channels Yes│  z85  φb≈0.068   spher 0.68   tubes + Z channels
   │  ICE_z95  channels Yes│  z95  φb≈0.150   spher 0.62   skeletal layer
   └───────────────────────┘
   ocean        T ≈ -1.8 °C  ── warm, brine-rich, connected vertical channels
```

## Per-slice parameters and their provenance

| run_id | depth z/H | T (°C) | S (ppt) | φ_brine | φ_gas (void) | sphericity | channels |
|--------|-----------|--------|---------|---------|--------------|-----------|----------|
| ICE_z05 | 0.05 | -19.1 | 7.0 | 0.022 | 0.020 | 0.85 | No |
| ICE_z15 | 0.15 | -17.3 | 5.5 | 0.019 | 0.015 | 0.84 | No |
| ICE_z25 | 0.25 | -15.4 | 4.8 | 0.018 | 0.012 | 0.83 | No |
| ICE_z35 | 0.35 | -13.6 | 4.5 | 0.019 | 0.010 | 0.82 | No |
| ICE_z45 | 0.45 | -11.8 | 4.3 | 0.020 | 0.010 | 0.80 | No |
| ICE_z55 | 0.55 | -10.0 | 4.3 | 0.023 | 0.010 | 0.78 | No |
| ICE_z65 | 0.65 |  -8.2 | 4.5 | 0.030 | 0.010 | 0.75 | No |
| ICE_z75 | 0.75 |  -6.3 | 5.0 | 0.041 | 0.012 | 0.72 | No |
| ICE_z85 | 0.85 |  -4.5 | 6.0 | 0.068* | 0.015 | 0.68 | Yes (0.018) |
| ICE_z95 | 0.95 |  -2.7 | 8.0 | 0.150* | 0.020 | 0.62 | Yes (0.060) |

\* For the two channel slices the total brine fraction is split between meshed
pockets (`VoF_incl_sphere`) and the percolating vertical network
(`channel_vof_target`): z85 = 0.050 + 0.018, z95 = 0.090 + 0.060.

### Temperature profile
Linear conductive winter profile, cold top to ocean-pinned bottom:
`T(z) = T_top + (T_bot - T_top)·(z/H)`, `T_top = -20 °C`, `T_bot = -1.8 °C`
(seawater freezing point at S ≈ 34). Values above are slice mid-points.

### Salinity profile
Classic FY-ice **C-shape** (Cox & Weeks 1983; Eicken 1992): elevated at the top
(initial entrapment) and bottom (skeletal layer), with a ~4.3 ppt interior
minimum. Bulk-salinity values, not brine salinity.

### Brine volume fraction — Frankenstein & Garner (1967)
`φb = S · (-49.185/T + 0.532) / 1000`, valid -0.5 ≥ T ≥ -22.9 °C, T in °C,
S in ppt. This is the standard low-order phase relation; Cox & Weeks (1983) /
Leppäranta & Manninen (1988) refine it and add the gas term used below. Note the
sharp rise toward the warm base — the physical driver of the soft skeletal layer.

### Gas (air) inclusions — `VoF_void_sphere`
Drained/entrapped air, 1–2 % by volume, slightly higher near the surface where
bubbles are retained (Cox & Weeks; Light et al. 2003 micro-CT). Modelled as true
voids (not meshed), distinct from the brine soft phase.

### Percolation — the "rule of fives" (Golden, Ackley & Lytle, *Science* 1998)
Brine becomes connected and fluid-permeable above φb ≈ 0.05 (≈ -5 °C, S ≈ 5).
The CSV crosses this between z75 (0.041, isolated) and z85 (0.068, connected), so
`generate_channels` switches **on** at z85/z95 — a microstructural transition, not
a tuning knob.

### Inclusion morphology — X-ray micro-CT / tomography
- **Shape vs depth:** cold ice holds near-spherical/short-ellipsoidal pockets
  (sphericity ~0.85); toward the warm base inclusions merge into vertical
  tubes/lamellae (sphericity ~0.6). Trend from Light et al. (2003), Cole et al.,
  Maus et al. (2015, 2021), Lieb-Lappen et al. (2017), Crabeck et al. (2016).
- **Orientation:** columnar ice has vertical c-axis brine layers ⇒
  `Growth_Direction=Z`, with `Growth_Concentration` rising 0.40→0.70 down-column
  as the columnar texture sharpens.
- **Size:** `r_avg` 0.030→0.045 and channel radii 0.020–0.025 are in *model
  units* (RVE edge `L=0.50`). Classical effective moduli are scale-invariant, so
  the absolute unit is free; map `L=0.50` to a physical RVE of a few mm
  (sub-mm pockets, ~1 mm channels — micro-CT range) when reporting.

### Pure-ice matrix
`E_matrix(T) = 9.36e9 + 0.0129e9·(-2 - T)` (ice stiffens ~13 MPa/°C as it cools,
matching the existing `T_m*` rows), `ν = 0.33`. Brine is the established K/G
soft-solid: `K = 2.2 GPa`, `G = 0.44 MPa` (ν ≈ 0.4999 near-fluid plateau) — see
the `fluid-cavity-vs-soft-liquid` note for why this beats `*FLUID CAVITY` in a
periodic RVE.

## Running it

```bash
# generate all decks (full_tensor=Yes => 6 modes per slice; ~60 solves)
SPAX_SEED=0 SPAX_MESH_ORDER=2 python3 Spatium_Standalone.py \
    parametric_sea_ice_column.csv out_column/

# solve each Job-ICE_zNN-{utx,uty,utz,ss12,ss13,ss23}.odb with Abaqus, then:
abaqus python Spatium_PostProcess.py parametric_sea_ice_column.csv out_column/ results_column.csv

# anisotropy / transverse-isotropy profile vs depth
python3 Spatium_PostProcess.py analyze rve-study results_column.csv
```

`results_column.csv` then has `E_x,E_y,E_z,G_xy,G_xz,G_yz,E_anisotropy,
E_z_over_xy` per slice — i.e. the depth profile. Add replicate seeds
(`SPAX_SEED=1,2,…`) per slice for confidence intervals; `Kappa=0` here because
the bending size-effect was shown to be non-microstructural (no MCST length
scale — see `bending-length-scale-artifact`), so first-order anisotropy is the
meaningful output.

## Other studies this toolkit supports (proposed)

1. **Temperature-gradient column (this CSV).** Effective `E(z), G(z),
   E_z/E_xy(z)` through the sheet; feeds a graded / laminated-plate macro model
   of the whole ice cover.
2. **Percolation transition sweep.** Refine φb across 0.03–0.10 with channels
   on/off to locate the mechanical-anisotropy onset at the rule-of-fives
   threshold (pairs naturally with fluid permeability).
3. **Morphology study from CT.** At fixed φb, sweep sphericity (0.6→0.9) and
   `Growth_Concentration` to isolate how pocket shape/orientation — the directly
   tomography-measured quantities — control `E_z/E_xy`.
4. **FY vs multi-year ice.** C-shape vs desalinated low-salinity profile;
   raise `VoF_void_sphere` for gas-rich summer/drained ice (`Is_Porous=Hybrid`).
5. **Full transverse-isotropy tensor of columnar ice.** `full_tensor=Yes` 6×6
   `C_ij` (`elasticity` subcommand) vs depth/temperature.
6. **Seasonal sweep.** Re-run the column at several `T_top` (winter→spring) to
   trace how the warm-base soft layer thickens as the sheet warms.
