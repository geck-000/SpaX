# Fluid-cavity brine model — proof of concept

A minimal, **non-periodic** test to get Abaqus `*FLUID CAVITY` working for a brine
pocket before wiring it into the periodic RVE batch. One centred interior sphere
in a cube, so the cavity surface is naturally **closed** — which is the single
biggest reason fluid cavities fail inside a periodic RVE (pockets get cut by the
cell boundary and the surface is no longer watertight).

## Why bother

The production `Inclusion_Type=Liquid` model turns the brine into a soft
isotropic solid via K/G → E,ν, giving **E ≈ 0.13 GPa, ν ≈ 0.49**. Because E is so
tiny, the high ν barely couples — the inclusion behaves almost like a void, and
the homogenised `ν_eff` stays flat at ~0.33 instead of rising with brine content
(the old kernel reports showed ν climbing to ~0.37). A fluid cavity restores the
**incompressible volumetric coupling** (brine resists volume change, carries zero
shear) for the right physical reason. This POC checks that the mechanics work and
that `ν_eff` actually moves.

## Files

| File | Run with | Purpose |
|------|----------|---------|
| `poc_gen.py` | `python3` (+ gmsh) | meshes the cube+sphere, writes both decks |
| `poc_fluidcavity.inp` | Abaqus | sphere = `*FLUID CAVITY` (compressible brine) |
| `poc_void.inp` | Abaqus | identical mesh, sphere = open pore (baseline) |
| `poc_extract.py` | `abaqus python` | prints E_eff, ν_eff, cavity pressure/volume |

## Run (on CSC)

```bash
python3 poc_gen.py                       # regenerate decks (only if you change params)
abaqus job=poc_void       interactive
abaqus job=poc_fluidcavity interactive
abaqus python poc_extract.py poc_void.odb
abaqus python poc_extract.py poc_fluidcavity.odb
```

## What to look for

- **Convergence.** The fluid-cavity job should converge in one increment (linear).
  If it dies, the message tells you which failure mode you hit:
  - *"cavity is not closed" / negative cavity volume* → surface orientation or an
    open surface (not an issue here; the sphere is interior).
  - *unconstrained DOF on the reference node* → the `CAVREF, 1, 3` pin is missing.
- **ν_eff(fluid cavity) > ν_eff(void).** The brine resists lateral contraction, so
  Poisson's ratio rises toward the matrix value. Raise `K_brine` in `poc_gen.py`
  toward incompressible and ν_eff should climb further. This is the effect the
  soft-solid Liquid model cannot produce.
- **PCAV > 0** under tension: the cavity actually pressurised (volume tried to grow,
  the compressible fluid pushed back).

## Key implementation details (the parts that usually trip people up)

1. **Closed surface.** Only works because the sphere is fully interior. In the
   periodic RVE you must either keep brine pockets off the cell boundary or cap
   boundary-cut pockets — see "Path to production" below.
2. **Surface normal.** `*Surface, type=ELEMENT` on the *matrix* element faces: the
   element's outward face normal points **into** the cavity (matrix is outside the
   sphere), which is exactly what `*FLUID CAVITY` requires. No flip needed.
3. **One cavity per pocket.** Each pocket is hydraulically isolated → its own
   closed surface + its own reference node + its own `*FLUID CAVITY`. Never lump
   pockets into one cavity.
4. **Compressible hydraulic fluid.** `*Fluid Behavior` + `*Fluid Density` +
   `*Fluid Bulk Modulus`. Omitting the bulk modulus makes the fluid fully
   incompressible and can volumetrically lock the cell.
5. **Reference node** is a free node (not in any element); pin its displacement
   DOFs (`CAVREF, 1, 3`) — only its pressure DOF matters.

## Path to production (periodic RVE)

The batch pipeline already has the hard parts: voids are meshed as unmeshed empty
cavities (`Spatium_GmshPeriodic.py`, "Void volumes get NO physical group"), and
the matrix-side interface surface tags are already collected (`incl_surface_tags`).
To add an opt-in `Inclusion_Type=FluidCavity` mode:

1. In `Spatium_GmshPeriodic.py`: treat brine pockets like voids (don't mesh the
   interior), but tag **each pocket separately** so its wall facets can be emitted
   as one closed surface. Reject (or flag) pockets that touch the cell boundary.
2. In `Spatium_Standalone.py`: per pocket, write `*Surface, type=ELEMENT`
   (matrix faces) + a cavity ref node + `*Fluid Cavity`/`*Fluid Behavior`. Run it
   **alongside** the existing Liquid path, not replacing it.
3. Periodicity caveat: a pocket split across the boundary is two open half-cavities
   that physically share one pressure. Simplest robust route — constrain brine
   pockets to be strictly interior for this mode; let air voids stay boundary-cut.
