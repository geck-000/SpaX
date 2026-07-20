# viz/ — RVE visualization

Tools to turn meshes and solved models into publication figures.

| File | Run with | Role |
|------|----------|------|
| `render_rve.py` | `python3` (PyVista) | Offscreen renders. `micro` mode: brine channels (blue) vs pockets (orange) in a faint ice cube from a `Job-*.inp`. `field` mode: von Mises on the warped, clipped RVE from a `.vtk`. Okabe-Ito colourblind-safe palette. |
| `odb_to_vtk.py` | `abaqus python` | Exports a solved `.odb` to legacy VTK (displacement vectors, von Mises, material tag) for ParaView or for `render_rve.py field`. |

**Usage**

```bash
# microstructure figure straight from a generated deck
python3 render_rve.py micro Job-RVE_a-utx.inp rve.png

# stress field: first export the ODB (needs Abaqus), then render
abaqus python odb_to_vtk.py Job-RVE_a-utx.odb rve.vtk
python3 render_rve.py field rve.vtk field.png
```

Requires `pyvista` (`pip install pyvista`); the ODB export requires an Abaqus
installation.
