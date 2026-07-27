# viz/ — RVE visualization

Tools to turn meshes and solved models into publication figures.

| File | Run with | Role |
|------|----------|------|
| `render_rve.py` | `python3` (PyVista) | Offscreen renders. `micro` mode: brine channels (blue) vs pockets (orange) in a faint ice cube from a `Job-*.inp`. `mesh` / `meshcut` modes: the tetrahedral mesh itself, intact or with one corner octant of whole elements removed. `field` mode: von Mises on the warped, clipped RVE from a `.vtk`. Okabe-Ito colourblind-safe palette. |
| `odb_to_vtk.py` | `abaqus python` | Exports a solved `.odb` to legacy VTK (displacement vectors, von Mises, material tag) for ParaView or for `render_rve.py field`. |

**Usage**

```bash
# microstructure figure straight from a generated deck
python3 render_rve.py micro Job-RVE_a-utx.inp rve.png

# the mesh itself, for the periodic-BC figure of the paper (Fig. 7):
# intact cell, then the same cell cut open. --nolegend/--noaxes drop the
# pyvista overlays so the labels can be set in TikZ; output is auto-trimmed.
python3 render_rve.py mesh    out_colseeds/Job-CSEED_z95_s1-utx.inp rve_mesh.png    20 12 --nolegend --noaxes
python3 render_rve.py meshcut out_colseeds/Job-CSEED_z95_s1-utx.inp rve_meshcut.png 20 12 --nolegend --noaxes

# stress field: first export the ODB (needs Abaqus), then render
abaqus python odb_to_vtk.py Job-RVE_a-utx.odb rve.vtk
python3 render_rve.py field rve.vtk field.png
```

Requires `pyvista` (`pip install pyvista`); the ODB export requires an Abaqus
installation.
