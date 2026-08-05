"""Export an Abaqus ODB to a legacy VTK UnstructuredGrid for ParaView.

Only Abaqus can read a .odb, so run this under `abaqus python odb_to_vtk.py
<odb> <out.vtk> [step] [frame]`. It writes an ASCII legacy .vtk that ParaView
opens directly, carrying:
  * POINT_DATA  U           (displacement vectors -> Filters > Warp By Vector)
  * CELL_DATA   mises       (von Mises stress, element-averaged)
  * CELL_DATA   material    (0 = matrix, 1 = inclusion/brine -> Threshold/colour)

Tetrahedra only (C3D4 / C3D4H, VTK type 10), which is what the RVE meshes use.
"""
import sys
from odbAccess import openOdb


def main():
    odb_path = sys.argv[1]
    out_path = sys.argv[2]
    step_name = sys.argv[3] if len(sys.argv) > 3 else None
    frame_idx = int(sys.argv[4]) if len(sys.argv) > 4 else -1

    odb = openOdb(path=odb_path, readOnly=True)
    step = odb.steps[step_name] if step_name else odb.steps[list(odb.steps.keys())[0]]
    frame = step.frames[frame_idx]

    # Single solid instance (exclude PBC/assembly helper instances).
    inst = None
    for name, ins in odb.rootAssembly.instances.items():
        if 'PBC' not in name.upper() and name.upper() != 'ASSEMBLY':
            inst = ins
            break
    if inst is None:
        inst = list(odb.rootAssembly.instances.values())[0]

    # ---- points (node label -> contiguous index) ----
    labels, coords, idx_of = [], [], {}
    for n in inst.nodes:
        idx_of[n.label] = len(labels)
        labels.append(n.label)
        c = list(n.coordinates) + [0.0, 0.0, 0.0]
        coords.append((c[0], c[1], c[2]))

    # ---- tetra cells + element label -> index ----
    cells, elem_labels, eidx_of = [], [], {}
    for e in inst.elements:
        if len(e.connectivity) != 4:
            continue
        eidx_of[e.label] = len(cells)
        cells.append([idx_of[c] for c in e.connectivity])
        elem_labels.append(e.label)

    # ---- material tag per element from element sets ----
    material = [0] * len(cells)
    for sname, sset in inst.elementSets.items():
        is_incl = ('SPHERE' in sname.upper() or 'INCL' in sname.upper()
                   or 'BRINE' in sname.upper())
        if not is_incl:
            continue
        for e in sset.elements:
            if e.label in eidx_of:
                material[eidx_of[e.label]] = 1

    # ---- displacement U at nodes ----
    U = [(0.0, 0.0, 0.0)] * len(labels)
    try:
        ufield = frame.fieldOutputs['U'].getSubset(region=inst)
        for v in ufield.values:
            if v.nodeLabel in idx_of:
                d = list(v.data) + [0.0, 0.0, 0.0]
                U[idx_of[v.nodeLabel]] = (d[0], d[1], d[2])
    except Exception as ex:
        sys.stderr.write('U export skipped: %s\n' % ex)

    # ---- von Mises per element (average over integration points) ----
    mises = [0.0] * len(cells)
    cnt = [0] * len(cells)
    try:
        sfield = frame.fieldOutputs['S'].getSubset(region=inst)
        for v in sfield.values:
            if v.elementLabel in eidx_of:
                j = eidx_of[v.elementLabel]
                mises[j] += float(v.mises)
                cnt[j] += 1
        for j in range(len(mises)):
            if cnt[j]:
                mises[j] /= cnt[j]
    except Exception as ex:
        sys.stderr.write('S export skipped: %s\n' % ex)

    odb.close()

    # ---- write legacy VTK ASCII ----
    with open(out_path, 'w') as f:
        f.write('# vtk DataFile Version 3.0\n')
        f.write('SpaX RVE %s frame %d\n' % (odb_path, frame_idx))
        f.write('ASCII\nDATASET UNSTRUCTURED_GRID\n')
        f.write('POINTS %d float\n' % len(coords))
        for x, y, z in coords:
            f.write('%g %g %g\n' % (x, y, z))
        f.write('CELLS %d %d\n' % (len(cells), 5 * len(cells)))
        for c in cells:
            f.write('4 %d %d %d %d\n' % (c[0], c[1], c[2], c[3]))
        f.write('CELL_TYPES %d\n' % len(cells))
        for _ in cells:
            f.write('10\n')
        f.write('POINT_DATA %d\n' % len(coords))
        f.write('VECTORS U float\n')
        for u in U:
            f.write('%g %g %g\n' % (u[0], u[1], u[2]))
        f.write('CELL_DATA %d\n' % len(cells))
        f.write('SCALARS mises float 1\nLOOKUP_TABLE default\n')
        for m in mises:
            f.write('%g\n' % m)
        f.write('SCALARS material int 1\nLOOKUP_TABLE default\n')
        for m in material:
            f.write('%d\n' % m)
    sys.stderr.write('wrote %s: %d nodes, %d tets (%d inclusion)\n'
                     % (out_path, len(coords), len(cells), sum(material)))


if __name__ == '__main__':
    main()
