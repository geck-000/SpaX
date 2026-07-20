# -*- coding: utf-8 -*-
"""
SpaX_GmshPeriodic.py — Gmsh-Based Periodic Mesh for RVE
=============================================================

Generates a periodic tetrahedral mesh using Gmsh, with guaranteed
matching nodes on opposite faces. Eliminates the PBC_Surface + Tie
approach entirely.

Pipeline:
  1. Receive sphere array from SpaX_Kernel (positions, radii)
  2. Build Gmsh geometry: RVE box + spheres (with periodic copies)
  3. Boolean fragment to create matrix/inclusion volumes
  4. Apply periodic mesh constraints on opposite faces
  5. Generate mesh and export as .inp
  6. Import into Abaqus as orphan mesh parts
  7. Apply direct node-to-node PBC equations

Requirements:
  - Gmsh installed on the system (gmsh executable in PATH)
  - This script runs OUTSIDE Abaqus (system Python 3 with gmsh package)
  - The generated .inp is then imported into Abaqus

Usage:
  From system Python:
    python SpaX_GmshPeriodic.py

  From Abaqus (import the generated mesh):
    import SpaX_GmshPeriodic_Import as GI
    GI.import_and_apply_pbc(inp_path, L=0.2, Disp=0.002)
"""

import numpy as np
import os
import sys
import math

try:
    import gmsh
    HAS_GMSH = True
except ImportError:
    HAS_GMSH = False


# =====================================================================
# Geometry construction
# =====================================================================

def _add_periodic_sphere_copies(cx, cy, cz, r, L, tol=1e-10):
    """
    For a sphere at (cx,cy,cz) with radius r in an RVE [0,L]^3,
    determine which periodic copies are needed and return their centres.
    
    A copy is needed in direction +d if the sphere crosses the d=0 face,
    and in -d if it crosses the d=L face. Corner/edge copies are also
    generated for multi-axis crossings.
    """
    from itertools import product
    
    shifts_x = [0]
    shifts_y = [0]
    shifts_z = [0]
    
    if cx - r < tol:
        shifts_x.append(L)
    if cx + r > L - tol:
        shifts_x.append(-L)
    if cy - r < tol:
        shifts_y.append(L)
    if cy + r > L - tol:
        shifts_y.append(-L)
    if cz - r < tol:
        shifts_z.append(L)
    if cz + r > L - tol:
        shifts_z.append(-L)
    
    copies = []
    for sx, sy, sz in product(shifts_x, shifts_y, shifts_z):
        if sx == 0 and sy == 0 and sz == 0:
            continue  # skip original
        copies.append((cx + sx, cy + sy, cz + sz))

    return copies


def _add_periodic_cylinder_copies(cx, cy, r, L, tol=1e-10):
    """Periodic XY copies for a vertical (Z) channel of radius r at (cx,cy).
    A channel spans the full RVE height, so it only needs XY tiling (never a
    Z copy): a copy in +x if it crosses x=0, in -x if it crosses x=L, likewise
    in y, plus the corner combinations."""
    from itertools import product
    shifts_x = [0]
    shifts_y = [0]
    if cx - r < tol:
        shifts_x.append(L)
    if cx + r > L - tol:
        shifts_x.append(-L)
    if cy - r < tol:
        shifts_y.append(L)
    if cy + r > L - tol:
        shifts_y.append(-L)
    return [(cx + sx, cy + sy) for sx, sy in product(shifts_x, shifts_y)
            if not (sx == 0 and sy == 0)]


def _match_periodic_surfaces(L, eps=None):
    """
    Match each surface on a min face (x=0, y=0 or z=0) to its unique
    translate on the opposite max face, using Gmsh's spatial bounding-box
    index. This is the canonical, robust approach from Gmsh's own
    periodic.py demo and correctly handles faces that the Boolean fragment
    has split into several sub-surfaces.

    For every sub-surface on a min face (the MASTER) we translate its full
    3D bounding box by +L along the axis and query getEntitiesInBoundingBox
    to find the matching sub-surface on the max face (the SLAVE). Master is
    always the min face so the periodic constraint graph is consistent across
    all three axes (shared edges/corners compose correctly).

    Returns list of (master_tag, slave_tag, [tx, ty, tz]).
    """
    if eps is None:
        eps = 1e-4 * L
    pairs = []
    used_slaves = set()

    for axis, trans in [(0, [L, 0, 0]), (1, [0, L, 0]), (2, [0, 0, L])]:
        # Slab that selects all surfaces lying on the min face of this axis.
        lo = [-eps, -eps, -eps]
        hi = [L + eps, L + eps, L + eps]
        hi[axis] = eps  # restrict to coord in [-eps, eps] on this axis
        smin = gmsh.model.getEntitiesInBoundingBox(
            lo[0], lo[1], lo[2], hi[0], hi[1], hi[2], 2)

        for dm, s in smin:
            bb = gmsh.model.getBoundingBox(dm, s)
            # Bounding box of the expected translate on the max face.
            tb = [bb[0] + trans[0], bb[1] + trans[1], bb[2] + trans[2],
                  bb[3] + trans[0], bb[4] + trans[1], bb[5] + trans[2]]
            cand = gmsh.model.getEntitiesInBoundingBox(
                tb[0] - eps, tb[1] - eps, tb[2] - eps,
                tb[3] + eps, tb[4] + eps, tb[5] + eps, 2)
            for dm2, s2 in cand:
                if s2 in used_slaves:
                    continue
                bb2 = gmsh.model.getBoundingBox(dm2, s2)
                if all(abs(bb2[k] - tb[k]) < eps for k in range(6)):
                    pairs.append((s, s2, list(trans)))
                    used_slaves.add(s2)
                    break

    return pairs


def _periodic_curve_pairs(L, eps=1e-6):
    """Single-axis periodic curve correspondence for the RVE boundary.

    A curve lying on one or more *max* faces is slaved, along the LOWEST such
    axis, to its -L translate. This builds a consistent master tree rooted at
    the all-min curves, so a cube-edge curve (shared by two faces) is never
    double-assigned. Matching is by centre-of-mass translate (unambiguous for
    distinct curves, unlike a bounding box which two arcs can share).

    This is the foundation of strict periodicity here: we drive 1D periodicity
    ourselves with these pairs (gmsh.model.mesh.setPeriodic(1, ...)) and then
    copy the 2D face meshes manually (see _manual_periodic_surface_copy).

    Why not Gmsh's own setPeriodic(2)? Its surface-derived curve correspondence
    mis-matches cube-edge curves once a face is split into many sub-surfaces by
    inclusion cuts (it maps an edge to the wrong parallel edge). That internal
    map drives the 2D copy and cannot be overridden through the API, so the copy
    fails with "curve X has N vertices, whereas correspondant has M". Driving 1D
    ourselves + a manual 2D copy sidesteps it entirely.

    Returns list of (master_curve, slave_curve, [tx, ty, tz]).
    """
    info = {}
    for dim, c in gmsh.model.getEntities(1):
        bb = gmsh.model.getBoundingBox(1, c)
        com = gmsh.model.occ.getCenterOfMass(1, c)
        onmax = [a for a in range(3)
                 if abs(bb[a] - L) < eps and abs(bb[a + 3] - L) < eps]
        info[c] = (com, onmax)

    def key(p):
        return (round(p[0], 7), round(p[1], 7), round(p[2], 7))

    idx = {key(v[0]): c for c, v in info.items()}
    pairs = []
    for c, (com, onmax) in info.items():
        if not onmax:
            continue
        axis = min(onmax)
        trans = [0.0, 0.0, 0.0]
        trans[axis] = L
        m = idx.get(key([com[0] - trans[0], com[1] - trans[1],
                         com[2] - trans[2]]))
        if m is not None and m != c:
            pairs.append((m, c, trans))
    return pairs


def _manual_periodic_surface_copy(L, face_pairs):
    """Replace each slave (max-face) sub-surface mesh with the exact translate
    of its master (min-face) sub-surface mesh, reusing the already-periodic
    boundary-curve nodes.

    Prerequisite: the 1D mesh must already be periodic (set the pairs from
    _periodic_curve_pairs via setPeriodic(1) and run generate(2) first) so the
    slave's bounding-curve nodes are exact translates of the master's. We then
    create new interior nodes for the slave as translated copies of the master
    interior nodes, map the master boundary nodes onto the existing slave
    boundary nodes, clear the slave's old surface mesh, and re-add the master
    triangulation (winding reversed for the opposite outward normal). After this
    the opposite RVE faces are conforming node-for-node, so the 3D mesh and the
    exported periodic node pairs are strictly periodic.

    Returns (n_pairs, n_missing_boundary_nodes); n_missing must be 0.
    """
    def key(p):
        return (round(p[0], 7), round(p[1], 7), round(p[2], 7))

    all_nt, _, _ = gmsh.model.mesh.getNodes()
    next_tag = int(max(all_nt)) + 1 if len(all_nt) else 1
    max_et = 0
    _, etags_all, _ = gmsh.model.mesh.getElements()
    for a in etags_all:
        if len(a):
            max_et = max(max_et, int(max(a)))
    next_elem = max_et + 1

    total_missing = 0
    for master, slave, trans in face_pairs:
        # Index slave boundary-curve nodes by coordinate.
        snode = {}
        for d, cv in gmsh.model.getBoundary([(2, slave)], oriented=False):
            nt, co, _ = gmsh.model.mesh.getNodes(1, cv, includeBoundary=True)
            for i, t in enumerate(nt):
                snode[key(co[3 * i:3 * i + 3])] = int(t)

        # Master interior nodes -> new translated slave interior nodes.
        mnt, mco, _ = gmsh.model.mesh.getNodes(2, master, includeBoundary=False)
        nodemap = {}
        ntags = []
        ncoords = []
        for i, t in enumerate(mnt):
            p = mco[3 * i:3 * i + 3]
            nodemap[int(t)] = next_tag
            ntags.append(next_tag)
            ncoords.extend([p[0] + trans[0], p[1] + trans[1], p[2] + trans[2]])
            next_tag += 1

        # Master boundary nodes -> existing slave boundary nodes.
        _dbg = os.environ.get('SPAX_TILT_DEBUG')
        for d, cv in gmsh.model.getBoundary([(2, master)], oriented=False):
            nt, co, _ = gmsh.model.mesh.getNodes(1, cv, includeBoundary=True)
            for i, t in enumerate(nt):
                if int(t) in nodemap:
                    continue
                p = co[3 * i:3 * i + 3]
                k = key([p[0] + trans[0], p[1] + trans[1], p[2] + trans[2]])
                if k in snode:
                    nodemap[int(t)] = snode[k]
                else:
                    total_missing += 1
                    if _dbg:
                        print("    [tiltdbg] MISS master face {} curve {}: node {} "
                              "at ({:.4f},{:.4f},{:.4f}) -> translate ({:.4f},{:.4f},"
                              "{:.4f}) not on slave {}".format(
                                  master, cv, int(t), p[0], p[1], p[2],
                                  p[0] + trans[0], p[1] + trans[1], p[2] + trans[2],
                                  slave))

        # Master triangles -> slave, winding reversed (opposite outward normal).
        etypes, etags, enodes = gmsh.model.mesh.getElements(2, master)
        gmsh.model.mesh.clear([(2, slave)])
        if ntags:
            gmsh.model.mesh.addNodes(2, slave, ntags, ncoords)
        new_etags = []
        new_enodes = []
        for ti, et in enumerate(etypes):
            nn = gmsh.model.mesh.getElementProperties(et)[3]
            etl = []
            enl = []
            for e in range(len(etags[ti])):
                tri = [int(enodes[ti][e * nn + k]) for k in range(nn)]
                if any(x not in nodemap for x in tri):
                    continue          # skip unmapped (reported via total_missing)
                mapped = [nodemap[x] for x in tri][::-1]
                etl.append(next_elem)
                next_elem += 1
                enl.extend(mapped)
            new_etags.append(etl)
            new_enodes.append(enl)
        gmsh.model.mesh.addElements(2, slave, etypes, new_etags, new_enodes)

    return len(face_pairs), total_missing


def _manual_periodic_curve_copy(L, curve_pairs, eps=1e-6):
    """Replace each slave boundary curve's 1D mesh with the exact translate of
    its master curve's mesh, reusing the shared end-point (vertex) nodes.

    Prerequisite: generate(1) has meshed all curves independently. For each
    (master, slave) pair from _periodic_curve_pairs we create new interior
    nodes on the slave as translated copies of the master interior nodes, map
    the master end nodes onto the existing slave end nodes (cube-corner /
    arc-junction point nodes, matched by translated coordinate), clear the
    slave's old line mesh and re-add the master's line elements. Afterwards
    every opposite boundary curve is an exact translate, so generate(2) meshes
    each face against a periodic boundary WITHOUT any Gmsh setPeriodic call —
    sidestepping Gmsh's internal 1D periodic node-copy, which intermittently
    fails ("Can't find periodic counterpart of mesh edge ... on curve ...") and
    was the dominant cause of exhausted mesh retries (worse with finer meshes).

    Ordering matters: a cube-edge curve lies on two max faces and is slaved,
    via the single-axis rule, to another edge that is itself a slave (a chain
    edge(x=L,y=L) -> edge(x=0,y=L) -> root edge(x=0,y=0)). We must finalise a
    master before copying any slave that reads from it, so we process pairs in
    order of increasing slave "depth" = number of max faces the slave lies on
    (1 before 2). Reading the master with getNodes then returns its already
    translated nodes, composing the chain correctly.

    Returns (n_pairs, n_missing_end_nodes); n_missing must be 0.
    """
    def key(p):
        return (round(p[0], 7), round(p[1], 7), round(p[2], 7))

    def depth(curve):
        bb = gmsh.model.getBoundingBox(1, curve)
        return sum(1 for a in range(3)
                   if abs(bb[a] - L) < eps and abs(bb[a + 3] - L) < eps)

    ordered = sorted(curve_pairs, key=lambda pr: depth(pr[1]))

    all_nt, _, _ = gmsh.model.mesh.getNodes()
    next_tag = int(max(all_nt)) + 1 if len(all_nt) else 1
    max_et = 0
    _, etags_all, _ = gmsh.model.mesh.getElements()
    for a in etags_all:
        if len(a):
            max_et = max(max_et, int(max(a)))
    next_elem = max_et + 1

    total_missing = 0
    for master, slave, trans in ordered:
        # Index slave end-point (vertex) nodes by coordinate.
        send = {}
        for d, pt in gmsh.model.getBoundary([(1, slave)], oriented=False):
            nt, co, _ = gmsh.model.mesh.getNodes(0, pt)
            for i, t in enumerate(nt):
                send[key(co[3 * i:3 * i + 3])] = int(t)

        # Master interior nodes -> new translated slave interior nodes.
        mnt, mco, _ = gmsh.model.mesh.getNodes(1, master, includeBoundary=False)
        nodemap = {}
        ntags = []
        ncoords = []
        for i, t in enumerate(mnt):
            p = mco[3 * i:3 * i + 3]
            nodemap[int(t)] = next_tag
            ntags.append(next_tag)
            ncoords.extend([p[0] + trans[0], p[1] + trans[1], p[2] + trans[2]])
            next_tag += 1

        # Master end nodes -> existing slave end nodes (translated coordinate).
        for d, pt in gmsh.model.getBoundary([(1, master)], oriented=False):
            nt, co, _ = gmsh.model.mesh.getNodes(0, pt)
            for i, t in enumerate(nt):
                if int(t) in nodemap:
                    continue
                p = co[3 * i:3 * i + 3]
                k = key([p[0] + trans[0], p[1] + trans[1], p[2] + trans[2]])
                if k in send:
                    nodemap[int(t)] = send[k]
                else:
                    total_missing += 1

        # Master line elements -> slave (same node order).
        etypes, etags, enodes = gmsh.model.mesh.getElements(1, master)
        gmsh.model.mesh.clear([(1, slave)])
        if ntags:
            # Curve nodes MUST carry a parametric coordinate (u) on the slave
            # curve, or the following generate(2) cannot recover boundary edges
            # ("No parametric coordinate on node ... classified on curve" -> the
            # 2D Delaunay spins forever). The new nodes lie exactly on the slave
            # curve (translated copies), so getParametrization gives the correct
            # u in the slave's own parametrization.
            uparams = gmsh.model.getParametrization(1, slave, ncoords)
            gmsh.model.mesh.addNodes(1, slave, ntags, ncoords, uparams)
        new_etags = []
        new_enodes = []
        for ti, et in enumerate(etypes):
            nn = gmsh.model.mesh.getElementProperties(et)[3]
            etl = []
            enl = []
            for e in range(len(etags[ti])):
                ln = [int(enodes[ti][e * nn + k]) for k in range(nn)]
                if any(x not in nodemap for x in ln):
                    continue
                mapped = [nodemap[x] for x in ln]
                etl.append(next_elem)
                next_elem += 1
                enl.extend(mapped)
            new_etags.append(etl)
            new_enodes.append(enl)
        gmsh.model.mesh.addElements(1, slave, etypes, new_etags, new_enodes)

    return len(curve_pairs), total_missing


# =====================================================================
# Wavy (mean-vertical) channels -- preserves cubic periodicity
# =====================================================================

def _add_wavy_channel(cx, cy, r, L, z_ext, amp, az_cos, az_sin, npts=25, zshift=0.0):
    """'Wavy' mean-vertical brine channel built as a boolean UNION of short
    cylinder segments following the centreline (NOT occ.addPipe -- that produces
    malformed swept solids in this OCC build, an off-centre blob even for a
    straight vertical sweep). The centreline leans along azimuth (az_cos, az_sin)
    by off(z) = amp*(1 - cos(2*pi*z/L)). This waveform is L-periodic AND has BOTH
    off(z)=0 and off'(z)=0 at z=0 and z=L, so the channel is exactly VERTICAL where
    it crosses the two z-faces: the cut there is a plain CIRCLE (not a tilted
    ellipse), identical on both faces. That matters because Gmsh's setPeriodic(1)
    aligns circle cuts fine (as for straight channels) but cannot seam-align a
    closed tilted-ellipse pair. The channel still leans up to atan(2*pi*amp/L) off
    vertical in the interior (max at z=L/4, 3L/4), bulging laterally by up to 2*amp.
    The segments span beyond both faces (z_ext) so the channel is continuous across
    z=0/z=L.

    CRITICAL for periodicity: the cylinder-segment JOINTS are sampled on an
    L-PERIODIC z-grid (spacing dz = L/M, half-shifted so no joint lands exactly on
    a z-face). Then the segment crossing z=0 and the segment crossing z=L are exact
    translates (endpoints differ by exactly L, off(z) is L-periodic), so the two
    z-face cuts have identical geometry AND identical curve decompositions -- the
    prerequisite for the manual periodic mesh copy. A uniform grid over the whole
    span does NOT guarantee this (joints fall at different wave phases near the two
    faces) and breaks CAD-level z-periodicity. `npts` sets the segment density
    (segments per period M = max(8, npts//2)). Returns the fused 3D tag(s)."""
    k = 2.0 * math.pi / L
    M = max(8, npts // 2)                     # segments per period L
    dz = L / M                                # L-periodic grid spacing
    # half-shifted grid covering [-z_ext, L + z_ext]; z=0 and z=L are never grid
    # points and the grid maps to itself under a +L shift.
    i_lo = int(math.floor((-z_ext) / dz - 0.5))
    i_hi = int(math.ceil((L + z_ext) / dz - 0.5))
    pts = []
    for i in range(i_lo, i_hi + 1):
        z = (i + 0.5) * dz + zshift
        off = amp * (1.0 - math.cos(k * z))  # vertical (circle cut) at z=0 and z=L
        pts.append((cx + off * az_cos, cy + off * az_sin, z))
    segs = []
    for j in range(len(pts) - 1):
        x0, y0, z0 = pts[j]
        x1, y1, z1 = pts[j + 1]
        segs.append((3, gmsh.model.occ.addCylinder(
            x0, y0, z0, x1 - x0, y1 - y0, z1 - z0, r)))
    if len(segs) == 1:
        return [segs[0][1]]
    fused, _ = gmsh.model.occ.fuse([segs[0]], segs[1:],
                                   removeObject=True, removeTool=True)
    return [t for (d, t) in fused if d == 3]


# =====================================================================
# Main mesh generation
# =====================================================================

def generate_periodic_mesh(sphere_array, L, L_mesh, output_dir,
                           Is_Porous='Composite', run_id='gmsh_rve',
                           mesh_order=1, optimise=True,
                           VoF_void_sphere=0.0, VoF_incl_sphere=0.0,
                           Inclusion_Type='Solid', gap_balls=None):
    """
    Generate a periodic RVE mesh using Gmsh.
    
    Parameters
    ----------
    sphere_array : np.ndarray, shape (N, 10)
        Sphere data from SpaX_Kernel. Columns:
        [cx, cy, cz, rx, ry, rz, rot1, rot2, rot3, sphericity]
    L : float
        RVE side length.
    L_mesh : float
        Target mesh element size.
    output_dir : str
        Directory for output files.
    Is_Porous : str
        'Composite' (matrix + inclusions), 'Porous' (matrix + voids),
        or 'Hybrid' (matrix + voids + inclusions).
    run_id : str
        Identifier for output filenames.
    mesh_order : int
        1 for linear tets (C3D4), 2 for quadratic (C3D10).
    optimise : bool
        Run Gmsh mesh optimisation.
    VoF_void_sphere : float
        Target void volume fraction (Hybrid mode only).
    VoF_incl_sphere : float
        Target inclusion volume fraction (Hybrid mode only).
    
    Returns
    -------
    result : dict
        mesh_file: path to .inp file
        n_nodes: total node count
        n_elements_matrix: element count in matrix
        n_elements_sphere: element count in inclusions (0 if porous)
        n_face_pairs: number of periodic face pairs
        periodicity: dict with matching stats per axis
    """
    if not HAS_GMSH:
        raise RuntimeError("Gmsh not available")
    
    gmsh.initialize()
    gmsh.clear()  # Clear any leftover state from previous failed attempts
    gmsh.option.setNumber("General.Terminal", 1)
    
    # Geometry tolerances for the Boolean fragment. Kept tight (1e-6): for valid
    # non-overlapping packings the fragment produces opposite-face vertices that
    # match within setPeriodic's point tolerance, so a larger value is not needed
    # and would risk merging genuine small features at high volume fraction.
    # (The non-conforming periodic-face problem is fixed by the mesh-size field
    # settings above, not by this tolerance.)
    gmsh.option.setNumber("Geometry.Tolerance", 1e-6)
    gmsh.option.setNumber("Geometry.ToleranceBoolean", 1e-6)
    gmsh.option.setNumber("Geometry.OCCFixDegenerated", 1)
    gmsh.option.setNumber("Geometry.OCCFixSmallEdges", 1)
    gmsh.option.setNumber("Geometry.OCCFixSmallFaces", 1)
    gmsh.option.setNumber("Geometry.OCCSewFaces", 1)
    
    # Optional 2D mesher override (SPAX_MESH_ALGO2D, e.g. 6=Frontal-Delaunay).
    # Left at the Gmsh default unless set — forcing it changed element counts
    # without reducing face slivers (those are handled by cap rejection in the
    # packer for the quadratic path).
    _algo2d = os.environ.get('SPAX_MESH_ALGO2D', '')
    if _algo2d:
        gmsh.option.setNumber("Mesh.Algorithm", int(_algo2d))

    # Global mesh size bounds will be set by the adaptive field below
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", L_mesh * 3.0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", L_mesh * 0.15)
    
    # NOTE: order-2 (C3D10) is NOT requested here. Setting Mesh.ElementOrder=2
    # before the custom periodic workflow makes generate(2)/(3) emit quadratic
    # surface/volume meshes whose mid-edge nodes the manual master->slave face
    # copy does not reproduce, leaving the boundary inconsistent and crashing
    # the 3D mesher. Instead we mesh LINEARLY through the whole periodic
    # workflow and promote to quadratic with setOrder(2) AFTER generate(3): the
    # boundary faces are planar, so the added mid-edge nodes land at edge
    # midpoints of already-periodic linear edges and are themselves exact
    # periodic translates (the coordinate-based pair matcher then pairs them).

    # ROOT-CAUSE FIX for non-conforming periodic faces.
    # By default (=1) Gmsh extends/interpolates element sizes inward from each
    # surface's boundary curves. It does this INDEPENDENTLY on a min face and
    # its max face, so corresponding curves get different 1D node counts (the
    # "curve X has N vertices, whereas correspondant has M" error) and the
    # periodic mesh copy that guarantees conforming external faces fails.
    # Forcing the element size to come purely from the background field — which
    # is symmetric and therefore identical on opposite faces — makes the master
    # and slave curves mesh with identical node counts, so setPeriodic can copy
    # the master face mesh onto the slave exactly. This is the key to periodicity.
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    
    print("\n" + "=" * 70)
    print("GMSH PERIODIC MESH GENERATION")
    print("=" * 70)
    print("  L={}, L_mesh={}, order={}, mode={}".format(L, L_mesh, mesh_order, Is_Porous))
    print("  Spheres: {}".format(len(sphere_array)))
    
    # ---- Build geometry ----
    box = gmsh.model.occ.addBox(0, 0, 0, L, L, L)
    
    # Identify unique parent spheres (skip periodic copies already in array)
    # A parent is the entry with centre in [0,L]^3
    parents = []
    tol_inside = 1e-10
    seen = set()
    
    for i in range(len(sphere_array)):
        cx, cy, cz = sphere_array[i, 0:3]
        # Channel sentinel: a vertical (Z) cylinder is stored in the sphere
        # array as a tall entry (rz >= ~L/2) tagged with sphericity ~0.01. Its
        # cross-section radius is the XY semi-axis, NOT max(rx,ry,rz) (which
        # would be the full half-height) — using the bounding sphere here built
        # a box-filling sphere and crashed the Boolean kernel.
        is_channel = (sphere_array[i, 9] <= 0.02 and
                      sphere_array[i, 5] >= 0.4 * L)
        if is_channel:
            r = sphere_array[i, 3]
        else:
            r = max(sphere_array[i, 3], sphere_array[i, 4], sphere_array[i, 5])

        # Check if this is a parent (centre inside or on boundary of RVE)
        if (-tol_inside <= cx <= L + tol_inside and
            -tol_inside <= cy <= L + tol_inside and
            -tol_inside <= cz <= L + tol_inside):

            # Avoid duplicates (periodic copies already in the array)
            key = (round(cx, 8), round(cy, 8), round(cz, 8), round(r, 8))
            if key not in seen:
                seen.add(key)
                parents.append((cx, cy, cz, r, i, is_channel))

    print("  Parent spheres: {}".format(len(parents)))

    # Create spheres and their periodic copies
    sphere_tags = []  # (3, tag) for all sphere volumes
    sphere_to_parent = {}  # tag -> parent_index
    # Channels extend past the z-faces so the Boolean fragment cuts them
    # cleanly (a cylinder ending exactly on a face would graze it and confuse
    # OCC's orientation); the box clips the overhang back to [0,L].
    z_ext = 0.1 * L

    # Optional channel tilt: a mean-vertical 'wavy' channel that leans by up to
    # SPAX_CHANNEL_TILT_DEG off Z (azimuth SPAX_CHANNEL_TILT_AZIMUTH_DEG) and
    # returns to the same XY on both z-faces (position+tangent) -> periodicity
    # preserved. amp=0 keeps the original straight cylinder (zero regression).
    _tilt_deg = float(os.environ.get('SPAX_CHANNEL_TILT_DEG', '0') or 0)
    _chan_amp = (L * math.tan(math.radians(_tilt_deg)) / (2.0 * math.pi)
                 if _tilt_deg > 1e-6 else 0.0)
    # Azimuth of the lean. If SPAX_CHANNEL_TILT_AZIMUTH_DEG is set, every channel
    # leans that way (in-plane biased); otherwise each channel gets its OWN azimuth
    # (deterministic hash of its centre) so the inclination has NO net in-plane
    # direction -- the clean way to test whether tilt DILUTES the vertical
    # anisotropy without introducing an artificial x/y anisotropy. (Wavy channels
    # are kept off the x/y faces, so per-channel azimuth never breaks periodicity.)
    _az_env = os.environ.get('SPAX_CHANNEL_TILT_AZIMUTH_DEG')
    _az_fixed = math.radians(float(_az_env)) if _az_env not in (None, '') else None

    def _chan_azimuth(cx, cy):
        if _az_fixed is not None:
            return math.cos(_az_fixed), math.sin(_az_fixed)
        h = math.sin(cx * 127.1 + cy * 311.7) * 43758.5453
        ang = 2.0 * math.pi * (h - math.floor(h))
        return math.cos(ang), math.sin(ang)

    if _chan_amp > 0.0:
        print("  Wavy channels: max tilt {:.1f} deg, amp {:.4f}, azimuth {}"
              .format(_tilt_deg, _chan_amp,
                      "{:.0f} deg".format(math.degrees(_az_fixed))
                      if _az_fixed is not None else "per-channel random"))

    for cx, cy, cz, r, parent_idx, is_channel in parents:
        if is_channel:
            if _chan_amp > 0.0:
                # Wavy channels stay clear of the x/y faces (packing margin += amp)
                # but they DO cross the z-faces. Build a union of cylinder segments
                # spanning a full period BEYOND each z-face (big_ext) so the channel
                # is continuous across z=0 and z=L -- the box clips it and, because
                # the waveform is L-periodic, the two faces get IDENTICAL cut-curve
                # decompositions that pair for the periodic mesh copy.
                big_ext = L + z_ext
                _acos, _asin = _chan_azimuth(cx, cy)
                for t in _add_wavy_channel(cx, cy, r, L, big_ext, _chan_amp,
                                           _acos, _asin, npts=25):
                    sphere_tags.append((3, t)); sphere_to_parent[t] = parent_idx
                continue
            # Vertical cylinder spanning the full height (Z), XY-periodic only.
            cyl = gmsh.model.occ.addCylinder(
                cx, cy, -z_ext, 0.0, 0.0, L + 2.0 * z_ext, r)
            sphere_tags.append((3, cyl))
            sphere_to_parent[cyl] = parent_idx
            for ccx, ccy in _add_periodic_cylinder_copies(cx, cy, r, L):
                cc = gmsh.model.occ.addCylinder(
                    ccx, ccy, -z_ext, 0.0, 0.0, L + 2.0 * z_ext, r)
                sphere_tags.append((3, cc))
                sphere_to_parent[cc] = parent_idx
            continue
        # Main sphere (clipped to RVE later via fragment)
        s = gmsh.model.occ.addSphere(cx, cy, cz, r)
        sphere_tags.append((3, s))
        sphere_to_parent[s] = parent_idx

        # Periodic copies
        copies = _add_periodic_sphere_copies(cx, cy, cz, r, L)
        for ccx, ccy, ccz in copies:
            sc = gmsh.model.occ.addSphere(ccx, ccy, ccz, r)
            sphere_tags.append((3, sc))
            sphere_to_parent[sc] = parent_idx
    
    print("  Total sphere volumes (incl copies): {}".format(len(sphere_tags)))
    
    # ---- Boolean fragment ----
    print("  Running Boolean fragment...")
    if sphere_tags:
        out, outmap = gmsh.model.occ.fragment([(3, box)], sphere_tags)
    else:
        out = [(3, box)]
        outmap = [[(3, box)]]
    gmsh.model.occ.synchronize()
    
    # Post-fragment healing: remove small edges and merge near-coincident
    # vertices that cause asymmetric curve counts on periodic faces.
    try:
        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()
    except:
        pass  # removeAllDuplicates not available in older Gmsh versions
    
    print("  Fragment produced {} volumes".format(len(out)))
    
    # ---- Identify matrix vs sphere volumes using outmap ----
    # outmap[0] = volumes that came from the box (ALL inside [0,L]^3)
    # outmap[1..N] = volumes that came from each sphere
    # Matrix = box_children NOT in any sphere_children
    # Sphere_inside = box_children AND sphere_children
    # Outside = sphere_children NOT in box_children
    
    box_children = set()
    for dim, tag in outmap[0]:
        if dim == 3:
            box_children.add(tag)
    
    sphere_children = set()
    for i in range(len(sphere_tags)):
        for dim, tag in outmap[i + 1]:
            if dim == 3:
                sphere_children.add(tag)
    
    matrix_vols = sorted(list(box_children - sphere_children))
    sphere_vols = sorted(list(box_children & sphere_children))
    outside_vols = sorted(list(sphere_children - box_children))
    
    print("  Matrix volumes: {} tags={}".format(len(matrix_vols), matrix_vols))
    
    if len(matrix_vols) == 0:
        gmsh.finalize()
        raise RuntimeError("Boolean fragment produced 0 matrix volumes. "
                           "Check sphere positions (overlapping or outside RVE).")
    print("  Sphere volumes: {} tags={}".format(len(sphere_vols), sphere_vols))
    print("  Outside volumes: {} (excluded)".format(len(outside_vols)))

    # ---- Remove OUTSIDE protrusions from the geometry ----
    # The part of an inclusion poking past a box face is redundant: the periodic
    # copy already supplies that material on the opposite face. Leaving it in the
    # model means the 2D mesher meshes its outer surface, which sits right
    # against the box cut-disc of the SAME inclusion -- and the 3D mesher then
    # reports "overlapping facets" on those two coincident surfaces and crashes
    # (seed-dependent: it bites whenever an inclusion crosses a face deeply
    # enough to leave a sizeable protrusion). recursive=True drops the
    # protrusion's now-unused surfaces/curves but keeps the cut-disc, which is
    # shared with the retained inside body. This is a geometry fix, independent
    # of mesh size, so it cannot be worked around by refinement.
    if outside_vols:
        try:
            gmsh.model.occ.remove([(3, t) for t in outside_vols], recursive=True)
            gmsh.model.occ.synchronize()
            print("  Removed {} outside protrusion volume(s) from geometry".format(
                len(outside_vols)))
        except Exception as _e:
            print("  WARN: could not remove outside volumes ({}); "
                  "meshing may hit overlapping facets".format(_e))

    # ---- Hybrid mode: classify sphere volumes as void or inclusion ----
    void_vols = []
    incl_vols = []
    
    if Is_Porous == 'Hybrid' and sphere_vols:
        # Classify by random assignment (matching kernel behaviour)
        # Group sphere volumes by parent sphere using centroid proximity
        import random
        
        V_RVE = L ** 3
        vol_data = []  # (vol_tag, parent_sphere_idx, volume, centroid)
        
        for vtag in sphere_vols:
            try:
                com = gmsh.model.occ.getCenterOfMass(3, vtag)
                vol = gmsh.model.occ.getMass(3, vtag)
            except:
                com = (0, 0, 0)
                vol = 0
            
            # Find closest parent sphere
            best_parent = 0
            best_dist = 1e30
            for pi, (scx, scy, scz, sr, _pidx, p_is_channel) in enumerate(parents):
                if p_is_channel:
                    # Channel match on XY only (cylinder spans full Z); compare
                    # against the nearest XY periodic image.
                    dmin = 1e30
                    for ccx, ccy in ([(scx, scy)]
                                     + _add_periodic_cylinder_copies(scx, scy, sr, L)):
                        d = math.sqrt((com[0]-ccx)**2 + (com[1]-ccy)**2)
                        if d < dmin:
                            dmin = d
                    if dmin < best_dist:
                        best_dist = dmin
                        best_parent = pi
                    continue
                d = math.sqrt((com[0]-scx)**2 + (com[1]-scy)**2 + (com[2]-scz)**2)
                if d < best_dist:
                    best_dist = d
                    best_parent = pi
                # Also check periodic copies
                for ccx, ccy, ccz in _add_periodic_sphere_copies(scx, scy, scz, sr, L):
                    d2 = math.sqrt((com[0]-ccx)**2 + (com[1]-ccy)**2 + (com[2]-ccz)**2)
                    if d2 < best_dist:
                        best_dist = d2
                        best_parent = pi
            
            vol_data.append((vtag, best_parent, vol, com))
        
        # Group by parent
        parent_groups = {}
        for vtag, pidx, vol, com in vol_data:
            if pidx not in parent_groups:
                parent_groups[pidx] = []
            parent_groups[pidx].append((vtag, vol, com))
        
        # Random assignment: shuffle parents, fill void quota first
        parent_list = list(parent_groups.keys())
        random.shuffle(parent_list)
        
        void_vof_accum = 0.0
        for pidx in parent_list:
            parent_vol = sum(v for _, v, _ in parent_groups[pidx])
            parent_vof = parent_vol / V_RVE
            
            if void_vof_accum < VoF_void_sphere:
                # Assign to void (air — unmeshed)
                for vtag, _, _ in parent_groups[pidx]:
                    void_vols.append(vtag)
                void_vof_accum += parent_vof
            else:
                # Assign to inclusion (MESHED — both Solid and Liquid modes)
                # Liquid mode only changes element type and material, not mesh topology
                for vtag, _, _ in parent_groups[pidx]:
                    incl_vols.append(vtag)
        
        print("  Hybrid classification:")
        print("    Void volumes (air): {} (VoF ~ {:.2%})".format(len(void_vols), void_vof_accum))
        print("    Inclusion volumes: {} ({})".format(
            len(incl_vols), 'liquid (K/G)' if Inclusion_Type == 'Liquid' else 'solid'))
    
    elif Is_Porous == 'Porous':
        void_vols = list(sphere_vols)
        incl_vols = []
    else:
        # Composite (both Solid and Liquid): all spheres are meshed inclusions
        void_vols = []
        incl_vols = list(sphere_vols)
    
    # ---- Physical groups for material assignment ----
    if matrix_vols:
        gmsh.model.addPhysicalGroup(3, matrix_vols, tag=1)
        gmsh.model.setPhysicalName(3, 1, "Matrix_Only")
    
    if incl_vols:
        gmsh.model.addPhysicalGroup(3, incl_vols, tag=2)
        gmsh.model.setPhysicalName(3, 2, "Sphere_Only")
    
    # Void volumes get NO physical group -> not meshed (empty cavities)
    if void_vols:
        print("  Void volumes (not meshed): {}".format(void_vols))
    
    # ---- Set mesh size (ADAPTIVE) ----
    # Fine mesh near inclusion surfaces, coarse in matrix interior.
    # Uses Gmsh Distance + Threshold fields for smooth grading.
    
    # Mesh size parameters
    # lc_fine sets the surface element size at inclusion boundaries. Finer
    # lc_fine resolves smaller spheres, which lets the packer use a lower
    # min-radius floor and reach higher VoF (env-tunable for sweeps).
    lc_fine = L_mesh * float(os.environ.get('SPAX_LC_FINE_MULT', '0.4'))
    lc_coarse = L_mesh * 2.5    # coarse size far from inclusions
    dist_min = 0.0              # start fine at surface
    dist_max = L * 0.15         # transition to coarse over 15% of L
    
    # Ensure global bounds
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc_coarse)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc_fine * 0.5)
    
    # Collect all inclusion surface tags (spheres + voids)
    incl_surface_tags = []
    all_incl_vols = list(sphere_vols)
    if void_vols:
        all_incl_vols.extend(void_vols)
    
    for dim, tag in out:
        if dim == 3 and tag in all_incl_vols:
            boundary = gmsh.model.getBoundary([(3, tag)], oriented=False)
            for bdim, btag in boundary:
                if bdim == 2:
                    incl_surface_tags.append(btag)
    
    if incl_surface_tags:
        # Distance field: distance from inclusion surfaces
        f_dist = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(f_dist, "SurfacesList", incl_surface_tags)
        gmsh.model.mesh.field.setNumber(f_dist, "Sampling", 100)
        
        # Threshold field: maps distance to mesh size
        # lc_fine at distance <= dist_min, lc_coarse at distance >= dist_max
        f_thresh = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(f_thresh, "InField", f_dist)
        gmsh.model.mesh.field.setNumber(f_thresh, "SizeMin", lc_fine)
        gmsh.model.mesh.field.setNumber(f_thresh, "SizeMax", lc_coarse)
        gmsh.model.mesh.field.setNumber(f_thresh, "DistMin", dist_min)
        gmsh.model.mesh.field.setNumber(f_thresh, "DistMax", dist_max)
        
        # Mesh-in-gap refinement: a local size ball at each narrow sphere-sphere
        # sliver the packer flagged. Forcing a small element size inside the gap
        # lets the mesher resolve a sub-element matrix sliver (which would
        # otherwise produce overlapping facets) WITHOUT the packer having to
        # widen the gap, so the pack keeps its volume fraction. Each ball:
        # SizeIn at the sliver centre, grading to lc_coarse over Thickness.
        gap_fields = []
        if gap_balls:
            for (gx, gy, gz, gsize) in gap_balls:
                f_ball = gmsh.model.mesh.field.add("Ball")
                gmsh.model.mesh.field.setNumber(f_ball, "XCenter", gx)
                gmsh.model.mesh.field.setNumber(f_ball, "YCenter", gy)
                gmsh.model.mesh.field.setNumber(f_ball, "ZCenter", gz)
                # Radius of the fully-refined core: the sliver's lateral reach
                # (a few element sizes) so the thin matrix region is covered.
                r_in = max(4.0 * gsize, 0.5 * lc_fine)
                gmsh.model.mesh.field.setNumber(f_ball, "Radius", r_in)
                gmsh.model.mesh.field.setNumber(f_ball, "Thickness", 3.0 * r_in)
                gmsh.model.mesh.field.setNumber(f_ball, "VIn", gsize)
                gmsh.model.mesh.field.setNumber(f_ball, "VOut", lc_coarse)
                gap_fields.append(f_ball)
            # allow the global floor to reach the smallest requested ball size
            min_ball = min(b[3] for b in gap_balls)
            gmsh.option.setNumber("Mesh.CharacteristicLengthMin",
                                  min(lc_fine * 0.5, min_ball * 0.9))
            print("  Mesh-in-gap: {} refinement ball(s), smallest size "
                  "{:.5f}".format(len(gap_balls), min_ball))

        # Also need fine mesh on RVE boundary faces for periodic matching
        boundary_tags = []
        cube_boundary = gmsh.model.getBoundary([(3, matrix_vols[0])], oriented=False)
        for bdim, btag in cube_boundary:
            if bdim == 2:
                boundary_tags.append(btag)

        if boundary_tags:
            # Distance field from RVE boundaries
            f_bdist = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(f_bdist, "SurfacesList", boundary_tags)
            gmsh.model.mesh.field.setNumber(f_bdist, "Sampling", 100)
            
            # Moderate size on boundaries for good periodic matching
            lc_boundary = L_mesh * 1.0
            f_bthresh = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(f_bthresh, "InField", f_bdist)
            gmsh.model.mesh.field.setNumber(f_bthresh, "SizeMin", lc_boundary)
            gmsh.model.mesh.field.setNumber(f_bthresh, "SizeMax", lc_coarse)
            gmsh.model.mesh.field.setNumber(f_bthresh, "DistMin", 0.0)
            gmsh.model.mesh.field.setNumber(f_bthresh, "DistMax", L * 0.05)
            
            # Minimum of all fields (surface grading + boundary + gap balls)
            f_min = gmsh.model.mesh.field.add("Min")
            gmsh.model.mesh.field.setNumbers(f_min, "FieldsList",
                                             [f_thresh, f_bthresh] + gap_fields)
            gmsh.model.mesh.field.setAsBackgroundMesh(f_min)
        elif gap_fields:
            f_min = gmsh.model.mesh.field.add("Min")
            gmsh.model.mesh.field.setNumbers(f_min, "FieldsList",
                                             [f_thresh] + gap_fields)
            gmsh.model.mesh.field.setAsBackgroundMesh(f_min)
        else:
            gmsh.model.mesh.field.setAsBackgroundMesh(f_thresh)
        
        print("  Adaptive mesh: lc_fine={:.5f}, lc_coarse={:.5f}, grade_dist={:.4f}".format(
            lc_fine, lc_coarse, dist_max))
    else:
        # No inclusions — uniform mesh (pure matrix)
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), L_mesh)
        print("  Uniform mesh: lc={:.5f}".format(L_mesh))
    
    # ---- DEBUG: map specific surface tags to their nearest inclusion ----
    _dbg = os.environ.get('SPAX_DEBUG_SURF', '')
    if _dbg:
        import math as _m
        targets = [int(t) for t in _dbg.split(',') if t.strip()]
        for st in targets:
            try:
                bb = gmsh.model.getBoundingBox(2, st)
            except Exception as _e:
                print("  [dbg] surface {} not found ({})".format(st, _e))
                continue
            sx = 0.5*(bb[0]+bb[3]); sy = 0.5*(bb[1]+bb[4]); sz = 0.5*(bb[2]+bb[5])
            ext = (bb[3]-bb[0], bb[4]-bb[1], bb[5]-bb[2])
            best = None
            for (px, py, pz, pr, pidx, pch) in parents:
                dx = sx-px; dx -= L*round(dx/L)
                dy = sy-py; dy -= L*round(dy/L)
                dz = sz-pz; dz -= L*round(dz/L)
                d = _m.sqrt(dx*dx+dy*dy+dz*dz)
                if best is None or d < best[0]:
                    best = (d, px, py, pz, pr, pidx, pch)
            print("  [dbg] surf {}: ctr=({:.3f},{:.3f},{:.3f}) ext=({:.3f},{:.3f},{:.3f}) "
                  "-> parent idx={} ch={} ctr=({:.3f},{:.3f},{:.3f}) r={:.3f} dist={:.4f}".format(
                      st, sx, sy, sz, ext[0], ext[1], ext[2],
                      best[5], best[6], best[1], best[2], best[3], best[4], best[0]))

    # ---- Apply periodic constraints (strict 1D + 2D, then 3D) ----
    # Conforming external faces are MANDATORY for a valid periodic RVE: the mesh
    # on each min face must be an exact translate of its max face. Gmsh's own
    # setPeriodic(2) is unreliable here — once inclusion cuts split a face into
    # many sub-surfaces it mis-matches cube-edge curves (mapping an edge to the
    # wrong parallel edge) and the 2D copy fails with "curve X has N vertices,
    # whereas correspondant has M". That internal correspondence drives the copy
    # and cannot be overridden through the API. So we enforce periodicity
    # ourselves in three stages:
    #   1. setPeriodic(1) on every boundary curve, using a single-axis master
    #      tree (_periodic_curve_pairs) so 1D meshes of opposite curves are
    #      exact translates and cube edges are never double-assigned. Gmsh owns
    #      the 1D periodic node-copy (handling parametric coords and loop
    #      closure correctly); when opposite-face curve decompositions differ it
    #      raises a *catchable* "Can't find periodic counterpart" — the caller's
    #      retry loop re-packs, and Fix-1 subprocess isolation makes that safe.
    #   2. generate(2): mesh all surfaces. Boundary curves are now periodic;
    #      face interiors are still independent.
    #   3. _manual_periodic_surface_copy: stamp each master face's 2D mesh onto
    #      its slave as an exact translate, reusing the periodic boundary nodes.
    # The result is conforming node-for-node on every opposite-face pair.
    #
    # (A fully manual 1D copy was tried to remove the last setPeriodic
    # dependency, but per-curve stamping can leave a slave face boundary that
    # is not a closed loop when OCC splits opposite faces differently — it then
    # fails on every retry instead of intermittently. setPeriodic(1) + retries
    # is strictly more robust. _manual_periodic_curve_copy is kept unused for
    # reference.)
    print("  Matching periodic face pairs...")
    face_pairs = _match_periodic_surfaces(L)
    print("  Matched {} surface pairs".format(len(face_pairs)))

    # Stage 1: periodic 1D curves. Works for straight AND wavy channels: the wavy
    # channel's 1-cos waveform is vertical where it crosses the z-faces, so its cut
    # there is a plain CIRCLE (like a straight channel), which setPeriodic(1) aligns
    # fine -- no closed-tilted-ellipse seam problem.
    affine_id = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    curve_pairs = _periodic_curve_pairs(L)
    for master_c, slave_c, trans in curve_pairs:
        aff = list(affine_id)
        aff[3], aff[7], aff[11] = trans[0], trans[1], trans[2]
        gmsh.model.mesh.setPeriodic(1, [slave_c], [master_c], aff)
    print("  Periodic curve pairs set: {}".format(len(curve_pairs)))

    # Stage 2: mesh 1D + 2D (master faces meshed with the size field; slave
    # boundary curves are periodic copies of their masters).
    gmsh.option.setNumber("Mesh.ToleranceInitialDelaunay", 1e-8)
    print("  Generating surface mesh (2D)...")
    gmsh.model.mesh.generate(2)

    # Stage 3: manual periodic copy of every slave face mesh.
    print("  Copying master face meshes onto slaves...")
    n_copied, n_missing = _manual_periodic_surface_copy(L, face_pairs)
    if n_missing:
        gmsh.finalize()
        raise RuntimeError(
            "Manual periodic copy left {} slave boundary node(s) unmatched — "
            "the 1D periodic mesh is inconsistent (retry packing).".format(
                n_missing))
    print("  Copied {} face pairs (boundary nodes all matched)".format(n_copied))

    # ---- Generate 3D mesh on the now-periodic boundary ----
    print("  Generating 3D mesh...")
    gmsh.model.mesh.generate(3)

    if optimise:
        print("  Optimising mesh...")
        gmsh.model.mesh.optimize("", force=True)

    # Promote linear -> quadratic AFTER the periodic linear mesh is complete.
    # Straight-sided (mid-nodes at edge midpoints, the default) so the boundary
    # mid-edge nodes stay exact periodic translates. Interior inclusion surfaces
    # are not on the periodic boundary, so leaving them straight-sided does not
    # affect periodicity (it is a minor geometric-accuracy trade we accept for
    # strict periodicity). Done before all node harvesting below so the new
    # mid-edge nodes are exported and PBC-paired.
    if mesh_order == 2:
        # A few skewed sliver tets near curved inclusion surfaces are
        # positive-VOLUME as linear C3D4 but go negative-JACOBIAN at the
        # quadratic Gauss points once midpoint nodes are added (Abaqus then
        # aborts: "volume zero/small/negative"). Fix the QUALITY of the LINEAR
        # mesh first (Netgen de-slivers by moving interior nodes; like the
        # default optimiser it respects the fixed periodic boundary, so the
        # 0-skipped pairing is preserved), THEN promote straight-sided so the
        # planar-boundary mid-nodes stay exact periodic translates. We do NOT
        # use optimize("HighOrder"): it slides boundary mid-nodes in-plane and
        # breaks periodicity.
        try:
            # Aggressive: optimise every element below 0.5 quality, repeat so
            # edge-swaps/relocations cascade through clustered slivers.
            gmsh.option.setNumber("Mesh.OptimizeThreshold", 0.5)
            for _ in range(int(os.environ.get('SPAX_OPT_PASSES', '3'))):
                gmsh.model.mesh.optimize("Netgen")
        except Exception as _e:
            print("    [Netgen optimise] skipped: {}".format(_e))
        print("  Promoting to quadratic (C3D10) via setOrder(2)...")
        # STRAIGHT-SIDED second order: place every mid-edge node at the edge
        # midpoint instead of projecting it onto the CAD surface. Two reasons:
        # (1) projecting inclusion-surface mid-nodes onto the spheres curves the
        #     elements and INVERTS the coarse ones on high curvature (minSICN<0);
        # (2) straight outer-boundary mid-nodes stay exact periodic translates.
        # The geometric cost (faceted inclusion surfaces) is small at this element
        # size and is the same approximation the linear mesh already makes.
        gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)
        gmsh.model.mesh.setOrder(2)

        # Validity GATE: any straight-sided C3D10 built on a too-distorted linear
        # tet has a negative Jacobian at its Gauss points and Abaqus aborts on
        # input ("volume zero/small/negative"). These slivers are STOCHASTIC
        # (a function of the random packing), so rather than chase them with
        # global refinement we reject the whole mesh and let the parent retry
        # loop RE-PACK (with escalated cap/sliver rejection) until a packing
        # meshes cleanly. minSICN < 0 == inverted; require a small positive
        # margin so Abaqus' stricter check also passes. Tunable via SPAX_MIN_SICN.
        try:
            etags3 = gmsh.model.mesh.getElements(3)[1]
            etags3 = np.concatenate(etags3) if len(etags3) else np.array([])
            q = gmsh.model.mesh.getElementQualities(etags3.astype('uint64'), "minSICN")
            min_sicn = float(np.min(q)) if len(q) else 1.0
            n_bad = int(np.sum(np.asarray(q) < float(os.environ.get('SPAX_MIN_SICN', '0.01'))))
            print("  Quadratic quality: min(minSICN)={:.4f}, {} element(s) below tol".format(
                min_sicn, n_bad))
            if n_bad > 0:
                raise RuntimeError(
                    "quadratic mesh has {} inverted/degenerate C3D10 element(s) "
                    "(min minSICN={:.4f}) — re-pack".format(n_bad, min_sicn))
        except RuntimeError:
            raise
        except Exception as _e:
            print("    [Jacobian gate] quality query unavailable: {}".format(_e))

    # ---- Get mesh statistics ----
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    n_nodes = len(node_tags)
    coords = node_coords.reshape(-1, 3)
    
    # Count elements per physical group (skip removed outside volumes)
    n_elem_matrix = 0
    n_elem_sphere = 0
    outside_set = set(outside_vols + void_vols)
    
    for dim, tag in out:
        if dim != 3:
            continue
        if tag in outside_set:
            continue
        try:
            etypes, etags, _ = gmsh.model.mesh.getElements(dim, tag)
        except:
            continue
        n_elem = sum(len(t) for t in etags)
        if tag in matrix_vols:
            n_elem_matrix += n_elem
        elif tag in sphere_vols:
            n_elem_sphere += n_elem
    
    print("\n  Mesh statistics:")
    print("    Nodes:            {}".format(n_nodes))
    print("    Matrix elements:  {}".format(n_elem_matrix))
    print("    Sphere elements:  {}".format(n_elem_sphere))
    print("    Total elements:   {}".format(n_elem_matrix + n_elem_sphere))
    
    # ---- Verify periodicity ----
    tol_match = 1e-10
    periodicity = {}
    
    for axis_name, axis_idx in [('X', 0), ('Y', 1), ('Z', 2)]:
        neg_nodes = []
        pos_nodes = []
        for i in range(len(coords)):
            if abs(coords[i][axis_idx]) < tol_match:
                neg_nodes.append(coords[i])
            if abs(coords[i][axis_idx] - L) < tol_match:
                pos_nodes.append(coords[i])
        
        ip = [j for j in range(3) if j != axis_idx]
        matched = 0
        for cn in neg_nodes:
            for cp in pos_nodes:
                if (abs(cn[ip[0]] - cp[ip[0]]) < tol_match and
                    abs(cn[ip[1]] - cp[ip[1]]) < tol_match):
                    matched += 1
                    break
        
        status = 'PERFECT' if matched == len(neg_nodes) else 'INCOMPLETE ({}/{})'.format(
            matched, len(neg_nodes))
        periodicity[axis_name] = {
            'neg': len(neg_nodes),
            'pos': len(pos_nodes),
            'matched': matched,
            'status': status
        }
        print("    {}: {} neg, {} pos, {} matched — {}".format(
            axis_name, len(neg_nodes), len(pos_nodes), matched, status))
    
    # ---- Export as custom .inp (only inside nodes and elements) ----
    # Gmsh's default writer includes ALL nodes even with SaveAll=0.
    # We write our own .inp with only nodes belonging to inside elements.
    inp_path = os.path.abspath(os.path.join(output_dir, '{}_gmsh.inp'.format(run_id)))
    
    # Collect nodes and elements from inside volumes only
    inside_node_tags = set()
    inside_elements = []  # (elem_tag, elem_type_name, node_tags)
    
    matrix_elem_tags = set()
    sphere_elem_tags = set()
    
    for vol_tag in matrix_vols:
        try:
            etypes, etags, enodes = gmsh.model.mesh.getElements(3, vol_tag)
            for ti, et in enumerate(etypes):
                name, _, _, nnodes, _, _ = gmsh.model.mesh.getElementProperties(et)
                for ei in range(len(etags[ti])):
                    elem_tag = int(etags[ti][ei])
                    en = [int(enodes[ti][ei * nnodes + k]) for k in range(nnodes)]
                    inside_elements.append((elem_tag, name, en))
                    inside_node_tags.update(en)
                    matrix_elem_tags.add(elem_tag)
        except:
            pass
    
    for vol_tag in incl_vols:
        try:
            etypes, etags, enodes = gmsh.model.mesh.getElements(3, vol_tag)
            for ti, et in enumerate(etypes):
                name, _, _, nnodes, _, _ = gmsh.model.mesh.getElementProperties(et)
                for ei in range(len(etags[ti])):
                    elem_tag = int(etags[ti][ei])
                    en = [int(enodes[ti][ei * nnodes + k]) for k in range(nnodes)]
                    inside_elements.append((elem_tag, name, en))
                    inside_node_tags.update(en)
                    sphere_elem_tags.add(elem_tag)
        except:
            pass
    
    # Get coordinates for inside nodes only
    all_node_tags, all_node_coords, _ = gmsh.model.mesh.getNodes()
    all_coords = all_node_coords.reshape(-1, 3)
    node_coord_map = {}
    for i, ntag in enumerate(all_node_tags):
        ntag_int = int(ntag)
        if ntag_int in inside_node_tags:
            node_coord_map[ntag_int] = all_coords[i]
    
    # Remap node labels to sequential 1..N
    old_to_new = {}
    new_label = 1
    for old_tag in sorted(inside_node_tags):
        old_to_new[old_tag] = new_label
        new_label += 1
    
    # Determine element type string for Abaqus
    # Gmsh "Tetrahedron 4" -> C3D4, "Tetrahedron 10" -> C3D10
    # Liquid mode: use C3D4H (hybrid linear tet) to handle near-incompressibility.
    # Same node count as C3D4, just adds internal pressure DOF.
    if Inclusion_Type == 'Liquid':
        elem_type_map = {
            'Tetrahedron 4': 'C3D4H',
            'Tetrahedron 10': 'C3D10H',
        }
        print("  Element type: C3D4H (hybrid, for near-incompressible inclusions)")
    else:
        elem_type_map = {
            'Tetrahedron 4': 'C3D4',
            'Tetrahedron 10': 'C3D10',
        }
    
    # Write .inp — flat format (no *PART wrapper) for better Abaqus compatibility
    with open(inp_path, 'w') as f:
        f.write('*HEADING\n')
        f.write('Gmsh periodic RVE mesh\n')
        f.write('*NODE, NSET=AllNodes\n')
        for old_tag in sorted(inside_node_tags):
            nl = old_to_new[old_tag]
            c = node_coord_map[old_tag]
            f.write('{}, {}, {}, {}\n'.format(nl, c[0], c[1], c[2]))
        
        # Write matrix elements
        abq_type = None
        for _, gtype, _ in inside_elements:
            if gtype in elem_type_map:
                abq_type = elem_type_map[gtype]
                break
        if abq_type is None:
            abq_type = 'C3D4'
        
        f.write('*ELEMENT, TYPE={}, ELSET=Matrix_Only\n'.format(abq_type))
        el = 1
        matrix_start = el
        for etag, gtype, enodes in inside_elements:
            if etag in matrix_elem_tags:
                conn = ', '.join([str(old_to_new[n]) for n in enodes])
                f.write('{}, {}\n'.format(el, conn))
                el += 1
        matrix_end = el - 1
        
        sphere_start = 0
        sphere_end = 0
        if sphere_elem_tags:
            f.write('*ELEMENT, TYPE={}, ELSET=Sphere_Only\n'.format(abq_type))
            sphere_start = el
            for etag, gtype, enodes in inside_elements:
                if etag in sphere_elem_tags:
                    conn = ', '.join([str(old_to_new[n]) for n in enodes])
                    f.write('{}, {}\n'.format(el, conn))
                    el += 1
            sphere_end = el - 1
        
        # Write explicit ELSET definitions using GENERATE for reliability
        f.write('*ELSET, ELSET=Matrix_Only, GENERATE\n')
        f.write('{}, {}, 1\n'.format(matrix_start, matrix_end))
        
        if sphere_start > 0:
            f.write('*ELSET, ELSET=Sphere_Only, GENERATE\n')
            f.write('{}, {}, 1\n'.format(sphere_start, sphere_end))
    
    print("\n  Written custom .inp: {}".format(inp_path))
    print("    Inside nodes: {}, elements: {} (matrix: {}-{}, sphere: {}-{})".format(
        len(inside_node_tags), el - 1,
        matrix_start, matrix_end, sphere_start, sphere_end))
    
    # ---- Export periodic node pairs (using remapped labels) ----
    match_path = os.path.abspath(os.path.join(output_dir, '{}_periodic_pairs.csv'.format(run_id)))
    
    print("  Nodes in inside volumes: {}".format(len(inside_node_tags)))
    
    # Build coordinate lookup with NEW labels
    new_node_coords = {}
    for old_tag in inside_node_tags:
        new_node_coords[old_to_new[old_tag]] = node_coord_map[old_tag]
    
    with open(match_path, 'w') as f:
        f.write("axis,neg_node,pos_node,y_or_ip0,z_or_ip1\n")
        
        n_pairs_written = 0
        n_pairs_skipped = 0
        
        for axis_name, axis_idx in [('X', 0), ('Y', 1), ('Z', 2)]:
            ip = [j for j in range(3) if j != axis_idx]
            
            neg_data = []
            pos_data = []
            for new_tag, c in new_node_coords.items():
                if abs(c[axis_idx]) < tol_match:
                    neg_data.append((new_tag, c[ip[0]], c[ip[1]]))
                if abs(c[axis_idx] - L) < tol_match:
                    pos_data.append((new_tag, c[ip[0]], c[ip[1]]))
            
            for ntag, ny, nz in neg_data:
                matched = False
                for ptag, py, pz in pos_data:
                    if abs(ny - py) < tol_match and abs(nz - pz) < tol_match:
                        f.write("{},{},{},{},{}\n".format(
                            axis_name, ntag, ptag, ny, nz))
                        n_pairs_written += 1
                        matched = True
                        break
                if not matched:
                    n_pairs_skipped += 1
    
    print("  Periodic pairs written: {}, skipped: {}".format(
        n_pairs_written, n_pairs_skipped))
    print("  Periodic pairs: {}".format(match_path))
    
    gmsh.finalize()
    
    result = {
        'mesh_file': inp_path,
        'match_file': match_path,
        'n_nodes': n_nodes,
        'n_elements_matrix': n_elem_matrix,
        'n_elements_sphere': n_elem_sphere,
        'n_face_pairs': len(face_pairs),
        'n_pairs_written': n_pairs_written,
        'n_pairs_skipped': n_pairs_skipped,
        'periodicity': periodicity,
        'matrix_vols': matrix_vols,
        'sphere_vols': sphere_vols,
        'matrix_label_range': (matrix_start, matrix_end),
        'sphere_label_range': (sphere_start, sphere_end),
        'Inclusion_Type': Inclusion_Type,
    }
    
    print("\n" + "=" * 70)
    print("GMSH PERIODIC MESH COMPLETE")
    print("=" * 70)
    
    return result


# =====================================================================
# Abaqus import script (runs inside Abaqus Python)
# =====================================================================

def write_abaqus_import_script(result, L, L_mesh, Disp, 
                                E_matrix, nu_matrix,
                                E_inclusion=0, nu_inclusion=0.3,
                                Is_Porous='Composite',
                                Mode='Uniaxial Tension X',
                                output_dir='.'):
    """
    Write a Python script that can be run inside Abaqus CAE to:
    1. Import the Gmsh .inp mesh
    2. Assign materials
    3. Apply direct PBC equations using the matched node pairs
    4. Set up loading and submit
    """
    inp_path = result['mesh_file']
    match_path = result['match_file']
    
    script_path = os.path.abspath(os.path.join(output_dir, 'import_gmsh_mesh.py'))
    
    with open(script_path, 'w') as f:
        # Get element label ranges
        m_start, m_end = result['matrix_label_range']
        s_start, s_end = result['sphere_label_range']
        
        f.write('# -*- coding: utf-8 -*-\n')
        f.write('# Auto-generated script to import Gmsh periodic mesh into Abaqus.\n')
        f.write('\n')
        f.write('from abaqus import mdb, session\n')
        f.write('from abaqusConstants import *\n')
        f.write('import csv\n')
        f.write('import os\n')
        f.write('\n')
        f.write('# ---- Parameters ----\n')
        f.write('L = {}\n'.format(L))
        f.write('L_mesh = {}\n'.format(L_mesh))
        f.write('Disp = {}\n'.format(Disp))
        f.write('E_matrix = {}\n'.format(E_matrix))
        f.write('nu_matrix = {}\n'.format(nu_matrix))
        f.write('E_inclusion = {}\n'.format(E_inclusion))
        f.write('nu_inclusion = {}\n'.format(nu_inclusion))
        f.write("Is_Porous = '{}'\n".format(Is_Porous))
        f.write("Mode = '{}'\n".format(Mode))
        f.write("inp_path = r'{}'\n".format(inp_path))
        f.write("match_path = r'{}'\n".format(match_path))
        f.write('matrix_label_start = {}\n'.format(m_start))
        f.write('matrix_label_end = {}\n'.format(m_end))
        f.write('sphere_label_start = {}\n'.format(s_start))
        f.write('sphere_label_end = {}\n'.format(s_end))
        f.write('''
# ---- Clean up ----
try:
    del mdb.models['Model-1']
except:
    pass
mdb.Model(name='Model-1')
m = mdb.models['Model-1']

# ---- Import mesh ----
print("Importing Gmsh mesh...")
m.PartFromInputFile(inputFileName=inp_path)

# List imported parts and their element sets
print("Imported parts:")
for pname in m.parts.keys():
    p = m.parts[pname]
    print("  {} : {} nodes, {} elements".format(pname, len(p.nodes), len(p.elements)))
    print("    Element sets: {}".format(list(p.sets.keys())))

# ---- Create assembly ----
a = m.rootAssembly
for pname in m.parts.keys():
    a.Instance(name=pname + '-1', part=m.parts[pname], dependent=ON)

# ---- Assign materials using element label ranges ----
# The .inp file has matrix elements 1..matrix_elem_end,
# sphere elements sphere_elem_start..sphere_elem_end
m.Material(name='Mat_Matrix')
m.materials['Mat_Matrix'].Elastic(table=((E_matrix, nu_matrix),))
m.HomogeneousSolidSection(name='Sec_Matrix', material='Mat_Matrix', thickness=None)

if Is_Porous in ('Composite', 'Hybrid') and E_inclusion > 0:
    m.Material(name='Mat_Inclusion')
    m.materials['Mat_Inclusion'].Elastic(table=((E_inclusion, nu_inclusion),))
    m.HomogeneousSolidSection(name='Sec_Inclusion', material='Mat_Inclusion', thickness=None)

part_name = list(m.parts.keys())[0]
p = m.parts[part_name]
print("  Part '{}': {} nodes, {} elements".format(part_name, len(p.nodes), len(p.elements)))

# Create Matrix_Only set from element labels
matrix_labels = list(range(matrix_label_start, matrix_label_end + 1))
matrix_elems = p.elements.sequenceFromLabels(matrix_labels)
p.Set(elements=matrix_elems, name='Matrix_Only')
p.SectionAssignment(region=p.sets['Matrix_Only'], sectionName='Sec_Matrix',
    offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
    thicknessAssignment=FROM_SECTION)
print("  Assigned Sec_Matrix to Matrix_Only ({} elements)".format(len(matrix_labels)))

# Create Sphere_Only set from element labels
if (Is_Porous in ('Composite', 'Hybrid')) and E_inclusion > 0 and sphere_label_end >= sphere_label_start and sphere_label_start > 0:
    sphere_labels = list(range(sphere_label_start, sphere_label_end + 1))
    sphere_elems = p.elements.sequenceFromLabels(sphere_labels)
    p.Set(elements=sphere_elems, name='Sphere_Only')
    p.SectionAssignment(region=p.sets['Sphere_Only'], sectionName='Sec_Inclusion',
        offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
        thicknessAssignment=FROM_SECTION)
    print("  Assigned Sec_Inclusion to Sphere_Only ({} elements)".format(len(sphere_labels)))

# ---- Create Reference Points ----
rp1 = a.ReferencePoint(point=(L/2., L/2., L/2.))
a.Set(name='RP-1', referencePoints=(a.referencePoints[rp1.id],))
rp2 = a.ReferencePoint(point=(L/2., L/2., L/2.))
a.Set(name='RP-2', referencePoints=(a.referencePoints[rp2.id],))
rp3 = a.ReferencePoint(point=(L/2., L/2., L/2.))
a.Set(name='RP-3', referencePoints=(a.referencePoints[rp3.id],))
print("Created RP-1, RP-2, RP-3")

# ---- Read periodic node pairs ----
print("Reading periodic pairs from {}...".format(match_path))
pairs = {'X': [], 'Y': [], 'Z': []}
with open(match_path, 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        axis = row['axis']
        neg = int(row['neg_node'])
        pos = int(row['pos_node'])
        pairs[axis].append((neg, pos))

for ax in ['X', 'Y', 'Z']:
    print("  {}: {} pairs".format(ax, len(pairs[ax])))

# ---- Find which instance each node belongs to ----
node_to_inst = {}
valid_labels = set()
for inst_name in a.instances.keys():
    inst = a.instances[inst_name]
    for node in inst.nodes:
        node_to_inst[node.label] = inst_name
        valid_labels.add(node.label)

print("  Total nodes in model: {}".format(len(valid_labels)))

# ---- Apply PBC equations ----
print("Applying direct PBC equations...")
axis_map = {'X': ('RP-1', 1), 'Y': ('RP-2', 2), 'Z': ('RP-3', 3)}
n_eq = 0
n_skipped = 0
constrained_nodes = set()  # track nodes already in equations

for axis in ['X', 'Y', 'Z']:
    rp_name, rp_dof = axis_map[axis]
    
    for idx, (neg_label, pos_label) in enumerate(pairs[axis]):
        # Skip if either node doesn't exist in the model
        if neg_label not in valid_labels or pos_label not in valid_labels:
            n_skipped += 1
            continue
        
        # Skip if either node is already constrained (edge/corner node
        # already handled by a previous face pair)
        if neg_label in constrained_nodes or pos_label in constrained_nodes:
            n_skipped += 1
            continue
        
        neg_inst = node_to_inst.get(neg_label)
        pos_inst = node_to_inst.get(pos_label)
        
        if neg_inst is None or pos_inst is None:
            n_skipped += 1
            continue
        
        neg_set = 'DN_{}N{}'.format(axis, idx)
        pos_set = 'DN_{}P{}'.format(axis, idx)
        
        neg_seq = a.instances[neg_inst].nodes.sequenceFromLabels((neg_label,))
        pos_seq = a.instances[pos_inst].nodes.sequenceFromLabels((pos_label,))
        
        a.Set(name=neg_set, nodes=neg_seq)
        a.Set(name=pos_set, nodes=pos_seq)
        
        for dof in [1, 2, 3]:
            eq_name = 'PBC_{}_{}d{}'.format(axis, idx, dof)
            if dof == rp_dof:
                m.Equation(name=eq_name, terms=(
                    ( 1.0, pos_set, dof),
                    (-1.0, neg_set, dof),
                    (-1.0, rp_name, dof)))
            else:
                m.Equation(name=eq_name, terms=(
                    ( 1.0, pos_set, dof),
                    (-1.0, neg_set, dof)))
            n_eq += 1
        
        constrained_nodes.add(neg_label)
        constrained_nodes.add(pos_label)

print("  Total equations: {}".format(n_eq))
print("  Skipped pairs (duplicate/missing): {}".format(n_skipped))
print("  Constrained nodes: {}".format(len(constrained_nodes)))

# Track which labels were skipped for shear mode reference
skip_labels = set()
for axis in ['X', 'Y', 'Z']:
    for idx, (neg_label, pos_label) in enumerate(pairs[axis]):
        if neg_label in constrained_nodes and pos_label in constrained_nodes:
            pass  # both constrained = this pair was processed
        elif neg_label not in valid_labels or pos_label not in valid_labels:
            skip_labels.add(neg_label)
            skip_labels.add(pos_label)

# ---- Pin centre node for symmetric displacement field ----
best_node = None
best_dist = 1e30
cx, cy, cz = L/2., L/2., L/2.
inst_name = list(a.instances.keys())[0]
for node in a.instances[inst_name].nodes:
    x, y, z = node.coordinates
    d = (x - cx)**2 + (y - cy)**2 + (z - cz)**2
    if d < best_dist:
        best_dist = d
        best_node = node

if best_node:
    seq = a.instances[inst_name].nodes.sequenceFromLabels((best_node.label,))
    a.Set(name='Fix_Ref_Centre', nodes=seq)
    m.DisplacementBC(name='Fix_Ref_Centre', createStepName='Initial',
        region=a.sets['Fix_Ref_Centre'],
        u1=0.0, u2=0.0, u3=0.0,
        ur1=UNSET, ur2=UNSET, ur3=UNSET)
    print("  Pinned centre node {} at ({},{},{})".format(
        best_node.label,
        best_node.coordinates[0], best_node.coordinates[1], best_node.coordinates[2]))

# ---- Step and loading ----
m.StaticStep(name='Step-1', previous='Initial',
    maxNumInc=100000, initialInc=0.1,
    minInc=1e-10, maxInc=0.1, nlgeom=OFF)

m.TabularAmplitude(
    name='LoadRamp', timeSpan=STEP, smooth=SOLVER_DEFAULT,
    data=((0.0, 0.0), (0.1, 0.05), (0.3, 0.2), (0.6, 0.5), (1.0, 1.0)))

m.FieldOutputRequest(name='F-Output-1', createStepName='Step-1',
    variables=('S', 'E', 'U', 'RF', 'EVOL'))
m.HistoryOutputRequest(name='H-Output-1', createStepName='Step-1',
    variables=('ALLSE', 'ALLWK'))

# No general contact needed - single conformal part with shared nodes at interfaces
print("  Single conformal part - no contact needed")

# Apply BCs based on Mode
mode_parts = Mode.strip().split()
direction = mode_parts[-1].upper() if len(mode_parts) > 2 else 'X'

if Mode.startswith('Uniaxial Tension'):
    bc_map = {
        'X': [('RP-1', Disp, UNSET, UNSET, 'LoadRamp'),
              ('RP-2', UNSET, UNSET, UNSET, UNSET),
              ('RP-3', UNSET, UNSET, UNSET, UNSET)],
        'Y': [('RP-1', UNSET, UNSET, UNSET, UNSET),
              ('RP-2', UNSET, Disp, UNSET, 'LoadRamp'),
              ('RP-3', UNSET, UNSET, UNSET, UNSET)],
        'Z': [('RP-1', UNSET, UNSET, UNSET, UNSET),
              ('RP-2', UNSET, UNSET, UNSET, UNSET),
              ('RP-3', UNSET, UNSET, Disp, 'LoadRamp')],
    }
    for rp_n, u1, u2, u3, amp in bc_map[direction]:
        m.DisplacementBC(name=rp_n, createStepName='Step-1',
            region=a.sets[rp_n], u1=u1, u2=u2, u3=u3,
            ur1=UNSET, ur2=UNSET, ur3=UNSET,
            amplitude=amp if amp != UNSET else UNSET,
            fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    print("  Applied Uniaxial Tension {}".format(direction))

elif Mode.startswith('Biaxial Tension'):
    plane = mode_parts[-1].upper() if mode_parts[-1].upper() in ('XY','XZ','YZ') else 'YZ'
    bc_map = {
        'XY': [('RP-1', Disp, UNSET, UNSET, UNSET),
               ('RP-2', UNSET, Disp, UNSET, UNSET),
               ('RP-3', UNSET, UNSET, UNSET, UNSET)],
        'XZ': [('RP-1', Disp, UNSET, UNSET, UNSET),
               ('RP-2', UNSET, UNSET, UNSET, UNSET),
               ('RP-3', UNSET, UNSET, Disp, UNSET)],
        'YZ': [('RP-1', UNSET, UNSET, UNSET, UNSET),
               ('RP-2', UNSET, Disp, UNSET, UNSET),
               ('RP-3', UNSET, UNSET, Disp, UNSET)],
    }
    for rp_n, u1, u2, u3, amp in bc_map[plane]:
        m.DisplacementBC(name=rp_n, createStepName='Step-1',
            region=a.sets[rp_n], u1=u1, u2=u2, u3=u3,
            ur1=UNSET, ur2=UNSET, ur3=UNSET,
            amplitude=amp if amp != UNSET else UNSET,
            fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    print("  Applied Biaxial Tension {}".format(plane))

elif Mode.startswith('Confined Compression'):
    bc_map = {
        'X': [('RP-1', -Disp, 0.0, 0.0, UNSET),
              ('RP-2', 0.0, 0.0, 0.0, UNSET),
              ('RP-3', 0.0, 0.0, 0.0, UNSET)],
        'Y': [('RP-1', 0.0, 0.0, 0.0, UNSET),
              ('RP-2', 0.0, -Disp, 0.0, UNSET),
              ('RP-3', 0.0, 0.0, 0.0, UNSET)],
        'Z': [('RP-1', 0.0, 0.0, 0.0, UNSET),
              ('RP-2', 0.0, 0.0, 0.0, UNSET),
              ('RP-3', 0.0, 0.0, -Disp, UNSET)],
    }
    for rp_n, u1, u2, u3, amp in bc_map[direction]:
        m.DisplacementBC(name=rp_n, createStepName='Step-1',
            region=a.sets[rp_n], u1=u1, u2=u2, u3=u3,
            ur1=UNSET, ur2=UNSET, ur3=UNSET,
            amplitude=amp if amp != UNSET else UNSET,
            fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    print("  Applied Confined Compression {}".format(direction))

elif Mode.startswith('Simple Shear'):
    shear_comp = mode_parts[-1].upper() if mode_parts[-1].upper() in ('S12','S13','S23') else 'S13'
    
    # Create RP-4 for shear
    rp4 = a.ReferencePoint(point=(L/2., L/2., L/2.))
    a.Set(name='RP-4', referencePoints=(a.referencePoints[rp4.id],))
    
    # Fix RP-1/2/3 to zero
    for rp_n in ('RP-1', 'RP-2', 'RP-3'):
        m.DisplacementBC(name=rp_n, createStepName='Step-1',
            region=a.sets[rp_n], u1=0.0, u2=0.0, u3=0.0,
            ur1=UNSET, ur2=UNSET, ur3=UNSET,
            amplitude=UNSET, fixed=OFF,
            distributionType=UNIFORM, fieldName='', localCsys=None)
    
    # Shear PBC: add RP-4 coupling to appropriate face pair
    # S13: DOF 1 jump across Z-faces (RP-4 DOF 1)
    # S23: DOF 2 jump across Z-faces (RP-4 DOF 2)
    # S12: DOF 2 jump across X-faces (RP-4 DOF 2)
    
    shear_config = {
        'S13': {'face_axis': 'Z', 'shear_dof': 1, 'rp4_bc': (Disp, UNSET, 0.0)},
        'S23': {'face_axis': 'Z', 'shear_dof': 2, 'rp4_bc': (UNSET, Disp, 0.0)},
        'S12': {'face_axis': 'X', 'shear_dof': 2, 'rp4_bc': (UNSET, Disp, 0.0)},
    }
    cfg = shear_config[shear_comp]
    face_axis = cfg['face_axis']
    shear_dof = cfg['shear_dof']
    
    # Add shear equations: for each pair on the shear face,
    # add RP-4 coupling on the shear DOF
    # The existing PBC equations have: u_pos_dof - u_neg_dof = 0 (tangential)
    # For shear: replace with: u_pos_dof - u_neg_dof - RP4_dof = 0
    
    # We need to delete the existing tangential equation and replace it
    face_prefix = face_axis + '0'  # e.g. 'Z0'
    n_shear_eq = 0
    
    for axis_label in ['X', 'Y', 'Z']:
        for idx in range(len(pairs.get(axis_label, []))):
            eq_name = 'PBC_{}_{}d{}'.format(axis_label, idx, shear_dof)
            
            # Only modify equations on the shear face axis
            if axis_label != face_axis:
                continue
            
            neg_set = 'DN_{}N{}'.format(axis_label, idx)
            pos_set = 'DN_{}P{}'.format(axis_label, idx)
            
            # Skip if sets don't exist (pair was skipped during initial PBC setup)
            if neg_set not in a.sets.keys() or pos_set not in a.sets.keys():
                continue
            
            # Delete existing equation for this DOF
            if eq_name in m.constraints.keys():
                del m.constraints[eq_name]
            
            # Create new equation with RP-4
            m.Equation(name=eq_name, terms=(
                ( 1.0, pos_set, shear_dof),
                (-1.0, neg_set, shear_dof),
                (-1.0, 'RP-4', shear_dof)))
            n_shear_eq += 1
    
    # RP-4 BC
    u1, u2, u3 = cfg['rp4_bc']
    m.DisplacementBC(name='RP-4', createStepName='Step-1',
        region=a.sets['RP-4'], u1=u1, u2=u2, u3=u3,
        ur1=UNSET, ur2=UNSET, ur3=UNSET,
        amplitude=UNSET, fixed=OFF,
        distributionType=UNIFORM, fieldName='', localCsys=None)
    
    print("  Applied Simple Shear {} ({} shear equations)".format(shear_comp, n_shear_eq))

elif Mode.startswith('Bending'):
    raise RuntimeError(
        "Bending is not supported by this legacy CAE import script; "
        "generate bending jobs with SpaX_Standalone.py (write_complete_inp).")

else:
    print("  WARNING: Unknown mode '{}'".format(Mode))

print("\\nModel setup complete - {}".format(Mode))
print("Submit with: mdb.Job(name='Job-gmsh', model='Model-1').submit()")
''')
    
    print("  Abaqus import script: {}".format(script_path))
    return script_path


# =====================================================================
# Sphere array I/O (for Abaqus <-> system Python data exchange)
# =====================================================================

def export_sphere_array(sphere_array, filepath):
    """Export sphere array to CSV (call from Abaqus or system Python)."""
    np.savetxt(filepath, sphere_array, delimiter=',',
               header='cx,cy,cz,rx,ry,rz,rot1,rot2,rot3,sphericity',
               comments='')
    print("  Exported {} spheres to {}".format(len(sphere_array), filepath))


def load_sphere_array(filepath):
    """Load sphere array from CSV."""
    data = np.loadtxt(filepath, delimiter=',', skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    print("  Loaded {} spheres from {}".format(len(data), filepath))
    return data


# =====================================================================
# Standalone usage
# =====================================================================

if __name__ == '__main__':
    import sys

    # ---- Isolated subprocess mesh mode ----
    # Run one mesh job from a pickled payload and write the result back. The
    # parent (SpaX_Standalone) invokes this in a child process so that a
    # native crash on degenerate geometry (SIGSEGV/SIGABRT inside Gmsh's
    # C++ mesher — uncatchable from Python) is contained here and surfaces to
    # the parent as a non-zero exit code it can catch and retry, instead of
    # taking down the whole parametric run.
    if len(sys.argv) >= 4 and sys.argv[1] == '--_subprocess_mesh':
        import pickle
        in_pkl, out_pkl = sys.argv[2], sys.argv[3]
        with open(in_pkl, 'rb') as _f:
            _p = pickle.load(_f)
        _res = generate_periodic_mesh(
            sphere_array=_p['sphere_array'], L=_p['L'], L_mesh=_p['L_mesh'],
            output_dir=_p['output_dir'], Is_Porous=_p['Is_Porous'],
            run_id=_p['run_id'], VoF_void_sphere=_p['VoF_void_sphere'],
            VoF_incl_sphere=_p['VoF_incl_sphere'],
            Inclusion_Type=_p['Inclusion_Type'],
            gap_balls=_p.get('gap_balls'),
            mesh_order=_p.get('mesh_order', 1))
        with open(out_pkl, 'wb') as _f:
            pickle.dump(_res, _f)
        sys.exit(0)

    import argparse

    parser = argparse.ArgumentParser(description='Generate periodic RVE mesh with Gmsh')
    parser.add_argument('--spheres', type=str, default='',
                        help='Path to sphere array CSV (from Abaqus export)')
    parser.add_argument('--L', type=float, default=0.2, help='RVE side length')
    parser.add_argument('--L_mesh', type=float, default=0.008, help='Element size')
    parser.add_argument('--Disp', type=float, default=0.002, help='Applied displacement')
    parser.add_argument('--E_matrix', type=float, default=9e9)
    parser.add_argument('--nu_matrix', type=float, default=0.3)
    parser.add_argument('--E_inclusion', type=float, default=1e9)
    parser.add_argument('--nu_inclusion', type=float, default=0.35)
    parser.add_argument('--mode', type=str, default='Composite',
                        help='Composite, Porous, or Hybrid')
    parser.add_argument('--loading', type=str, default='Uniaxial Tension X',
                        help='Loading mode: "Uniaxial Tension X/Y/Z", "Simple Shear S12/S13/S23", '
                             '"Biaxial Tension XY/XZ/YZ", "Confined Compression X/Y/Z"')
    parser.add_argument('--output', type=str, default='.',
                        help='Output directory')
    parser.add_argument('--run_id', type=str, default='gmsh_rve')
    parser.add_argument('--VoF_void', type=float, default=0.0,
                        help='Void volume fraction (Hybrid mode)')
    parser.add_argument('--VoF_incl', type=float, default=0.0,
                        help='Inclusion volume fraction (Hybrid mode)')
    
    args = parser.parse_args()
    
    # Load or generate sphere array
    if args.spheres and os.path.isfile(args.spheres):
        sphere_array = load_sphere_array(args.spheres)
    else:
        print("No sphere CSV provided — generating random test spheres")
        np.random.seed(42)
        n_spheres = 10
        sphere_array = np.zeros((n_spheres, 10))
        for i in range(n_spheres):
            r = np.random.uniform(0.015, 0.035)
            cx = np.random.uniform(0, args.L)
            cy = np.random.uniform(0, args.L)
            cz = np.random.uniform(0, args.L)
            sphere_array[i] = [cx, cy, cz, r, r, r, 0, 0, 0, 1.0]
        
        # Save for reference
        csv_path = os.path.join(args.output, '{}_spheres.csv'.format(args.run_id))
        export_sphere_array(sphere_array, csv_path)
    
    # Generate mesh
    result = generate_periodic_mesh(
        sphere_array=sphere_array,
        L=args.L, L_mesh=args.L_mesh,
        output_dir=args.output,
        Is_Porous=args.mode,
        run_id=args.run_id,
        VoF_void_sphere=args.VoF_void,
        VoF_incl_sphere=args.VoF_incl
    )
    
    # Write Abaqus import script
    script = write_abaqus_import_script(
        result=result,
        L=args.L, L_mesh=args.L_mesh, Disp=args.Disp,
        E_matrix=args.E_matrix, nu_matrix=args.nu_matrix,
        E_inclusion=args.E_inclusion, nu_inclusion=args.nu_inclusion,
        Is_Porous=args.mode,
        Mode=args.loading,
        output_dir=args.output
    )
    
    print("\n" + "=" * 70)
    print("WORKFLOW:")
    print("  1. In Abaqus CAE:")
    print("     execfile(r'{}')".format(script))
    print("  2. Submit job:")
    print("     mdb.Job(name='Job-gmsh', model='Model-1').submit()")
    print("=" * 70)
