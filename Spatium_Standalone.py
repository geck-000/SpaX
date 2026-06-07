#!/usr/bin/env python3
"""
Spatium_Standalone.py
=====================
Generate complete Abaqus .inp files for RVE homogenisation
WITHOUT Abaqus CAE. Requires only: Python 3, numpy, gmsh.

    pip install numpy gmsh

Usage:
    python Spatium_Standalone.py parametric_sea_ice_v2.csv /path/to/output

    # Or single RVE:
    python Spatium_Standalone.py --single --L 0.55 --VoF 0.10 --r_avg 0.05 ...

This generates solver-ready .inp files (UTX, SS13, Bending) that can be
submitted directly to Abaqus solver on CSC:
    abaqus job=Job-SI_cold_r1-utx input=Job-SI_cold_r1-utx cpus=8
"""

from __future__ import print_function
import numpy as np
import os
import csv
import math
import sys

# =====================================================================
# OCTREE for fast neighbour queries (from Spatium_Kernel)
# =====================================================================

class OctreeNode:
    def __init__(self, x_min, x_max, y_min, y_max, z_min, z_max,
                 capacity=4, depth=0, max_depth=10):
        self.x_min, self.x_max = x_min, x_max
        self.y_min, self.y_max = y_min, y_max
        self.z_min, self.z_max = z_min, z_max
        self.capacity = capacity
        self.depth = depth
        self.max_depth = max_depth
        self.spheres = []
        self.children = []

    def subdivide(self):
        xm = 0.5*(self.x_min+self.x_max)
        ym = 0.5*(self.y_min+self.y_max)
        zm = 0.5*(self.z_min+self.z_max)
        for x0,x1 in [(self.x_min,xm),(xm,self.x_max)]:
            for y0,y1 in [(self.y_min,ym),(ym,self.y_max)]:
                for z0,z1 in [(self.z_min,zm),(zm,self.z_max)]:
                    self.children.append(OctreeNode(x0,x1,y0,y1,z0,z1,
                        self.capacity,self.depth+1,self.max_depth))

    def insert(self, center, radius):
        cx,cy,cz = center
        if not (self.x_min<=cx<=self.x_max and self.y_min<=cy<=self.y_max
                and self.z_min<=cz<=self.z_max):
            return False
        if len(self.spheres) < self.capacity or self.depth == self.max_depth:
            self.spheres.append((center, radius))
            return True
        if not self.children:
            self.subdivide()
        for child in self.children:
            if child.insert(center, radius):
                return True
        return False

    def query(self, center, radius):
        cx,cy,cz = center
        if (cx+radius < self.x_min or cx-radius > self.x_max or
            cy+radius < self.y_min or cy-radius > self.y_max or
            cz+radius < self.z_min or cz-radius > self.z_max):
            return []
        result = list(self.spheres)
        for child in self.children:
            result.extend(child.query(center, radius))
        return result


class Octree:
    def __init__(self, L, capacity=4, max_depth=10):
        self.root = OctreeNode(-0.2*L, 1.2*L, -0.2*L, 1.2*L, -0.2*L, 1.2*L,
                               capacity, max_depth=max_depth)
    def insert(self, center, radius):
        self.root.insert(center, radius)
    def query(self, center, radius):
        return self.root.query(center, radius)


# =====================================================================
# CHANNEL GENERATION (from Spatium_Kernel)
# =====================================================================

def generate_channels(L, channel_vof_target, r_channel_avg, r_channel_std,
                      min_distance, max_iterations, octree, L_mesh=0.0):
    """Generate vertical cylindrical channels through Z-direction."""
    channel_vof = 0.0
    channels = []  # list of (cx, cy, radius)
    
    for iteration in range(max_iterations):
        if channel_vof >= channel_vof_target:
            break
        
        cx = L * np.random.rand()
        cy = L * np.random.rand()
        radius = max(np.random.normal(r_channel_avg, r_channel_std), L_mesh, 0.001)
        if radius <= 0:
            continue
        
        # Check overlap with existing channels (periodic in XY)
        overlap = False
        for ex, ey, er in channels:
            dx = cx - ex; dx -= L * round(dx / L)
            dy = cy - ey; dy -= L * round(dy / L)
            if math.sqrt(dx*dx + dy*dy) < radius + er + min_distance:
                overlap = True; break
        if overlap: continue
        
        # Check overlap with spheres via octree
        # Query along entire Z-column
        query_r = radius + min_distance
        candidates = octree.query((cx, cy, L/2), max(query_r, L))
        sphere_overlap = False
        for (sc, sr) in candidates:
            dx = cx - sc[0]; dy = cy - sc[1]
            dist_xy = math.sqrt(dx*dx + dy*dy)
            if dist_xy < radius + sr + min_distance:
                sphere_overlap = True; break
        if sphere_overlap: continue
        
        # Boundary portion check: reject channels that create thin slivers
        # near XY boundaries (cap must be substantial)
        eff_r = radius + min_distance / 2.0
        bad_boundary = False
        for coord in [cx, cy]:
            for face_dist in [coord, L - coord]:
                if face_dist < eff_r:
                    # Channel crosses this face
                    cap_depth = eff_r - face_dist
                    body_depth = 2 * eff_r - cap_depth
                    # Reject if cap or body is too thin
                    if cap_depth < min_distance * 0.5 or body_depth < min_distance:
                        bad_boundary = True
                        break
            if bad_boundary: break
        if bad_boundary: continue
        
        # Accept — add periodic copies if needed
        from itertools import product
        x_shifts = [0]; y_shifts = [0]
        if cx - eff_r < 0: x_shifts.append(L)
        if cx + eff_r > L: x_shifts.append(-L)
        if cy - eff_r < 0: y_shifts.append(L)
        if cy + eff_r > L: y_shifts.append(-L)
        
        for dx, dy in product(x_shifts, y_shifts):
            xn, yn = cx + dx, cy + dy
            if not any(abs(xn-e[0])<1e-10 and abs(yn-e[1])<1e-10 for e in channels):
                channels.append((xn, yn, radius))
        
        channel_vof += math.pi * radius**2 * L / L**3
    
    print("  Channels: {} placed, VoF = {:.4f} (target {:.4f})".format(
        len(channels), channel_vof, channel_vof_target))
    
    return np.array(channels) if channels else np.empty((0, 3))


# =====================================================================
# RSA SPHERE PACKING with Octree acceleration
# =====================================================================

def generate_sphere_packing(L, r_avg, r_std, VoF_target, min_distance,
                             max_iterations, sphericity_avg=0.85,
                             sphericity_std=0.1, growth_direction='Random',
                             growth_concentration=0.0, min_radius=None,
                             sliver_gap=0.0):
    """
    Random Sequential Adsorption sphere/ellipsoid packing.
    Returns numpy array shape (N, 10):
        [cx, cy, cz, rx, ry, rz, rot1, rot2, rot3, sphericity]

    min_radius : float or None
        Hard floor on the adaptive radius. Packing stops shrinking below this.
        Pass a mesh-resolvable size (e.g. ~0.5*L_mesh) so the adaptive schedule
        cannot fill the remaining VoF with thousands of sub-element-size
        inclusions that blow up the Gmsh entity count and wreck periodic
        meshing. If None, falls back to 10% of r_avg (geometry-only runs).
    sliver_gap : float
        Mesh-resolvable feature floor (typically ~half a surface element,
        i.e. ~0.5*lc_fine). Any boundary cap, inter-inclusion gap, or
        edge/corner grazing thinner than this is geometrically degenerate —
        the surface mesher produces invalid elements there and Gmsh's C++
        kernel can crash. Rejecting such placements at the packing stage
        prevents the degeneracy at its source (cheaper and more reliable than
        catching the downstream mesher crash). 0 disables the extra guard
        (legacy behaviour, thresholds tied to min_distance only).
    """
    V_RVE = L ** 3
    spheres = []  # list of (cx, cy, cz, rx, ry, rz)
    octree = Octree(L, capacity=8, max_depth=12)
    current_vof = 0.0
    n_consecutive_fails = 0
    current_r = r_avg
    # Floor on the adaptive radius: never shrink below a mesh-resolvable size.
    r_floor = min_radius if (min_radius and min_radius > 0) else r_avg * 0.10
    # Effective surface separation: a feature (cap, inter-sphere gap, edge
    # grazing) below this is an unmeshable sliver. max() so it never weakens
    # the user's physical min_distance.
    sep = max(min_distance, sliver_gap)
    
    for iteration in range(max_iterations):
        if current_vof >= VoF_target:
            break
        
        # Random centre
        cx = np.random.uniform(0, L)
        cy = np.random.uniform(0, L)
        cz = np.random.uniform(0, L)
        
        # Random radius with std (from current_r, not original r_avg)
        r = max(np.random.normal(current_r, min(r_std, current_r * 0.3)),
                current_r * 0.3)
        
        # Boundary handling: prevent shallow face intersections that cause
        # asymmetric Gmsh periodic topology. A sphere must either:
        # (a) NOT intersect any face (centre > r + margin from face), or
        # (b) Intersect DEEPLY (centre < r - margin from face)
        # The "danger zone" (shallow grazing) is eliminated. Width is the
        # mesh-resolvable sep so a grazing cap can never be a sub-element sliver.
        margin = max(min_distance * 2.0, sep)  # minimum intersection depth
        coords = [cx, cy, cz]
        for i in range(3):
            v = coords[i]
            dist_lo = v            # distance to face at 0
            dist_hi = L - v        # distance to face at L

            # Check face at 0
            if abs(dist_lo - r) < margin:
                # In danger zone — push to cross deeply
                coords[i] = r * 0.5
            # Check face at L
            if abs(dist_hi - r) < margin:
                coords[i] = L - r * 0.5
        cx, cy, cz = coords
        
        # Random sphericity
        sph = np.clip(np.random.normal(sphericity_avg, sphericity_std), 0.5, 1.0)
        
        # Compute semi-axes from equivalent radius and sphericity
        r_eq = r
        if sph < 0.999:
            r_long = r_eq / (sph ** (1.0/3.0))
            r_short = r_long * sph
        else:
            r_long = r_eq
            r_short = r_eq
        
        # Assign semi-axes using Von Mises-Fisher orientation distribution
        # Maps growth_concentration to vMF concentration parameter kappa
        if growth_direction != 'Random' and growth_concentration > 0.01:
            # Sample direction biased toward preferred axis via vMF-like distribution
            kappa_vmf = growth_concentration * 30.0  # 0.6 -> kappa=18, 0.8 -> 24
            
            # Preferred axis unit vector
            if growth_direction == 'Z': mu = np.array([0., 0., 1.])
            elif growth_direction == 'X': mu = np.array([1., 0., 0.])
            elif growth_direction == 'Y': mu = np.array([0., 1., 0.])
            else: mu = np.array([0., 0., 1.])
            
            # Sample direction: mix preferred + random, weighted by kappa
            rand_dir = np.random.randn(3)
            rand_dir /= np.linalg.norm(rand_dir)
            d = kappa_vmf * mu + rand_dir
            d /= np.linalg.norm(d)
            
            # Assign semi-axes: r_long along d, r_short perpendicular
            # Project onto axes: the axis closest to d gets r_long
            abs_d = np.abs(d)
            dominant = np.argmax(abs_d)
            axes = [r_short, r_short, r_short]
            axes[dominant] = r_long
            rx, ry, rz = axes
        else:
            # Random orientation
            axes = [r_long, r_short, r_short]
            np.random.shuffle(axes)
            rx, ry, rz = axes
        
        # Sliver rejection (geometric, mesh-aware). Reject placements whose
        # cut geometry would be thinner than one resolvable element (sep):
        #   - a thin spherical cap where the sphere grazes a face,
        #   - a thin body left inside the RVE,
        #   - a shallow grazing of a cube EDGE or CORNER (the sphere surface
        #     skims the edge/corner without crossing it deeply), which leaves a
        #     wedge/needle the surface mesher cannot triangulate.
        r_check = max(rx, ry, rz)
        reject_portion = False
        cc = [cx, cy, cz]
        for i, v in enumerate(cc):
            for face_dist in [v, L - v]:  # distance to each face pair
                if face_dist < r_check:  # sphere crosses this face
                    cap_height = r_check - face_dist
                    if cap_height < sep:
                        reject_portion = True
                        break
                    # Also check the body (portion inside RVE)
                    body_height = 2 * r_check - cap_height
                    if body_height < sep:
                        reject_portion = True
                        break
            if reject_portion:
                break
        # Edge / corner grazing: distance from centre to each cube edge line
        # (two axes pinned to a face value) and each corner point. A sphere is
        # a sliver risk when its surface lies within sep of the feature but it
        # does not engulf the feature deeply (|d - r| < sep).
        if not reject_portion and sliver_gap > 0.0:
            # Edges: pick two distinct axes, each pinned to 0 or L.
            for a in range(3):
                for b in range(a + 1, 3):
                    k = 3 - a - b  # the free axis
                    for fa in (0.0, L):
                        for fb in (0.0, L):
                            d_edge = math.hypot(cc[a] - fa, cc[b] - fb)
                            if abs(d_edge - r_check) < sep:
                                reject_portion = True
                                break
                        if reject_portion: break
                    if reject_portion: break
                if reject_portion: break
        if not reject_portion and sliver_gap > 0.0:
            for fx in (0.0, L):
                for fy in (0.0, L):
                    for fz in (0.0, L):
                        d_corner = math.sqrt((cx - fx) ** 2 + (cy - fy) ** 2
                                             + (cz - fz) ** 2)
                        if abs(d_corner - r_check) < sep:
                            reject_portion = True
                            break
                    if reject_portion: break
                if reject_portion: break
        if reject_portion:
            n_consecutive_fails += 1
            continue
        
        # Check overlap using octree (O(log N) instead of O(N)). The required
        # surface separation is `sep`, so two inclusions are never closer than
        # one resolvable element — eliminating the thin matrix wedge between
        # near-tangent spheres that the volume mesher chokes on.
        overlap = False
        r_check = max(rx, ry, rz)
        query_r = r_check + r_avg * 2.0 + sep

        for dx in [0, L, -L]:
            for dy in [0, L, -L]:
                for dz in [0, L, -L]:
                    candidates = octree.query((cx+dx, cy+dy, cz+dz), query_r)
                    for (sc, sr) in candidates:
                        dist = math.sqrt((cx+dx-sc[0])**2 + (cy+dy-sc[1])**2 + (cz+dz-sc[2])**2)
                        if dist < r_check + sr + sep:
                            overlap = True; break
                    if overlap: break
                if overlap: break
            if overlap: break
        
        if not overlap:
            spheres.append((cx, cy, cz, rx, ry, rz))
            octree.insert((cx, cy, cz), r_check)
            # Insert periodic images into octree too
            for dx in [L, -L, 0]:
                for dy in [L, -L, 0]:
                    for dz in [L, -L, 0]:
                        if dx == 0 and dy == 0 and dz == 0: continue
                        octree.insert((cx+dx, cy+dy, cz+dz), r_check)
            vol = 4.0/3.0 * math.pi * rx * ry * rz
            current_vof += vol / V_RVE
            n_consecutive_fails = 0
        else:
            n_consecutive_fails += 1
            # Adaptive radius reduction (matching kernel behaviour).
            # Try harder at the current radius (25 consecutive misses) before
            # shrinking, and shrink gently (8% per step) down to a 10% floor.
            # RSA jams well below the target VoF; a gentler schedule keeps
            # packing small spheres into the remaining gaps instead of bailing
            # out at ~60% of target, while large-radius placements are exhausted
            # first for a more realistic size distribution.
            if n_consecutive_fails >= 25:
                current_r *= 0.92
                n_consecutive_fails = 0
                print("    [Adaptive] r_avg reduced to {:.5f} after {} placements (VoF={:.1%})".format(
                    current_r, len(spheres), current_vof))
                if current_r < r_floor:
                    print("    [Adaptive] Minimum radius {:.5f} reached, stopping "
                          "(keeps inclusions mesh-resolvable)".format(r_floor))
                    break
    
    # Build sphere array
    N = len(spheres)
    Sphere_array = np.zeros((N, 10))
    for i, (cx, cy, cz, rx, ry, rz) in enumerate(spheres):
        sph = min(rx, ry, rz) / max(rx, ry, rz) if max(rx, ry, rz) > 0 else 1.0
        Sphere_array[i] = [cx, cy, cz, rx, ry, rz, 0, 0, 0, sph]
    
    print("  Packing: {} spheres, VoF = {:.4f} (target {:.4f}, {:.1%} achieved)".format(
        N, current_vof, VoF_target, current_vof/VoF_target if VoF_target > 0 else 1.0))
    
    return Sphere_array, octree


# =====================================================================
# LESICAR EQ.14 INTEGRAL CONSTRAINTS (from Spatium_SecondOrder)
# =====================================================================

def compute_lesicar_constraints(nodes, elements, L, bending_plane='xz'):
    """
    Compute Lesicar Eq.14 integral constraints for the bending .inp.
    
    These constrain the mean fluctuation on each face to zero, removing
    spurious modes that soften the bending response.
    
    Returns list of constraints, each: (dof, [(node_label, weight), ...], rp_E_coeff, rp_K_coeff)
    """
    tol = L * 0.005
    
    # Identify face nodes
    face_nodes = {('x',0):[], ('x',1):[], ('y',0):[], ('y',1):[], ('z',0):[], ('z',1):[]}
    for label, (x, y, z) in nodes.items():
        if x < tol: face_nodes[('x',0)].append(label)
        if x > L-tol: face_nodes[('x',1)].append(label)
        if y < tol: face_nodes[('y',0)].append(label)
        if y > L-tol: face_nodes[('y',1)].append(label)
        if z < tol: face_nodes[('z',0)].append(label)
        if z > L-tol: face_nodes[('z',1)].append(label)
    
    # Compute tributary areas: for each face node, sum up 1/3 of each
    # adjacent face triangle's area
    face_areas = {}  # (face_key) -> {node_label: area}
    
    # C3D4 face definitions (0-indexed): 4 faces per tet
    tet_faces = [(0,1,2), (0,1,3), (1,2,3), (0,2,3)]
    
    for face_key, fnodes in face_nodes.items():
        fnode_set = set(fnodes)
        areas = {n: 0.0 for n in fnodes}
        
        for label, (etype, elset, conn) in elements.items():
            if len(conn) < 4: continue
            # Check each element face
            for fi in tet_faces:
                face_node_labels = [conn[i] for i in fi if i < len(conn)]
                # Check if all 3 face nodes are on this boundary
                if all(n in fnode_set for n in face_node_labels):
                    # Compute triangle area
                    p0 = np.array(nodes[face_node_labels[0]])
                    p1 = np.array(nodes[face_node_labels[1]])
                    p2 = np.array(nodes[face_node_labels[2]])
                    area = 0.5 * np.linalg.norm(np.cross(p1-p0, p2-p0))
                    # Distribute equally to 3 nodes
                    for n in face_node_labels:
                        areas[n] += area / 3.0
        
        face_areas[face_key] = areas
    
    # Build constraints: for each negative face, for each DOF,
    # sum(A_i * u_i) should equal the prescribed displacement integral
    constraints = []
    
    for axis_name, axis_idx in [('x', 0), ('y', 1), ('z', 2)]:
        neg_face = (axis_name, 0)
        areas = face_areas.get(neg_face, {})
        fnodes = face_nodes.get(neg_face, [])
        
        if not fnodes or not areas:
            continue
        
        # Zeroth-moment: sum A_i * r_dof = 0 for each DOF
        for dof in [1, 2, 3]:
            node_weights = [(n, areas[n]) for n in fnodes if areas.get(n, 0) > 1e-20]
            if not node_weights:
                continue
            
            # Compute RHS: sum A_i * u_prescribed
            # u_prescribed depends on bending_plane and node coordinates
            rp_E_coeff = 0.0
            rp_K_coeff = 0.0
            
            for n, a in node_weights:
                x = nodes[n]
                x_c = [x[0]-L/2, x[1]-L/2, x[2]-L/2]
                
                if bending_plane == 'xz':
                    if dof == 1:
                        rp_E_coeff += a * x_c[0]  # membrane: x1
                        rp_K_coeff += a * (-x_c[0]*x_c[2])  # bending: -x1*x3
                    elif dof == 3:
                        rp_K_coeff += a * 0.5*x_c[0]**2  # bending: x1^2/2
                elif bending_plane == 'yz':
                    if dof == 2:
                        rp_E_coeff += a * x_c[1]
                        rp_K_coeff += a * (-x_c[1]*x_c[2])
                    elif dof == 3:
                        rp_K_coeff += a * 0.5*x_c[1]**2
                elif bending_plane == 'xy':
                    if dof == 1:
                        rp_E_coeff += a * x_c[0]
                        rp_K_coeff += a * (-x_c[0]*x_c[1])
                    elif dof == 2:
                        rp_K_coeff += a * 0.5*x_c[0]**2
            
            constraints.append({
                'face': neg_face,
                'dof': dof,
                'nodes': node_weights,
                'rp_E_coeff': rp_E_coeff,
                'rp_K_coeff': rp_K_coeff,
            })
    
    return constraints


# =====================================================================
# COMPLETE .INP WRITER
# =====================================================================

def write_complete_inp(gmsh_inp_path, pairs_csv_path, output_inp_path,
                       L, E_matrix, nu_matrix, E_incl, nu_incl,
                       Is_Porous, mode, disp, Inclusion_Type='Solid',
                       matrix_label_range=(1, 1), sphere_label_range=(0, 0),
                       bending_plane='xz', kappa=0.0,
                       nlgeom='OFF'):
    """
    Read the partial Gmsh .inp (nodes + elements) and the periodic pairs CSV,
    then write a complete solver-ready Abaqus .inp file.
    
    Parameters
    ----------
    gmsh_inp_path : str
        Path to the Gmsh-generated .inp (nodes + elements only)
    pairs_csv_path : str
        Path to the periodic pairs CSV (axis, neg_node, pos_node, ...)
    output_inp_path : str
        Path for the complete output .inp
    L : float
        RVE side length
    E_matrix, nu_matrix : float
        Matrix material properties
    E_incl, nu_incl : float
        Inclusion material properties (already converted from K/G if Liquid)
    Is_Porous : str
        'Porous', 'Composite', or 'Hybrid'
    mode : str
        'utx' (uniaxial tension X), 'ss13' (simple shear), 'bend' (bending)
    disp : float
        Applied displacement
    Inclusion_Type : str
        'Solid' or 'Liquid'
    matrix_label_range : tuple
        (start, end) element labels for matrix
    sphere_label_range : tuple
        (start, end) element labels for inclusions
    """
    
    # ---- Read Gmsh .inp for nodes and elements ----
    nodes = {}       # label -> (x, y, z)
    elements = {}    # label -> (type, elset, [node_labels])
    current_elset = None
    current_type = None
    section = None
    
    with open(gmsh_inp_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('**'):
                continue
            
            if line.startswith('*NODE'):
                section = 'node'
                continue
            elif line.startswith('*ELEMENT'):
                section = 'element'
                # Parse TYPE and ELSET
                parts = line.split(',')
                for p in parts:
                    p = p.strip()
                    if p.upper().startswith('TYPE='):
                        current_type = p.split('=')[1].strip()
                    if p.upper().startswith('ELSET='):
                        current_elset = p.split('=')[1].strip()
                continue
            elif line.startswith('*'):
                section = None
                continue
            
            if section == 'node':
                parts = line.split(',')
                label = int(parts[0])
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                nodes[label] = (x, y, z)
            
            elif section == 'element':
                parts = [int(p.strip()) for p in line.split(',')]
                label = parts[0]
                conn = parts[1:]
                elements[label] = (current_type, current_elset, conn)
    
    n_nodes = len(nodes)
    n_elements = len(elements)
    print("    Read {} nodes, {} elements from Gmsh .inp".format(n_nodes, n_elements))
    
    # ---- Read periodic pairs ----
    pairs = {'X': [], 'Y': [], 'Z': []}
    with open(pairs_csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            axis = row['axis']
            neg = int(row['neg_node'])
            pos = int(row['pos_node'])
            pairs[axis].append((neg, pos))
    
    n_pairs = sum(len(v) for v in pairs.values())
    print("    Read {} periodic pairs (X:{}, Y:{}, Z:{})".format(
        n_pairs, len(pairs['X']), len(pairs['Y']), len(pairs['Z'])))
    
    # ---- Find centre node (closest to L/2, L/2, L/2) ----
    centre = np.array([L/2, L/2, L/2])
    best_node = None
    best_dist = 1e30
    tol = L * 0.01
    for label, (x, y, z) in nodes.items():
        # Exclude boundary nodes
        if x < tol or x > L-tol or y < tol or y > L-tol or z < tol or z > L-tol:
            continue
        d = math.sqrt((x-centre[0])**2 + (y-centre[1])**2 + (z-centre[2])**2)
        if d < best_dist:
            best_dist = d
            best_node = label
    
    print("    Centre node: {} at dist {:.4f} from centre".format(best_node, best_dist))
    
    # ---- Reference point nodes (assembly-level) ----
    rp_base = n_nodes + 100  # offset to avoid conflicts
    rp_nodes = {
        'RP-1': (rp_base + 1, L, L/2, L/2),
        'RP-2': (rp_base + 2, L/2, L, L/2),
        'RP-3': (rp_base + 3, L/2, L/2, L),
        'RP-4': (rp_base + 4, L/2, L/2, L/2),
    }
    
    # ---- Determine element type ----
    elem_type = 'C3D4H' if Inclusion_Type == 'Liquid' else 'C3D4'
    
    # ---- Determine loading ----
    # step_name is consumed at the *Step card below; step_desc only labels the
    # header comment. The actual PBC equations and step BCs for every mode are
    # written further down (see the assembly/loading sections), so this block
    # just needs to name the step for all supported modes.
    step_names = {
        'utx': 'Step-UTX',  'uty': 'Step-UTY',  'utz': 'Step-UTZ',
        'ss12': 'Step-SS12', 'ss13': 'Step-SS13', 'ss23': 'Step-SS23',
    }
    step_descs = {
        'utx': 'Uniaxial Tension along X',
        'uty': 'Uniaxial Tension along Y',
        'utz': 'Uniaxial Tension along Z',
        'ss12': 'Simple Shear S12',
        'ss13': 'Simple Shear S13',
        'ss23': 'Simple Shear S23',
    }
    if mode in step_names:
        step_name = step_names[mode]
        step_desc = step_descs[mode]
    elif mode == 'bend':
        step_name = 'Step-Bending'
        step_desc = 'Second-Order Bending ({})'.format(bending_plane)
    else:
        raise ValueError(
            "Unsupported mode: {}. Use one of utx/uty/utz/ss12/ss13/ss23 or bend.".format(mode))
    
    # ---- Bending-specific setup ----
    # For bending mode, we need coordinate-dependent PBC equations
    # instead of the standard RP-based ones
    is_bending = (mode == 'bend')
    
    if is_bending:
        # Compute coordinate-dependent coefficients for each pair
        # RP_E: membrane strain (free, DOF 1)
        # RP_K: curvature (prescribed, DOF 1 = Kappa)
        rp_base_bend = n_nodes + 200
        rp_E_id = rp_base_bend + 1
        rp_K_id = rp_base_bend + 2
        
        rp_nodes_bend = {
            'RP_E': (rp_E_id, 2.5*L, 2.0*L, 2.0*L),
            'RP_K': (rp_K_id, 2.0*L, 2.0*L, 2.0*L),
        }
    
    # ---- Write complete .inp ----
    inst_name = 'PART-1-1'
    
    with open(output_inp_path, 'w') as f:
        # ============================================================
        # HEADER
        # ============================================================
        f.write('*Heading\n')
        f.write('** Spatium RVE - {} (standalone)\n'.format(step_desc))
        f.write('*Preprint, echo=NO, model=NO, history=NO, contact=NO\n')
        f.write('**\n')
        
        # ============================================================
        # PART (nodes, elements, sections — no materials here)
        # ============================================================
        f.write('*Part, name=PART-1\n')
        
        f.write('*Node\n')
        for label in sorted(nodes.keys()):
            x, y, z = nodes[label]
            f.write('{}, {}, {}, {}\n'.format(label, x, y, z))
        
        m_start, m_end = matrix_label_range
        f.write('*Element, type={}, elset=Matrix_Only\n'.format(elem_type))
        for label in sorted(elements.keys()):
            etype, elset, conn = elements[label]
            if elset == 'Matrix_Only' or (m_start <= label <= m_end):
                conn_str = ', '.join(str(n) for n in conn)
                f.write('{}, {}\n'.format(label, conn_str))
        
        s_start, s_end = sphere_label_range
        has_spheres = s_start > 0 and s_end >= s_start
        if has_spheres:
            f.write('*Element, type={}, elset=Sphere_Only\n'.format(elem_type))
            for label in sorted(elements.keys()):
                etype, elset, conn = elements[label]
                if elset == 'Sphere_Only' or (s_start <= label <= s_end):
                    conn_str = ', '.join(str(n) for n in conn)
                    f.write('{}, {}\n'.format(label, conn_str))
        
        f.write('*Elset, elset=Matrix_Only, generate\n')
        f.write('{}, {}, 1\n'.format(m_start, m_end))
        if has_spheres:
            f.write('*Elset, elset=Sphere_Only, generate\n')
            f.write('{}, {}, 1\n'.format(s_start, s_end))
        
        f.write('*Solid Section, elset=Matrix_Only, material=Mat_Matrix\n')
        f.write(',\n')
        if has_spheres:
            f.write('*Solid Section, elset=Sphere_Only, material=Mat_Inclusion\n')
            f.write(',\n')
        
        f.write('*End Part\n')
        f.write('**\n')
        
        # ============================================================
        # ASSEMBLY (instance, RPs, node sets, equations)
        # ============================================================
        f.write('*Assembly, name=Assembly\n')
        f.write('**\n')
        f.write('*Instance, name={}, part=PART-1\n'.format(inst_name))
        f.write('*End Instance\n')
        f.write('**\n')
        
        # Assembly-level reference point nodes
        # Always create RP-1 through RP-4
        for rp_name, (rp_id, rx, ry, rz) in rp_nodes.items():
            f.write('*Node\n')
            f.write('{}, {}, {}, {}\n'.format(rp_id, rx, ry, rz))
            f.write('*Nset, nset={}\n'.format(rp_name))
            f.write('{},\n'.format(rp_id))
        
        if is_bending:
            for rp_name, (rp_id, rx, ry, rz) in rp_nodes_bend.items():
                f.write('*Node\n')
                f.write('{}, {}, {}, {}\n'.format(rp_id, rx, ry, rz))
                f.write('*Nset, nset={}\n'.format(rp_name))
                f.write('{},\n'.format(rp_id))
        
        # Centre node set (for pinning)
        f.write('*Nset, nset=Fix_Ref_Centre, instance={}\n'.format(inst_name))
        f.write('{},\n'.format(best_node))
        f.write('**\n')
        
        # ---- PBC EQUATIONS ----
        # Kernel pattern: each face pair has ONLY its normal DOF through
        # the face RP. Other DOFs get 2-term equations (zero jump).
        # For shear modes, RP-4 couples the tangential shear DOF.
        
        # Determine shear config
        shear_config = {
            'ss12': {'face': 'X', 'shear_dof': 2},
            'ss13': {'face': 'Z', 'shear_dof': 1},
            'ss23': {'face': 'Z', 'shear_dof': 2},
        }
        is_shear = mode in shear_config
        
        # Normal DOF for each face
        normal_dof = {'X': 1, 'Y': 2, 'Z': 3}
        face_rp = {'X': 'RP-1', 'Y': 'RP-2', 'Z': 'RP-3'}
        
        if is_bending:
            f.write('** Second-Order Bending PBCs\n')
            constrained = set()
            # (node, dof) already eliminated as the leading/dependent term of an
            # *Equation. Abaqus permits a given (node,dof) to lead only ONE
            # equation; the bending PBCs and the Lesicar integral constraints
            # below must not both claim the same one (else: over-constraint).
            used_dep = set()
            eq_count = 0

            for axis in ['X', 'Y', 'Z']:
                for neg, pos in pairs[axis]:
                    if neg in constrained or pos in constrained:
                        continue
                    xn = nodes[neg]; xp = nodes[pos]
                    x1p, x2p, x3p = xp[0]-L/2, xp[1]-L/2, xp[2]-L/2
                    x1n, x2n, x3n = xn[0]-L/2, xn[1]-L/2, xn[2]-L/2
                    
                    if bending_plane == 'xz':
                        coeffs = {
                            1: (x1p-x1n, -(x1p*x3p-x1n*x3n)),
                            2: (0.0, 0.0),
                            3: (0.0, 0.5*(x1p**2-x1n**2)),
                        }
                    elif bending_plane == 'yz':
                        coeffs = {
                            1: (0.0, 0.0),
                            2: (x2p-x2n, -(x2p*x3p-x2n*x3n)),
                            3: (0.0, 0.5*(x2p**2-x2n**2)),
                        }
                    else:  # xy
                        coeffs = {
                            1: (x1p-x1n, -(x1p*x2p-x1n*x2n)),
                            2: (0.0, 0.5*(x1p**2-x1n**2)),
                            3: (0.0, 0.0),
                        }
                    
                    for dof in [1, 2, 3]:
                        c_mem, c_bend = coeffs[dof]
                        n_terms = 2
                        if abs(c_mem) > 1e-15: n_terms += 1
                        if abs(c_bend) > 1e-15: n_terms += 1
                        
                        f.write('*Equation\n{}\n'.format(n_terms))
                        f.write('{}.{}, {}, 1.\n'.format(inst_name, pos, dof))
                        f.write('{}.{}, {}, -1.\n'.format(inst_name, neg, dof))
                        if abs(c_mem) > 1e-15:
                            f.write('RP_E, 1, {}\n'.format(-c_mem))
                        if abs(c_bend) > 1e-15:
                            f.write('RP_K, 1, {}\n'.format(-c_bend))
                        used_dep.add((pos, dof))
                        eq_count += 1
                    constrained.add(neg); constrained.add(pos)
            print("    Bending equations: {}".format(eq_count))
            
            # Lesicar Eq.14 integral constraints (always included)
            f.write('** Lesicar Eq.14 integral constraints\n')
            lesicar_constraints = compute_lesicar_constraints(
                nodes, elements, L, bending_plane)
            n_lesicar = 0
            n_lesicar_skipped = 0
            for con in lesicar_constraints:
                nw = con['nodes']
                if len(nw) < 2: continue
                cdof = con['dof']

                # Abaqus eliminates the FIRST listed node of the equation. The
                # raw node order puts the same node first for every constraint,
                # which would eliminate that (node,dof) repeatedly. Rotate the
                # term list so it LEADS with a node not yet used as a dependent
                # DOF (by a bending PBC eq or an earlier Lesicar eq). An integral
                # constraint spans hundreds of face nodes, so a free one
                # essentially always exists; if none does, the constraint is
                # redundant with those already written, so skip it.
                lead_idx = next((k for k, (lab, _w) in enumerate(nw)
                                 if (lab, cdof) not in used_dep), None)
                if lead_idx is None:
                    n_lesicar_skipped += 1
                    continue
                nw = [nw[lead_idx]] + nw[:lead_idx] + nw[lead_idx + 1:]
                used_dep.add((nw[0][0], cdof))

                # Count terms: node terms + RP_E + RP_K
                n_terms = len(nw)
                if abs(con['rp_E_coeff']) > 1e-20: n_terms += 1
                if abs(con['rp_K_coeff']) > 1e-20: n_terms += 1

                f.write('*Equation\n{}\n'.format(n_terms))
                for node_label, weight in nw:
                    f.write('{}.{}, {}, {}\n'.format(
                        inst_name, node_label, cdof, weight))
                if abs(con['rp_E_coeff']) > 1e-20:
                    f.write('RP_E, 1, {}\n'.format(-con['rp_E_coeff']))
                if abs(con['rp_K_coeff']) > 1e-20:
                    f.write('RP_K, 1, {}\n'.format(-con['rp_K_coeff']))
                n_lesicar += 1
            
            print("    Lesicar constraints: {} ({} redundant skipped)".format(
                n_lesicar, n_lesicar_skipped))
        
        else:
            f.write('** Periodic Boundary Conditions\n')
            eq_count = 0
            n_3term = 0
            n_dropped = 0

            # Abaqus eliminates the FIRST listed DOF of every *Equation (it
            # becomes the dependent DOF). A given (node, dof) may therefore be
            # the leading term of only ONE equation; reusing it in a second
            # equation is an over-constraint that the input processor rejects
            # ("*EQUATION ... node ... used as dependent more than once").
            # Edge nodes are shared by 2 faces and corner nodes by 3, so Gmsh's
            # face-pairwise matching lists them as `pos` in two/three axes. We
            # therefore pick, per DOF, whichever endpoint is still free as the
            # dependent term (flipping the equation's sign if we must lead with
            # `neg` instead of `pos`); if BOTH endpoints are already eliminated
            # the relation is implied by the ones already written, so we drop it.
            used_dep = set()  # (node, dof) already eliminated as a dependent DOF

            for axis in ['X', 'Y', 'Z']:
                n_dof = normal_dof[axis]
                rp = face_rp[axis]

                for neg, pos in pairs[axis]:
                    for dof in [1, 2, 3]:
                        # Check if this DOF gets RP coupling
                        use_rp = None

                        if dof == n_dof:
                            # Normal DOF always through face RP
                            use_rp = rp
                        elif is_shear:
                            cfg = shear_config[mode]
                            if axis == cfg['face'] and dof == cfg['shear_dof']:
                                # Shear DOF through RP-4
                                use_rp = 'RP-4'

                        # Choose the dependent (leading) node: prefer `pos`, fall
                        # back to `neg`, drop if both are already eliminated.
                        if (pos, dof) not in used_dep:
                            dep, indep, sign = pos, neg, 1.0
                        elif (neg, dof) not in used_dep:
                            dep, indep, sign = neg, pos, -1.0
                        else:
                            n_dropped += 1
                            continue
                        used_dep.add((dep, dof))

                        if use_rp:
                            # dep - indep - (sign)*RP = 0  (sign flips with dep)
                            f.write('*Equation\n3\n')
                            f.write('{}.{}, {}, {:g}\n'.format(inst_name, dep, dof, sign))
                            f.write('{}.{}, {}, {:g}\n'.format(inst_name, indep, dof, -sign))
                            f.write('{}, {}, {:g}\n'.format(use_rp, dof, -sign))
                            n_3term += 1
                        else:
                            # 2-term: zero displacement jump
                            f.write('*Equation\n2\n')
                            f.write('{}.{}, {}, {:g}\n'.format(inst_name, dep, dof, sign))
                            f.write('{}.{}, {}, {:g}\n'.format(inst_name, indep, dof, -sign))

                        eq_count += 1

            print("    PBC equations: {} ({} 3-term, {} 2-term, {} redundant dropped)".format(
                eq_count, n_3term, eq_count - n_3term, n_dropped))
        
        f.write('*End Assembly\n')
        f.write('**\n')
        
        # ============================================================
        # MATERIALS (model level, after assembly)
        # ============================================================
        f.write('*Material, name=Mat_Matrix\n')
        f.write('*Elastic\n')
        f.write('{}, {}\n'.format(E_matrix, nu_matrix))
        if has_spheres:
            f.write('*Material, name=Mat_Inclusion\n')
            f.write('*Elastic\n')
            f.write('{}, {}\n'.format(E_incl, nu_incl))
        f.write('**\n')
        
        # ============================================================
        # INITIAL BOUNDARY CONDITIONS (model level, before step)
        # ============================================================
        f.write('** BOUNDARY CONDITIONS\n')
        f.write('*Boundary\n')
        f.write('Fix_Ref_Centre, 1, 1\n')
        f.write('Fix_Ref_Centre, 2, 2\n')
        f.write('Fix_Ref_Centre, 3, 3\n')
        
        if is_bending:
            f.write('RP_E, 2, 6\n')
            f.write('RP_K, 2, 6\n')
        f.write('**\n')
        
        # Amplitude
        f.write('*Amplitude, name=LoadRamp, time=STEP TIME\n')
        f.write('0., 0., 0.01, 0.001, 0.025, 0.005, 0.05, 0.02\n')
        f.write('0.1, 0.05, 0.2, 0.1, 0.36, 0.2, 0.52, 0.4\n')
        f.write('0.68, 0.6, 0.84, 0.8, 1., 1.\n')
        
        # ============================================================
        # STEP
        # ============================================================
        nlgeom_str = 'YES' if nlgeom.upper() in ('ON', 'YES', 'TRUE') else 'NO'
        f.write('*Step, name={}, nlgeom={}, inc=100000\n'.format(step_name, nlgeom_str))
        f.write('*Static\n')
        f.write('0.1, 1., 1e-10, 0.1\n')
        f.write('**\n')
        
        # Step BCs for each loading mode
        f.write('** LOADING\n')
        
        if mode == 'utx':
            f.write('*Boundary, amplitude=LoadRamp\n')
            f.write('RP-1, 1, 1, {}\n'.format(disp))
        elif mode == 'uty':
            f.write('*Boundary, amplitude=LoadRamp\n')
            f.write('RP-2, 2, 2, {}\n'.format(disp))
        elif mode == 'utz':
            f.write('*Boundary, amplitude=LoadRamp\n')
            f.write('RP-3, 3, 3, {}\n'.format(disp))
        elif mode in ('ss12', 'ss13', 'ss23'):
            # Fix all 3 face RPs
            f.write('*Boundary\n')
            f.write('RP-1, 1, 1\nRP-1, 2, 2\nRP-1, 3, 3\n')
            f.write('*Boundary\n')
            f.write('RP-2, 1, 1\nRP-2, 2, 2\nRP-2, 3, 3\n')
            f.write('*Boundary\n')
            f.write('RP-3, 1, 1\nRP-3, 2, 2\nRP-3, 3, 3\n')
            # Drive shear through RP-4
            cfg = shear_config[mode]
            f.write('*Boundary, amplitude=LoadRamp\n')
            f.write('RP-4, {}, {}, {}\n'.format(cfg['shear_dof'], cfg['shear_dof'], disp))
            # Fix other RP-4 DOFs
            for d in [1, 2, 3]:
                if d != cfg['shear_dof']:
                    f.write('RP-4, {}, {}\n'.format(d, d))
        elif mode == 'bend':
            f.write('*Boundary, amplitude=LoadRamp\n')
            f.write('RP_K, 1, 1, {}\n'.format(kappa))
        
        f.write('**\n')
        
        # Output
        f.write('*Output, field, frequency=1\n')
        f.write('*Node Output\nU, RF\n')
        f.write('*Element Output\n')
        f.write('S, E, LE, IVOL\n' if not is_bending else 'S, E, EVOL, COORD\n')
        f.write('*Output, history, frequency=1\n')
        if mode == 'bend':
            f.write('*Node Output, nset=RP_K\n')
        elif is_shear:
            f.write('*Node Output, nset=RP-4\n')
        else:
            f.write('*Node Output, nset={}\n'.format(
                {'utx': 'RP-1', 'uty': 'RP-2', 'utz': 'RP-3'}[mode]))
        f.write('U, RF\n')
        f.write('*End Step\n')
    
    print("    Written: {}".format(output_inp_path))
    return output_inp_path


# =====================================================================
# BATCH PIPELINE
# =====================================================================

def mesh_in_subprocess(sphere_array, L, L_mesh, output_dir, mode, run_id,
                       VoF_void, VoF_incl, Inclusion_Type):
    """Run generate_periodic_mesh in an isolated child process.

    Gmsh's C++ mesher can SIGSEGV/SIGABRT on degenerate inclusion geometry
    (a tiny sphere cap / sub-element sliver that MeshAdapt cannot repair).
    That native crash is uncatchable from Python and would kill the whole
    parametric run. Running the mesh job in a child process contains the
    crash: it becomes a non-zero exit code here, which we raise as a normal
    exception so the caller's retry loop can re-pack and try again.

    Output mesh/match/.inp files are written to output_dir by the child and
    persist on disk; only the small result dict is returned via a pickle.
    Works identically with the system Gmsh on any machine (laptop or HPC).
    """
    import pickle, subprocess, tempfile
    script_dir = os.path.dirname(os.path.abspath(__file__))
    gmsh_module = os.path.join(script_dir, 'Spatium_GmshPeriodic.py')
    payload = dict(sphere_array=sphere_array, L=L, L_mesh=L_mesh,
                   output_dir=output_dir, Is_Porous=mode, run_id=run_id,
                   VoF_void_sphere=VoF_void, VoF_incl_sphere=VoF_incl,
                   Inclusion_Type=Inclusion_Type)
    tmpd = tempfile.mkdtemp(prefix='spax_mesh_')
    in_pkl = os.path.join(tmpd, 'in.pkl')
    out_pkl = os.path.join(tmpd, 'out.pkl')
    try:
        with open(in_pkl, 'wb') as f:
            pickle.dump(payload, f)
        # Hard wall-clock cap: a degenerate boundary can make Gmsh's 2D/3D
        # mesher spin indefinitely (edge recovery loop) rather than crash.
        # Timing out turns that into a catchable failure the retry loop handles.
        mesh_timeout = float(os.environ.get('SPAX_MESH_TIMEOUT', '900'))
        try:
            proc = subprocess.run(
                [sys.executable, gmsh_module, '--_subprocess_mesh',
                 in_pkl, out_pkl], timeout=mesh_timeout)
        except subprocess.TimeoutExpired:
            raise RuntimeError(
                "mesh subprocess exceeded {:.0f}s (degenerate geometry — "
                "stuck in mesher)".format(mesh_timeout))
        if proc.returncode != 0:
            raise RuntimeError(
                "mesh subprocess crashed (exit {}) — degenerate geometry".format(
                    proc.returncode))
        if not os.path.exists(out_pkl):
            raise RuntimeError("mesh subprocess produced no result")
        with open(out_pkl, 'rb') as f:
            return pickle.load(f)
    finally:
        import shutil
        shutil.rmtree(tmpd, ignore_errors=True)


def process_csv(csv_path, output_dir):
    """
    Read a parametric study CSV and generate complete .inp files
    for each RVE (UTX + SS13 load cases).
    """
    import Spatium_GmshPeriodic as SG
    
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    print("=" * 70)
    print("SPATIUM STANDALONE .INP GENERATOR")
    print("=" * 70)
    print("  CSV: {}".format(csv_path))
    print("  Output: {}".format(output_dir))
    print("  RVEs: {}".format(len(rows)))
    
    os.makedirs(output_dir, exist_ok=True)
    
    for idx, params in enumerate(rows):
        run_id = params['run_id']
        print("\n" + "=" * 70)
        print("[{}/{}] {}".format(idx + 1, len(rows), run_id))
        print("=" * 70)
        
        # Parse parameters
        L = float(params['L'])
        L_mesh = float(params['L_mesh'])
        r_avg = float(params['r_avg'])
        r_std = float(params.get('r_std', r_avg * 0.2) or r_avg * 0.2)
        VoF = float(params['VoF_sphere'])
        E_matrix = float(params['E_matrix'])
        nu_matrix = float(params['nu_matrix'])
        Is_Porous = params['Is_Porous']
        Disp = float(params['Disp'])
        min_dist = float(params.get('min_distance', 0.002) or 0.002)
        max_iter = int(params.get('max_iterations', 200000) or 200000)
        sph_avg = float(params.get('sphericity_avg', 0.85) or 0.85)
        sph_std = float(params.get('sphericity_std', 0.1) or 0.1)
        growth_dir = params.get('Growth_Direction', 'Random')
        growth_conc = float(params.get('Growth_Concentration', 0.0) or 0.0)
        
        Inclusion_Type = params.get('Inclusion_Type', 'Solid')
        K_incl = float(params.get('K_inclusion', 0) or 0)
        G_incl = float(params.get('G_inclusion', 0) or 0)
        VoF_void = float(params.get('VoF_void_sphere', 0) or 0)
        VoF_incl = float(params.get('VoF_incl_sphere', 0) or 0)
        
        # Compute inclusion material
        if Inclusion_Type == 'Liquid' and K_incl > 0 and G_incl > 0:
            E_incl = 9.0 * K_incl * G_incl / (3.0 * K_incl + G_incl)
            nu_incl = (3.0 * K_incl - 2.0 * G_incl) / (2.0 * (3.0 * K_incl + G_incl))
            print("  Liquid: K={:.3e}, G={:.3e} -> E={:.3e}, nu={:.6f}".format(
                K_incl, G_incl, E_incl, nu_incl))
        else:
            E_incl = float(params.get('E_sphere_inclusion', E_matrix) or E_matrix)
            nu_incl = float(params.get('nu_sphere_inclusion', nu_matrix) or nu_matrix)
        
        # Determine Gmsh mode
        if VoF_void > 0 and VoF_incl > 0:
            gmsh_mode = 'Hybrid'
        elif Is_Porous == 'Porous':
            gmsh_mode = 'Porous'
        else:
            gmsh_mode = 'Composite'
        
        # Step 1: Generate sphere packing
        # Floor the adaptive radius at a mesh-resolvable size. Below ~0.75*L_mesh
        # the packing emits spheres/caps too small for the surface mesher to
        # resolve, which produces degenerate facets that crash Gmsh's 3D mesher
        # (a C-level core dump, not a catchable exception). A sweep over seeds
        # showed 0.6*L_mesh segfaults intermittently while 0.75*L_mesh meshes
        # reliably; we floor here for robustness (VoF ~0.15-0.17 on the sea-ice
        # cases). Lower it only together with stronger thin-cap/sliver rejection.
        # Couple the floor to the inclusion mesh size: a sphere is cleanly
        # meshable when its radius spans a few surface elements (lc_fine).
        # min_radius = FLOOR_MULT * lc_fine, with lc_fine = LC_FINE_MULT*L_mesh.
        lc_fine_mult = float(os.environ.get('SPAX_LC_FINE_MULT', '0.4'))
        floor_mult = float(os.environ.get('SPAX_FLOOR_MULT', str(0.75 / lc_fine_mult)))
        r_floor = max(r_avg * 0.10, floor_mult * lc_fine_mult * L_mesh)
        # Keep the inter-inclusion gap mesh-resolvable to avoid matrix slivers
        # (a sub-lc_fine gap makes degenerate facets that crash the 3D mesher).
        gap_mult = float(os.environ.get('SPAX_GAP_MULT', '0.0'))
        if gap_mult > 0.0:
            min_dist = max(min_dist, gap_mult * lc_fine_mult * L_mesh)
        # Geometric sliver rejection, applied as an ESCALATION rather than
        # globally. Refusing thin cut-geometry (face cap, inter-sphere gap,
        # edge/corner grazing) at the packing stage removes the degeneracy that
        # crashes Gmsh's mesher — but it also costs packing density (~6% VoF at
        # full strength). So we keep it OFF for the first attempts (full VoF for
        # the seeds that mesh fine) and ramp it in only once a packing keeps
        # failing, trading density for meshability solely to rescue an otherwise
        # skipped row. A/B over 10 seeds: full-on turned 9/10 -> 10/10 by
        # rescuing one pathological seed but cost every seed ~6% VoF; escalation
        # captures the rescue without taxing the seeds that never needed it.
        #
        # SPAX_SLIVER_MULT = max strength (units of lc_fine, default 0.5).
        # SPAX_SLIVER_START = failed attempts before escalation begins (def 2).
        sliver_mult = float(os.environ.get('SPAX_SLIVER_MULT', '0.5'))
        sliver_start = int(os.environ.get('SPAX_SLIVER_START', '2'))
        sliver_gap_max = sliver_mult * lc_fine_mult * L_mesh
        max_retries = int(os.environ.get('SPAX_MAX_RETRIES', '6'))

        def sliver_for_attempt(n):
            # n = 0-based attempt index. 0 until sliver_start, then ramp
            # linearly to sliver_gap_max on the final attempt.
            if sliver_gap_max <= 0.0 or n < sliver_start:
                return 0.0
            span = max(1, (max_retries - 1) - sliver_start)
            return min(1.0, (n - sliver_start + 1) / float(span)) * sliver_gap_max

        sliver_gap = sliver_for_attempt(0)
        print("  Step 1: Sphere packing...")
        Sphere_array, octree = generate_sphere_packing(
            L=L, r_avg=r_avg, r_std=r_std,
            VoF_target=VoF, min_distance=min_dist,
            max_iterations=max_iter,
            sphericity_avg=sph_avg, sphericity_std=sph_std,
            growth_direction=growth_dir,
            growth_concentration=growth_conc,
            min_radius=r_floor, sliver_gap=sliver_gap)

        # Step 1b: Generate channels if requested
        gen_channels = params.get('generate_channels', 'No').strip().lower() in ('yes', 'true', '1')
        channel_vof_target = float(params.get('channel_vof_target', 0) or 0)
        channel_array = np.empty((0, 3))
        
        if gen_channels and channel_vof_target > 0.001:
            r_ch_avg = float(params.get('r_channel_avg', r_avg * 0.5) or r_avg * 0.5)
            r_ch_std = float(params.get('r_channel_std', r_ch_avg * 0.2) or r_ch_avg * 0.2)
            channel_array = generate_channels(
                L=L, channel_vof_target=channel_vof_target,
                r_channel_avg=r_ch_avg, r_channel_std=r_ch_std,
                min_distance=min_dist, max_iterations=max_iter,
                octree=octree, L_mesh=L_mesh)
            
            # Convert channels to sphere array entries (cylinders as tall ellipsoids)
            # GmshPeriodic handles cylinders via the sphere array with rz >> rx,ry
            for i in range(channel_array.shape[0]):
                ch_x, ch_y, ch_r = channel_array[i]
                # Represent as very tall ellipsoid spanning full Z
                ch_entry = np.array([[ch_x, ch_y, L/2, ch_r, ch_r, L/2, 0, 0, 0, 0.01]])
                Sphere_array = np.vstack((Sphere_array, ch_entry))
        
        # Step 2: Gmsh periodic mesh (with retry on mesh failure)
        print("  Step 2: Gmsh mesh...")
        
        # Each retry re-packs and re-meshes in an isolated subprocess, so a
        # crash or an intermittent periodic-matching failure costs only one
        # attempt and never kills the run. Retries are cheap and safe, so we
        # allow several (max_retries, set above) to push the success rate up,
        # with sliver rejection escalating per attempt (sliver_for_attempt).
        result = None
        for attempt in range(max_retries):
            try:
                result = mesh_in_subprocess(
                    sphere_array=Sphere_array,
                    L=L, L_mesh=L_mesh,
                    output_dir=output_dir,
                    mode=gmsh_mode,
                    run_id=run_id,
                    VoF_void=VoF_void,
                    VoF_incl=VoF_incl,
                    Inclusion_Type=Inclusion_Type)
                
                # Check for empty mesh
                n_total = result.get('n_elements_matrix', 0) + result.get('n_elements_sphere', 0)
                if n_total == 0:
                    raise RuntimeError("Empty mesh (0 elements)")

                # Strict-periodicity gate: a non-zero skipped count means some
                # boundary nodes have no periodic partner in a meshed volume —
                # i.e. material on one face maps to a void on the opposite face
                # (a grazing-void asymmetry). Those nodes would get no PBC
                # equation, so the RVE is not strictly periodic. Treat it as a
                # soft failure: re-pack with escalated sliver rejection (which
                # rejects the grazing void) rather than emit a defective RVE.
                n_skip = result.get('n_pairs_skipped', 0)
                if n_skip > 0:
                    raise RuntimeError(
                        "{} periodic pair(s) unmatched — void/boundary not "
                        "strictly periodic".format(n_skip))

                break  # success
            except Exception as e:
                print("    Mesh attempt {}/{} failed: {}".format(
                    attempt + 1, max_retries, str(e)[:80]))
                if attempt < max_retries - 1:
                    # Escalate sliver rejection for the next packing: stays 0
                    # for the first SPAX_SLIVER_START retries (full VoF), then
                    # ramps up to rescue a stubborn seed at a density cost.
                    sliver_gap = sliver_for_attempt(attempt + 1)
                    if sliver_gap > 0.0:
                        print("    Retrying with new packing "
                              "(sliver rejection on, gap={:.5f})...".format(
                                  sliver_gap))
                    else:
                        print("    Retrying with new packing...")
                    np.random.seed(np.random.randint(0, 100000))
                    Sphere_array, octree = generate_sphere_packing(
                        L=L, r_avg=r_avg, r_std=r_std,
                        VoF_target=VoF, min_distance=min_dist * 1.2,
                        max_iterations=max_iter,
                        sphericity_avg=sph_avg, sphericity_std=sph_std,
                        growth_direction=growth_dir,
                        growth_concentration=growth_conc,
                        min_radius=r_floor, sliver_gap=sliver_gap)
                    result = None
                else:
                    print("    SKIPPING {} after {} failures".format(run_id, max_retries))
                    result = None
        
        if result is None:
            continue
        
        gmsh_inp = result['mesh_file']
        pairs_csv = result['match_file']
        m_range = result['matrix_label_range']
        s_range = result['sphere_label_range']
        
        # Validate mesh is not empty
        n_total = result.get('n_elements_matrix', 0) + result.get('n_elements_sphere', 0)
        if n_total == 0:
            print("  ERROR: Empty mesh (0 elements). Skipping {}".format(run_id))
            continue
        
        # Step 3: Write .inp files (all in output_dir, flat)
        print("  Step 3: Writing .inp files...")
        
        # Parse loading modes from CSV
        Mode = params.get('Mode', 'Uniaxial Tension X')
        Mode2 = params.get('Mode2', '')
        Kappa_val = float(params.get('Kappa', 0.0) or 0.0)
        bp = params.get('Bending_Plane', 'xz')
        full_tensor = params.get('full_tensor', 'No').strip().lower() in ('yes', 'true', '1')
        nlgeom_flag = params.get('nlgeom_flag', 'OFF')
        
        # Common kwargs for all write_complete_inp calls
        common = dict(
            L=L, E_matrix=E_matrix, nu_matrix=nu_matrix,
            E_incl=E_incl, nu_incl=nu_incl,
            Is_Porous=Is_Porous, Inclusion_Type=Inclusion_Type,
            matrix_label_range=m_range, sphere_label_range=s_range,
            nlgeom=nlgeom_flag)
        
        def mode_short(m):
            m = m.lower()
            if 'uniaxial' in m and 'x' in m: return 'utx'
            if 'uniaxial' in m and 'y' in m: return 'uty'
            if 'uniaxial' in m and 'z' in m: return 'utz'
            if 'shear' in m and '12' in m: return 'ss12'
            if 'shear' in m and '13' in m: return 'ss13'
            if 'shear' in m and '23' in m: return 'ss23'
            if 'uniaxial' in m: return 'utx'
            if 'shear' in m: return 'ss13'
            return 'utx'
        
        if full_tensor:
            all_modes = ['utx', 'uty', 'utz', 'ss12', 'ss13', 'ss23']
            print("    Full tensor mode: generating {} load cases".format(len(all_modes)))
            for ms in all_modes:
                path = os.path.join(output_dir, 'Job-{}-{}.inp'.format(run_id, ms))
                write_complete_inp(gmsh_inp, pairs_csv, path,
                    mode=ms, disp=Disp, **common)
        else:
            ms = mode_short(Mode)
            path1 = os.path.join(output_dir, 'Job-{}-{}.inp'.format(run_id, ms))
            write_complete_inp(gmsh_inp, pairs_csv, path1,
                mode=ms, disp=Disp, **common)
            
            if Mode2:
                Disp2 = float(params.get('Disp2', Disp) or Disp)
                ms2 = mode_short(Mode2)
                path2 = os.path.join(output_dir, 'Job-{}-{}.inp'.format(run_id, ms2))
                write_complete_inp(gmsh_inp, pairs_csv, path2,
                    mode=ms2, disp=Disp2, **common)
        
        # Bending (always if Kappa > 0)
        if Kappa_val > 0:
            path_ben = os.path.join(output_dir, 'Job-{}-ben.inp'.format(run_id))
            write_complete_inp(gmsh_inp, pairs_csv, path_ben,
                mode='bend', disp=0, bending_plane=bp, kappa=Kappa_val, **common)

        # Remove the Gmsh mesh .inp and periodic-pairs .csv now that they have
        # been folded into the solver-ready Job-*.inp files. The output_dir is
        # meant to hold only the Job .inp files; these per-run intermediates
        # ({run_id}_gmsh.inp, {run_id}_periodic_pairs.csv) would otherwise be
        # left behind alongside them.
        for intermediate in (gmsh_inp, pairs_csv):
            try:
                if intermediate and os.path.exists(intermediate):
                    os.remove(intermediate)
            except OSError as e:
                print("    Warning: could not remove intermediate {}: {}".format(
                    intermediate, e))

        print("  Done: {} ({} nodes, {} elements)".format(
            run_id, result.get('n_nodes', '?'), 
            result.get('n_elements_matrix', 0) + result.get('n_elements_sphere', 0)))
    
    # Final sweep: drop any mesh intermediates left by runs that were skipped
    # after the mesh was already written (e.g. failing the strict-periodicity
    # gate). Only the two known Gmsh patterns are removed, so Job-*.inp files
    # and any unrelated user files in output_dir are untouched.
    import glob
    leftovers = (glob.glob(os.path.join(output_dir, '*_gmsh.inp'))
                 + glob.glob(os.path.join(output_dir, '*_periodic_pairs.csv')))
    for f in leftovers:
        try:
            os.remove(f)
        except OSError as e:
            print("Warning: could not remove leftover {}: {}".format(f, e))

    print("\n" + "=" * 70)
    print("GENERATION COMPLETE: {} RVEs".format(len(rows)))
    print("Output: {}".format(output_dir))
    print("=" * 70)


# =====================================================================
# CLI ENTRY POINT
# =====================================================================

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python Spatium_Standalone.py <csv_file> <output_dir>")
        print("  Generates complete .inp files for Abaqus solver (no CAE needed)")
        print("  Requires: pip install numpy gmsh")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    output_dir = sys.argv[2]
    
    process_csv(csv_file, output_dir)
