#!/usr/bin/env python3
"""
SpaX_Standalone.py
=====================
Generate complete Abaqus .inp files for RVE homogenisation
WITHOUT Abaqus CAE. Requires only: Python 3, numpy, gmsh.

    pip install numpy gmsh

Usage:
    python SpaX_Standalone.py parametric_sea_ice_v2.csv /path/to/output

    # Or single RVE:
    python SpaX_Standalone.py --single --L 0.55 --VoF 0.10 --r_avg 0.05 ...

This generates solver-ready .inp files (UTX, SS13, Bending) that can be
submitted directly to the Abaqus solver:
    abaqus job=Job-SI_cold_r1-utx input=Job-SI_cold_r1-utx cpus=8
"""

from __future__ import print_function
import numpy as np
import os
import csv
import math
import sys

# =====================================================================
# OCTREE for fast neighbour queries (from SpaX_Kernel)
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
# CHANNEL GENERATION (from SpaX_Kernel)
# =====================================================================

def _channel_sphere_cap(cx, cy, octree, L, sep):
    """Max channel radius at (cx,cy) before a vertical channel touches any
    inclusion. ELLIPSOID-AWARE: a Z channel only meets an inclusion's XY
    cross-section, whose extent toward the channel is the *XY* radial
    r_xy(u)=1/sqrt((ux/rx)^2+(uy/ry)^2) -- NOT the bounding sphere max(rx,ry,rz).
    For Z-elongated inclusions the long axis rz is irrelevant to a Z channel, so
    using the bounding sphere would push channels needlessly far away. Periodic
    inclusion images are already in the octree. Returns inf if no inclusion."""
    cap = float('inf')
    for (sc, sax) in octree.query((cx, cy, L / 2.0), L):
        dx = cx - sc[0]; dy = cy - sc[1]
        dxy = math.sqrt(dx*dx + dy*dy)
        if dxy < 1e-12:
            return 0.0
        ux, uy = dx / dxy, dy / dxy
        r_xy = 1.0 / math.sqrt((ux/sax[0])**2 + (uy/sax[1])**2)
        c = dxy - r_xy - sep
        if c < cap:
            cap = c
    return cap


def _densify_channels(prim, L, sep, octree, vof_target, current_vof, A,
                      max_rounds=300, max_step=0.05):
    """Grow + perturb densification for channels (2D analogue of the sphere
    pass). Grow each channel radius into its free XY headroom -- clearance to
    other channels (min-image), to inclusions (ellipsoid-aware XY), and to the
    XY-face danger zone -- then nudge it off its tightest neighbour. Channels
    span the full Z height, so only XY matters."""
    n = len(prim)
    if n == 0 or current_vof >= vof_target:
        return prim, current_vof
    margin = sep
    stagnant = 0
    for _ in range(max_rounds):
        if current_vof >= vof_target:
            break
        vof0 = current_vof
        order = list(range(n)); np.random.shuffle(order)
        grew = False
        for i in order:
            cx, cy, R = prim[i]
            head = max_step * R
            for j in range(n):                       # other channels (min-image)
                if j == i:
                    continue
                ex, ey, er = prim[j]
                dx = cx - ex; dx -= L * round(dx / L)
                dy = cy - ey; dy -= L * round(dy / L)
                g = math.sqrt(dx*dx + dy*dy) - er - R - sep
                if g < head:
                    head = g
            scap = _channel_sphere_cap(cx, cy, octree, L, sep) - R   # inclusions
            if scap < head:
                head = scap
            for c in (cx, cy):                       # XY-face danger zone
                for d in (c, L - c):
                    if R < d:
                        h = (d - margin) - R
                        if h < head:
                            head = h
            if head <= 1e-12:
                continue
            Rn = R + head
            current_vof += math.pi * (Rn*Rn - R*R) / A
            prim[i] = (cx, cy, Rn)
            grew = True
            if current_vof >= vof_target:
                break
        if current_vof >= vof_target:
            break
        moved = False
        for i in order:
            cx, cy, R = prim[i]
            wx = wy = 0.0; worst = float('inf')
            for j in range(n):
                if j == i:
                    continue
                ex, ey, er = prim[j]
                dx = cx - ex; dx -= L * round(dx / L)
                dy = cy - ey; dy -= L * round(dy / L)
                clr = math.sqrt(dx*dx + dy*dy) - er - R
                if clr < worst:
                    worst = clr; wx, wy = dx, dy
            norm = math.sqrt(wx*wx + wy*wy)
            if norm < 1e-12:
                continue
            step = 0.10 * R
            nx = cx + step * wx / norm
            ny = cy + step * wy / norm
            if not (0.0 <= nx <= L and 0.0 <= ny <= L):
                continue
            if (abs(nx - R) < margin or abs((L - nx) - R) < margin or
                    abs(ny - R) < margin or abs((L - ny) - R) < margin):
                continue
            ok = True
            for j in range(n):
                if j == i:
                    continue
                ex, ey, er = prim[j]
                dx = nx - ex; dx -= L * round(dx / L)
                dy = ny - ey; dy -= L * round(dy / L)
                if math.sqrt(dx*dx + dy*dy) < R + er + sep:
                    ok = False; break
            if ok and _channel_sphere_cap(nx, ny, octree, L, sep) >= R:
                prim[i] = (nx, ny, R); moved = True
        if current_vof - vof0 < 1e-6:
            stagnant += 1
            if stagnant >= 6:
                break
        else:
            stagnant = 0
        if not grew and not moved:
            break
    return prim, current_vof


def generate_channels(L, channel_vof_target, r_channel_avg, r_channel_std,
                      min_distance, max_iterations, octree, L_mesh=0.0,
                      densify=True):
    """Generate vertical (Z) cylindrical channels: RSA in XY with adaptive
    radius reduction, then growth+perturbation densification, then periodic XY
    copies. Channel<->inclusion clearance is ellipsoid-aware (see
    `_channel_sphere_cap`)."""
    A = L * L                       # channel VoF = pi r^2 L / L^3 = pi r^2 / L^2
    sep = min_distance
    r_floor = max(L_mesh, 0.001)
    # Wavy-channel tilt (SPAX_CHANNEL_TILT_DEG): a channel leaning by up to `amp`
    # off Z must not cross a vertical face (that would break XY periodicity), so
    # inflate the face danger-zone by amp -> channels stay clear of x/y faces.
    _tilt = float(os.environ.get('SPAX_CHANNEL_TILT_DEG', '0') or 0)
    _amp = (L * math.tan(math.radians(_tilt)) / (2.0 * math.pi)) if _tilt > 1e-6 else 0.0
    margin = sep + _amp
    prim = []                       # primary channels (cx, cy, r), centres in [0,L)
    channel_vof = 0.0
    current_r = max(r_channel_avg, r_floor)
    fails = 0

    for _ in range(max_iterations):
        if channel_vof >= channel_vof_target:
            break
        cx = L * np.random.rand()
        cy = L * np.random.rand()
        radius = max(np.random.normal(current_r, min(r_channel_std, current_r*0.3)),
                     current_r * 0.3, r_floor)

        ok = True
        for ex, ey, er in prim:                 # channel-channel (periodic XY)
            dx = cx - ex; dx -= L * round(dx / L)
            dy = cy - ey; dy -= L * round(dy / L)
            if math.sqrt(dx*dx + dy*dy) < radius + er + sep:
                ok = False; break
        if ok and _channel_sphere_cap(cx, cy, octree, L, sep) < radius:
            ok = False                          # inclusion clearance (XY)
        if ok:                                  # XY-face danger zone
            for c in (cx, cy):
                if abs(c - radius) < margin or abs((L - c) - radius) < margin:
                    ok = False; break
        if not ok:
            fails += 1
            if fails >= 25:
                current_r *= 0.92
                fails = 0
                if current_r < r_floor:
                    break
            continue

        prim.append((cx, cy, radius))
        channel_vof += math.pi * radius**2 / A
        fails = 0

    if densify and channel_vof < channel_vof_target:
        vof_before = channel_vof
        prim, channel_vof = _densify_channels(
            prim, L, sep, octree, channel_vof_target, channel_vof, A)
        print("    [Densify] channel VoF {:.4f} -> {:.4f}".format(vof_before, channel_vof))

    # Expand primaries to periodic XY copies (preserved return contract).
    from itertools import product
    channels = []
    for cx, cy, radius in prim:
        eff_r = radius + sep / 2.0
        x_shifts = [0]; y_shifts = [0]
        if cx - eff_r < 0: x_shifts.append(L)
        if cx + eff_r > L: x_shifts.append(-L)
        if cy - eff_r < 0: y_shifts.append(L)
        if cy + eff_r > L: y_shifts.append(-L)
        for dx, dy in product(x_shifts, y_shifts):
            xn, yn = cx + dx, cy + dy
            if not any(abs(xn-e[0]) < 1e-10 and abs(yn-e[1]) < 1e-10 for e in channels):
                channels.append((xn, yn, radius))

    print("  Channels: {} placed ({} with copies), VoF = {:.4f} (target {:.4f})".format(
        len(prim), len(channels), channel_vof, channel_vof_target))

    return np.array(channels) if channels else np.empty((0, 3))


# =====================================================================
# RSA SPHERE PACKING with Octree acceleration
# =====================================================================

def _ellipsoid_radial(rx, ry, rz, ux, uy, uz):
    """Distance from an axis-aligned ellipsoid centre to its surface along the
    unit direction (ux,uy,uz). Exact for the centre-line interaction that
    dominates packing: r(u) = 1 / sqrt((ux/rx)^2 + (uy/ry)^2 + (uz/rz)^2)."""
    return 1.0 / math.sqrt((ux/rx)**2 + (uy/ry)**2 + (uz/rz)**2)


def _sphere_channel_clear(cx, cy, rx, ry, channels, L, sep):
    """True if an inclusion centred at (cx,cy) with XY equatorial semi-axes
    (rx,ry) clears every vertical (Z) channel by at least `sep`.

    Channels span the full RVE height, so the binding constraint is the
    inclusion's *widest* XY footprint (its equatorial ellipse) -- this is the
    mirror of `_channel_sphere_cap`, enforced from the sphere side so spheres
    packed AFTER the channels avoid them. ELLIPSOID-AWARE: the inclusion extent
    toward a channel along the in-plane direction u is r_xy(u) =
    1/sqrt((ux/rx)^2+(uy/ry)^2), not the bounding circle max(rx,ry). `channels`
    are primaries (centre in [0,L)); periodic images handled by min-image."""
    if not channels:
        return True
    for (chx, chy, R) in channels:
        dx = cx - chx; dx -= L * round(dx / L)
        dy = cy - chy; dy -= L * round(dy / L)
        dxy = math.sqrt(dx*dx + dy*dy)
        if dxy < 1e-12:
            return False
        ux, uy = dx / dxy, dy / dxy
        r_xy = 1.0 / math.sqrt((ux/rx)**2 + (uy/ry)**2)
        if dxy < r_xy + R + sep:
            return False
    return True


# ---------------------------------------------------------------------------
# True (off-axis) ellipsoid-ellipsoid surface gap via GJK support mapping.
#
# The packing/densify hot loop uses the centre-line radial extent
# (`_ellipsoid_radial`), which is exact only along the inter-centre line. Two
# tilted axis-aligned ellipsoids can approach closer OFF that line, leaving a
# matrix sliver thinner than the centre-line gap (measured ~0.78x worst case).
# GJK on the Minkowski difference gives the exact minimum surface distance for
# strictly-convex bodies; the ellipsoid support map is closed-form. Validated
# against a brute-force surface search to machine precision at ~0.3 ms/pair, so
# it is affordable for a one-shot repair over the (broad-phase pruned) near
# pairs -- but NOT for the inner grow loop, hence repair runs post-densify.
# ---------------------------------------------------------------------------

def _gjk_support(cx, cy, cz, A2, B2, dx, dy, dz):
    """Support point of the Minkowski difference (ellipsoid A at origin minus
    ellipsoid B at (cx,cy,cz)) in direction d. A2,B2 are (rx^2,ry^2,rz^2).
    Ellipsoid support in dir d: c + A2 d / sqrt(d.A2 d)."""
    ax, ay, az = A2[0]*dx, A2[1]*dy, A2[2]*dz
    da = math.sqrt(dx*ax + dy*ay + dz*az) or 1e-300
    bx, by, bz = B2[0]*dx, B2[1]*dy, B2[2]*dz
    db = math.sqrt(dx*bx + dy*by + dz*bz) or 1e-300
    return (ax/da - (cx - bx/db), ay/da - (cy - by/db), az/da - (cz - bz/db))


def _cp_seg(a, b):
    abx, aby, abz = b[0]-a[0], b[1]-a[1], b[2]-a[2]
    t = -(a[0]*abx + a[1]*aby + a[2]*abz) / (abx*abx + aby*aby + abz*abz + 1e-300)
    if t <= 0.0:
        return a, (a,)
    if t >= 1.0:
        return b, (b,)
    return (a[0]+t*abx, a[1]+t*aby, a[2]+t*abz), (a, b)


def _cp_tri(a, b, c):
    """Closest point on triangle abc to the origin (Ericson, Real-Time
    Collision Detection)."""
    abx, aby, abz = b[0]-a[0], b[1]-a[1], b[2]-a[2]
    acx, acy, acz = c[0]-a[0], c[1]-a[1], c[2]-a[2]
    d1 = -(abx*a[0] + aby*a[1] + abz*a[2]); d2 = -(acx*a[0] + acy*a[1] + acz*a[2])
    if d1 <= 0 and d2 <= 0:
        return a, (a,)
    d3 = -(abx*b[0] + aby*b[1] + abz*b[2]); d4 = -(acx*b[0] + acy*b[1] + acz*b[2])
    if d3 >= 0 and d4 <= d3:
        return b, (b,)
    vc = d1*d4 - d3*d2
    if vc <= 0 and d1 >= 0 and d3 <= 0:
        v = d1/(d1-d3); return (a[0]+v*abx, a[1]+v*aby, a[2]+v*abz), (a, b)
    d5 = -(abx*c[0] + aby*c[1] + abz*c[2]); d6 = -(acx*c[0] + acy*c[1] + acz*c[2])
    if d6 >= 0 and d5 <= d6:
        return c, (c,)
    vb = d5*d2 - d1*d6
    if vb <= 0 and d2 >= 0 and d6 <= 0:
        w = d2/(d2-d6); return (a[0]+w*acx, a[1]+w*acy, a[2]+w*acz), (a, c)
    va = d3*d6 - d5*d4
    if va <= 0 and (d4-d3) >= 0 and (d5-d6) >= 0:
        w = (d4-d3)/((d4-d3)+(d5-d6))
        return (b[0]+w*(c[0]-b[0]), b[1]+w*(c[1]-b[1]), b[2]+w*(c[2]-b[2])), (b, c)
    den = 1.0/(va+vb+vc); v = vb*den; w = vc*den
    return (a[0]+abx*v+acx*w, a[1]+aby*v+acy*w, a[2]+abz*v+acz*w), (a, b, c)


def _cp_simplex(s):
    """Closest point to origin on a 1-4 vertex simplex; returns (point, reduced
    simplex) -- the GJK distance sub-algorithm."""
    if len(s) == 1:
        return s[0], s
    if len(s) == 2:
        return _cp_seg(s[0], s[1])
    if len(s) == 3:
        return _cp_tri(s[0], s[1], s[2])
    best = None; bp = None; bs = None
    for tri in ((0,1,2), (0,1,3), (0,2,3), (1,2,3)):
        p, sub = _cp_tri(s[tri[0]], s[tri[1]], s[tri[2]])
        d = p[0]*p[0] + p[1]*p[1] + p[2]*p[2]
        if best is None or d < best:
            best = d; bp = p; bs = sub
    return bp, list(bs)


def _gjk_gap_vec(cx, cy, cz, A2, B2, max_iter=30):
    """As _gjk_gap, but also returns the closest point on the Minkowski
    difference (the separation vector). Returns (gap, closest) with closest=None
    when the ellipsoids intersect."""
    dx, dy, dz = -cx, -cy, -cz
    if dx == 0 and dy == 0 and dz == 0:
        dx = 1.0
    s = [_gjk_support(cx, cy, cz, A2, B2, dx, dy, dz)]
    closest = s[0]
    for _ in range(max_iter):
        dx, dy, dz = -closest[0], -closest[1], -closest[2]
        if dx*dx + dy*dy + dz*dz < 1e-20:
            return 0.0, None
        p = _gjk_support(cx, cy, cz, A2, B2, dx, dy, dz)
        if (p[0]*dx + p[1]*dy + p[2]*dz) - (closest[0]*dx + closest[1]*dy + closest[2]*dz) < 1e-10:
            break
        s.append(p)
        closest, s = _cp_simplex(s)
        s = list(s)
        if closest[0]**2 + closest[1]**2 + closest[2]**2 < 1e-20:
            return 0.0, None
    return math.sqrt(closest[0]**2 + closest[1]**2 + closest[2]**2), closest


def _gjk_gap(cx, cy, cz, A2, B2, max_iter=30):
    """Minimum surface-to-surface distance between two axis-aligned ellipsoids
    (A at origin, B at (cx,cy,cz)). 0.0 if they intersect."""
    return _gjk_gap_vec(cx, cy, cz, A2, B2, max_iter)[0]


def _ellipsoid_gap(A, B, L):
    """True min surface gap between inclusions A,B (each (cx,cy,cz,rx,ry,rz)),
    using the nearest periodic image of B."""
    dx = B[0]-A[0]; dx -= L*round(dx/L)
    dy = B[1]-A[1]; dy -= L*round(dy/L)
    dz = B[2]-A[2]; dz -= L*round(dz/L)
    return _gjk_gap(dx, dy, dz, (A[3]**2, A[4]**2, A[5]**2), (B[3]**2, B[4]**2, B[5]**2))


def _point_ellipse_dist(px, py, a, b):
    """Distance from an external point (px,py) to the boundary of an axis-aligned
    ellipse (semi-axes a,b centred at origin). Newton iteration on the boundary
    angle from a robust initial guess; exact to ~1e-12 in a few steps. Used for
    the TRUE channel<->inclusion gap: a vertical channel meets the inclusion's
    XY equatorial ellipse, and the off-axis nearest point on that ellipse is
    closer than the centre-line radial proxy."""
    px, py = abs(px), abs(py)            # first-quadrant symmetry
    if px < 1e-15 and py < 1e-15:
        return min(a, b)
    t = math.atan2(py * a, px * b)       # good initial parameter
    for _ in range(12):
        ct, st = math.cos(t), math.sin(t)
        ex, ey = a * ct, b * st          # boundary point
        # derivative of distance^2 wrt t = 0  <=>  (E-P).E' = 0
        exd, eyd = -a * st, b * ct
        f = (ex - px) * exd + (ey - py) * eyd
        fp = (exd*exd + (ex - px)*(-a*ct)) + (eyd*eyd + (ey - py)*(-b*st))
        if abs(fp) < 1e-18:
            break
        t -= f / fp
        t = min(math.pi/2.0, max(0.0, t))
    ct, st = math.cos(t), math.sin(t)
    return math.hypot(a*ct - px, b*st - py)


def _channel_inclusion_gap(ch, A, L):
    """True surface gap between a vertical channel ch=(chx,chy,R) and inclusion
    A=(cx,cy,cz,rx,ry,rz): the channel sees A's XY equatorial ellipse (rx,ry),
    so it is the point-to-ellipse distance from the channel axis minus R, over
    the nearest periodic image."""
    dx = ch[0]-A[0]; dx -= L*round(dx/L)
    dy = ch[1]-A[1]; dy -= L*round(dy/L)
    return _point_ellipse_dist(dx, dy, A[3], A[4]) - ch[2]


def _channel_inclusion_gap_mid(ch, A, L):
    """True channel<->inclusion gap AND the world-space midpoint of the thin
    matrix sliver, for mesh refinement. The sliver is thinnest at the inclusion
    equator (z=cz): the channel axis sees the inclusion's XY (rx,ry) ellipse, so
    the nearest ellipse point and the facing channel-surface point bound the
    sliver. Returns (gap, (mx,my,mz))."""
    dx = ch[0]-A[0]; dx -= L*round(dx/L)     # channel relative to inclusion
    dy = ch[1]-A[1]; dy -= L*round(dy/L)
    a, b = A[3], A[4]
    px, py = abs(dx), abs(dy)
    sgx = 1.0 if dx >= 0 else -1.0
    sgy = 1.0 if dy >= 0 else -1.0
    if px < 1e-15 and py < 1e-15:
        return _channel_inclusion_gap(ch, A, L), (A[0], A[1], A[2])
    t = math.atan2(py * a, px * b)
    for _ in range(12):
        ct, st = math.cos(t), math.sin(t)
        ex, ey = a * ct, b * st
        exd, eyd = -a * st, b * ct
        f = (ex - px) * exd + (ey - py) * eyd
        fp = (exd*exd + (ex - px)*(-a*ct)) + (eyd*eyd + (ey - py)*(-b*st))
        if abs(fp) < 1e-18:
            break
        t -= f / fp
        t = min(math.pi/2.0, max(0.0, t))
    ct, st = math.cos(t), math.sin(t)
    # nearest ellipse point (inclusion-centred, signs restored)
    enx, eny = sgx * a * ct, sgy * b * st
    # facing channel-surface point: from channel axis toward the ellipse point
    vx, vy = enx - dx, eny - dy
    vn = math.hypot(vx, vy) or 1e-300
    cnx, cny = dx + ch[2] * vx / vn, dy + ch[2] * vy / vn
    mx = A[0] + 0.5 * (enx + cnx)
    my = A[1] + 0.5 * (eny + cny)
    return (vn - ch[2]), (mx, my, A[2])


def _ellipsoid_gap_mid(A, B, L):
    """True min surface gap between inclusions A,B AND the world-space midpoint
    of the closest-approach segment, over the nearest periodic image of B.
    Returns (gap, (mx,my,mz)). The midpoint is the centre of the thin matrix
    sliver, used to place a local mesh-refinement ball. Witness points come from
    the converged GJK separation direction (exact for the smooth single-point
    contact of two ellipsoids)."""
    dx = B[0]-A[0]; dx -= L*round(dx/L)
    dy = B[1]-A[1]; dy -= L*round(dy/L)
    dz = B[2]-A[2]; dz -= L*round(dz/L)
    A2 = (A[3]**2, A[4]**2, A[5]**2)
    B2 = (B[3]**2, B[4]**2, B[5]**2)
    g, closest = _gjk_gap_vec(dx, dy, dz, A2, B2)
    if closest is None:
        return g, None
    # dir = -closest is the separation direction (from A's surface toward B).
    ex, ey, ez = -closest[0], -closest[1], -closest[2]
    na = math.sqrt(A2[0]*ex*ex + A2[1]*ey*ey + A2[2]*ez*ez) or 1e-300
    nb = math.sqrt(B2[0]*ex*ex + B2[1]*ey*ey + B2[2]*ez*ez) or 1e-300
    # closest point on A (A-centred frame) and on B (B at (dx,dy,dz))
    ax, ay, az = A2[0]*ex/na, A2[1]*ey/na, A2[2]*ez/na
    bx, by, bz = dx - B2[0]*ex/nb, dy - B2[1]*ey/nb, dz - B2[2]*ez/nb
    mx = A[0] + 0.5*(ax + bx)
    my = A[1] + 0.5*(ay + by)
    mz = A[2] + 0.5*(az + bz)
    return g, (mx, my, mz)


def _collect_gap_balls(spheres, L, lc_fine, channels=None,
                       thresh_mult=1.0, resolve=0.5, broad=2.5,
                       channel_thresh_mult=1.5):
    """Locate narrow matrix slivers and return mesh-refinement balls for the
    mesher: list of (x, y, z, size). A pair whose TRUE gap falls below the
    threshold gets a ball at the sliver centre forcing the local element size to
    ~resolve*gap, so the mesher puts ~1/resolve elements across the gap and
    resolves it -- letting the packer keep the inclusions TIGHT (full VoF)
    instead of widening the gap (which costs VoF).

    Two sliver types are refined:
      * sphere-sphere (TRUE GJK gap < thresh_mult*lc_fine) -- the off-axis
        slivers that crash the mesher;
      * channel<->inclusion (gap < channel_thresh_mult*lc_fine) -- a vertical
        channel sees the inclusion equator head-on, so its gap == the XY
        centre-line clearance, which the packer holds at exactly one element.
        An exactly-one-element matrix sheet meshes unreliably (overlapping
        facets), so it is refined too.

    Balls near a periodic face are duplicated across it so the size field stays
    periodic (required for matched periodic surface meshes). Returns world-space
    (x,y,z,size); size is floored so the element count stays bounded."""
    n = len(spheres)
    if n < 2 or lc_fine <= 0.0:
        return []
    rad = _ellipsoid_radial
    thresh = thresh_mult * lc_fine
    reach = broad * thresh
    size_floor = 0.18 * lc_fine          # don't refine below this (cost cap)
    balls = []
    for i in range(n):
        xi, yi, zi, rxi, ryi, rzi = spheres[i]
        for j in range(i+1, n):
            xj, yj, zj, rxj, ryj, rzj = spheres[j]
            dx = xi-xj; dx -= L*round(dx/L)
            dy = yi-yj; dy -= L*round(dy/L)
            dz = zi-zj; dz -= L*round(dz/L)
            D = math.sqrt(dx*dx + dy*dy + dz*dz)
            if D < 1e-12:
                continue
            ux, uy, uz = dx/D, dy/D, dz/D
            clg = (D - rad(rxi, ryi, rzi, ux, uy, uz)
                   - rad(rxj, ryj, rzj, ux, uy, uz))
            if clg >= reach:                 # broad-phase prune
                continue
            g, mid = _ellipsoid_gap_mid(spheres[i], spheres[j], L)
            if mid is None or g >= thresh:
                continue
            size = max(size_floor, resolve * g)
            mx = mid[0] % L; my = mid[1] % L; mz = mid[2] % L
            balls.append((mx, my, mz, size))
    # channel<->inclusion slivers
    cthresh = channel_thresh_mult * lc_fine
    for ch in (channels or []):
        for i in range(n):
            g, mid = _channel_inclusion_gap_mid(ch, spheres[i], L)
            if g >= cthresh:
                continue
            size = max(size_floor, resolve * max(g, 0.0))
            balls.append((mid[0] % L, mid[1] % L, mid[2] % L, size))
    if not balls:
        return balls
    # Periodic duplication: a ball within `margin` of a face needs a twin on the
    # opposite face so matching periodic surfaces see the same size field.
    margin = 2.5 * lc_fine
    out = []
    for (x, y, z, sz) in balls:
        shifts_x = [0.0] + ([L] if x < margin else []) + ([-L] if x > L-margin else [])
        shifts_y = [0.0] + ([L] if y < margin else []) + ([-L] if y > L-margin else [])
        shifts_z = [0.0] + ([L] if z < margin else []) + ([-L] if z > L-margin else [])
        for sx in shifts_x:
            for sy in shifts_y:
                for sz_ in shifts_z:
                    out.append((x+sx, y+sy, z+sz_, sz))
    return out


def _repair_offaxis_slivers(spheres, L, gap_target, V_RVE, current_vof,
                            broad=2.5, max_passes=8, channels=None,
                            channel_gap=None):
    """Post-densify repair of off-axis matrix slivers.

    The grow pass keeps inclusion surfaces `sep` apart along the inter-centre
    line, but two tilted ellipsoids can approach closer OFF that line, leaving a
    sub-element sliver the surface mesher chokes on (degenerate facets -> Gmsh
    crash). This pass measures the TRUE surface gap (`_ellipsoid_gap`, GJK) for
    each near pair and, where it falls below `gap_target` -- a mesh-resolvable
    ELEMENT floor (~0.6*lc_fine, NOT the larger centre-line `sep`) -- SHRINKS the
    larger inclusion (uniform 3-axis scale, preserving sphericity and
    orientation) by bisection until the true gap reaches `gap_target`.

    Targeting an absolute element floor rather than `sep` matters: densify packs
    many pairs to a centre-line gap of exactly `sep`, whose true gap is always a
    bit under `sep`, so demanding true>=sep would shrink most of the pack (~15%
    VoF). The mesher only needs a resolvable feature, so we lift just the genuine
    sub-element OUTLIERS and leave the rest, costing almost no volume fraction.

    Shrinking an inclusion only ever WIDENS its gaps, so each ordered sweep is
    monotone and introduces no new violations; a few passes clear the pack. When
    `channels` is given, inclusion<->channel slivers are repaired the same way
    (the channel is fixed; only the inclusion shrinks), using the true
    point-to-ellipse gap. Returns (spheres, current_vof)."""
    n = len(spheres)
    if n < 1 or gap_target <= 0.0:
        return spheres, current_vof
    rad = _ellipsoid_radial
    reach = broad * gap_target
    chans = channels or []
    # Channels need a FULL element of matrix, not the 0.6*lc_fine sphere floor:
    # a thin channel<->inclusion gap runs the full RVE height, so the mesher
    # tiles a tall matrix sheet that collapses to overlapping facets unless it
    # is at least one element thick. Default to ~1.7x the sphere floor.
    cgap = channel_gap if (channel_gap and channel_gap > 0) else gap_target * (1.0/0.6)
    n_fixed = 0
    for _ in range(max_passes):
        viol = []                               # (gap, i, j)  j=-1-k => channel k
        for i in range(n):
            xi, yi, zi, rxi, ryi, rzi = spheres[i]
            for j in range(i+1, n):
                xj, yj, zj, rxj, ryj, rzj = spheres[j]
                dx = xi-xj; dx -= L*round(dx/L)
                dy = yi-yj; dy -= L*round(dy/L)
                dz = zi-zj; dz -= L*round(dz/L)
                D = math.sqrt(dx*dx + dy*dy + dz*dz)
                if D < 1e-12:
                    continue
                ux, uy, uz = dx/D, dy/D, dz/D
                clg = (D - rad(rxi, ryi, rzi, ux, uy, uz)
                       - rad(rxj, ryj, rzj, ux, uy, uz))
                if clg < reach:                 # broad-phase prune
                    g = _ellipsoid_gap(spheres[i], spheres[j], L)
                    if g < gap_target - 1e-9:
                        viol.append((g, i, j))
            for k, ch in enumerate(chans):      # inclusion <-> channel
                g = _channel_inclusion_gap(ch, spheres[i], L)
                if g < cgap - 1e-9:
                    viol.append((g, i, -1 - k))
        if not viol:
            break
        viol.sort()
        touched = set()
        for _g, i, j in viol:
            if i in touched or j in touched:
                continue                        # geometry changed; next pass
            if j < 0:                           # channel: shrink the inclusion
                bi = i
                ch = chans[-1 - j]
                cx, cy, cz, rx, ry, rz = spheres[bi]
                lo, hi = 0.3, 1.0
                for _ in range(24):
                    mid = 0.5*(lo+hi)
                    if _channel_inclusion_gap(
                            ch, (cx, cy, cz, rx*mid, ry*mid, rz*mid), L) >= cgap:
                        hi = mid
                    else:
                        lo = mid
            else:
                bi = i if max(spheres[i][3:6]) >= max(spheres[j][3:6]) else j
                other = spheres[j] if bi == i else spheres[i]
                cx, cy, cz, rx, ry, rz = spheres[bi]
                lo, hi = 0.3, 1.0               # scale on bi's semi-axes
                for _ in range(24):
                    mid = 0.5*(lo+hi)
                    trial = (cx, cy, cz, rx*mid, ry*mid, rz*mid)
                    if _ellipsoid_gap(trial, other, L) >= gap_target:
                        hi = mid
                    else:
                        lo = mid
            s = hi
            cx, cy, cz, rx, ry, rz = spheres[bi]
            current_vof += (4.0/3.0*math.pi*(rx*ry*rz)*(s*s*s - 1.0)) / V_RVE
            spheres[bi] = (cx, cy, cz, rx*s, ry*s, rz*s)
            touched.add(i)
            if j >= 0:
                touched.add(j)
            n_fixed += 1
    if n_fixed:
        print("    [Off-axis] shrank {} inclusion(s) to clear sub-element "
              "slivers (sphere gap >= {:.5f}, channel gap >= {:.5f})".format(
                  n_fixed, gap_target, cgap))
    return spheres, current_vof


def _densify_packing(spheres, L, sep, min_distance, VoF_target,
                     current_vof, V_RVE, max_rounds=400, max_step=0.05,
                     channels=None, z_bias=0.0, grow_axis=2):
    """Ellipsoid-aware growth + perturbation densification of a jammed RSA pack.

    Pure RSA jams well below the target volume fraction (~60-65% of target for
    the inclusion sizes used here) because it never rearranges. Inclusions are
    axis-aligned ellipsoids (semi-axes rx,ry,rz along x,y,z; no rotation), so
    overlap, headroom and boundary clearance are measured with the *radial*
    centre-line extent r(u) (see `_ellipsoid_radial`) rather than the bounding
    sphere -- elongated inclusions can interlock and reach target instead of
    being held a bounding-sphere apart.

    Each round: (1) GROW pass scales every inclusion's three semi-axes by a
    single factor (so sphericity AND orientation are preserved) into the radial
    clearance to its tightest neighbour and up to the per-axis boundary
    danger-zone cap; (2) PERTURB pass nudges each inclusion off the neighbour
    with least radial clearance to open room for the next grow pass
    (Jodrey-Tory). Growth is sequential with immediate write-back, so each
    inclusion is sized against the current geometry of all others -- closing the
    full clearance keeps every pair >= sep apart. Rounds repeat until the target
    VoF is reached or progress stalls. O(N^2) per round; N is small (~100-300).

    z_bias : float, default 0.0 (isotropic)
        Anisotropic-growth strength. With z_bias=w>0 the grow pass elongates
        each inclusion preferentially along `grow_axis` (Z by default) as it
        fills headroom, instead of scaling all three axes equally -- mimicking
        the vertically elongated brine channels of sea ice. The biased axis
        grows as s**(1+w) and the other two as s**(1-w/2), so the volume gain is
        still s**3 for the same overall factor s (the distribution changes, not
        the budget). Because the radial extent is then non-linear in s, the
        per-inclusion factor is found by bisection against the same neighbour /
        boundary / channel constraints. w>0 LOWERS the final sphericity (a
        deliberate modelling choice); w=0 keeps the original uniform behaviour.
    grow_axis : int, default 2
        Preferred elongation axis when z_bias>0 (0=x, 1=y, 2=z).
    """
    n = len(spheres)
    if n == 0 or current_vof >= VoF_target:
        return spheres, current_vof
    margin = max(min_distance * 2.0, sep)
    rad = _ellipsoid_radial

    stagnant = 0
    for _ in range(max_rounds):
        if current_vof >= VoF_target:
            break
        vof_round_start = current_vof
        order = list(range(n))
        np.random.shuffle(order)
        # --- grow pass ---
        grew = False
        for i in order:
            cx, cy, cz, rx, ry, rz = spheres[i]
            if z_bias <= 0.0:
                # Uniform growth: radial extent scales linearly with the single
                # factor s, so the binding factor is a closed-form min over
                # neighbour / boundary / channel constraints.
                s_max = 1.0 + max_step           # cap growth rate per round
                for j in range(n):
                    if j == i:
                        continue
                    jx, jy, jz, jrx, jry, jrz = spheres[j]
                    dx = cx - jx; dx -= L * round(dx / L)
                    dy = cy - jy; dy -= L * round(dy / L)
                    dz = cz - jz; dz -= L * round(dz / L)
                    D = math.sqrt(dx*dx + dy*dy + dz*dz)
                    if D < 1e-12:
                        s_max = 1.0; break
                    ux, uy, uz = dx/D, dy/D, dz/D
                    avail = D - rad(jrx, jry, jrz, ux, uy, uz) - sep   # room for I
                    if avail <= 0.0:
                        s_max = 1.0; break
                    s_j = avail / rad(rx, ry, rz, ux, uy, uz)
                    if s_j < s_max:
                        s_max = s_j
                if s_max <= 1.0 + 1e-9:
                    continue
                # Per-axis boundary cap: the extent toward a face is the semi-axis
                # on that face's axis. A face the inclusion does not cross caps
                # growth at `margin` short of it; a crossed face only deepens it.
                for c, r_ax in ((cx, rx), (cy, ry), (cz, rz)):
                    for d in (c, L - c):
                        if r_ax < d:
                            cap = (d - margin) / r_ax
                            if cap < s_max:
                                s_max = cap
                if channels:
                    # Channel cap: the grown equatorial ellipse must stay `sep`
                    # clear of every vertical channel (XY-only -- a Z channel
                    # spans the full height).
                    for (chx, chy, R) in channels:
                        dx = cx - chx; dx -= L * round(dx / L)
                        dy = cy - chy; dy -= L * round(dy / L)
                        dxy = math.sqrt(dx*dx + dy*dy)
                        if dxy < 1e-12:
                            s_max = 1.0; break
                        ux, uy = dx / dxy, dy / dxy
                        avail = dxy - R - sep
                        if avail <= 0.0:
                            s_max = 1.0; break
                        s_c = avail / (1.0 / math.sqrt((ux/rx)**2 + (uy/ry)**2))
                        if s_c < s_max:
                            s_max = s_c
                if s_max <= 1.0 + 1e-9:
                    continue
                nrx, nry, nrz = rx*s_max, ry*s_max, rz*s_max
            else:
                # Anisotropic growth: elongate `grow_axis` as s**(1+w), the other
                # two as s**(1-w/2) (volume gain still s**3). The radial extent is
                # non-linear in s, so bisect the largest s that still satisfies
                # every constraint.
                w = z_bias
                ex = [1.0 - 0.5*w, 1.0 - 0.5*w, 1.0 - 0.5*w]
                ex[grow_axis] = 1.0 + w

                def _axes(s, _rx=rx, _ry=ry, _rz=rz, _ex=ex):
                    return _rx*s**_ex[0], _ry*s**_ex[1], _rz*s**_ex[2]

                def _fits(s):
                    ax, ay, az = _axes(s)
                    for j in range(n):
                        if j == i:
                            continue
                        jx, jy, jz, jrx, jry, jrz = spheres[j]
                        dx = cx - jx; dx -= L * round(dx / L)
                        dy = cy - jy; dy -= L * round(dy / L)
                        dz = cz - jz; dz -= L * round(dz / L)
                        D = math.sqrt(dx*dx + dy*dy + dz*dz)
                        if D < 1e-12:
                            return False
                        ux, uy, uz = dx/D, dy/D, dz/D
                        if (D - rad(jrx, jry, jrz, ux, uy, uz) - sep
                                < rad(ax, ay, az, ux, uy, uz)):
                            return False
                    for c, r_o, R in ((cx, rx, ax), (cy, ry, ay), (cz, rz, az)):
                        for d in (c, L - c):
                            if r_o < d and R > d - margin:
                                return False
                    if channels:
                        for (chx, chy, Rc) in channels:
                            dx = cx - chx; dx -= L * round(dx / L)
                            dy = cy - chy; dy -= L * round(dy / L)
                            dxy = math.sqrt(dx*dx + dy*dy)
                            if dxy < 1e-12:
                                return False
                            ux, uy = dx / dxy, dy / dxy
                            avail = dxy - Rc - sep
                            if avail <= 0.0:
                                return False
                            if 1.0 / math.sqrt((ux/ax)**2 + (uy/ay)**2) > avail:
                                return False
                    return True

                hi = 1.0 + max_step
                if _fits(hi):
                    s = hi
                else:
                    if not _fits(1.0 + 1e-6):
                        continue                 # no room to grow at all
                    lo = 1.0
                    for _ in range(20):
                        mid = 0.5*(lo + hi)
                        if _fits(mid):
                            lo = mid
                        else:
                            hi = mid
                    s = lo
                if s <= 1.0 + 1e-9:
                    continue
                nrx, nry, nrz = _axes(s)
            current_vof += (4.0/3.0*math.pi*(nrx*nry*nrz - rx*ry*rz)) / V_RVE
            spheres[i] = (cx, cy, cz, nrx, nry, nrz)
            grew = True
            if current_vof >= VoF_target:
                break
        if current_vof >= VoF_target:
            break
        # --- perturb pass ---
        moved = False
        for i in order:
            cx, cy, cz, rx, ry, rz = spheres[i]
            wx = wy = wz = 0.0
            worst = float('inf')             # least radial clearance neighbour
            for j in range(n):
                if j == i:
                    continue
                jx, jy, jz, jrx, jry, jrz = spheres[j]
                dx = cx - jx; dx -= L * round(dx / L)
                dy = cy - jy; dy -= L * round(dy / L)
                dz = cz - jz; dz -= L * round(dz / L)
                D = math.sqrt(dx*dx + dy*dy + dz*dz)
                if D < 1e-12:
                    continue
                ux, uy, uz = dx/D, dy/D, dz/D
                clr = D - rad(rx, ry, rz, ux, uy, uz) - rad(jrx, jry, jrz, ux, uy, uz)
                if clr < worst:
                    worst = clr; wx, wy, wz = dx, dy, dz
            norm = math.sqrt(wx*wx + wy*wy + wz*wz)
            if norm < 1e-12:
                continue
            step = 0.10 * max(rx, ry, rz)
            nx = cx + step * wx / norm
            ny = cy + step * wy / norm
            nz = cz + step * wz / norm
            if not (0.0 <= nx <= L and 0.0 <= ny <= L and 0.0 <= nz <= L):
                continue
            bad = False                       # per-axis face danger-zone
            for c, r_ax in ((nx, rx), (ny, ry), (nz, rz)):
                if abs(c - r_ax) < margin or abs((L - c) - r_ax) < margin:
                    bad = True; break
            if bad:
                continue
            ok = True                          # radial separation to neighbours
            for j in range(n):
                if j == i:
                    continue
                jx, jy, jz, jrx, jry, jrz = spheres[j]
                dx = nx - jx; dx -= L * round(dx / L)
                dy = ny - jy; dy -= L * round(dy / L)
                dz = nz - jz; dz -= L * round(dz / L)
                D = math.sqrt(dx*dx + dy*dy + dz*dz)
                if D < 1e-12:
                    ok = False; break
                ux, uy, uz = dx/D, dy/D, dz/D
                if D < rad(rx, ry, rz, ux, uy, uz) + rad(jrx, jry, jrz, ux, uy, uz) + sep:
                    ok = False; break
            if ok and channels and not _sphere_channel_clear(nx, ny, rx, ry, channels, L, sep):
                ok = False
            if ok:
                spheres[i] = (nx, ny, nz, rx, ry, rz)
                moved = True
        if current_vof - vof_round_start < 1e-5:
            stagnant += 1
            if stagnant >= 6:
                break
        else:
            stagnant = 0
        if not grew and not moved:
            break
    return spheres, current_vof


def generate_sphere_packing(L, r_avg, r_std, VoF_target, min_distance,
                             max_iterations, sphericity_avg=0.85,
                             sphericity_std=0.1, growth_direction='Random',
                             growth_concentration=0.0, min_radius=None,
                             sliver_gap=0.0, densify=True, channels=None,
                             offaxis_floor=None, offaxis_channel_floor=None,
                             z_bias=0.0):
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
    offaxis_floor : float or None
        Mesh-resolvable TRUE-gap floor (~0.6*lc_fine) for the post-densify
        off-axis sliver repair. After densification, any inclusion pair whose
        exact (GJK) surface gap falls below this is shrunk just enough to clear
        it -- catching the off-axis slivers the centre-line grow cap misses.
        None falls back to `sep` (strict; shrinks more, costs more VoF).
    channels : list of (chx, chy, R) or None
        Vertical (Z) channels already placed (primaries, centre in [0,L)). When
        given, inclusions are packed and grown to clear every channel's XY
        footprint by `sep` (ellipsoid-aware, mirrors `_channel_sphere_cap`). The
        channel-first ordering lets channels claim their XY space before the
        spheres densify and saturate the plane; passing them here keeps the
        spheres out of the channels. None for runs without channels.
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
        
        # Random sphericity. `sphericity` here is the ratio of the two equal
        # semi-axes to the distinct one, so s<1 is a PROLATE needle (one long
        # axis, two short) and s>1 an OBLATE plate (one short axis, two long).
        # The ceiling used to sit at 1.0, which made plates unreachable: every
        # inclusion was a needle however the deck was written. Columnar sea ice
        # carries much of its basal brine in LAYERS between ice platelets, which
        # are oblate with the short axis horizontal, so the range must extend
        # above 1 for that morphology to be expressible at all.
        # Upper limit is generous rather than physical: high-aspect plates are
        # what the sliver rejection below and the mesher decide on, not this
        # clip. Basal brine layers are reported at aspect ratios of order ten.
        sph = np.clip(np.random.normal(sphericity_avg, sphericity_std), 0.3, 12.0)

        # Compute semi-axes from equivalent radius and sphericity. The algebra
        # is unchanged for s<1, so prolate packings are bit-identical to before;
        # for s>1 it simply returns r_long < r_short, and the assignment below
        # then puts the SHORT axis along the growth direction, giving a plate
        # whose plane contains the other two axes.
        # Volume-preserving: the ellipsoid must have the volume of a sphere of
        # radius r_eq whatever its shape, so that r_avg means the size the deck
        # asks for and the inclusion COUNT is not a function of sphericity.
        #
        #   distinct axis a = r_eq s^(-2/3),  pair axis b = r_eq s^(1/3)
        #   a b^2 = r_eq^3  for every s,      and b/a = s  (the aspect ratio)
        #
        # The previous form (a = r_eq s^(-1/3), b = a s) gave V = V_sphere * s:
        # 40% under-volume at the s=0.6 of the warm slices and, once s>1 became
        # reachable, eight times over-volume at s=8. The packer iterates on the
        # volume FRACTION so the achieved VoF was still met, and effective
        # properties are scale-invariant at fixed fraction and shape, so earlier
        # prolate results are unaffected; what was wrong is the size and hence
        # the number of inclusions per cell, which is what sets the realisation
        # scatter -- and, for plates, the morphology outright.
        r_eq = r
        if abs(sph - 1.0) > 1e-3:
            r_long = r_eq * (sph ** (-2.0/3.0))   # distinct axis
            r_short = r_eq * (sph ** (1.0/3.0))   # the two equal axes
        else:
            r_long = r_eq
            r_short = r_eq
        
        # Assign semi-axes using Von Mises-Fisher orientation distribution
        # Maps growth_concentration to vMF concentration parameter kappa
        if growth_direction != 'Random' and growth_concentration > 0.01:
            # Orientation is drawn from a von Mises-Fisher distribution about the
            # preferred axis and then snapped to the nearest Cartesian axis,
            # because the mesher builds axis-aligned ellipsoids only. The snap
            # means the fabric is represented as a MIXTURE of three aligned
            # populations whose weights the concentration controls; it is not a
            # continuous orientation distribution, which is a limitation of the
            # geometry kernel rather than of the sampling.
            kappa_vmf = growth_concentration * 30.0

            if growth_direction == 'Z': mu = np.array([0., 0., 1.])
            elif growth_direction == 'X': mu = np.array([1., 0., 0.])
            elif growth_direction == 'Y': mu = np.array([0., 1., 0.])
            else: mu = np.array([0., 0., 1.])

            # Wood's exact vMF sampler, in the form that closes for p=3.
            # The previous draw, d = kappa*mu + unit_random, was NOT a vMF
            # sample: with kappa >= 3 the mu term dominates a unit vector
            # absolutely, so every draw landed within arcsin(1/kappa) of the
            # preferred axis and the snap below always returned it. Measured
            # over 200k draws the aligned fraction was 1.0000 at every
            # concentration from 0.1 to 0.9 -- the parameter did nothing.
            u = np.random.rand()
            w = 1.0 + (1.0 / kappa_vmf) * np.log(
                u + (1.0 - u) * np.exp(-2.0 * kappa_vmf))
            # a unit vector in the plane orthogonal to mu
            tmp = np.array([1.0, 0.0, 0.0])
            if abs(np.dot(tmp, mu)) > 0.9:
                tmp = np.array([0.0, 1.0, 0.0])
            e1 = np.cross(mu, tmp); e1 /= np.linalg.norm(e1)
            e2 = np.cross(mu, e1)
            phi = 2.0 * np.pi * np.random.rand()
            v = np.cos(phi) * e1 + np.sin(phi) * e2
            d = w * mu + np.sqrt(max(0.0, 1.0 - w * w)) * v
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
        # The cap left at a face is governed by the semi-axis NORMAL to that
        # face, not by the bounding radius. Testing every face against
        # max(rx,ry,rz) overestimates the cap height on the two short axes of a
        # prolate inclusion -- a genuine sliver reads as comfortably thick and
        # is accepted, and the mesher then fails to build a consistent periodic
        # boundary. (Harmless while inclusions were meshed as bounding spheres,
        # since then the bounding radius *was* the geometry.)
        r_check = max(rx, ry, rz)
        r_axis = (rx, ry, rz)
        reject_portion = False
        cc = [cx, cy, cz]
        for i, v in enumerate(cc):
            r_ax = r_axis[i]
            for face_dist in [v, L - v]:  # distance to each face pair
                if face_dist < r_ax:  # inclusion crosses this face
                    cap_height = r_ax - face_dist
                    if cap_height < sep:
                        reject_portion = True
                        break
                    # Also check the body (portion inside RVE)
                    body_height = 2 * r_ax - cap_height
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

        # The octree payload is each inclusion's semi-axes (rx,ry,rz); RSA's
        # placement test uses the bounding sphere max(axes) (unchanged), while
        # channel placement reads the axes for an ellipsoid-aware XY test.
        for dx in [0, L, -L]:
            for dy in [0, L, -L]:
                for dz in [0, L, -L]:
                    candidates = octree.query((cx+dx, cy+dy, cz+dz), query_r)
                    for (sc, sax) in candidates:
                        sr = max(sax)
                        dist = math.sqrt((cx+dx-sc[0])**2 + (cy+dy-sc[1])**2 + (cz+dz-sc[2])**2)
                        if dist < r_check + sr + sep:
                            overlap = True; break
                    if overlap: break
                if overlap: break
            if overlap: break

        # Channel clearance: an inclusion packed after the channels must keep its
        # XY equatorial footprint `sep` clear of every vertical channel (mirror
        # of the channel<->inclusion test, enforced here from the sphere side).
        if not overlap and channels and not _sphere_channel_clear(
                cx, cy, rx, ry, channels, L, sep):
            overlap = True

        if not overlap:
            spheres.append((cx, cy, cz, rx, ry, rz))
            octree.insert((cx, cy, cz), (rx, ry, rz))
            # Insert periodic images into octree too
            for dx in [L, -L, 0]:
                for dy in [L, -L, 0]:
                    for dz in [L, -L, 0]:
                        if dx == 0 and dy == 0 and dz == 0: continue
                        octree.insert((cx+dx, cy+dy, cz+dz), (rx, ry, rz))
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
    
    # Densify: grow the jammed RSA packing into the empty headroom to approach
    # the target VoF (RSA alone jams ~35-40% short). The octree is then rebuilt
    # from the grown radii so downstream consumers (channel placement, meshing)
    # see the true inclusion sizes rather than the pre-growth ones.
    if densify and current_vof < VoF_target:
        vof_before = current_vof
        grow_axis = {'X': 0, 'Y': 1, 'Z': 2}.get(
            str(growth_direction).strip().upper()[:1], 2)
        spheres, current_vof = _densify_packing(
            spheres, L, sep, min_distance, VoF_target, current_vof, V_RVE,
            channels=channels, z_bias=z_bias, grow_axis=grow_axis)
        # Off-axis sliver repair: the grow pass separates surfaces along the
        # centre line, but tilted ellipsoids can sit closer off it; this shrinks
        # the few offending pairs to a true (GJK) gap of `sep` so the mesher
        # never sees a sub-element matrix sliver. On by default; SPAX_OFFAXIS=0
        # restores the centre-line-only behaviour.
        if os.environ.get('SPAX_OFFAXIS', '1') != '0':
            gap_target = offaxis_floor if (offaxis_floor and offaxis_floor > 0) else sep
            spheres, current_vof = _repair_offaxis_slivers(
                spheres, L, gap_target, V_RVE, current_vof, channels=channels,
                channel_gap=offaxis_channel_floor)
        octree = Octree(L, capacity=8, max_depth=12)
        for (cx, cy, cz, rx, ry, rz) in spheres:
            octree.insert((cx, cy, cz), (rx, ry, rz))
            for dx in [L, -L, 0]:
                for dy in [L, -L, 0]:
                    for dz in [L, -L, 0]:
                        if dx == 0 and dy == 0 and dz == 0:
                            continue
                        octree.insert((cx+dx, cy+dy, cz+dz), (rx, ry, rz))
        print("    [Densify] VoF {:.4f} -> {:.4f} by growing {} inclusions".format(
            vof_before, current_vof, len(spheres)))

    # Land on the target fraction rather than the first value past it.
    #
    # The grow/densify rounds stop as soon as the running VoF crosses the
    # target, so the packing keeps whatever the last growth step added -- an
    # overshoot of several percent, which then propagates into every modulus
    # through the knockdown law. A single uniform scaling of the semi-axes
    # removes it exactly, since VoF goes as the cube of the scale factor.
    #
    # Only ever applied downwards. Shrinking increases every pair clearance and
    # every distance to a face, so it cannot create an overlap or a sliver that
    # the placement tests already accepted; growing could do both, so an
    # undershooting pack (the packer stalled and could not reach the target) is
    # left alone and reported as it stands.
    if VoF_target > 0 and current_vof > VoF_target and spheres:
        f = (VoF_target / current_vof) ** (1.0 / 3.0)
        spheres = [(cx, cy, cz, rx * f, ry * f, rz * f)
                   for (cx, cy, cz, rx, ry, rz) in spheres]
        print("    [Trim] VoF {:.4f} -> {:.4f} by scaling semi-axes x{:.4f}".format(
            current_vof, VoF_target, f))
        current_vof = VoF_target

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
# LESICAR EQ.14 INTEGRAL CONSTRAINTS (from SpaX_SecondOrder)
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
                
                if bending_plane == 'torsion':
                    # twist about z: u_x = -alpha*z*y, u_y = +alpha*z*x
                    if dof == 1:
                        rp_K_coeff += a * (-x_c[2]*x_c[1])
                    elif dof == 2:
                        rp_K_coeff += a * (x_c[2]*x_c[0])
                elif bending_plane == 'xz':
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


def bending_pbc_coeffs(xp, xn, L, bending_plane):
    """Coordinate-dependent RP coupling coefficients for one periodic pair.

    Returns {dof: (c_mem, c_bend)} so that the second-order periodicity relation
    on each pair is  u_pos[dof] - u_neg[dof] = c_mem*RP_E + c_bend*RP_K , where
    RP_E is the macroscopic membrane strain and RP_K the prescribed curvature.
    Coordinates are taken relative to the RVE centre (L/2).

    `bending_plane='torsion'` selects the torsion mode instead, for which the
    macroscopic field is a rigid twist about z at rate alpha,

        u_x = -alpha * z * y,   u_y = +alpha * z * x,   u_z = 0,

    so RP_K carries alpha and there is no membrane term. Torsion is the
    canonical couple-stress probe and, unlike bending, imposes no plate-like
    kinematic on the cube -- which is the point of having it: the cube-versus-
    plate extraction bias that dominates the bending control should be absent.
    The field is not periodic, but its DIFFERENCE across a face pair is, which
    is the same property the bending modes rely on.
    """
    x1p, x2p, x3p = xp[0]-L/2, xp[1]-L/2, xp[2]-L/2
    x1n, x2n, x3n = xn[0]-L/2, xn[1]-L/2, xn[2]-L/2
    if bending_plane == 'torsion':
        return {
            1: (0.0, -(x3p*x2p - x3n*x2n)),
            2: (0.0,  (x3p*x1p - x3n*x1n)),
            3: (0.0, 0.0),
        }
    if bending_plane == 'xz':
        return {
            1: (x1p-x1n, -(x1p*x3p-x1n*x3n)),
            2: (0.0, 0.0),
            3: (0.0, 0.5*(x1p**2-x1n**2)),
        }
    elif bending_plane == 'yz':
        return {
            1: (0.0, 0.0),
            2: (x2p-x2n, -(x2p*x3p-x2n*x3n)),
            3: (0.0, 0.5*(x2p**2-x2n**2)),
        }
    else:  # xy
        return {
            1: (x1p-x1n, -(x1p*x2p-x1n*x2n)),
            2: (0.0, 0.5*(x1p**2-x1n**2)),
            3: (0.0, 0.0),
        }


def build_bending_pbc_equations(pairs, nodes, L, bending_plane, tol=1e-15):
    """Build the second-order (bending) periodic *Equation set.

    Each periodic node pair contributes, per DOF, the relation
        u_pos[dof] - u_neg[dof] - c_mem*RP_E - c_bend*RP_K = 0.

    Abaqus eliminates the FIRST-listed (node,dof) of every *Equation, so a given
    (node,dof) may lead only one equation. Edge/corner nodes are periodic slaves
    on 2-3 faces, so a single node can be `pos` in more than one axis; we
    therefore pick, per DOF, whichever endpoint is still free as the eliminated
    term (flipping the equation sign when we must lead with `neg`). When BOTH
    endpoints are already eliminated the relation closes a loop of pairs; since
    the macroscopic bending field is a single-valued (path-independent) function
    of position, that loop relation is implied exactly by the equations already
    written and is dropped.

    Returns (equations, used_dep, n_dropped) where each equation is a list of
    terms (is_node, name, dof, coeff): is_node True -> `name` is a node label,
    is_node False -> `name` is 'RP_E'/'RP_K'.
    """
    used_dep = set()   # (node, dof) already eliminated as a dependent DOF
    equations = []
    n_dropped = 0

    for axis in ['X', 'Y', 'Z']:
        for neg, pos in pairs[axis]:
            coeffs = bending_pbc_coeffs(nodes[pos], nodes[neg], L, bending_plane)
            for dof in [1, 2, 3]:
                c_mem, c_bend = coeffs[dof]

                if (pos, dof) not in used_dep:
                    dep, indep, s = pos, neg, 1.0
                elif (neg, dof) not in used_dep:
                    dep, indep, s = neg, pos, -1.0
                else:
                    n_dropped += 1
                    continue

                terms = [(True, dep, dof, 1.0), (True, indep, dof, -1.0)]
                if abs(c_mem) > tol:
                    terms.append((False, 'RP_E', 1, -c_mem * s))
                if abs(c_bend) > tol:
                    terms.append((False, 'RP_K', 1, -c_bend * s))

                used_dep.add((dep, dof))
                equations.append(terms)

    return equations, used_dep, n_dropped


def build_lesicar_equations(lesicar_constraints, used_dep, nodes, L, ftol_factor=0.01):
    """Assign a dependent (leading) node to each Lesicar Eq.14 integral constraint.

    Two hazards must be avoided, both of which Abaqus rejects as the dependent
    DOF "already eliminated by another equation":

    1. Over-constraint: a (node, dof) already eliminated (by a bending PBC eq or
       an earlier Lesicar eq) must not lead again -> skip nodes in `used_dep`.
    2. Dependency cycle: the leading node is expressed in terms of EVERY other
       node on its face. If that lead is an edge/corner node, its cross-axis
       periodic partner sits on the SAME face and is itself a dependent whose
       PBC equation references the lead -> a cycle. Leading instead with an
       INTERIOR-OF-FACE node (on exactly one cube face) guarantees its only
       periodic partner is on the OPPOSITE face, never in this term list, so no
       equation references the lead back and the dependent graph stays acyclic.

    Mutates `used_dep`. Returns (equations, n_skipped); each equation is a list
    of terms (is_node, name, dof, coeff) as in build_bending_pbc_equations.
    """
    ftol = L * ftol_factor

    def faces_on(lab):
        x, y, z = nodes[lab]
        return ((x < ftol) + (x > L - ftol)
                + (y < ftol) + (y > L - ftol)
                + (z < ftol) + (z > L - ftol))

    equations = []
    n_skipped = 0
    for con in lesicar_constraints:
        nw = con['nodes']
        if len(nw) < 2:
            continue
        cdof = con['dof']

        # Prefer an interior-of-face free node; fall back to any free node;
        # skip if none (constraint is then redundant with those written).
        lead_idx = next((k for k, (lab, _w) in enumerate(nw)
                         if (lab, cdof) not in used_dep and faces_on(lab) == 1), None)
        if lead_idx is None:
            lead_idx = next((k for k, (lab, _w) in enumerate(nw)
                             if (lab, cdof) not in used_dep), None)
        if lead_idx is None:
            n_skipped += 1
            continue
        nw = [nw[lead_idx]] + nw[:lead_idx] + nw[lead_idx + 1:]
        used_dep.add((nw[0][0], cdof))

        terms = [(True, lab, cdof, w) for lab, w in nw]
        if abs(con['rp_E_coeff']) > 1e-20:
            terms.append((False, 'RP_E', 1, -con['rp_E_coeff']))
        if abs(con['rp_K_coeff']) > 1e-20:
            terms.append((False, 'RP_K', 1, -con['rp_K_coeff']))
        equations.append(terms)

    return equations, n_skipped


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

    # Gmsh -> Abaqus C3D10 mid-edge node remap. Gmsh's 10-node tet orders its
    # last two mid-edge nodes as [..., m03, m23, m13] whereas Abaqus C3D10 wants
    # [..., m03, m13, m23] (nodes 9 and 10 swapped). Without this every element
    # reads as a near-planar, zero/negative-volume tet and the input processor
    # aborts. Linear (4-node) elements are untouched.
    _n_remapped = 0
    for label, (etype, elset, conn) in elements.items():
        if len(conn) == 10:
            conn[8], conn[9] = conn[9], conn[8]
            _n_remapped += 1
    if _n_remapped:
        print("    Remapped {} C3D10 elements to Abaqus mid-node order".format(_n_remapped))
    
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
    
    # ---- Determine element type (per phase) ----
    # Order from the mesh actually read: SPAX_MESH_ORDER=2 emits quadratic
    # 10-node tets (Gmsh "Tetrahedron 10"), locking-free in bending; else linear.
    # HYBRID (the 'H' suffix, an internal pressure DOF) is only needed for a
    # (near-)INCOMPRESSIBLE phase -- it prevents volumetric locking but ADDS a
    # variable per element and so enlarges/slows the solve. The brine inclusions
    # (nu~0.49) need it; the compressible MATRIX (nu~0.33) does NOT, and it is
    # the majority of elements, so making the matrix non-hybrid shrinks the
    # system with no accuracy loss. Key the suffix on each phase's nu; force
    # all-hybrid (legacy) with SPAX_FORCE_HYBRID=1.
    _max_conn = max((len(c) for (_t, _e, c) in elements.values()), default=4)
    _base = 'C3D10' if _max_conn >= 10 else 'C3D4'
    _force_hyb = os.environ.get('SPAX_FORCE_HYBRID', '0') == '1'
    _NU_HYB = float(os.environ.get('SPAX_HYBRID_NU', '0.45'))
    def _etype(nu):
        return _base + ('H' if (_force_hyb or nu >= _NU_HYB) else '')
    elem_type_matrix = _etype(nu_matrix)
    elem_type_sphere = _etype(nu_incl)
    print("    Element types: matrix={} (nu={:.3f}), inclusion={} (nu={:.3f})".format(
        elem_type_matrix, nu_matrix, elem_type_sphere, nu_incl))
    
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
    elif mode == 'tors':
        step_name = 'Step-Torsion'
        step_desc = 'Second-Order Torsion (twist about z)'
    else:
        raise ValueError(
            "Unsupported mode: {}. Use one of utx/uty/utz/ss12/ss13/ss23, "
            "bend or tors.".format(mode))

    # ---- Second-order setup ----
    # Bending and torsion share all of the machinery below and differ only in
    # the prescribed macroscopic field, which enters through the coefficient
    # functions keyed on `so_plane`. In torsion RP_K carries the twist rate
    # alpha and the membrane reference point RP_E is unused.
    is_torsion = (mode == 'tors')
    is_bending = (mode == 'bend') or is_torsion
    so_plane = 'torsion' if is_torsion else bending_plane

    if is_bending:
        # Compute coordinate-dependent coefficients for each pair
        # RP_E: membrane strain (free, DOF 1)
        # RP_K: curvature (prescribed, DOF 1 = Kappa)
        rp_base_bend = n_nodes + 200
        rp_E_id = rp_base_bend + 1
        rp_K_id = rp_base_bend + 2
        
        # Torsion has no membrane term, so every RP_E coefficient is zero and
        # the node would appear in no equation. Abaqus then treats it as
        # inactive and rejects the boundary condition on it, so it is simply
        # not created in that mode.
        rp_nodes_bend = {
            'RP_K': (rp_K_id, 2.0*L, 2.0*L, 2.0*L),
        }
        if not is_torsion:
            rp_nodes_bend['RP_E'] = (rp_E_id, 2.5*L, 2.0*L, 2.0*L)
    
    # ---- Write complete .inp ----
    inst_name = 'PART-1-1'
    
    with open(output_inp_path, 'w') as f:
        # ============================================================
        # HEADER
        # ============================================================
        f.write('*Heading\n')
        f.write('** SpaX RVE - {} (standalone)\n'.format(step_desc))
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
        f.write('*Element, type={}, elset=Matrix_Only\n'.format(elem_type_matrix))
        for label in sorted(elements.keys()):
            etype, elset, conn = elements[label]
            if elset == 'Matrix_Only' or (m_start <= label <= m_end):
                conn_str = ', '.join(str(n) for n in conn)
                f.write('{}, {}\n'.format(label, conn_str))
        
        s_start, s_end = sphere_label_range
        has_spheres = s_start > 0 and s_end >= s_start
        if has_spheres:
            f.write('*Element, type={}, elset=Sphere_Only\n'.format(elem_type_sphere))
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
            f.write('** Second-Order {} PBCs\n'.format(
                'Torsion' if is_torsion else 'Bending'))
            # Build the periodic *Equation set. `used_dep` (the (node,dof) pairs
            # eliminated as dependent DOFs) is shared with the Lesicar integral
            # constraints below so the two never claim the same dependent DOF.
            bend_eqs, used_dep, n_dropped = build_bending_pbc_equations(
                pairs, nodes, L, so_plane)
            for terms in bend_eqs:
                f.write('*Equation\n{}\n'.format(len(terms)))
                for is_node, name, dof, coeff in terms:
                    # .12g keeps coordinate-derived RP coefficients essentially
                    # exact while printing 1.0/-1.0 cleanly as 1/-1.
                    if is_node:
                        f.write('{}.{}, {}, {:.12g}\n'.format(inst_name, name, dof, coeff))
                    else:
                        f.write('{}, {}, {:.12g}\n'.format(name, dof, coeff))
            print("    Bending equations: {} ({} redundant dropped)".format(
                len(bend_eqs), n_dropped))
            
            # Lesicar Eq.14 integral constraints (always included). The lead
            # selection (interior-of-face, no over-constraint, no dependency
            # cycle) lives in build_lesicar_equations and is unit-tested.
            f.write('** Lesicar Eq.14 integral constraints\n')
            lesicar_constraints = compute_lesicar_constraints(
                nodes, elements, L, so_plane)
            lesicar_eqs, n_lesicar_skipped = build_lesicar_equations(
                lesicar_constraints, used_dep, nodes, L)
            for terms in lesicar_eqs:
                f.write('*Equation\n{}\n'.format(len(terms)))
                for is_node, name, dof, coeff in terms:
                    if is_node:
                        f.write('{}.{}, {}, {}\n'.format(inst_name, name, dof, coeff))
                    else:
                        f.write('{}, {}, {:.12g}\n'.format(name, dof, coeff))
            print("    Lesicar constraints: {} ({} redundant skipped)".format(
                len(lesicar_eqs), n_lesicar_skipped))
        
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
            # RP_E (membrane strain) and RP_K (curvature) are plain assembly
            # nodes coupled to the model only through DOF 1. Pin their remaining
            # translational DOFs (2,3) so they are not left as free, unconstrained
            # DOFs (which would make the stiffness matrix singular). Only DOFs 1-3
            # exist on a plain node, so constraining 2-3 (not 2-6) avoids spurious
            # "boundary condition on inactive dof" warnings.
            if not is_torsion:
                f.write('RP_E, 2, 3\n')
            f.write('RP_K, 2, 3\n')
        f.write('**\n')
        
        # Amplitude -- LINEAR (amplitude == step time). The post-processor uses
        # eps_macro = frameValue * eng_strain, i.e. it treats step time as the
        # strain axis. That proxy is only correct when the ramp is linear; a
        # back-loaded ramp makes frameValue overstate the strain in the 10-40%
        # fit window and halves the extracted E_eff/G_eff.
        f.write('*Amplitude, name=LoadRamp, time=STEP TIME\n')
        f.write('0., 0., 1., 1.\n')
        
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
        elif mode in ('bend', 'tors'):
            # In torsion RP_K carries the twist rate alpha (radians per unit
            # length); `kappa` is reused as its magnitude, so a torsion deck
            # needs no new column.
            f.write('*Boundary, amplitude=LoadRamp\n')
            f.write('RP_K, 1, 1, {}\n'.format(kappa))
        
        f.write('**\n')
        
        # Output
        f.write('*Output, field, frequency=1\n')
        f.write('*Node Output\nU, RF\n')
        f.write('*Element Output\n')
        # EVOL (whole-element volume) is what the post-processor volume-averages
        # over; it must be requested for every mode. (Was IVOL for first-order,
        # which the post-processor doesn't read -> KeyError 'EVOL'.)
        f.write('S, E, LE, EVOL\n' if not is_bending else 'S, E, EVOL, COORD\n')
        f.write('*Output, history, frequency=1\n')
        if mode in ('bend', 'tors'):
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

def build_slabs(params, L):
    """Cell-spanning brine layers from the deck row, or [] if none are asked for.

    `slab_vof` is the brine fraction carried by the layers; the pocket phase
    keeps whatever `VoF_incl_sphere` still asks for, so a deck can run layers
    only, pockets only, or both. Each layer is a box of thickness t pierced by
    ice bridges of total area fraction b, and the brine it carries is
    t*(1-b)*L^2, which is what fixes t:

        n_slabs * t * (1 - b) / L = slab_vof

    Layers are spread evenly through the cell on the slab normal. Even spacing
    rather than random placement because in series the transverse modulus
    depends only on the total layer thickness and the bridge fraction, not on
    where the layers sit, and even spacing keeps them clear of each other and of
    the faces without a rejection loop.
    """
    import SpaX_GmshPeriodic as _gp

    n_slabs = int(float(params.get('n_slabs', 0) or 0))
    slab_vof = float(params.get('slab_vof', 0.0) or 0.0)
    if n_slabs <= 0 or slab_vof <= 0:
        return []

    b = float(params.get('bridge_fraction', 0.0) or 0.0)
    n_bridges = int(float(params.get('n_bridges', 4) or 4))
    axis = {'x': 0, 'y': 1, 'z': 2}.get(
        str(params.get('slab_axis', 'x')).strip().lower(), 0)

    t = slab_vof * L / (n_slabs * max(1.0 - b, 1e-6))
    pitch = L / n_slabs
    if t >= pitch:
        raise ValueError(
            'slab_vof %.3f needs thickness %.4f per layer but the pitch is only '
            '%.4f: use more layers or a lower fraction' % (slab_vof, t, pitch))

    # Correlation of bridge positions between layers. Zero reproduces the
    # independent placement, under which the drained modulus falls as n^-1.14
    # with cell size and the cell does not homogenise; one stacks the bridges
    # so the load path is straight. Default zero so existing decks are
    # unchanged, and set explicitly where it is being tested.
    corr = float(params.get('bridge_correlation', 0.0) or 0.0)

    per_layer, b_real = _gp.place_bridges_layers(
        L, b, n_bridges, n_slabs, correlation=corr, seed=7919)

    slabs = []
    for k in range(n_slabs):
        slabs.append(dict(origin=(k + 0.5) * pitch - 0.5 * t,
                          thickness=t, axis=axis, bridges=per_layer[k]))
    print("    [Slabs] {} layer(s) normal to {}, t={:.4f} ({:.1f}% of L), "
          "b={:.4f} over {} bridge(s), correlation {:.2f}".format(
              n_slabs, 'xyz'[axis], t, 100.0 * n_slabs * t / L, b_real,
              n_bridges, corr))
    return slabs


def mesh_in_subprocess(sphere_array, L, L_mesh, output_dir, mode, run_id,
                       VoF_void, VoF_incl, Inclusion_Type, gap_balls=None,
                       slabs=None):
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
    gmsh_module = os.path.join(script_dir, 'SpaX_GmshPeriodic.py')
    payload = dict(sphere_array=sphere_array, L=L, L_mesh=L_mesh,
                   output_dir=output_dir, Is_Porous=mode, run_id=run_id,
                   VoF_void_sphere=VoF_void, VoF_incl_sphere=VoF_incl,
                   Inclusion_Type=Inclusion_Type, gap_balls=gap_balls or [],
                   slabs=slabs or [],
                   # SPAX_MESH_ORDER=2 -> quadratic C3D10 (locking-free in
                   # bending). Default 1 (C3D4) keeps existing runs unchanged.
                   mesh_order=int(os.environ.get('SPAX_MESH_ORDER', '1')))
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


def _mode_short(m):
    m = m.lower()
    if 'shear' in m:
        if '12' in m: return 'ss12'
        if '13' in m: return 'ss13'
        if '23' in m: return 'ss23'
        return 'ss13'
    # Uniaxial: pick the axis from an explicit x/y/z token. NB the word
    # "uniaxial" itself contains an 'x', so strip the mode keywords before
    # scanning for the direction (otherwise every uniaxial mode reads as X).
    d = m.replace('uniaxial', ' ').replace('tension', ' ').replace('compression', ' ')
    toks = d.split()
    for t in reversed(toks):
        if t in ('x', 'y', 'z'): return 'ut' + t
    if 'z' in d: return 'utz'
    if 'y' in d: return 'uty'
    if 'x' in d: return 'utx'
    return 'utx'


def _expected_job_paths(params, output_dir):
    """Job-*.inp file(s) a row is expected to produce, mirroring the write
    logic in _generate_one_row. Used by the resume check."""
    run_id = params['run_id']
    full_tensor = params.get('full_tensor', 'No').strip().lower() in ('yes', 'true', '1')
    if full_tensor:
        modes = ['utx', 'uty', 'utz', 'ss12', 'ss13', 'ss23']
    else:
        modes = [_mode_short(params.get('Mode', 'Uniaxial Tension X'))]
        if params.get('Mode2', ''):
            modes.append(_mode_short(params['Mode2']))
    if float(params.get('Kappa', 0.0) or 0.0) > 0:
        modes.append('ben')
    seen = []  # de-dup, preserve order (Mode/Mode2 may map to the same short name)
    for ms in modes:
        if ms not in seen:
            seen.append(ms)
    return [os.path.join(output_dir, 'Job-{}-{}.inp'.format(run_id, ms)) for ms in seen]


def _resume_enabled():
    return os.environ.get('SPAX_RESUME', '').strip().lower() in ('1', 'yes', 'true', 'on')


def _generate_one_row(task):
    """Generate one RVE's .inp file(s). Independent per row, so this runs
    in a worker process and packing+meshing parallelise across CPUs.
    Returns (run_id, status) where status is 'ok', 'exists' or 'skipped'."""
    idx, n_rows, params, output_dir = task
    # Independent RNG per worker: forked children inherit the parent's
    # numpy state, which would make parallel rows pack identically. Reseed
    # from SPAX_SEED (reproducible) or fresh OS entropy (default).
    _seed = os.environ.get('SPAX_SEED')
    if _seed:
        np.random.seed((int(_seed) + idx * 2654435761) % (2 ** 32))
    else:
        np.random.seed(None)
    run_id = params['run_id']
    print("\n" + "=" * 70)
    print("[{}/{}] {}".format(idx + 1, n_rows, run_id))
    print("=" * 70)

    # Resume (opt-in via SPAX_RESUME): if this row's expected Job-*.inp already
    # exist and are non-empty, skip the expensive pack+mesh. Lets a timed-out or
    # partial generation job pick up where it left off; a fresh run (flag unset)
    # still regenerates everything.
    if _resume_enabled():
        expected = _expected_job_paths(params, output_dir)
        if expected and all(os.path.exists(p) and os.path.getsize(p) > 0
                            for p in expected):
            print("  RESUME: {} Job file(s) already present -> skip".format(
                len(expected)))
            return (run_id, 'exists')

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
    # Enabled by default (1.0 -> one fine element): the ellipsoid-aware
    # densification pass packs inclusion *surfaces* `min_dist` apart along the
    # inter-centre line, but two ellipsoids can sit closer off that line, so the
    # matrix slivers between grown inclusions collapse below element size unless
    # the gap floor is a full element. A real gmsh run on the sea-ice RVE
    # confirmed 0.5 (half element) still crashed the mesher on the first two
    # attempts at 100% target VoF, while 1.0 meshes first-try and still hits
    # ~100% (densification recovers VoF, so the larger gap costs no final volume
    # fraction). Set 0 to restore the legacy min_distance-only behaviour.
    gap_mult = float(os.environ.get('SPAX_GAP_MULT', '1.0'))
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

    # Quadratic (C3D10) needs a THICKER minimum boundary cap than linear: a thin
    # face triangle that meshes fine and is positive-volume as a linear C3D4
    # becomes a negative-Jacobian quadratic element (Abaqus then aborts). The
    # escalation below never rescues this because Gmsh succeeds at attempt 0 (the
    # sliver only fails downstream in Abaqus, so the retry loop is never entered).
    # So for order-2 we reject grazing caps from attempt 0 with a constant floor
    # of SPAX_SLIVER_MULT_Q * lc_fine (default 1.0, vs the linear 0.5). Costs a
    # little VoF; eliminates the face slivers at the DEFAULT element size.
    _mesh_order = int(os.environ.get('SPAX_MESH_ORDER', '1'))
    sliver_floor_q = (float(os.environ.get('SPAX_SLIVER_MULT_Q', '1.0'))
                      * lc_fine_mult * L_mesh) if _mesh_order == 2 else 0.0

    def sliver_for_attempt(n):
        # n = 0-based attempt index. Quadratic carries a constant cap-rejection
        # floor from attempt 0; on top of that the linear escalation ramps from
        # sliver_start to sliver_gap_max on the final attempt.
        ramp = 0.0
        if sliver_gap_max > 0.0 and n >= sliver_start:
            span = max(1, (max_retries - 1) - sliver_start)
            ramp = min(1.0, (n - sliver_start + 1) / float(span)) * sliver_gap_max
        return max(ramp, sliver_floor_q)

    # Channel options resolved once; the packing itself is channel-FIRST.
    gen_channels = params.get('generate_channels', 'No').strip().lower() in ('yes', 'true', '1')
    channel_vof_target = float(params.get('channel_vof_target', 0) or 0)
    do_channels = gen_channels and channel_vof_target > 0.001
    r_ch_avg = float(params.get('r_channel_avg', r_avg * 0.5) or r_avg * 0.5)
    r_ch_std = float(params.get('r_channel_std', r_ch_avg * 0.2) or r_ch_avg * 0.2)

    # True-gap floor for the off-axis sliver repair: a mesh-resolvable fraction
    # of the fine element size (lc_fine = lc_fine_mult*L_mesh). The repair lifts
    # only the genuine sub-element outliers to this floor, so it costs almost no
    # VoF. SPAX_OFFAXIS_FRAC tunes it (0 -> fall back to the centre-line sep).
    offaxis_frac = float(os.environ.get('SPAX_OFFAXIS_FRAC', '0.6'))
    offaxis_floor = offaxis_frac * lc_fine_mult * L_mesh if offaxis_frac > 0 else None
    # Channel<->inclusion slivers run the FULL RVE height, so the mesher tiles a
    # tall matrix sheet that collapses to overlapping facets unless it is at
    # least one element thick. Hold these gaps to a full lc_fine (vs the
    # 0.6*lc_fine sphere-sphere floor). SPAX_OFFAXIS_CHANNEL_FRAC tunes it.
    offaxis_channel_frac = float(os.environ.get('SPAX_OFFAXIS_CHANNEL_FRAC', '1.0'))
    offaxis_channel_floor = (offaxis_channel_frac * lc_fine_mult * L_mesh
                             if offaxis_channel_frac > 0 else None)

    # Dense (esp. channel/hybrid) packs let two near-tangent ellipsoids tilt
    # into a sub-element OFF-AXIS sliver (worst sphere-sphere true gap
    # ~0.88*lc_fine at the default sep -> Gmsh "overlapping facets" crash).
    # Rather than WIDEN the gap to a full element (which caps VoF -- widening
    # the densify sep over-widens every gap to fix the worst outlier, costing
    # ~5 VoF points), we keep the inclusions TIGHT and instead REFINE THE MESH
    # IN THE GAP: each sub-element sliver gets a local size ball (below) so the
    # mesher resolves it. This recovers the VoF the gap-widening sacrificed.
    # SPAX_CHANNEL_SEP still widens the channel-side sep if ever needed (default
    # 1.0 = no extra widening; channel<->inclusion gaps have no off-axis penalty
    # and are already mesh-safe at sep>=lc_fine).
    chan_sep = float(os.environ.get('SPAX_CHANNEL_SEP', '1.0'))
    if do_channels and chan_sep > 1.0:
        min_dist = max(min_dist, chan_sep * lc_fine_mult * L_mesh)

    # Mesh-in-gap refinement: collect local size balls at narrow sphere-sphere
    # slivers so the mesher resolves tight gaps instead of the packer widening
    # them. On by default; SPAX_GAP_REFINE=0 disables (legacy widen-only path).
    gap_refine = os.environ.get('SPAX_GAP_REFINE', '1') != '0'
    gap_resolve = float(os.environ.get('SPAX_GAP_RESOLVE', '0.5'))
    lc_fine = lc_fine_mult * L_mesh

    # Anisotropic densification: grow inclusions preferentially along the
    # Growth_Direction axis (physical brine channels elongate vertically) rather
    # than scaling all axes equally. Opt-in (default 0 = isotropic, unchanged);
    # w>0 lowers final sphericity by design. SPAX_ZGROW_BIAS sets the strength.
    z_bias = float(os.environ.get('SPAX_ZGROW_BIAS', '0.0'))

    def _pack_row(sliver_gap, md_scale=1.0):
        """Pack one RVE's inclusions + channels (channel-FIRST). Channels are
        placed on an empty octree so they claim their XY space before the
        spheres densify and saturate the plane; the channel primaries are then
        handed to `generate_sphere_packing`, which packs and grows the spheres
        to clear them (ellipsoid-aware, both directions). Returns the combined
        sphere array (inclusions + channels-as-tall-ellipsoids). Re-run per
        attempt so a re-pack keeps its channels (they used to be added once,
        outside the retry loop, and were silently lost on every retry).
        `md_scale` widens the min-distance on retries (the old retry used 1.2)."""
        # --- Packing decoupling for controlled mesh-convergence studies ---
        # SPAX_LOAD_PACKING=<dir>: reuse a frozen sphere array (<run_id>.npy),
        # skipping the packer entirely, so the SAME geometry can be meshed at
        # several L_mesh (the packer's r_floor/sliver floors scale with L_mesh
        # and would otherwise change the packing). gap_balls are recomputed for
        # the current lc_fine. Tall channel ellipsoids (rz >= L/2) are excluded
        # from the gap-ball pass, as in the normal path.
        _load_dir = os.environ.get('SPAX_LOAD_PACKING', '')
        if _load_dir:
            _pf = os.path.join(_load_dir, run_id + '.npy')
            sphere_array = np.load(_pf)
            print("  [Packing] loaded {} entries from {} (packer skipped)".format(
                len(sphere_array), _pf))
            gap_balls = []
            if gap_refine:
                incl = [tuple(row[:6]) for row in sphere_array
                        if row[5] < 0.49 * L]
                gap_balls = _collect_gap_balls(incl, L, lc_fine, channels=None,
                                               resolve=gap_resolve)
            return sphere_array, None, gap_balls
        md = min_dist * md_scale
        channel_array = np.empty((0, 3))
        channel_prims = None
        if do_channels:
            print("  Step 1: Channel packing (channel-first)...")
            empty_oct = Octree(L, capacity=8, max_depth=12)
            channel_array = generate_channels(
                L=L, channel_vof_target=channel_vof_target,
                r_channel_avg=r_ch_avg, r_channel_std=r_ch_std,
                min_distance=md, max_iterations=max_iter,
                octree=empty_oct, L_mesh=L_mesh)
            # Primaries (centre in [0,L)) feed the sphere-side clearance; the
            # helper's min-image reconstructs the periodic neighbours.
            channel_prims = [(c[0], c[1], c[2]) for c in channel_array
                             if 0.0 <= c[0] < L and 0.0 <= c[1] < L]

        print("  Step 1b: Sphere packing...")
        sphere_array, oct_ = generate_sphere_packing(
            L=L, r_avg=r_avg, r_std=r_std,
            VoF_target=VoF, min_distance=md,
            max_iterations=max_iter,
            sphericity_avg=sph_avg, sphericity_std=sph_std,
            growth_direction=growth_dir,
            growth_concentration=growth_conc,
            min_radius=r_floor, sliver_gap=sliver_gap,
            channels=channel_prims, offaxis_floor=offaxis_floor,
            offaxis_channel_floor=offaxis_channel_floor, z_bias=z_bias)

        # Mesh-in-gap balls from the INCLUSION pack (before channels are
        # appended): channels are vertical cylinders meshed head-on and never
        # need gap refinement, so only sphere-sphere slivers are collected.
        gap_balls = []
        if gap_refine:
            incl = [tuple(row[:6]) for row in sphere_array]
            gap_balls = _collect_gap_balls(incl, L, lc_fine, channels=channel_prims,
                                           resolve=gap_resolve)
            if gap_balls:
                print("    [Mesh-in-gap] {} refinement ball(s) at narrow "
                      "sphere-sphere slivers".format(len(gap_balls)))

        # Append channels as tall ellipsoids (rz = L/2 spans full Z); the
        # GmshPeriodic path treats these like any other inclusion.
        for i in range(channel_array.shape[0]):
            ch_x, ch_y, ch_r = channel_array[i]
            ch_entry = np.array([[ch_x, ch_y, L/2, ch_r, ch_r, L/2, 0, 0, 0, 0.01]])
            sphere_array = np.vstack((sphere_array, ch_entry))
        # Freeze this packing (full array incl. channels) for reuse across
        # meshes: SPAX_SAVE_PACKING=<dir> -> <run_id>.npy.
        _save_dir = os.environ.get('SPAX_SAVE_PACKING', '')
        if _save_dir:
            try:
                os.makedirs(_save_dir)
            except OSError:
                pass
            _pf = os.path.join(_save_dir, run_id + '.npy')
            np.save(_pf, sphere_array)
            print("  [Packing] saved {} entries -> {}".format(len(sphere_array), _pf))
        return sphere_array, oct_, gap_balls

    sliver_gap = sliver_for_attempt(0)
    Sphere_array, octree, gap_balls = _pack_row(sliver_gap)

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
                Inclusion_Type=Inclusion_Type,
                gap_balls=gap_balls,
                slabs=build_slabs(params, L))

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
                Sphere_array, octree, gap_balls = _pack_row(sliver_gap, md_scale=1.2)
                result = None
            else:
                print("    SKIPPING {} after {} failures".format(run_id, max_retries))
                result = None

    if result is None:
        return (run_id, 'skipped')

    gmsh_inp = result['mesh_file']
    pairs_csv = result['match_file']
    m_range = result['matrix_label_range']
    s_range = result['sphere_label_range']

    # Validate mesh is not empty
    n_total = result.get('n_elements_matrix', 0) + result.get('n_elements_sphere', 0)
    if n_total == 0:
        print("  ERROR: Empty mesh (0 elements). Skipping {}".format(run_id))
        return (run_id, 'skipped')

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

    if full_tensor:
        all_modes = ['utx', 'uty', 'utz', 'ss12', 'ss13', 'ss23']
        print("    Full tensor mode: generating {} load cases".format(len(all_modes)))
        for ms in all_modes:
            path = os.path.join(output_dir, 'Job-{}-{}.inp'.format(run_id, ms))
            write_complete_inp(gmsh_inp, pairs_csv, path,
                mode=ms, disp=Disp, **common)
    else:
        ms = _mode_short(Mode)
        path1 = os.path.join(output_dir, 'Job-{}-{}.inp'.format(run_id, ms))
        write_complete_inp(gmsh_inp, pairs_csv, path1,
            mode=ms, disp=Disp, **common)

        if Mode2:
            Disp2 = float(params.get('Disp2', Disp) or Disp)
            ms2 = _mode_short(Mode2)
            path2 = os.path.join(output_dir, 'Job-{}-{}.inp'.format(run_id, ms2))
            write_complete_inp(gmsh_inp, pairs_csv, path2,
                mode=ms2, disp=Disp2, **common)

    # Second-order probe (always if Kappa > 0). Bending_Plane selects which:
    # 'xz'/'yz'/'xy' give the bending modes, 'torsion' the twist probe. The two
    # are mutually exclusive per deck, so no existing deck changes behaviour.
    if Kappa_val > 0:
        if str(bp).strip().lower() == 'torsion':
            path_tor = os.path.join(output_dir, 'Job-{}-tor.inp'.format(run_id))
            write_complete_inp(gmsh_inp, pairs_csv, path_tor,
                mode='tors', disp=0, bending_plane='torsion',
                kappa=Kappa_val, **common)
        else:
            path_ben = os.path.join(output_dir, 'Job-{}-ben.inp'.format(run_id))
            write_complete_inp(gmsh_inp, pairs_csv, path_ben,
                mode='bend', disp=0, bending_plane=bp, kappa=Kappa_val,
                **common)

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
    return (run_id, 'ok')


def _gen_workers(n_rows):
    """Worker count: SPAX_GEN_WORKERS, else SLURM_CPUS_PER_TASK, else all cores;
    never more than the number of rows."""
    env = os.environ.get('SPAX_GEN_WORKERS')
    if env:
        n = int(env)
    else:
        n = int(os.environ.get('SLURM_CPUS_PER_TASK') or 0) or (os.cpu_count() or 1)
    return max(1, min(n, n_rows))


def process_csv(csv_path, output_dir):
    """
    Read a parametric study CSV and generate complete .inp files
    for each RVE (UTX + SS13 load cases).
    """
    
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    print("=" * 70)
    print("SPAX STANDALONE .INP GENERATOR")
    print("=" * 70)
    print("  CSV: {}".format(csv_path))
    print("  Output: {}".format(output_dir))
    print("  RVEs: {}".format(len(rows)))
    
    os.makedirs(output_dir, exist_ok=True)
    
    n_rows = len(rows)
    workers = _gen_workers(n_rows)
    tasks = [(i, n_rows, p, output_dir) for i, p in enumerate(rows)]

    if workers <= 1:
        results = [_generate_one_row(t) for t in tasks]
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        print("  Parallel generation: {} workers over {} rows".format(workers, n_rows))
        results = []
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(_generate_one_row, t): t for t in tasks}
            for fut in as_completed(futs):
                try:
                    results.append(fut.result())
                except Exception as e:
                    rid = futs[fut][2].get('run_id', '?')
                    print("  [{}] worker crashed: {}".format(rid, e))
                    results.append((rid, 'skipped'))

    n_ok = sum(1 for _, s in results if s == 'ok')
    n_exists = sum(1 for _, s in results if s == 'exists')
    n_skip = n_rows - n_ok - n_exists
    if n_exists:
        print("\n  Generated {} new, {} already present (resume), {} skipped "
              "(of {} rows)".format(n_ok, n_exists, n_skip, n_rows))
    else:
        print("\n  Generated {}/{} RVEs ({} skipped)".format(
            n_ok, n_rows, n_skip))

    
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
        print("Usage: python SpaX_Standalone.py <csv_file> <output_dir>")
        print("  Generates complete .inp files for Abaqus solver (no CAE needed)")
        print("  Requires: pip install numpy gmsh")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    output_dir = sys.argv[2]
    
    process_csv(csv_file, output_dir)
