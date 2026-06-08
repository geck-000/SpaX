"""
Verify build_bending_pbc_equations on a synthetic periodic cube whose mesh has
genuine shared edge/corner nodes (the case the old per-pair skip broke).

Checks, for every bending_plane:
  1. No (node,dof) is eliminated (led) more than once  -> no Abaqus over-constraint.
  2. Every written *Equation is exactly satisfied by the prescribed macro
     bending field  -> coefficients and sign-flips are correct.
  3. Per DOF, the written equations span the SAME connected-component partition
     as the full intended pair graph  -> only redundant cycle-closing equations
     were dropped, never a real (bridge) periodicity constraint.
  4. n_dropped > 0  -> the redundant-cycle drop path is actually exercised.

Run:  python3 tests/test_bending_pbc.py     (needs only numpy)
"""
from __future__ import print_function
import os
import sys

# Import the helpers from the repo root regardless of where this is run from.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Spatium_Standalone import build_bending_pbc_equations, build_lesicar_equations

L = 2.0
M = 4              # divisions per edge -> (M+1) points per edge
EPS = L * 1e-9


def build_cube_boundary():
    """All grid nodes with >=1 coord on a face. Returns (nodes, pairs)."""
    step = L / M
    coords = [round(i * step, 9) for i in range(M + 1)]
    nodes = {}            # label -> (x,y,z)
    lookup = {}           # (x,y,z) rounded -> label
    label = 1

    def key(x, y, z):
        return (round(x, 6), round(y, 6), round(z, 6))

    for x in coords:
        for y in coords:
            for z in coords:
                on_face = (x in (coords[0], coords[-1]) or
                           y in (coords[0], coords[-1]) or
                           z in (coords[0], coords[-1]))
                if not on_face:
                    continue
                nodes[label] = (x, y, z)
                lookup[key(x, y, z)] = label
                label += 1

    pairs = {'X': [], 'Y': [], 'Z': []}
    lo, hi = coords[0], coords[-1]
    for lab, (x, y, z) in nodes.items():
        if abs(x - lo) < EPS:
            pairs['X'].append((lab, lookup[key(hi, y, z)]))   # (neg, pos)
        if abs(y - lo) < EPS:
            pairs['Y'].append((lab, lookup[key(x, hi, z)]))
        if abs(z - lo) < EPS:
            pairs['Z'].append((lab, lookup[key(x, y, hi)]))
    return nodes, pairs


def macro_disp(x, y, z, a, b, plane):
    """Prescribed macro field; centred at L/2. RP_E=a, RP_K=b."""
    x1, x2, x3 = x - L / 2, y - L / 2, z - L / 2
    if plane == 'xz':
        return [a * x1 - b * (x1 * x3), 0.0, b * 0.5 * x1 ** 2]
    if plane == 'yz':
        return [0.0, a * x2 - b * (x2 * x3), b * 0.5 * x2 ** 2]
    # xy
    return [a * x1 - b * (x1 * x2), b * 0.5 * x1 ** 2, 0.0]


class UF(object):
    def __init__(self): self.p = {}
    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]; x = self.p[x]
        return x
    def union(self, a, b):
        self.p[self.find(a)] = self.find(b)


def run_plane(plane):
    nodes, pairs = build_cube_boundary()
    eqs, used_dep, n_dropped = build_bending_pbc_equations(pairs, nodes, L, plane)

    # ---- 1. no (node,dof) led twice ----
    dep_keys = [(eq[0][1], eq[0][2]) for eq in eqs]      # first term is the dep
    assert len(dep_keys) == len(set(dep_keys)), \
        "{}: a (node,dof) leads more than one equation".format(plane)
    assert len(eqs) == len(used_dep), \
        "{}: eq count {} != used_dep {}".format(plane, len(eqs), len(used_dep))

    # ---- 2. macro field satisfies every equation ----
    a, b = 0.013, 0.027
    disp = {lab: macro_disp(x, y, z, a, b, plane) for lab, (x, y, z) in nodes.items()}
    max_res = 0.0
    for eq in eqs:
        total = 0.0
        for is_node, name, dof, coeff in eq:
            if is_node:
                total += coeff * disp[name][dof - 1]
            else:                       # RP term, always dof 1
                total += coeff * (a if name == 'RP_E' else b)
        max_res = max(max_res, abs(total))
    assert max_res < 1e-9, "{}: macro field residual {:.3e}".format(plane, max_res)

    # ---- 3. coverage: written graph spans the full pair graph per dof ----
    all_nodes = set(nodes.keys())

    def partition(edges):
        uf = UF()
        for n in all_nodes:
            uf.find(n)
        for n1, n2 in edges:
            uf.union(n1, n2)
        comp = {}
        for n in all_nodes:
            comp.setdefault(uf.find(n), set()).add(n)
        return frozenset(frozenset(v) for v in comp.values())

    full_edges = [(neg, pos) for ax in ('X', 'Y', 'Z') for (neg, pos) in pairs[ax]]
    full_part = partition(full_edges)
    for dof in (1, 2, 3):
        written_edges = []
        for eq in eqs:
            nd = [name for (is_node, name, df, c) in eq if is_node and df == dof]
            if len(nd) == 2:
                written_edges.append((nd[0], nd[1]))
        assert partition(written_edges) == full_part, \
            "{} dof{}: written equations do not span the full pair graph " \
            "(a bridge constraint was dropped)".format(plane, dof)

    # ---- 4. redundant-cycle drops actually happened ----
    assert n_dropped > 0, "{}: expected some redundant drops, got 0".format(plane)

    print("  [{}] eqs={:<4d} dropped={:<3d} dep-unique=OK macro-res={:.1e} spans-full=OK".format(
        plane, len(eqs), n_dropped, max_res))


def synthetic_lesicar(nodes, L, tol):
    """Build Lesicar-style integral constraints (one per negative face per dof):
    every node on the face is a term. Mirrors compute_lesicar_constraints' shape
    so build_lesicar_equations can pick a leading node. RP coeffs are irrelevant
    to the lead choice, so set to 0."""
    faces = {'x': [], 'y': [], 'z': []}
    for lab, (x, y, z) in nodes.items():
        if x < tol: faces['x'].append(lab)
        if y < tol: faces['y'].append(lab)
        if z < tol: faces['z'].append(lab)
    cons = []
    for ax, labs in faces.items():
        for dof in (1, 2, 3):
            cons.append({'face': (ax, 0), 'dof': dof,
                         'nodes': [(l, 1.0) for l in labs],
                         'rp_E_coeff': 0.0, 'rp_K_coeff': 0.0})
    return cons


def has_dependency_cycle(equations):
    """equations: list of term-lists (is_node, name, dof, coeff). The first term
    is the eliminated (dependent) DOF. Returns a cycle (list) among dependent
    DOFs if one exists, else None. Abaqus cannot triangularize a cyclic set."""
    dep = {(t[0][1], t[0][2]) for t in equations}     # (name, dof) that are deps
    adj = {}
    for t in equations:
        d = (t[0][1], t[0][2])
        adj.setdefault(d, [])
        for is_node, name, dofv, _c in t[1:]:
            m = (name, dofv)
            if m in dep:
                adj[d].append(m)
    WHITE, GRAY, BLACK = 0, 1, 2
    color = {d: WHITE for d in adj}
    parent = {}
    for s in list(adj):
        if color[s] != WHITE:
            continue
        stack = [(s, iter(adj[s]))]
        color[s] = GRAY
        while stack:
            node, it = stack[-1]
            adv = False
            for m in it:
                cm = color.get(m, BLACK)
                if cm == WHITE:
                    parent[m] = node; color[m] = GRAY
                    stack.append((m, iter(adj[m]))); adv = True; break
                elif cm == GRAY:
                    cyc = [node]; x = node
                    while x != m and x in parent:
                        x = parent[x]; cyc.append(x)
                    return list(reversed(cyc))
            if not adv:
                color[node] = BLACK; stack.pop()
    return None


def run_plane_full(plane):
    """PBC + Lesicar together: the combined dependent set must be over-constraint
    free (no duplicate leads) AND acyclic (the bug where an edge node leads a
    Lesicar constraint and cycles with its cross-axis PBC partner)."""
    nodes, pairs = build_cube_boundary()
    pbc_eqs, used_dep, _ = build_bending_pbc_equations(pairs, nodes, L, plane)
    cons = synthetic_lesicar(nodes, L, EPS * 10)
    les_eqs, _ = build_lesicar_equations(cons, used_dep, nodes, L)

    # Convert PBC terms (node labels) to the same (is_node, name, dof, coeff)
    # shape the cycle checker expects; node names are plain labels here.
    all_eqs = pbc_eqs + les_eqs

    leads = [(t[0][1], t[0][2]) for t in all_eqs]
    assert len(leads) == len(set(leads)), \
        "{}: duplicate leading (node,dof) across PBC+Lesicar".format(plane)

    cyc = has_dependency_cycle(all_eqs)
    assert cyc is None, "{}: dependency cycle among eliminated DOFs: {}".format(plane, cyc)

    # Every Lesicar lead must be an interior-of-face node (on exactly one face).
    ftol = L * 0.01
    def faces_on(lab):
        x, y, z = nodes[lab]
        return ((x < ftol) + (x > L - ftol) + (y < ftol) + (y > L - ftol)
                + (z < ftol) + (z > L - ftol))
    for t in les_eqs:
        lab = t[0][1]
        assert faces_on(lab) == 1, \
            "{}: Lesicar led by non-interior node {} (faces_on={})".format(
                plane, lab, faces_on(lab))

    print("  [{}] pbc={:<4d} lesicar={:<2d} dep-unique=OK acyclic=OK interior-leads=OK".format(
        plane, len(pbc_eqs), len(les_eqs)))


def main():
    print("Testing build_bending_pbc_equations on a {}^3 periodic cube boundary".format(M + 1))
    for plane in ('xz', 'yz', 'xy'):
        run_plane(plane)
    print("\nTesting PBC + Lesicar combined (over-constraint + cycle freedom)")
    for plane in ('xz', 'yz', 'xy'):
        run_plane_full(plane)
    print("ALL BENDING PBC TESTS PASSED")


if __name__ == '__main__':
    main()
