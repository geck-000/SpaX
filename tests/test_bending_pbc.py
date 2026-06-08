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

from Spatium_Standalone import build_bending_pbc_equations

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


def main():
    print("Testing build_bending_pbc_equations on a {}^3 periodic cube boundary".format(M + 1))
    for plane in ('xz', 'yz', 'xy'):
        run_plane(plane)
    print("ALL BENDING PBC TESTS PASSED")


if __name__ == '__main__':
    main()
