"""Nodal-averaged B-bar for one phase of a CalculiX deck.

    python3 nodalbbar.py <in-ccx.inp> <out.inp> [--elset Sphere_Only]

Splits the phase's response in two so the bulk modulus never appears in an
element-local operator and cannot lock:

    K = K_dev   per tet, element type U5 (deviatoric only)
      + K_vol   per node, built here out of ccx primitives

The nodal volumetric strain is
    V_a     = sum over the phase's elements at node a of V_e/4
    theta_a = (1/V_a) sum over those elements of (V_e/4) div(u)|_e
            = sum over the 1-ring of  b_n^a . u_n
and its energy is (1/2) K V_a theta_a^2.

Rather than build a patch element, theta is carried as one extra dof on a
dummy node, tied by an *EQUATION and given energy by a grounded SPRING1 --
all existing ccx features, so no change to the matrix-structure code.

The scaling is folded into the equation so every spring is identical:
with  t_a = sqrt(K V_a) theta_a  the energy is exactly (1/2) t_a^2, i.e. unit
stiffness for all of them, and one *SPRING card covers the lot. Emitting one
card per node would otherwise mean tens of thousands of element sets.

    *EQUATION :  t_a  -  sum_n sqrt(K V_a) b_n^a . u_n  =  0

t_a is written first because ccx eliminates the leading term as the dependent
dof, and t_a is the one that must go.
"""
import sys
import numpy as np


def main():
    argv = sys.argv[1:]
    elset = 'Sphere_Only'
    if '--elset' in argv:
        i = argv.index('--elset')
        elset = argv[i + 1]
        del argv[i:i + 2]
    src, dst = argv[0], argv[1]
    lines = open(src).read().split('\n')

    def iskw(l):
        return l.startswith('*') and not l.startswith('**')

    # --- parse nodes, the phase's elements, and the material -------------
    co, els, maxnode = {}, [], 0
    mode = None
    for ln in lines:
        if iskw(ln):
            u = ln.upper().replace(' ', '')
            if u.startswith('*NODE') and 'OUTPUT' not in u and 'PRINT' not in u \
               and 'FILE' not in u:
                mode = 'n'
            elif u.startswith('*ELEMENT') and ('ELSET=' + elset).upper() in u:
                mode = 'e'
            else:
                mode = None
            continue
        f = [x.strip() for x in ln.split(',') if x.strip()]
        if mode == 'n' and len(f) >= 4:
            n = int(f[0])
            co[n] = (float(f[1]), float(f[2]), float(f[3]))
            maxnode = max(maxnode, n)
        elif mode == 'e' and len(f) >= 5:
            els.append([int(x) for x in f[1:5]])
    if not els:
        raise SystemExit('nodalbbar: no elements in ELSET=%s' % elset)

    E = nu = None
    for i, ln in enumerate(lines):
        if ln.strip().lower().startswith('*material, name=mat_inclusion'):
            for j in range(i + 1, min(i + 5, len(lines))):
                f = [x.strip() for x in lines[j].split(',') if x.strip()]
                if len(f) == 2 and not lines[j].startswith('*'):
                    E, nu = float(f[0]), float(f[1])
                    break
            break
    if E is None:
        # patch-test decks name their material something else; fall back to the
        # first *ELASTIC card in the deck.
        for i, ln in enumerate(lines):
            if ln.strip().upper().startswith('*ELASTIC'):
                f = [x.strip() for x in lines[i + 1].split(',') if x.strip()]
                if len(f) >= 2:
                    E, nu = float(f[0]), float(f[1])
                    break
    if E is None:
        raise SystemExit('nodalbbar: no elastic card found')
    K = E / (3.0 * (1.0 - 2.0 * nu))

    # --- nodal volumes and the averaged divergence operator --------------
    Va = {}
    b = {}                       # b[a][n] = 3-vector
    for e in els:
        p = np.array([co[n] for n in e])
        J = np.array([p[1] - p[0], p[2] - p[0], p[3] - p[0]]).T
        det = np.linalg.det(J)
        vol = abs(det) / 6.0
        if vol <= 0.0:
            continue
        Ji = np.linalg.inv(J)
        gr = np.zeros((4, 3))
        gr[1:, :] = Ji                      # rows of J^-1 are grad L1..L3
        gr[0, :] = -gr[1:, :].sum(axis=0)   # grads sum to zero
        w = vol / 4.0
        for a in e:
            Va[a] = Va.get(a, 0.0) + w
            ba = b.setdefault(a, {})
            for k, n in enumerate(e):
                ba[n] = ba.get(n, np.zeros(3)) + w * gr[k]
    for a in b:
        for n in b[a]:
            b[a][n] /= Va[a]

    # --- emit --------------------------------------------------------------
    nodes = sorted(Va)
    tnode = {a: maxnode + 1 + i for i, a in enumerate(nodes)}
    out, done_step = [], False
    for ln in lines:
        if iskw(ln):
            u = ln.upper().replace(' ', '')
            if u.startswith('*ELEMENT') and ('ELSET=' + elset).upper() in u:
                out.append('*USER ELEMENT,TYPE=U5,NODES=4,'
                           'INTEGRATIONPOINTS=1,MAXDOF=3')
                out.append(ln.replace('C3D4', 'U5').replace('c3d4', 'U5'))
                continue
            if u.startswith('*STEP') and not done_step:
                done_step = True
                out.append('** SPAX nodal-averaged B-bar: %d node patches'
                           % len(nodes))
                out.append('*NODE')
                for a in nodes:
                    x, y, z = co[a]
                    out.append('%d, %.12e, %.12e, %.12e'
                               % (tnode[a], x, y, z))
                out.append('*ELEMENT,TYPE=SPRING1,ELSET=BBARSPR')
                for i, a in enumerate(nodes):
                    out.append('%d, %d' % (10 ** 8 + i, tnode[a]))
                out.append('*SPRING,ELSET=BBARSPR')
                out.append('1')
                out.append('1.0')
                # Each theta node carries three dofs but only dof 1 is used by
                # the spring and its *EQUATION. Left free, dofs 2 and 3 are
                # unconstrained and the matrix is singular -- SPOOLES simply
                # stops after "Factoring the system of equations".
                out.append('*BOUNDARY')
                for a in nodes:
                    out.append('%d, 2, 3, 0.' % tnode[a])
                # Node set for reading theta back out. The printed dof 1 is
                # t_a = sqrt(K V_a) theta_a, so divide by that to recover the
                # nodal volumetric strain -- the field whose smoothness is the
                # displacement-method analogue of the pressure check.
                out.append('*NSET,NSET=BBARTHETA')
                tl = [tnode[a] for a in nodes]
                for i in range(0, len(tl), 8):
                    out.append(','.join(str(x) for x in tl[i:i + 8]))
                for a in nodes:
                    sc = (K * Va[a]) ** 0.5
                    terms = [(tnode[a], 1, 1.0)]
                    for n, v in b[a].items():
                        for d in range(3):
                            if abs(v[d]) > 0.0:
                                terms.append((n, d + 1, -sc * float(v[d])))
                    out.append('*EQUATION')
                    out.append(str(len(terms)))
                    for nd, dof, c in terms:
                        out.append('%d, %d, %.12e' % (nd, dof, c))
        if iskw(ln) and ln.strip().upper().startswith('*END STEP'):
            pass
        out.append(ln)

    # theta output, just before the step closes
    for i in range(len(out) - 1, -1, -1):
        if out[i].strip().upper().startswith('*END STEP'):
            out[i:i] = ['*NODE PRINT,NSET=BBARTHETA', 'U']
            break

    open(dst, 'w').write('\n'.join(out))
    with open(dst + '.theta', 'w') as fh:
        for a in nodes:
            fh.write('%d %d %.12e\n' % (a, tnode[a], (K * Va[a]) ** 0.5))
    ring = [len(b[a]) for a in nodes]
    print('nodalbbar: %s -> %s' % (src, dst))
    print('  %d elements retyped to U5 in ELSET=%s' % (len(els), elset))
    print('  %d node patches, 1-ring %d..%d nodes (mean %.1f)'
          % (len(nodes), min(ring), max(ring), sum(ring) / len(ring)))
    print('  K = %.4e from the Mat_Inclusion card' % K)
    print('  one *SPRING of unit stiffness; V_a folded into the *EQUATIONs')


if __name__ == '__main__':
    main()
