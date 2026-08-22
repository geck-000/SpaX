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

    # material actually assigned to the target elset, from its section card
    matname = None
    for ln in lines:
        u = ln.upper().replace(' ', '')
        if u.startswith('*SOLIDSECTION') and ('ELSET=' + elset).upper() in u:
            for fld in ln.split(','):
                if fld.strip().upper().startswith('MATERIAL='):
                    matname = fld.split('=', 1)[1].strip()
            break
    if matname is None:
        raise SystemExit('nodalbbar: no *SOLID SECTION for ELSET=%s' % elset)
    E = nu = None
    for i, ln in enumerate(lines):
        if ln.strip().upper().replace(' ', '').startswith(
                '*MATERIAL,NAME=' + matname.upper()):
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

    # --- emit ----------------------------------------------------------
    #
    # The volumetric term goes in as U6 PATCH ELEMENTS, not as *EQUATION +
    # SPRING1. That first design was verified correct in every part -- and in
    # ccx it still came out 1.55x too stiff on a confined block, because
    # adjacent patches overlap so each mesh dof appears as an independent term
    # in ~24 equations and ccx's MPC machinery does not survive that. One
    # patch alone reproduced an independent Python assembly exactly; the full
    # set did not. As elements the MPCs are gone.
    #
    # *USER ELEMENT fixes NODES= per TYPE, and ring sizes vary, so one type is
    # declared per distinct ring size.
    nodes = sorted(Va)
    ring = {a: [a] + sorted(n for n in b[a] if n != a) for a in nodes}
    bysize = {}
    for a in nodes:
        bysize.setdefault(len(ring[a]), []).append(a)
    alnum = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    if len(bysize) > len(alnum):
        raise SystemExit('nodalbbar: %d distinct ring sizes, more than the '
                         'available U6 type suffixes' % len(bysize))
    suffix = {sz: alnum[i] for i, sz in enumerate(sorted(bysize))}

    out, done_step = [], False
    eid = 10 ** 8
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
                out.append('** SPAX nodal-averaged B-bar: %d patch elements'
                           % len(nodes))
                for sz in sorted(bysize):
                    t = 'U6' + suffix[sz]
                    out.append('*USER ELEMENT,TYPE=%s,NODES=%d,'
                               'INTEGRATIONPOINTS=1,MAXDOF=3' % (t, sz))
                    out.append('*ELEMENT,TYPE=%s,ELSET=BBARP' % t)
                    for a in bysize[sz]:
                        eid += 1
                        r = ring[a]
                        row = [str(eid)] + [str(x) for x in r]
                        # ccx takes at most 16 fields per line
                        out.append(','.join(row[:16]))
                        for k in range(16, len(row), 16):
                            out.append(','.join(row[k:k + 16]))
                out.append('*SOLID SECTION,ELSET=BBARP,MATERIAL=%s' % matname)
        out.append(ln)

    open(dst, 'w').write('\n'.join(out))
    rsz = [len(ring[a]) for a in nodes]
    print('nodalbbar: %s -> %s' % (src, dst))
    print('  %d elements retyped to U5 in ELSET=%s' % (len(els), elset))
    print('  %d U6 patch elements, ring %d..%d nodes (mean %.1f), '
          '%d types' % (len(nodes), min(rsz), max(rsz),
                        sum(rsz) / len(rsz), len(bysize)))
    print('  K = %.4e from material %s' % (K, matname))
    print('  volumetric term delivered as elements -- no MPCs')


if __name__ == '__main__':
    main()
