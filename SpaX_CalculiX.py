#!/usr/bin/env python3
"""
SpaX_CalculiX.py -- solve the SpaX RVE decks with CalculiX (ccx) instead of
Abaqus, and read the results back as the same quantities SpaX_PostProcess
pulls out of an ODB.

The generator (SpaX_Standalone.py) writes an Abaqus deck. CalculiX reads a
dialect of the same language, so the deck does not have to be regenerated --
it has to be TRANSLATED, and the translation is small but not optional:

  *PART/*INSTANCE/*ASSEMBLY   CalculiX has no part-instance model. The mesh is
                              flattened to one global node/element numbering.
                              (The generator already emits globally unique
                              labels -- part nodes 1..N, reference points at
                              N+100 and beyond -- so flattening is a pure
                              deletion of the wrapper cards.)
  *EQUATION                   CalculiX resolves node SET names everywhere else,
                              but inside *EQUATION it wants a node NUMBER. Every
                              `RP-1, 1, -1` term is rewritten to `146, 1, -1`.
                              This is the one silent-failure trap in the whole
                              translation: left alone, ccx reports
                              "*ERROR reading *EQUATION" once per PBC equation
                              and stops, so it is loud rather than silent, but
                              a converter that dropped those terms instead
                              would produce an unconstrained cell that solves
                              happily and means nothing.
  C3D4H / C3D10H              CalculiX has no hybrid (mixed u/p) elements. The
                              H is stripped -- see the warning in strip_hybrid.
  *OUTPUT/*NODE OUTPUT/...    replaced by *NODE PRINT / *EL PRINT, which write
                              the .dat text file this module reads back.
  *AMPLITUDE, time=STEP TIME  ccx does not know the TIME parameter; step time
                              is its default, so dropping it preserves meaning.
  *STEP, name=, nlgeom=NO     ccx takes neither; NLGEOM is a bare flag.

What comes back is NOT a reduced result. The .dat carries per-element stress,
strain and volume exactly as the ODB does, so the volume-averaging and the
modulus fit are the same arithmetic on the same numbers -- this module hands
its frames to SpaX_PostProcess._reduce_* rather than re-deriving E, nu and G.

Usage
-----
    python3 SpaX_CalculiX.py convert <dir|deck.inp> [outdir]
    python3 SpaX_CalculiX.py solve   <dir> [--cpus N] [--jobs M]
    python3 SpaX_CalculiX.py all     <dir> [--cpus N] [--jobs M]

then post-process exactly as for Abaqus, except with plain python3 and no
licence:

    python3 SpaX_PostProcess.py params.csv <dir> results.csv

Converted decks are written next to the originals as `<job>-ccx.inp`, so
`Job-RVE_a-utx.inp` (Abaqus) and `Job-RVE_a-utx-ccx.inp` (CalculiX) coexist and
the post-processor can tell which solver produced a given result.
"""

from __future__ import print_function

import os
import re
import subprocess
import sys

import numpy as np


# =====================================================================
# DECK CONVERSION
# =====================================================================

class DeckError(ValueError):
    """A deck card this converter does not know how to translate.

    Raised rather than passed through. An unrecognised card that ccx also does
    not recognise costs a warning and a wrong answer; one that ccx recognises
    with different semantics costs a wrong answer with no warning at all.
    """


# Cards that exist only to describe the Abaqus part-instance hierarchy, which
# CalculiX does not have. Deleting them flattens the model; nothing else in the
# deck depends on them once the `PART-1-1.` prefixes are stripped off the
# equation terms.
_WRAPPER_CARDS = frozenset([
    '*PART', '*END PART', '*ASSEMBLY', '*END ASSEMBLY',
    '*INSTANCE', '*END INSTANCE', '*PREPRINT',
])

# Modes whose extraction never reads a strain field. Printing E for them roughly
# doubles a .dat that is already the largest artefact of the run: on a 3M-element
# cell the strain block alone is a couple of GB, and extract_second_order reads
# S, EVOL and the mesh only. Uniaxial modes DO need it (Poisson's ratio is a
# ratio of volume-averaged strains) and shear modes keep it so the full 6x6
# tensor route has both tensors available.
_NO_STRAIN_MODES = frozenset(['bend', 'tors'])


# CalculiX reads free-format input by copying each comma-separated field into a
# fixed 20-character buffer, and a longer field is TRUNCATED rather than
# rejected. Python's repr routinely exceeds that: a coordinate of
# 3.8163916471489756E-17 is 22 characters, truncates to "3.8163916471489756E-",
# and ccx stops with "*ERROR reading *NODE".
#
# The error is the lucky case. Truncation only fails loudly when it lands
# mid-exponent; a field that truncates to a still-valid number is accepted
# silently at the wrong value. An *Equation coefficient of 3.469446951953614e-18
# (21 characters) truncates to 3.469446951953614e-1 -- a periodicity constraint
# off by seventeen orders of magnitude, in a deck that solves without a murmur.
#
# So every numeric field is re-emitted at a width that cannot truncate rather
# than only the ones observed to break. 12 significant digits spans at most 19
# characters ("-1.23456789012e-308") and is far more precision than mesh
# coordinates or constraint coefficients carry meaning at.
_CCX_FIELD_WIDTH = 20
_NUM_DIGITS = 12

_INT_RE = re.compile(r'^[-+]?\d+$')
_FLOAT_RE = re.compile(r'^[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?$')


def _fmt_field(tok):
    """One data-line field, guaranteed to survive ccx's 20-character buffer.

    Integers (node and element labels, degrees of freedom) are passed through
    exactly -- reformatting a label would be a correctness bug, and no label in
    these decks comes close to the limit. Set names are passed through. Floats
    are re-rendered at 12 significant digits.
    """
    t = tok.strip()
    if not t or _INT_RE.match(t):
        return t
    if _FLOAT_RE.match(t):
        v = float(t.replace('D', 'E').replace('d', 'e'))
        s = '%.*g' % (_NUM_DIGITS, v)
        if len(s) > _CCX_FIELD_WIDTH:
            # Unreachable for finite doubles at 12 digits, but a wrong number
            # here is invisible, so check rather than assume.
            raise DeckError(
                "cannot render {} within {} characters".format(
                    t, _CCX_FIELD_WIDTH))
        return s
    if len(t) > _CCX_FIELD_WIDTH:
        raise DeckError(
            "field '{}' is {} characters; CalculiX truncates at {}".format(
                t, len(t), _CCX_FIELD_WIDTH))
    return t


def _fmt_data(line):
    """Re-emit a comma-separated data line with every field width-safe."""
    return ', '.join(_fmt_field(t) for t in line.split(','))


def _card(line):
    """The keyword of a card line, upper-cased, without its parameters."""
    return line.strip().split(',')[0].strip().upper()


def _params(line):
    """{PARAM: value} for a card line. Values keep their original case."""
    out = {}
    for part in line.strip().split(',')[1:]:
        if '=' in part:
            k, v = part.split('=', 1)
            out[k.strip().upper()] = v.strip()
        else:
            out[part.strip().upper()] = ''
    return out


def strip_hybrid(etype):
    """C3D10H -> C3D10. Returns (type, was_hybrid).

    CalculiX has no hybrid elements, so the mixed displacement/pressure
    formulation is simply lost. That matters only for a phase the generator
    chose H for, i.e. one with nu >= SPAX_HYBRID_NU (default 0.45) -- in these
    decks the brine, at nu ~= 0.49.

    A near-incompressible phase meshed with plain C3D4 volumetrically locks: the
    constant-strain tetrahedron cannot represent an isochoric deformation, so
    the phase comes back too stiff. C3D10 is far more forgiving. If a deck is
    going to CalculiX and it has a near-incompressible phase, generate it with
    SPAX_MESH_ORDER=2; the converter warns when it strips H off a linear
    element for exactly this reason.
    """
    e = etype.strip().upper()
    if e.endswith('H'):
        return e[:-1], True
    return e, False


def _resolve_node(token, nsets, where):
    """An *EQUATION term's node reference -> an integer node label.

    The generator writes either `PART-1-1.1234` (an instance-qualified mesh
    node) or `RP-1` (an assembly node set holding exactly one reference point).
    CalculiX wants the bare number in both cases.
    """
    t = token.strip()
    m = re.match(r'^[A-Za-z0-9_\-]+\.(\d+)$', t)
    if m:
        return int(m.group(1))
    if re.match(r'^\d+$', t):
        return int(t)
    labels = nsets.get(t.upper())
    if labels is None:
        raise DeckError(
            "{}: *EQUATION refers to '{}', which is neither a node number nor a "
            "node set defined in this deck".format(where, t))
    if len(labels) != 1:
        raise DeckError(
            "{}: *EQUATION refers to node set '{}', which holds {} nodes. "
            "CalculiX needs a single node number per term; only single-node "
            "reference-point sets can be translated.".format(
                where, t, len(labels)))
    return labels[0]


def _scan_nsets(lines):
    """{SETNAME: [node labels]} for every *Nset in the deck.

    A first pass, because *EQUATION may name a set defined further down and the
    emit pass has to resolve it to a number immediately.
    """
    nsets = {}
    name = None
    for raw in lines:
        s = raw.strip()
        if not s or s.startswith('**'):
            continue
        if s.startswith('*'):
            name = None
            if _card(s) == '*NSET':
                p = _params(s)
                if 'GENERATE' in p:
                    raise DeckError(
                        "*Nset with GENERATE is not translated (no deck the "
                        "generator writes uses it): " + s)
                name = p.get('NSET', '').upper()
                nsets.setdefault(name, [])
            continue
        if name is not None:
            for tok in s.split(','):
                tok = tok.strip()
                if tok:
                    nsets[name].append(int(tok))
    return nsets


def read_mesh(inp_path):
    """(nodes, elements) from an Abaqus or CalculiX deck.

    nodes    : {label: (x, y, z)}
    elements : {label: (type, elset, [node labels])}

    Used for element centroids in the second-order extraction. The ODB route
    takes centroids from the UNDEFORMED instance geometry, so reading them from
    the deck gives the same numbers rather than an approximation of them.
    """
    nodes = {}
    elements = {}
    section = None
    etype = None
    elset = None
    with open(inp_path, 'r') as f:
        for raw in f:
            s = raw.strip()
            if not s or s.startswith('**'):
                continue
            if s.startswith('*'):
                c = _card(s)
                if c == '*NODE':
                    section = 'node'
                elif c == '*ELEMENT':
                    section = 'element'
                    p = _params(s)
                    etype = p.get('TYPE', '')
                    elset = p.get('ELSET', '')
                else:
                    section = None
                continue
            if section == 'node':
                parts = s.split(',')
                nodes[int(parts[0])] = (float(parts[1]), float(parts[2]),
                                        float(parts[3]))
            elif section == 'element':
                parts = [int(p) for p in s.split(',') if p.strip()]
                elements[parts[0]] = (etype, elset, parts[1:])
    return nodes, elements


def convert_deck(src_path, dst_path, frd=False):
    """Translate one Abaqus deck into a CalculiX deck.

    Returns a metadata dict: mode, element types, the elsets and reference-point
    sets that were given output requests, and how many element blocks lost their
    hybrid formulation.
    """
    with open(src_path, 'r') as f:
        lines = f.read().splitlines()

    nsets = _scan_nsets(lines)
    src_name = os.path.basename(src_path)

    # Which mode this deck is, taken from the step name the generator wrote.
    # Only used to decide whether the strain field is worth printing.
    mode = None
    for s in lines:
        if _card(s) == '*STEP':
            nm = _params(s).get('NAME', '').lower()
            mode = {'step-utx': 'utx', 'step-uty': 'uty', 'step-utz': 'utz',
                    'step-ss12': 'ss12', 'step-ss13': 'ss13',
                    'step-ss23': 'ss23', 'step-bending': 'bend',
                    'step-torsion': 'tors'}.get(nm)
            break
    want_strain = mode not in _NO_STRAIN_MODES

    out = ['*HEADING',
           'SpaX RVE deck translated for CalculiX from {}'.format(src_name)]
    elsets = []            # element sets that get an *EL PRINT
    rp_sets = []           # single-node sets that get a *NODE PRINT
    etypes = {}            # elset -> translated element type
    n_hybrid = 0
    hybrid_linear = False

    section = None         # what the current data lines belong to
    eq_left = 0            # *EQUATION terms still to read
    i = 0
    n = len(lines)

    while i < n:
        raw = lines[i]
        s = raw.strip()
        i += 1

        if not s or s.startswith('**'):
            continue

        if s.startswith('*'):
            c = _card(s)
            p = _params(s)
            section = None

            if c in _WRAPPER_CARDS:
                continue

            if c == '*HEADING':
                continue                      # replaced by ours, above

            if c == '*NODE':
                out.append('*NODE')
                section = 'pass'
                continue

            if c == '*ELEMENT':
                et, was_h = strip_hybrid(p.get('TYPE', ''))
                es = p.get('ELSET', '')
                if was_h:
                    n_hybrid += 1
                    if not et.endswith('10'):
                        hybrid_linear = True
                etypes[es] = et
                if es not in elsets:
                    elsets.append(es)
                out.append('*ELEMENT, TYPE={}, ELSET={}'.format(et, es))
                section = 'pass'
                continue

            if c == '*ELSET':
                gen = ', GENERATE' if 'GENERATE' in p else ''
                es = p.get('ELSET', '')
                if es not in elsets:
                    elsets.append(es)
                out.append('*ELSET, ELSET={}{}'.format(es, gen))
                section = 'pass'
                continue

            if c == '*NSET':
                # `instance=` is meaningless once the assembly is flattened.
                name = p.get('NSET', '')
                out.append('*NSET, NSET={}'.format(name))
                if len(nsets.get(name.upper(), [])) == 1:
                    rp_sets.append(name)
                section = 'pass'
                continue

            if c == '*SOLID SECTION':
                out.append('*SOLID SECTION, ELSET={}, MATERIAL={}'.format(
                    p.get('ELSET', ''), p.get('MATERIAL', '')))
                # Abaqus wants a (here empty) thickness/data line; ccx does not
                # accept one for a solid section, so drop it.
                if i < n and lines[i].strip() in (',', ''):
                    i += 1
                section = None
                continue

            if c == '*EQUATION':
                out.append('*EQUATION')
                # Next data line is the term count; the terms follow one per
                # line. ccx accepts that layout (it does not require the
                # 4-per-line packing of the manual's examples).
                section = 'eqcount'
                continue

            if c == '*MATERIAL':
                out.append('*MATERIAL, NAME={}'.format(p.get('NAME', '')))
                section = None
                continue

            if c == '*ELASTIC':
                out.append('*ELASTIC')
                section = 'pass'
                continue

            if c == '*BOUNDARY':
                amp = p.get('AMPLITUDE')
                out.append('*BOUNDARY' + (', AMPLITUDE={}'.format(amp)
                                          if amp else ''))
                section = 'pass'
                continue

            if c == '*AMPLITUDE':
                # TIME=STEP TIME is ccx's default and ccx rejects the parameter
                # by name ("parameter not recognized"), so drop it. TOTAL TIME
                # it does understand, so keep that one if it ever appears.
                tp = p.get('TIME', '').upper()
                card = '*AMPLITUDE, NAME={}'.format(p.get('NAME', ''))
                if tp == 'TOTAL TIME':
                    card += ', TIME=TOTAL TIME'
                elif tp and tp != 'STEP TIME':
                    raise DeckError("unsupported *Amplitude TIME=" + tp)
                out.append(card)
                section = 'pass'
                continue

            if c == '*STEP':
                inc = p.get('INC')
                nlg = p.get('NLGEOM', '').upper()
                card = '*STEP'
                if inc:
                    card += ', INC={}'.format(inc)
                if nlg in ('YES', 'ON', 'TRUE'):
                    card += ', NLGEOM'
                out.append(card)
                section = None
                continue

            if c == '*STATIC':
                # ccx's default is SPOOLES, a direct solver whose memory grows
                # far faster than the model. A stock ccx is often linked
                # against nothing else (check with `ldd`: no PARDISO, no
                # PaStiX), and the campaign's production cells run to millions
                # of elements, where SPOOLES will not fit. SPAX_CCX_SOLVER
                # switches to one of ccx's built-in iterative solvers, which
                # need very little memory -- at the cost of convergence trouble
                # on exactly this kind of model, where a 70x phase contrast and
                # a near-incompressible phase make the system ill-conditioned.
                # Left unset the deck says nothing and ccx uses its default.
                solver = os.environ.get('SPAX_CCX_SOLVER', '').strip().upper()
                if solver:
                    allowed = ('SPOOLES', 'ITERATIVE CHOLESKY',
                               'ITERATIVE SCALING', 'PARDISO', 'PASTIX')
                    if solver not in allowed:
                        raise DeckError(
                            "SPAX_CCX_SOLVER={} is not one of {}".format(
                                solver, ', '.join(allowed)))
                    out.append('*STATIC, SOLVER={}'.format(solver))
                else:
                    out.append('*STATIC')
                section = 'pass'
                continue

            if c in ('*OUTPUT', '*NODE OUTPUT', '*ELEMENT OUTPUT'):
                # Abaqus output requests. ccx's equivalents are synthesised at
                # *END STEP instead, because what has to be requested depends on
                # the element and node SETS in the deck rather than on the field
                # names Abaqus uses.
                section = 'drop'
                continue

            if c == '*END STEP':
                out.extend(_output_requests(elsets, rp_sets, want_strain, frd))
                out.append('*END STEP')
                section = None
                continue

            raise DeckError("{}: unhandled card '{}'".format(src_name, s))

        # ---- data lines ----
        if section == 'drop' or section is None:
            continue

        if section == 'pass':
            out.append(_fmt_data(s))
            continue

        if section == 'eqcount':
            eq_left = int(s.split(',')[0])
            out.append(str(eq_left))
            section = 'eqterm'
            continue

        if section == 'eqterm':
            parts = [t.strip() for t in s.split(',')]
            if len(parts) < 3:
                raise DeckError("{}: malformed *Equation term '{}'".format(
                    src_name, s))
            node = _resolve_node(parts[0], nsets, src_name)
            out.append('{}, {}, {}'.format(node, _fmt_field(parts[1]),
                                           _fmt_field(parts[2])))
            eq_left -= 1
            if eq_left == 0:
                section = None
            continue

    with open(dst_path, 'w') as f:
        f.write('\n'.join(out) + '\n')

    return {
        'src': src_path,
        'dst': dst_path,
        'mode': mode,
        'elsets': elsets,
        'rp_sets': rp_sets,
        'etypes': etypes,
        'n_hybrid_stripped': n_hybrid,
        'hybrid_linear': hybrid_linear,
        'strain_printed': want_strain,
    }


def _output_requests(elsets, rp_sets, want_strain, frd):
    """The ccx result requests that stand in for the Abaqus *Output block.

    S, E and EVOL must share ONE *EL PRINT card per set. Splitting them across
    two cards on the same set makes ccx print the volume block twice under a
    single header -- harmless once you know, but it silently doubles any total
    taken over the block.
    """
    req = []
    for nm in rp_sets:
        # U carries the imposed macro deformation (the post-processor reads it
        # in preference to the step-time proxy); RF at the driven reference
        # point is the macroscopic stress resultant, and is the independent
        # check on the volume average.
        req.append('*NODE PRINT, NSET={}'.format(nm))
        req.append('U,RF')
    keys = 'S,E,EVOL' if want_strain else 'S,EVOL'
    for es in elsets:
        req.append('*EL PRINT, ELSET={}'.format(es))
        req.append(keys)
    if frd:
        # Only for looking at the cell in cgx/ParaView. Off by default: the .frd
        # is the same order of magnitude as the .dat and nothing here reads it.
        req.append('*NODE FILE')
        req.append('U')
        req.append('*EL FILE')
        req.append('S,E')
    return req


# =====================================================================
# SOLVING
# =====================================================================

def ccx_executable():
    return os.environ.get('SPAX_CCX', 'ccx')


def solve(inp_path, cpus=None, timeout=None, log_path=None):
    """Run CalculiX on a converted deck. Returns (ok, log_text).

    ccx is invoked with the job name (no extension) and writes every result
    next to it. `ok` means ccx reported a finished job AND left a .dat behind;
    a diverged step leaves the .dat empty or short, which the reader reports
    rather than averaging over nothing.
    """
    job = inp_path[:-4] if inp_path.lower().endswith('.inp') else inp_path
    workdir = os.path.dirname(os.path.abspath(job)) or '.'
    env = dict(os.environ)
    if cpus:
        env['OMP_NUM_THREADS'] = str(cpus)
        env['CCX_NPROC_STIFFNESS'] = str(cpus)
    try:
        proc = subprocess.Popen(
            [ccx_executable(), os.path.basename(job)],
            cwd=workdir, env=env,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            universal_newlines=True)
        log, _ = proc.communicate(timeout=timeout) if timeout else proc.communicate()
        rc = proc.returncode
    except OSError as e:
        return False, "could not run '{}': {}".format(ccx_executable(), e)

    if log_path:
        with open(log_path, 'w') as f:
            f.write(log)

    dat = job + '.dat'
    ok = (rc == 0 and 'Job finished' in log
          and os.path.isfile(dat) and os.path.getsize(dat) > 0)
    return ok, log


# =====================================================================
# READING ccx RESULTS (.dat)
# =====================================================================

# Each result block in a .dat opens with a header naming the quantity, the set
# and the time, e.g.
#   stresses (elem, integ.pnt.,sxx,syy,szz,sxy,sxz,syz) for set MATRIX_ONLY and time  0.1000000E+01
_BLOCK_RE = re.compile(
    r'^\s*(?P<what>[a-z][a-z .]*?)\s*(?:\([^)]*\))?\s*'
    r'for set\s+(?P<set>\S+)\s+and time\s+(?P<time>[-+0-9.eEdD]+)\s*$')

_STEP_RE = re.compile(r'^\s*S\s*T\s*E\s*P\s+(\d+)\s*$')

# Leading word of the header -> (frame key, numbers per data row).
# stresses/strains rows are (element, integration point, 6 components); volume
# rows are (element, volume); nodal rows are (node, 3 components).
_BLOCK_KIND = {
    'displacements': ('disp', 4),
    'forces': ('force', 4),
    'stresses': ('stress', 8),
    'strains': ('strain', 8),
    'volume': ('volume', 2),
}


def _numbers(buf, per_row, what, dat_path):
    """A block's buffered text -> an (nrow, per_row) float array."""
    if not buf:
        return np.zeros((0, per_row))
    arr = np.array(buf.split(), dtype=np.float64)
    if arr.size % per_row:
        raise ValueError(
            "{}: '{}' block has {} numbers, not a multiple of {} per row"
            .format(dat_path, what, arr.size, per_row))
    return arr.reshape(-1, per_row)


def _mean_per_element(labels, data):
    """Collapse the integration-point rows of an element to their mean.

    The same collapse SpaX_PostProcess._field_by_label applies to an ODB field:
    C3D10 has 4 integration points, and the whole-element volume in EVOL weights
    ONE row per element. Without this a quadratic cell over-counts every moment
    and every volume fourfold.
    """
    ulab, inv = np.unique(labels, return_inverse=True)
    if len(ulab) == len(labels):
        return labels.astype(np.int64), data
    sums = np.zeros((len(ulab), data.shape[1]), dtype=float)
    np.add.at(sums, inv, data)
    counts = np.bincount(inv).astype(float)
    return ulab.astype(np.int64), sums / counts[:, None]


def read_dat(dat_path):
    """Parse a ccx .dat into a list of frames, in solution order.

    Each frame is
        {'time': float,
         'disp' : {SET: {node: ndarray(3)}},
         'force': {SET: {node: ndarray(3)}},
         'stress'|'strain': {SET: (labels ndarray(n,), data ndarray(n, 6))},
         'volume': {SET: (labels ndarray(n,), vols ndarray(n,))}}

    with stresses and strains already collapsed to one row per element. A
    linear step gives one frame; ccx solves a linear static step in a single
    increment whatever increment size the *STATIC card asks for, which is the
    right answer and not an approximation of the ten-increment one -- see the
    single-increment branch in SpaX_PostProcess.extract_first_order.
    """
    if not os.path.isfile(dat_path):
        raise IOError("CalculiX .dat not found: " + dat_path)

    frames = []
    by_key = {}
    step = 1
    cur = None            # (frame dict, kind, setname, per_row)
    buf = []

    def flush():
        if cur is None:
            return
        frame, kind, setname, per_row = cur
        rows = _numbers(' '.join(buf), per_row, kind, dat_path)
        if kind in ('disp', 'force'):
            frame[kind][setname] = dict(
                (int(r[0]), r[1:4].copy()) for r in rows)
        elif kind == 'volume':
            frame[kind][setname] = (rows[:, 0].astype(np.int64),
                                    rows[:, 1].copy())
        else:
            lab, dat = _mean_per_element(rows[:, 0].astype(np.int64),
                                         rows[:, 2:8])
            frame[kind][setname] = (lab, dat)

    with open(dat_path, 'r') as f:
        for raw in f:
            line = raw.rstrip('\n')
            if not line.strip():
                continue
            m = _STEP_RE.match(line)
            if m:
                flush()
                cur = None
                buf = []
                step = int(m.group(1))
                continue
            if line.lstrip().upper().startswith('INCREMENT'):
                continue
            m = _BLOCK_RE.match(line)
            if m:
                flush()
                buf = []
                what = m.group('what').strip().split()[0].lower()
                kind_row = _BLOCK_KIND.get(what)
                if kind_row is None:
                    # A quantity nothing here reads (e.g. TOTALS output someone
                    # added by hand). Skip its data rather than misparsing it.
                    cur = None
                    continue
                kind, per_row = kind_row
                t = float(m.group('time').replace('D', 'E'))
                key = (step, round(t, 12))
                frame = by_key.get(key)
                if frame is None:
                    frame = {'time': t, 'step': step, 'disp': {}, 'force': {},
                             'stress': {}, 'strain': {}, 'volume': {}}
                    by_key[key] = frame
                    frames.append(frame)
                cur = (frame, kind, m.group('set'), per_row)
                continue
            if cur is not None:
                buf.append(line)
    flush()

    frames.sort(key=lambda fr: (fr['step'], fr['time']))
    return frames


def job_dat(odb_dir, run_id, short):
    """Path of the ccx .dat for one RVE and load case, or None."""
    p = os.path.join(odb_dir, 'Job-{}-{}-ccx.dat'.format(run_id, short))
    return p if os.path.isfile(p) else None


def job_inp(odb_dir, run_id, short):
    """Path of the ccx deck for one RVE and load case (needed for the mesh)."""
    return os.path.join(odb_dir, 'Job-{}-{}-ccx.inp'.format(run_id, short))


def _set_arrays(frame, kind):
    """Concatenate a frame's per-set element arrays into one (labels, data)."""
    labs, dats = [], []
    for setname in sorted(frame[kind].keys()):
        l, d = frame[kind][setname]
        labs.append(l)
        dats.append(d)
    if not labs:
        return np.zeros(0, dtype=np.int64), np.zeros((0, 6))
    return np.concatenate(labs), np.concatenate(dats, axis=0)


def _volume_map(frame):
    """{element label: volume} across every printed element set."""
    out = {}
    for setname in frame['volume']:
        lab, vol = frame['volume'][setname]
        out.update(zip(lab.tolist(), vol.tolist()))
    return out


def _rp_value(frame, kind, setname, dof_idx):
    """One reference point's U or RF component, or None if it was not printed.

    ccx upper-cases set names in the .dat, so 'RP-1' comes back as 'RP-1' but
    'Fix_Ref_Centre' as 'FIX_REF_CENTRE'; match case-insensitively.
    """
    if setname is None:
        return None
    want = setname.upper()
    for key in frame[kind]:
        if key.upper() == want:
            vals = frame[kind][key]
            if not vals:
                return None
            node = sorted(vals.keys())[0]
            return float(vals[node][dof_idx])
    return None


# =====================================================================
# EXTRACTION -- the same quantities SpaX_PostProcess pulls from an ODB
# =====================================================================

def extract_first_order(dat_path, s_comp, eng_strain, L):
    """E or G and nu from a ccx .dat.

    The twin of SpaX_PostProcess.extract_first_order, reading the same
    per-element stress/strain/volume from a different file. The accumulation
    below mirrors that function line for line; the fit that turns the frame
    series into a modulus is not repeated here -- it is imported, so the two
    routes cannot drift apart in the part that decides the number.
    """
    import SpaX_PostProcess as spp

    s_idx = spp.stress_index(s_comp)
    is_shear = s_comp in ('S12', 'S13', 'S23')
    V_RVE = L ** 3

    rp_name, rp_dof = spp.RP_FOR_SCOMP.get(s_comp, (None, None))
    t1t2 = {'S11': (1, 2), 'S22': (0, 2), 'S33': (0, 1)}.get(s_comp, (None, None))

    frames = read_dat(dat_path)
    if not frames:
        raise RuntimeError(
            "{} holds no result blocks (the step did not converge, or the deck "
            "carried no *EL PRINT request)".format(dat_path))

    stress_strain_data = []
    for frame in frames:
        s_lab, s_dat = _set_arrays(frame, 'stress')
        if len(s_lab) == 0:
            continue
        vol_map = _volume_map(frame)
        w = spp._aligned(vol_map, s_lab.tolist())

        total_stress_vol = float(np.dot(s_dat[:, s_idx], w))
        total_vol = float(w.sum())

        # Achieved soft-phase fraction, from the meshed inclusion set. Voids are
        # not meshed and so are invisible here; they show up in `porosity`.
        total_incl_vol = 0.0
        for setname in frame['volume']:
            if setname.upper() == 'SPHERE_ONLY':
                total_incl_vol = float(frame['volume'][setname][1].sum())

        axial = trans1 = trans2 = 0.0
        e_lab, e_dat = _set_arrays(frame, 'strain')
        if len(e_lab):
            erow = dict((int(l), i) for i, l in enumerate(e_lab.tolist()))
            idx = np.fromiter((erow.get(l, -1) for l in s_lab.tolist()),
                              dtype=int, count=len(s_lab))
            estr = np.zeros((len(s_lab), e_dat.shape[1]))
            valid = idx >= 0
            if valid.any():
                estr[valid] = e_dat[idx[valid]]
            axial = float(np.dot(estr[:, s_idx], w))
            if not is_shear and t1t2[0] is not None:
                trans1 = float(np.dot(estr[:, t1t2[0]], w))
                trans2 = float(np.dot(estr[:, t1t2[1]], w))

        rp_u = _rp_value(frame, 'disp', rp_name, rp_dof)
        eps_macro = (rp_u / L) if rp_u is not None else (frame['time'] * eng_strain)

        stress_strain_data.append({
            'sigma': total_stress_vol / V_RVE,
            'eps': eps_macro,
            'eps_axial_solid': (axial / total_vol) if total_vol > 0 and not is_shear else 0.0,
            'eps_trans1': (trans1 / total_vol) if total_vol > 0 and not is_shear else 0.0,
            'eps_trans2': (trans2 / total_vol) if total_vol > 0 and not is_shear else 0.0,
            'v_solid': total_vol,
            'v_incl': total_incl_vol,
            # Not used by the fit. The reaction at the driven reference point is
            # the macroscopic stress resultant, so RF/L^2 measures the same
            # sigma_bar the volume average produces, by a completely different
            # route: one integrates the element stresses, the other reads a
            # single constraint force. They agree to solver tolerance when the
            # periodic equations are right, and diverge when they are not --
            # which is the failure this translation is most exposed to.
            'sigma_rf': (_rp_value(frame, 'force', rp_name, rp_dof) or 0.0) / (L * L),
        })

    results = spp._reduce_first_order(stress_strain_data, V_RVE, is_shear)
    gap = _equilibrium_gap(stress_strain_data)
    if gap is not None:
        results['equilibrium_gap'] = gap
        if gap > _EQUILIBRIUM_TOL:
            print("    WARNING {}: volume-averaged stress and reference-point "
                  "reaction are {:.3%} apart".format(
                      os.path.basename(dat_path), gap))
    return results


# Relative gap between the two independent measurements of sigma_bar beyond
# which the solve should not be believed. A converged direct solve on a correct
# deck leaves ~1e-9; ccx's ITERATIVE CHOLESKY on the same deck leaves ~2e-3,
# and its E_eff is wrong by 0.15%. So this sits well above the direct-solver
# noise and well below anything that would pass unnoticed.
_EQUILIBRIUM_TOL = 1e-4


def _equilibrium_gap(stress_strain_data):
    """Worst relative disagreement between the two measurements of sigma_bar.

    Both are the macroscopic stress. The volume average integrates the element
    stresses over the cell; the reaction force is the constraint force at the
    reference point driving the periodic jump, and equals sigma_bar * L^2 by
    work conjugacy. Nothing ties them together except the model being right and
    the linear system being solved, so the gap detects two different failures:

      * a wrong constraint set -- an equation dropped in translation, a
        coefficient truncated, a dependent DOF claimed twice -- in a deck that
        nonetheless solves without complaint; and
      * an under-converged iterative solve, which returns a plausible modulus
        that is simply wrong. This is the one that matters in practice, because
        ccx's iterative solvers are the only ones that fit a large cell in
        memory and their tolerance is not adjustable from the deck.

    Returns None when the reaction was not printed for this mode.
    """
    worst = None
    for d in stress_strain_data:
        a, b = d['sigma'], d.get('sigma_rf', 0.0)
        scale = max(abs(a), abs(b))
        if scale <= 0.0 or not b:
            continue
        rel = abs(a - b) / scale
        worst = rel if worst is None else max(worst, rel)
    return worst


def extract_second_order(dat_path, L, Kappa, Bending_Plane, inp_path=None):
    """D_RVE and E_bending from a ccx bending .dat.

    The twin of SpaX_PostProcess.extract_second_order. Element centroids come
    from the deck rather than from a COORD field, which is what the ODB route
    does too -- it builds them from the undeformed instance geometry.
    """
    import SpaX_PostProcess as spp

    if inp_path is None:
        inp_path = dat_path[:-4] + '.inp'
    if not os.path.isfile(inp_path):
        raise IOError(
            "the bending extraction needs the deck for element centroids, but "
            "{} is missing".format(inp_path))

    if Bending_Plane == 'xz':
        s_idx, coord_idx = 0, 2
    elif Bending_Plane == 'yz':
        s_idx, coord_idx = 1, 2
    elif Bending_Plane == 'xy':
        s_idx, coord_idx = 0, 1
    else:
        raise ValueError("unknown Bending_Plane: " + str(Bending_Plane))
    coord_bar = L / 2.0

    nodes, elements = read_mesh(inp_path)
    cent = {}
    for label in elements:
        conn = elements[label][2]
        acc = 0.0
        for nid in conn:
            acc += nodes[nid][coord_idx]
        cent[label] = acc / len(conn)

    frames = read_dat(dat_path)
    if not frames:
        raise RuntimeError(
            "bending .dat {} holds no result blocks (analysis did not "
            "converge)".format(dat_path))
    frame = frames[-1]

    s_lab, s_dat = _set_arrays(frame, 'stress')
    vol_map = _volume_map(frame)
    labs = s_lab.tolist()
    w = spp._aligned(vol_map, labs)
    # An element with no centroid contributes nothing, matching the ODB route's
    # `continue` on a missing centroid.
    w = w * np.fromiter((1.0 if l in cent else 0.0 for l in labs),
                        dtype=float, count=len(labs))
    z_rel = spp._aligned(cent, labs) - coord_bar

    ncomp = min(6, s_dat.shape[1])
    sigma_vol = np.zeros(6)
    sigma_vol[:ncomp] = (s_dat[:, :ncomp] * w[:, None]).sum(axis=0)
    moment_vol = float(np.dot(s_dat[:, s_idx] * z_rel, w))
    total_vol = float(w.sum())

    rpk_u = _rp_value(frame, 'disp', 'RP_K', 0)
    kappa_actual = rpk_u if (rpk_u is not None and abs(rpk_u) > 1e-30) \
        else Kappa * frame['time']

    return spp._reduce_second_order(sigma_vol, moment_vol, total_vol,
                                    kappa_actual, L, s_idx)


def extract_averages(dat_path, L):
    """Volume-averaged stress and strain tensors at the last frame.

    The six numbers per tensor that the full 6x6 elasticity route needs, in the
    same Voigt order and with the same V_RVE = L^3 normalisation
    SpaX_PostProcess.extract_principals uses. The per-slice principal and SCF
    fields that function also writes are NOT produced here.
    """
    frames = read_dat(dat_path)
    if not frames:
        raise RuntimeError("no result blocks in " + dat_path)
    frame = frames[-1]
    V_RVE = L ** 3

    import SpaX_PostProcess as spp
    s_lab, s_dat = _set_arrays(frame, 'stress')
    vol_map = _volume_map(frame)
    w = spp._aligned(vol_map, s_lab.tolist())
    out = {'V_solid': float(w.sum())}
    keys = ['S11', 'S22', 'S33', 'S12', 'S13', 'S23']
    for k in range(6):
        out[keys[k]] = float(np.dot(s_dat[:, k], w)) / V_RVE

    e_lab, e_dat = _set_arrays(frame, 'strain')
    if len(e_lab):
        we = spp._aligned(vol_map, e_lab.tolist())
        keys = ['E11', 'E22', 'E33', 'E12', 'E13', 'E23']
        for k in range(6):
            out[keys[k]] = float(np.dot(e_dat[:, k], we)) / V_RVE
    return out


# =====================================================================
# CLI
# =====================================================================

def _decks_in(path):
    """Abaqus decks to convert: one file, or every Job-*.inp in a directory."""
    if os.path.isfile(path):
        return [path]
    out = []
    for name in sorted(os.listdir(path)):
        if name.endswith('-ccx.inp') or not name.endswith('.inp'):
            continue
        if name.startswith('Job-'):
            out.append(os.path.join(path, name))
    return out


def _ccx_name(src, outdir):
    base = os.path.basename(src)[:-4] + '-ccx.inp'
    return os.path.join(outdir or os.path.dirname(src) or '.', base)


def cmd_convert(args):
    src = args[0]
    outdir = args[1] if len(args) > 1 else None
    if outdir and not os.path.isdir(outdir):
        os.makedirs(outdir)
    decks = _decks_in(src)
    if not decks:
        print("no Job-*.inp decks found in {}".format(src))
        return 1
    warned = False
    for d in decks:
        dst = _ccx_name(d, outdir)
        meta = convert_deck(d, dst, frd=os.environ.get('SPAX_CCX_FRD') == '1')
        note = ''
        if meta['n_hybrid_stripped']:
            note = '  [hybrid dropped x{}]'.format(meta['n_hybrid_stripped'])
        print("  {} -> {}{}".format(os.path.basename(d),
                                    os.path.basename(dst), note))
        if meta['hybrid_linear'] and not warned:
            warned = True
            print("\n  WARNING: this deck asked for HYBRID LINEAR tetrahedra "
                  "(C3D4H), which\n"
                  "  CalculiX does not have. A phase with nu >= 0.45 meshed "
                  "with plain C3D4\n"
                  "  volumetrically locks and comes back too stiff. Regenerate "
                  "with\n"
                  "  SPAX_MESH_ORDER=2 before trusting a CalculiX modulus for "
                  "this cell.\n")
    print("converted {} deck(s)".format(len(decks)))
    return 0


def cmd_solve(args):
    path = args[0]
    cpus = None
    jobs = 1
    i = 1
    while i < len(args):
        if args[i] == '--cpus':
            cpus = int(args[i + 1]); i += 2
        elif args[i] == '--jobs':
            jobs = int(args[i + 1]); i += 2
        else:
            i += 1

    if os.path.isfile(path):
        decks = [path]
    else:
        decks = [os.path.join(path, n) for n in sorted(os.listdir(path))
                 if n.endswith('-ccx.inp')]
    if not decks:
        print("no *-ccx.inp decks found in {} -- run `convert` first".format(path))
        return 1

    print("solving {} deck(s) with {} ({} at a time, {} cpu(s) each)".format(
        len(decks), ccx_executable(), jobs, cpus or 1))

    failed = []
    if jobs <= 1:
        for d in decks:
            ok, _ = solve(d, cpus=cpus, log_path=d[:-4] + '.log')
            print("  {:<44} {}".format(os.path.basename(d),
                                       'ok' if ok else 'FAILED'))
            if not ok:
                failed.append(d)
    else:
        # Several ccx processes at once, each on `cpus` threads. ccx is a single
        # process per job, so this is the only parallelism available across load
        # cases.
        from multiprocessing.pool import ThreadPool
        pool = ThreadPool(jobs)
        results = pool.map(
            lambda d: (d, solve(d, cpus=cpus, log_path=d[:-4] + '.log')[0]),
            decks)
        pool.close()
        pool.join()
        for d, ok in results:
            print("  {:<44} {}".format(os.path.basename(d),
                                       'ok' if ok else 'FAILED'))
            if not ok:
                failed.append(d)

    if failed:
        print("\n{} job(s) failed; see the .log next to each deck".format(
            len(failed)))
        return 1
    print("all jobs finished")
    return 0


def main(argv):
    if len(argv) < 2 or argv[1] in ('-h', '--help', 'help'):
        print(__doc__)
        return 0
    cmd, args = argv[1], argv[2:]
    if cmd == 'convert':
        return cmd_convert(args)
    if cmd == 'solve':
        return cmd_solve(args)
    if cmd == 'all':
        rc = cmd_convert(args[:1])
        return rc or cmd_solve(args)
    print("unknown command: {}\n".format(cmd))
    print(__doc__)
    return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv))
