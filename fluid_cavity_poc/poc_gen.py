"""Single interior spherical cavity in a cube -> two Abaqus decks.

  poc_fluidcavity.inp : the sphere is an EMPTY cavity filled by an Abaqus
                        *FLUID CAVITY (compressible hydraulic fluid = brine).
  poc_void.inp        : identical mesh, sphere left as an open pore (no fluid).

Run both, extract E_eff / nu_eff (see poc_extract.py), and compare: the fluid
cavity should resist lateral contraction -> nu_eff noticeably higher than the
open void, approaching the matrix as K_brine grows. That is the incompressible
coupling the current soft-solid Liquid model (E~0.13 GPa, nu~0.49) cannot
reproduce, because its tiny E makes the high nu irrelevant.

Non-periodic on purpose: a single centred interior sphere keeps the cavity
surface naturally CLOSED (the #1 failure mode in a periodic RVE, where pockets
get cut by the boundary). Minimal RBM constraints leave the lateral faces free
so Poisson contraction is visible.

Needs: python3 + gmsh module. Abaqus is only needed to solve (run on CSC).
Usage:  python3 poc_gen.py [output_dir]
"""
import gmsh, numpy as np, os, sys, contextlib

# ---- parameters ----
L   = 1.0
r   = 0.30                      # interior sphere, r < L/2  -> VoF ~ 11.3%
ctr = (L/2, L/2, L/2)
lc  = 0.09                      # target element size
E_m, nu_m = 9.5e9, 0.33         # ice matrix
K_brine   = 2.2e9               # brine bulk modulus (compressibility)
rho_brine = 1000.0              # reference density (required, irrelevant in statics)
disp      = 0.01 * L            # 1% nominal axial strain
outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("poc")
box = gmsh.model.occ.addBox(0, 0, 0, L, L, L)
sph = gmsh.model.occ.addSphere(ctr[0], ctr[1], ctr[2], r)
# matrix = box MINUS sphere -> the sphere interior is an empty cavity
out, _ = gmsh.model.occ.cut([(3, box)], [(3, sph)], removeObject=True, removeTool=True)
gmsh.model.occ.synchronize()
mat_vol = out[0][1]

# Spherical wall = matrix boundary surface NOT lying on a cube face.
wall_surfs = []
for dim, s in gmsh.model.getBoundary([(3, mat_vol)], oriented=False):
    com = gmsh.model.occ.getCenterOfMass(2, s)
    on_face = any(abs(com[i]) < 1e-6 or abs(com[i] - L) < 1e-6 for i in range(3))
    if not on_face:
        wall_surfs.append(s)

gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc * 0.4)
with open(os.devnull, "w") as dn, contextlib.redirect_stdout(dn):
    gmsh.model.mesh.generate(3)

# ---- read mesh ----
ntags, ncoords, _ = gmsh.model.mesh.getNodes()
ncoords = ncoords.reshape(-1, 3)
coord = {int(t): ncoords[i] for i, t in enumerate(ntags)}

et, etags, enodes = gmsh.model.mesh.getElements(3, mat_vol)
tet_i = list(et).index(4)                          # type 4 = 4-node tet
tets = enodes[tet_i].reshape(-1, 4).astype(int)
tet_ids = etags[tet_i].astype(int)

wall_tris = []
for s in wall_surfs:
    tt, _, tn = gmsh.model.mesh.getElements(2, s)
    if 2 in tt:                                     # type 2 = 3-node tri
        ti = list(tt).index(2)
        wall_tris.append(tn[ti].reshape(-1, 3).astype(int))
wall_tris = np.vstack(wall_tris) if wall_tris else np.zeros((0, 3), int)

# Map face node-triple -> (elem id, Abaqus C3D4 face label).
#   S1:1-2-3  S2:1-4-2  S3:2-4-3  S4:3-4-1
FACES = [('S1', (0, 1, 2)), ('S2', (0, 3, 1)), ('S3', (1, 3, 2)), ('S4', (2, 3, 0))]
face_map = {}
for eid, tet in zip(tet_ids, tets):
    for lbl, idx in FACES:
        face_map[frozenset(int(tet[i]) for i in idx)] = (int(eid), lbl)

cav_faces = []
for tri in wall_tris:
    key = frozenset(int(x) for x in tri)
    if key in face_map:
        cav_faces.append(face_map[key])
assert len(cav_faces) == len(wall_tris), \
    "unmapped wall facets ({} of {})".format(len(cav_faces), len(wall_tris))
# The matrix element's outward face normal at the wall points INTO the sphere
# (the matrix is outside the cavity) -> exactly the orientation a *FLUID CAVITY
# wants. No flip needed.

# ---- boundary node sets ----
xs = {t: c[0] for t, c in coord.items()}
ys = {t: c[1] for t, c in coord.items()}
zs = {t: c[2] for t, c in coord.items()}
eps = 1e-6
XMIN = sorted(t for t in coord if abs(xs[t]) < eps)
XMAX = sorted(t for t in coord if abs(xs[t] - L) < eps)
def near(t, x, y, z): return abs(xs[t]-x) < eps and abs(ys[t]-y) < eps and abs(zs[t]-z) < eps
n000 = next(t for t in coord if near(t, 0, 0, 0))
n0L0 = next(t for t in coord if near(t, 0, L, 0))
n00L = next(t for t in coord if near(t, 0, 0, L))
ref_id = max(coord) + 1

print("matrix tets={}  wall tris={}  XMIN/XMAX={}/{}  cavity ref node={}".format(
    len(tets), len(wall_tris), len(XMIN), len(XMAX), ref_id))

# ---- writers ----
def write_nset(f, name, ids):
    f.write("*Nset, nset={}\n".format(name))
    for i in range(0, len(ids), 10):
        f.write(", ".join(str(x) for x in ids[i:i+10]) + ",\n")

def write_common(f, with_cavity):
    f.write("*Node\n")
    for t in sorted(coord):
        c = coord[t]
        f.write("{}, {:.10e}, {:.10e}, {:.10e}\n".format(t, c[0], c[1], c[2]))
    if with_cavity:
        f.write("{}, {:.10e}, {:.10e}, {:.10e}\n".format(ref_id, *ctr))
    f.write("*Element, type=C3D4\n")
    for eid, tet in zip(tet_ids, tets):
        f.write("{}, {}, {}, {}, {}\n".format(eid, tet[0], tet[1], tet[2], tet[3]))
    f.write("*Elset, elset=ALLMAT, generate\n{}, {}, 1\n".format(tet_ids.min(), tet_ids.max()))
    write_nset(f, "XMIN", XMIN)
    write_nset(f, "XMAX", XMAX)
    write_nset(f, "N000", [n000]); write_nset(f, "N0L0", [n0L0]); write_nset(f, "N00L", [n00L])
    if with_cavity:
        write_nset(f, "CAVREF", [ref_id])
        f.write("*Surface, type=ELEMENT, name=CAVWALL\n")
        for eid, lbl in cav_faces:
            f.write("{}, {}\n".format(eid, lbl))
    f.write("*Solid Section, elset=ALLMAT, material=ICE\n")
    f.write("*Material, name=ICE\n*Elastic\n{:.4e}, {}\n".format(E_m, nu_m))
    if with_cavity:
        f.write("*Fluid Behavior, name=BRINE\n")
        f.write("*Fluid Density\n{:.1f},\n".format(rho_brine))
        f.write("*Fluid Bulk Modulus\n{:.6e},\n".format(K_brine))
        f.write("*Fluid Cavity, name=BRINE_CAV, ref node=CAVREF, "
                "surface=CAVWALL, behavior=BRINE\n")
        f.write("\n")   # *Fluid Cavity needs an (empty) data line, else Abaqus
                        # errors "NO DATA CARDS PROVIDED FOR *FLUID CAVITY"

def write_step(f, with_cavity):
    f.write("*Step, name=Tension, nlgeom=NO\n*Static\n0.1, 1., 1e-08, 1.\n")
    f.write("*Boundary\n")
    f.write("XMIN, 1, 1, 0.0\n")
    f.write("N000, 2, 2, 0.0\nN000, 3, 3, 0.0\nN0L0, 3, 3, 0.0\nN00L, 2, 2, 0.0\n")
    if with_cavity:
        f.write("CAVREF, 1, 3, 0.0\n")             # pin the free cavity ref node
    f.write("*Boundary\nXMAX, 1, 1, {:.6e}\n".format(disp))
    f.write("*Output, field\n*Node Output\nU, RF\n")
    f.write("*Element Output, directions=YES\nS, E, EVOL\n")
    if with_cavity:
        # PCAV (cavity pressure) and CVOL (cavity volume) are nodal output on
        # the cavity reference node -- there is no *Fluid Cavity Output keyword.
        f.write("*Output, history\n*Node Output, nset=CAVREF\nPCAV, CVOL\n")
    f.write("*End Step\n")

for name, cav in [("fluidcavity", True), ("void", False)]:
    path = os.path.join(outdir, "poc_{}.inp".format(name))
    with open(path, "w") as f:
        f.write("*Heading\n Fluid-cavity POC ({}) L={} r={} K_brine={:.2e}\n".format(
            name, L, r, K_brine))
        write_common(f, cav)
        write_step(f, cav)
    print("wrote", path)

gmsh.finalize()
