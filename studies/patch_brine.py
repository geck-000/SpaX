"""Brine-modulus sensitivity: patch the inclusion *Elastic card on the two
base meshes (out_brine/) to make the modulus-sweep decks. Geometry is identical
across the sweep; only the soft-phase E,nu (derived from K,G) changes.

Variants (per microstructure iso & chan):
  G sweep (K=2.2 GPa fixed): G x{0.1,1,10,100,1000} of baseline G=4.4e5 Pa
  K sweep (G=4.4e5 fixed)  : K x{0.1,10} of baseline K=2.2e9 Pa
Writes the patched decks + rve_brine.csv (run_id list for post-processing) and
removes the unsuffixed baseline decks so the solver only runs the variants.
"""
import os, csv, glob

OUT = "out_brine"
MODES = ("utx", "utz")
MICROS = ("iso", "chan")
K0, G0 = 2.2e9, 440029.33528897085

def E_nu(K, G):
    E = 9.0 * K * G / (3.0 * K + G)
    nu = (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G))
    return E, nu

# (tag, K, G)
VARIANTS = [
    ("G0p1",  K0, G0 * 0.1),
    ("G1",    K0, G0),            # baseline
    ("G10",   K0, G0 * 10),
    ("G100",  K0, G0 * 100),
    ("G1000", K0, G0 * 1000),
    ("K0p1",  K0 * 0.1, G0),
    ("K10",   K0 * 10,  G0),
]

def patch_deck(src, dst, E, nu):
    with open(src) as f:
        lines = f.readlines()
    out, i = [], 0
    while i < len(lines):
        out.append(lines[i])
        if lines[i].strip() == "*Material, name=Mat_Inclusion":
            out.append(lines[i+1])                 # *Elastic
            out.append("{}, {}\n".format(repr(E), repr(nu)))  # patched data line
            i += 3                                  # skip orig *Elastic + data
            continue
        i += 1
    with open(dst, "w") as f:
        f.writelines(out)

def main():
    # read a baseline CSV row as the template for rve_brine.csv
    with open("rve_brine_base.csv") as f:
        base_rows = {r["run_id"]: r for r in csv.DictReader(f)}
    fieldnames = list(next(iter(base_rows.values())).keys())

    new_rows = []
    n = 0
    for micro in MICROS:
        tmpl = base_rows["BRINE_" + micro]
        for tag, K, G in VARIANTS:
            E, nu = E_nu(K, G)
            rid = "BRINE_{}_{}".format(micro, tag)
            for m in MODES:
                src = os.path.join(OUT, "Job-BRINE_{}-{}.inp".format(micro, m))
                dst = os.path.join(OUT, "Job-{}-{}.inp".format(rid, m))
                patch_deck(src, dst, E, nu)
                n += 1
            r = dict(tmpl); r["run_id"] = rid
            r["K_inclusion"] = repr(K); r["G_inclusion"] = repr(G)
            new_rows.append(r)
    # write run-id list for post-processing
    with open("rve_brine.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader(); w.writerows(new_rows)
    # remove unsuffixed baseline decks so solver runs only variants
    for micro in MICROS:
        for m in MODES:
            p = os.path.join(OUT, "Job-BRINE_{}-{}.inp".format(micro, m))
            if os.path.exists(p): os.remove(p)
    print("patched {} decks, {} variant RVEs -> rve_brine.csv".format(n, len(new_rows)))
    # sanity: show patched E,nu per variant
    for tag, K, G in VARIANTS:
        E, nu = E_nu(K, G)
        print("  {:6s} K={:.2e} G={:.2e} -> E={:.3e} nu={:.5f}".format(tag, K, G, E, nu))

if __name__ == "__main__":
    main()
