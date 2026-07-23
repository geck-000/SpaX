#!/usr/bin/env python3
"""Regenerate the study-based revision figures for main_rev.tex.

Reviewer asks (Section E): z/H on the vertical axis, larger fonts, a coordinate
frame consistent with Fig.2. Writes PDF + PNG for each of study_scfdepth,
ice_column_profiles and study_coltensor.
Colours: Okabe-Ito colourblind-safe (matching viz/render_rve.py).

Usage: run from the directory holding the result CSVs (as with the other
analyzers), optionally sending the figures elsewhere for upload:

    cd results && python3 ../analysis/make_rev_figs.py [--out DIR]
"""
import csv, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Shared Okabe-Ito palette, enlarged fonts, and the z/H-on-vertical depth axis
# (see figstyle.py); every depth figure in the paper draws from the same frame.
from figstyle import (BLUE, ORANGE, GREEN, VERM, PURPLE, SKY, BLACK,
                      depth_axis, orient_labels)
from figstyle import apply as _apply_style
_apply_style()

# CSVs are read by bare filename from the working directory, per the
# convention of the other analyzers; figures land alongside them unless
# --out redirects them.
_argv = sys.argv[1:]
RES = "."
OUT = os.path.expanduser(_argv[_argv.index("--out") + 1]) if "--out" in _argv else "."
os.makedirs(OUT, exist_ok=True)


def load(name):
    with open(os.path.join(RES, name)) as f:
        return list(csv.DictReader(f))


def zfrac(run_id):
    """ICE_z05 -> 0.05, ICE_z95 -> 0.95."""
    tok = run_id.lower().split("z")[-1]
    return int(tok) / 100.0


def save(fig, stem):
    for ext in ("pdf", "png"):
        p = os.path.join(OUT, f"{stem}.{ext}")
        fig.savefig(p, dpi=200)
        print("wrote", p)
    plt.close(fig)


# ===========================================================================
# 1) scfdepth : SCF percentiles vs z/H + normalized first-failure macro stress
# ===========================================================================
def fig_scfdepth():
    r = load("results_failure.csv")
    z = np.array([zfrac(x["run_id"]) for x in r])
    p50 = np.array([float(x["SCF_p50"]) for x in r])
    p90 = np.array([float(x["SCF_p90"]) for x in r])
    p99 = np.array([float(x["SCF_p99"]) for x in r])
    mc99 = np.array([float(x["MCnorm_p99"]) for x in r])

    # normalized first-failure macro stress: inversely proportional to the P99
    # criterion, scaled to 1 at the cold surface (relative load capacity)
    ff_scf = p99[0] / p99
    ff_mc = mc99[0] / mc99

    fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 6.2), sharey=True)

    a.fill_betweenx(z, p50, p99, color=BLUE, alpha=0.15,
                    label="P50-P99 spread")
    a.plot(p50, z, "--", color=BLUE, marker="o", label="P50")
    a.plot(p90, z, ":", color=GREEN, marker="s", label="P90")
    a.plot(p99, z, "-", color=VERM, marker="D", label="P99 (robust)")
    a.set_xlabel(r"Matrix SCF $=\sigma_1^{\max}/\bar\sigma_{11}$")
    depth_axis(a)
    a.set_title("(a)")
    a.legend(loc="upper right", framealpha=0.9)

    b.plot(ff_scf, z, "-", color=VERM, marker="D",
           label="max-principal (SCF)")
    b.plot(ff_mc, z, "-", color=PURPLE, marker="^",
           label="Mohr-Coulomb")
    b.set_xlabel("First-failure macro stress\n(normalized to surface)")
    depth_axis(b)
    orient_labels(b)
    b.set_title("(b)")
    b.legend(loc="center left", framealpha=0.9)
    b.set_xlim(0, 1.10)

    fig.tight_layout()
    save(fig, "study_scfdepth")


# ===========================================================================
# 2) ice_column_profiles : effective moduli + phase fractions vs z/H
# ===========================================================================
def fig_column():
    r = load("results_column.csv")
    z = np.array([zfrac(x["run_id"]) for x in r])
    Exy = np.array([float(x["E_x"]) for x in r]) / 1e9
    Ez = np.array([float(x["E_z"]) for x in r]) / 1e9
    brine = np.array([float(x["VoF_incl_sphere"]) for x in r]) * 100
    chan = np.array([float(x["channel_vof_target"]) for x in r]) * 100
    gas = np.array([float(x["VoF_void_sphere"]) for x in r]) * 100

    fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 6.2), sharey=True)

    a.plot(Exy, z, "-", color=BLUE, marker="o", label=r"$E_x=E_y$ (in-plane)")
    a.plot(Ez, z, "--", color=VERM, marker="D", label=r"$E_z$ (vertical)")
    a.set_xlabel(r"Effective Young's modulus (GPa)")
    depth_axis(a)
    a.set_title("(a)")
    a.legend(loc="upper left", framealpha=0.9)

    b.plot(brine, z, "-", color=ORANGE, marker="o", label="brine pockets")
    b.plot(chan, z, "-", color=BLUE, marker="s", label="brine channels")
    b.plot(gas, z, ":", color=BLACK, marker="^", label="gas voids")
    b.set_xlabel("Phase volume fraction (%)")
    depth_axis(b)
    orient_labels(b)
    b.set_title("(b)")
    b.legend(loc="upper right", framealpha=0.9)

    fig.tight_layout()
    save(fig, "ice_column_profiles")


# ===========================================================================
# 3) study_coltensor : full-tensor moduli + anisotropy ratios vs z/H
# ===========================================================================
def fig_coltensor():
    r = load("results_coltensor.csv")
    z = np.array([float(x["z"]) for x in r])
    Ex = np.array([float(x["E_x"]) for x in r]) / 1e9
    Ey = np.array([float(x["E_y"]) for x in r]) / 1e9
    Ez = np.array([float(x["E_z"]) for x in r]) / 1e9
    Er = np.array([float(x["E_ratio"]) for x in r])
    Gr = np.array([float(x["G_ratio"]) for x in r])

    fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 6.2), sharey=True)

    a.plot(Ex, z, "-", color=BLUE, marker="o", label=r"$E_x$")
    a.plot(Ey, z, "--", color=GREEN, marker="s", label=r"$E_y$")
    a.plot(Ez, z, "-.", color=VERM, marker="D", label=r"$E_z$")
    a.set_xlabel(r"Directional Young's modulus (GPa)")
    depth_axis(a)
    a.set_title("(a)")
    a.legend(loc="upper left", framealpha=0.9)

    b.axvline(1.0, color="0.5", lw=1.2, ls=":")
    b.plot(Er, z, "-", color=VERM, marker="D", label=r"$E_z/E_{xy}$")
    b.plot(Gr, z, "-", color=BLUE, marker="o", label=r"$G_{ax}/G_{xy}$")
    b.set_xlabel("Anisotropy ratio")
    depth_axis(b)
    orient_labels(b)
    b.set_title("(b)")
    b.legend(loc="upper right", framealpha=0.9)

    fig.tight_layout()
    save(fig, "study_coltensor")


if __name__ == "__main__":
    fig_scfdepth()
    fig_column()
    fig_coltensor()
    print("done")
