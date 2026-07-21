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

# ---- Okabe-Ito CVD-safe palette -------------------------------------------
BLUE, ORANGE, GREEN = "#0072B2", "#E69F00", "#009E73"
VERM, PURPLE, SKY, BLACK = "#D55E00", "#CC79A7", "#56B4E9", "#222222"

plt.rcParams.update({
    "font.size": 14, "axes.titlesize": 15, "axes.labelsize": 15,
    "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 12.5,
    "axes.linewidth": 1.0, "lines.linewidth": 2.2, "lines.markersize": 7,
    "figure.dpi": 120, "savefig.bbox": "tight", "axes.grid": True,
    "grid.alpha": 0.30, "grid.linewidth": 0.7,
})

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


def depth_axis(ax):
    ax.set_ylim(1.0, 0.0)                 # 0 (surface) at top, 1 (base) at bottom
    ax.set_ylabel(r"Normalized depth $z/H$")
    ax.axhspan(0.80, 1.0, color=SKY, alpha=0.10, zorder=0)  # warm bottom band


def orient_labels(ax):
    """Cold-surface / warm-bottom orientation, as ticks outside the right edge.

    Kept out of the plot area: as floating text it sat at bottom centre, where
    both the curves and the legend want to be.
    """
    tw = ax.twinx()
    tw.set_ylim(1.0, 0.0)
    tw.set_yticks([0.0, 1.0])
    tw.set_yticklabels(["surface (cold)", "bottom (warm)"], fontsize=11,
                       color="0.35")
    tw.tick_params(axis="y", length=0)
    tw.grid(False)
    for spine in tw.spines.values():
        spine.set_visible(False)


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
