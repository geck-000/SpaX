"""Figures for case study 3: the Novik Bay cantilever beams of Gogolaze et al. (2026).

Two figures.

  cantilever_knockdown.png   The constitutive step, which is the whole content of
                             the case study: the RVE knockdown E/E_matrix against
                             total soft fraction, split into the channel-free and
                             channel-bearing branches, with the three empirical
                             E(v_b) laws overlaid on the same axis. This is where
                             the morphology split at phi ~ 0.24 is visible, and
                             where the beam's own depth range is marked.

  cantilever_beam.png        The consequence: E(z) through the beam, the root
                             normal-stress distribution, and the deflection
                             profile, for each constitutive law against the
                             homogeneous reference.

Run:  python analysis/plot_cantilever.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle
from figstyle import BLACK, BLUE, GREEN, ORANGE, PURPLE, VERM
import matplotlib.pyplot as plt

from case_study_3_scoping import (B, E_EFF, E_HOMOG, F, GAS, H, L, MATRIX_FACTOR,
                                  E_MATRIX_0, brine_profile, composite_beam,
                                  consistent_subset, karulina, knockdown_curve,
                                  knockdown_exponent, load_database, rve_profile,
                                  vaudrey, weeks)

figstyle.apply()
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "results")

LAWS = [("Weeks & Assur", weeks, ORANGE),
        ("Vaudrey", vaudrey, PURPLE),
        ("Karulina", karulina, GREEN)]


def fig_knockdown(db):
    """The constitutive comparison, on the axis the RVE actually computes."""
    free = knockdown_curve(db, percolated=False)
    perc = knockdown_curve(db, percolated=True)

    fig, ax = plt.subplots(figsize=(8.2, 6.0))

    # beam's own depth range, so the reader sees where the law is being used
    z = np.linspace(0, H, 400)
    phi_beam = brine_profile(z) + GAS
    ax.axvspan(phi_beam.min(), phi_beam.max(), color=figstyle.SKY, alpha=0.13,
               zorder=0, label="range spanned by the beam")

    # Three populations, kept apart: the brine partition is a second axis, and
    # pooling it is what produced the earlier non-monotone E(z).
    keep = consistent_subset(db)
    other = db[db.channels & ~db.index.isin(keep.index)]

    def binned(sub):
        g = sub.groupby(sub.phi.round(3)).agg(
            n=("r", "size"), r=("r", "mean"), sd=("r", lambda s: s.std(ddof=0)))
        return g.sort_index()

    for sub, colour, marker, lbl in [
            (db[~db.channels], BLUE, "o", "channel-free"),
            (keep, VERM, "s", r"percolated, $f_{\mathrm{ch}}\approx0.4$ (used)"),
            (other, "0.55", "D", r"percolated, other $f_{\mathrm{ch}}$")]:
        g = binned(sub)
        w, t = g[g.n >= 3], g[g.n < 3]
        ax.errorbar(w.index, w.r, yerr=w.sd, fmt=marker, color=colour,
                    ms=6, lw=0, elinewidth=1.4, capsize=3, label=lbl, zorder=3)
        ax.plot(t.index, t.r, marker, color=colour, ms=5, mfc="white",
                mew=1.3, lw=0, zorder=3)

    # the fitted law actually used for the beam, over the range it is fitted on
    k, rms, _, p_lo, p_hi = knockdown_exponent(db)
    pf = np.linspace(p_lo, p_hi, 200)
    ax.plot(pf, (1 - pf) ** k, "-", color=BLACK, lw=2.6, zorder=5,
            label=r"fit $(1-\phi)^{%.2f}$ (%.1f%% RMS)" % (k, rms))

    # empirical laws, expressed as E/E_0 so they share the axis
    phi = np.linspace(0.02, 0.55, 300)
    vb = np.clip(phi - GAS, 1e-6, None)
    for lbl, law, colour in LAWS:
        r = law(vb) / law(np.array([1e-6]))[0]
        r = np.where(r > 0, r, np.nan)
        ax.plot(phi, r, "--", color=colour, lw=1.9, label=lbl + " (empirical)")

    ax.set_xlabel(r"Total soft fraction $\phi = \phi_{\mathrm{brine}} + \phi_{\mathrm{gas}}$")
    ax.set_ylabel(r"$E/E_{\mathrm{matrix}}$")
    ax.set_xlim(0.0, 0.55)
    ax.set_ylim(0.0, 1.02)
    ax.legend(loc="upper right", framealpha=0.93)
    ax.set_title("Micromechanical knockdown against the empirical $E(v_b)$ laws",
                 pad=10)
    ax.text(0.245, 0.055, "open symbols:\nfewer than 3 packings", fontsize=10.5,
            color="0.35", ha="left")
    fig.tight_layout()
    p = os.path.join(OUT, "cantilever_knockdown.png")
    fig.savefig(p)
    plt.close(fig)
    print("wrote", p)


def fig_beam(db):
    """E(z), root stress and deflection, per constitutive law."""
    k, _, _, _, _ = knockdown_exponent(db)
    z = np.linspace(0, H, 2001)
    zh = z / H
    vb = brine_profile(z)
    phi = vb + GAS

    profiles = [(lbl, law(vb), colour) for lbl, law, colour in LAWS]
    profiles.append(("SpaX-RVE", rve_profile(phi, k, MATRIX_FACTOR * E_MATRIX_0),
                     VERM))

    I_h = B * H ** 3 / 12
    s_h = F * L / I_h * (H / 2)

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 6.0))

    # ---- (a) modulus profile ----
    ax = axes[0]
    for lbl, E, colour in profiles:
        Ep = np.where(E > 0, E, np.nan)
        ax.plot(Ep / 1e9, zh, color=colour, label=lbl)
    ax.axvline(E_HOMOG / 1e9, color=BLACK, ls=":", lw=1.8, label="homogeneous (measured)")
    ax.axvline(E_EFF / 1e9, color=BLACK, ls="-.", lw=1.5, label=r"$E_{\mathrm{eff}}$ (root-corr.)")
    figstyle.depth_axis(ax)
    ax.set_xlabel(r"$E$ (GPa)")
    ax.set_xlim(-0.4, 4.0)
    ax.set_title("(a) modulus profile")
    ax.legend(loc="lower right", fontsize=10.5, framealpha=0.93)

    # ---- (b) root-section normal stress ----
    ax = axes[1]
    for lbl, E, colour in profiles:
        if np.any(E <= 0):
            continue                       # Vaudrey: construction fails
        _, z_n, _, sig, _ = composite_beam(E, z)
        ax.plot(sig / 1e3, zh, color=colour, label=lbl)
        # the neutral axis: where sigma crosses zero
        ax.plot([0.0], [z_n / H], "o", color=colour, ms=8,
                mfc="white", mew=2.0, zorder=6)
    sig_h = F * L / I_h * (H / 2 - z)
    ax.plot(sig_h / 1e3, zh, color=BLACK, ls=":", lw=1.9, label="homogeneous")
    ax.plot([0.0], [0.5], "o", color=BLACK, ms=8, mfc="white", mew=2.0, zorder=6)
    ax.axvline(0, color="0.5", lw=1.0)
    figstyle.depth_axis(ax, label=False)
    ax.set_xlabel(r"$\sigma$ at the root, $x=0$ (kPa)")
    ax.set_title("(b) normal stress; circles mark the neutral axis")
    ax.legend(loc="upper left", fontsize=10.5, framealpha=0.93)
    figstyle.orient_labels(ax)

    # ---- (c) deflection ----
    ax = axes[2]
    x = np.linspace(0, L, 200)
    for lbl, E, colour in profiles:
        if np.any(E <= 0):
            continue
        E_int, _, I_eff, _, _ = composite_beam(E, z)
        w = F * x ** 2 * (3 * L - x) / (6 * E_int * I_eff)
        ax.plot(x, w * 1e3, color=colour, label=lbl)
    w_h = F * x ** 2 * (3 * L - x) / (6 * E_HOMOG * I_h)
    ax.plot(x, w_h * 1e3, color=BLACK, ls=":", lw=1.9, label="homogeneous")
    ax.plot([L], [4.151], "*", color=BLACK, ms=17, zorder=6,
            label="measured tip (Beam 3)")
    ax.invert_yaxis()
    ax.set_xlabel(r"$x$ along the beam (m)")
    ax.set_ylabel(r"deflection $w$ (mm)")
    ax.set_title("(c) deflection: the level test")
    ax.legend(loc="lower left", fontsize=10.5, framealpha=0.93)

    fig.suptitle("Beam 3 of Gogolaze et al. (2026) under four constitutive laws "
                 "(Vaudrey omitted where $E<0$)", y=1.00, fontsize=15)
    fig.tight_layout()
    p = os.path.join(OUT, "cantilever_beam.png")
    fig.savefig(p)
    plt.close(fig)
    print("wrote", p)


if __name__ == "__main__":
    db = load_database()
    fig_knockdown(db)
    fig_beam(db)
