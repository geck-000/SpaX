"""Scoping calculation for case study 3: RVE-derived E(z) in a cantilever ice beam.

Reproduces the composite-beam construction of Gogolaze et al. (2026) and drives it
with four E(v_b) laws, three empirical and one from the SpaX RVE database.

Two things are established here, and they set the design of the study:

  1. The neutral-axis position and the whole normal-stress field are invariant to
     any constant rescaling of E(z).  The matrix calibration factor of
     Section sec:field therefore cancels identically from the stress comparison,
     which is what makes that comparison a genuine prediction.
  2. The tip deflection depends on the level of E(z) alone.  Deflection and
     stress are thus independent tests of the two halves of the model.

IMPORTANT (corrected 2026-08-06).  The brine volume in these beams never falls
below 64 per-mille, so *every* depth is above the percolation threshold and the
vertically connected channel network of Section sec:aniso is present throughout.
The knockdown must therefore be taken from the channel-BEARING population.  An
earlier version of this script pooled the channel-free packings and consequently
overestimated E(z), most badly at the base.

Run:  python analysis/case_study_3_scoping.py
"""

import glob
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RESULT_DIRS = [
    os.path.join(HERE, os.pardir, "results"),
    os.path.join(HERE, os.pardir, ".claude", "worktrees",
                 "skeletal-eringen-weibull", "results"),
]

# --- Beam 3 of Gogolaze et al. (2026), Tables 2 and 5 ------------------------
H, B, L, F = 0.32, 0.60, 1.96, 2002.0     # m, m, m (load point), N
E_HOMOG = 0.785e9                          # apparent modulus, their eq. (2)
E_EFF = 1.421e9                            # root-corrected, their eq. (19)
GAS = 0.02                                 # gas fraction carried by the RVE cells
MATRIX_FACTOR = 0.49                       # Section sec:field
E_MATRIX_0 = 9.36e9                        # defect-free matrix of the base cells

VB_POLY = (0.29315, -5.124, 85.977)        # their eq. (14), z in cm, v_b per-mille


def brine_profile(z_m):
    a, b, c = VB_POLY
    zc = z_m * 100.0
    return (a * zc**2 + b * zc + c) / 1000.0


# --- The empirical E(v_b) laws ----------------------------------------------
def weeks(v):
    """Weeks & Assur (1967), their eq. (3), with E_0 = 9.5 GPa."""
    return 9.5e9 * (1 - np.sqrt(v)) ** 4


def vaudrey(v):
    """Vaudrey (1977), their eq. (4).  Goes negative at high v_b."""
    return (5.31 - 0.436 * np.sqrt(v * 1000)) * 1e9


def karulina(v):
    """Karulina et al. (2019), their eq. (5)."""
    return 3.1031e9 * np.exp(-3.385 * np.sqrt(v))


# --- RVE database ------------------------------------------------------------
def load_database():
    """Pool every homogenisation run.

    Two things matter here and both were got wrong at first.  ``VoF_sphere``
    counts only the meshed pockets and EXCLUDES the channel network, so the
    physically meaningful abscissa is phi = VoF_sphere + channel_vof_target.
    And ``channels_frac``, the share of brine carried by the channels rather
    than the pockets, is a second independent axis: at phi = 0.24 a
    channel-dominated cell is half again stiffer than a pocket-dominated one.
    """
    rows = []
    for d in RESULT_DIRS:
        for p in glob.glob(os.path.join(d, "results_*.csv")):
            try:
                t = pd.read_csv(p)
            except Exception:
                continue
            if not {"E_matrix", "VoF_sphere", "E_x"} <= set(t.columns):
                continue
            chan = t.get("generate_channels", pd.Series(["No"] * len(t)))
            for c in ("channel_vof_target", "VoF_incl_sphere"):
                if c not in t:
                    t[c] = 0.0
            t = t.assign(channels=chan.astype(str).str.lower().isin(["yes", "true", "1"]),
                         source=os.path.basename(p))
            rows.append(t[["run_id", "source", "channels", "VoF_sphere",
                           "VoF_incl_sphere", "channel_vof_target", "E_matrix", "E_x"]])
    a = pd.concat(rows).dropna(subset=["E_x", "VoF_sphere"])
    a = a.drop_duplicates(subset=["run_id", "VoF_sphere", "E_x"])
    cvt = pd.to_numeric(a.channel_vof_target, errors="coerce").fillna(0.0)
    incl = pd.to_numeric(a.VoF_incl_sphere, errors="coerce").fillna(0.0)
    return a.assign(r=a.E_x / a.E_matrix,
                    phi=a.VoF_sphere + cvt,
                    channels_frac=(cvt / (cvt + incl).replace(0, np.nan)))


# Baseline brine partition: the production column's value above percolation, and
# the one studies/make_cantilever.py builds the new decks at.
CF_LO, CF_HI = 0.35, 0.45


def consistent_subset(db):
    """Channel-bearing cells at the baseline brine partition."""
    return db[db.channels & db.channels_frac.between(CF_LO, CF_HI)]


def knockdown_curve(db, percolated=True, phi_min=0.05):
    """Binned mean knockdown r(phi) for the requested population."""
    sub = db[db.channels & (db.phi >= phi_min)] if percolated else db[~db.channels]
    g = sub.groupby(sub.phi.round(3)).agg(
        n=("r", "size"), r=("r", "mean"), sd=("r", lambda s: s.std(ddof=0)))
    return g.sort_index()


# The beam spans phi = 0.084 to 0.242; the law is fitted over that range and a
# little beyond, and is not used outside it.  Fitting the full 0.06-0.52 span
# instead lets the high-phi plateau drag the exponent and leaves the fit 19%
# stiff at the base, which is where the bending lever arm is longest.
PHI_FIT = (0.060, 0.260)


def knockdown_exponent(db, phi_range=PHI_FIT):
    """Fit r = (1-phi)^k over the consistent-partition subset, in-range only.

    One parameter, one morphology.  Pooling every channel fraction instead
    scatters at 23% RMS and produces a non-monotone E(z).
    """
    sub = consistent_subset(db)
    sub = sub[sub.phi.between(*phi_range)]
    x, y = np.log(1 - sub.phi), np.log(sub.r)
    k = float((x * y).sum() / (x * x).sum())
    rms = float(np.sqrt(np.mean((((1 - sub.phi) ** k - sub.r) / sub.r) ** 2)) * 100)
    return k, rms, len(sub), float(sub.phi.min()), float(sub.phi.max())


def rve_profile(phi, k, E_matrix):
    """The fitted knockdown evaluated on the beam profile."""
    return E_matrix * (1.0 - phi) ** k


# --- Composite-beam construction, their eqs. (7)-(13) ------------------------
def composite_beam(E, z):
    """Return (E_int, z_n, I_eff, sigma(z), w_tip) for a modulus profile E(z)."""
    E_int = np.trapz(E, z) / H                              # eq. (7), per unit h
    z_n = H - np.trapz(E * (H - z), z) / np.trapz(E, z)     # eq. (8)
    y = z - z_n
    I_eff = B * np.trapz(E / E_int * y**2, z)               # eq. (9)
    sigma = F * L / I_eff * (E / E_int) * (-y)              # eq. (12), tension +ve
    w_tip = F * L**3 / (3 * E_int * I_eff)                  # eq. (10), Euler-Bernoulli
    return E_int, z_n, I_eff, sigma, w_tip


def main():
    db = load_database()
    free = db[~db.channels]
    k, rms, n_sub, p_lo, p_hi = knockdown_exponent(db)
    print(f"database: {len(db)} runs, {len(free)} channel-free "
          f"(phi <= {free.phi.max():.3f}), "
          f"{int(db.channels.sum())} channel-bearing "
          f"(phi <= {db[db.channels].phi.max():.3f})")
    print(f"  consistent-partition subset (channels_frac {CF_LO}-{CF_HI}): "
          f"n={n_sub}, phi {p_lo:.3f}-{p_hi:.3f}")
    print(f"  fitted knockdown (1-phi)^k: k = {k:.3f}, RMS {rms:.2f}%\n")

    z = np.linspace(0, H, 3201)
    vb = brine_profile(z)
    phi = vb + GAS
    print(f"beam: v_b top {vb[0]:.3f}, min {vb.min():.3f}, base {vb[-1]:.3f}")
    print(f"      every depth is percolated (min v_b {vb.min():.3f} >> ~0.05)\n")

    E_matrix = MATRIX_FACTOR * E_MATRIX_0
    laws = [
        ("Weeks & Assur", weeks(vb)),
        ("Vaudrey", vaudrey(vb)),
        ("Karulina", karulina(vb)),
        ("SpaX-RVE (percolated)", rve_profile(phi, k, E_matrix)),
    ]

    I_h = B * H**3 / 12
    s_h = F * L / I_h * (H / 2)
    hdr = (f"{'model':<24}{'E_int GPa':>10}{'z_n cm':>9}{'sig_top kPa':>13}"
           f"{'d vs homog':>12}{'sig_bot kPa':>13}{'w_tip mm':>10}")
    print(hdr)
    print("-" * len(hdr))
    for name, E in laws:
        E_int, z_n, _, sig, w = composite_beam(E, z)
        if np.any(E <= 0):
            print(f"{name:<24}{E_int/1e9:>10.3f}{z_n*100:>9.2f}"
                  f"{'--':>13}{'--':>12}{'--':>13}{'--':>10}  (E<0 below z~27cm)")
            continue
        print(f"{name:<24}{E_int/1e9:>10.3f}{z_n*100:>9.2f}{sig[0]/1e3:>13.1f}"
              f"{(sig[0]/s_h - 1)*100:>11.1f}%{sig[-1]/1e3:>13.1f}{w*1000:>10.3f}")
    print(f"{'Homogeneous (measured)':<24}{E_HOMOG/1e9:>10.3f}{H/2*100:>9.2f}"
          f"{s_h/1e3:>13.1f}{0.0:>11.1f}%{-s_h/1e3:>13.1f}{4.15:>10.3f}")
    print(f"{'their E_eff (root-corr)':<24}{E_EFF/1e9:>10.3f}")

    print("\n-- Invariance of stress and neutral axis to the matrix factor --")
    for f_ in [1.00, 0.49, 0.42]:
        E = rve_profile(phi, k, f_ * E_MATRIX_0)
        E_int, z_n, _, sig, w = composite_beam(E, z)
        print(f"  factor {f_:.2f} -> E_int {E_int/1e9:6.3f} GPa, z_n {z_n*100:6.3f} cm, "
              f"sig_top {sig[0]/1e3:7.2f} kPa, w_tip {w*1000:6.3f} mm")

    print("\n-- Sensitivity to the brine partition (what rve_gogo_chanfrac.csv settles) --")
    sub = consistent_subset(db)
    for lbl, kk in [("channel-free branch",
                     float((np.log(1 - free.phi) * np.log(free.r)).sum()
                           / (np.log(1 - free.phi) ** 2).sum())),
                    (f"consistent partition (k={k:.2f})", k)]:
        E = rve_profile(phi, kk, E_matrix)
        E_int, z_n, _, sig, _ = composite_beam(E, z)
        print(f"  {lbl:<34} E_int {E_int/1e9:5.2f} GPa, z_n {z_n*100:5.2f} cm, "
              f"sig_top {sig[0]/1e3:6.1f} kPa ({(sig[0]/s_h-1)*100:+5.1f}%)")
    print(f"  subset spans channels_frac {sub.channels_frac.min():.2f}"
          f"-{sub.channels_frac.max():.2f}; the sweep will widen this to 0.20-0.85")


if __name__ == "__main__":
    main()
