"""Offline analysis of the four sea-ice column studies.

Pure pandas / numpy / matplotlib -- NO Abaqus needed. Run on any laptop:

    python analyze_studies.py

Reads whichever of these exist and makes one PNG each:
    results_morph.csv -> study_morphology.png   (E_z/E_x over sphericity x growth)
    results_perc.csv  -> study_percolation.png  (anisotropy onset vs phi_b)
    results_mono.csv  -> study_monotonic.png    (monotonic E(z) vs Marchenko 2024)
    results_seas.csv  -> study_seasonal.png     (E(z) profiles, soft-layer thickening)

Each results CSV is produced on the Abaqus machine by SpaX_PostProcess.py and
carries at least: run_id, E_x, E_z (E_eff == E_x for these utx+utz runs).
"""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
GPa = 1e9

def num(s):
    return pd.to_numeric(s, errors="coerce")

def load(path):
    if not os.path.exists(path):
        print(f"[skip] {path} not found"); return None
    df = pd.read_csv(path)
    for c in ("E_x", "E_y", "E_z", "E_eff"):
        if c in df.columns:
            df[c] = num(df[c])
    return df

# ---- Marchenko (2024) vibrating-beam fit ----------------------------------
def marchenko(depth_from_top):
    Ebot, M, n = 1.67, 2.63, 0.5
    zeta = 1.0 - depth_from_top          # normalized height above base
    return Ebot * ((M - 1) * zeta**n + 1.0)

# ===========================================================================
def do_morphology():
    df = load("results_morph.csv")
    if df is None: return
    sph = df.run_id.str.extract(r"s(\d+)g")[0].astype(float) / 100
    grw = df.run_id.str.extract(r"g(\d+)")[0].astype(float) / 100
    Ex, Ez = df.E_x.values / GPa, df.E_z.values / GPa
    aniso = Ez / Ex
    df2 = pd.DataFrame({"sph": sph, "grw": grw, "Ex": Ex, "Ez": Ez, "aniso": aniso})
    P = df2.pivot_table(index="grw", columns="sph", values="aniso")
    E = df2.pivot_table(index="grw", columns="sph", values="Ex")

    fig, ax = plt.subplots(1, 2, figsize=(13, 5.4))
    for a, M, ttl, cb in ((ax[0], P, "(a) vertical anisotropy  E_z / E_x", "E_z/E_x"),
                          (ax[1], E, "(b) horizontal modulus  E_x (GPa)", "E_x (GPa)")):
        im = a.imshow(M.values, origin="lower", aspect="auto", cmap="viridis",
                      extent=[M.columns.min()-.05, M.columns.max()+.05,
                              M.index.min()-.15, M.index.max()+.15])
        a.set_xticks(M.columns); a.set_yticks(M.index)
        a.set_xlabel("sphericity (elongated <- -> round)")
        a.set_ylabel("growth concentration (orientation)")
        a.set_title(ttl)
        for iy, gy in enumerate(M.index):
            for ix, sx in enumerate(M.columns):
                v = M.values[iy, ix]
                if np.isfinite(v):
                    a.text(sx, gy, f"{v:.2f}", ha="center", va="center",
                           color="white" if v < np.nanmean(M.values) else "black",
                           fontsize=9)
        fig.colorbar(im, ax=a, fraction=0.046, label=cb)
    fig.tight_layout(); fig.savefig("study_morphology.png", dpi=160)
    print("wrote study_morphology.png")

# ===========================================================================
def do_percolation():
    df = load("results_perc.csv")
    if df is None: return
    phi = df.run_id.str.extract(r"p(\d+)(?:off|on)")[0].astype(float) / 1000
    chan = df.run_id.str.endswith("on")
    Ex, Ez = df.E_x.values / GPa, df.E_z.values / GPa
    d = pd.DataFrame({"phi": phi, "chan": chan.values,
                      "Ex": Ex, "Ez": Ez, "aniso": Ez / Ex}).sort_values("phi")
    off, on = d[~d.chan], d[d.chan]

    fig, ax = plt.subplots(1, 2, figsize=(13, 5.4))
    a = ax[0]
    a.plot(off.phi, off.Ex, "o-", color="tab:blue", label="E_x  pockets only")
    a.plot(off.phi, off.Ez, "^--", color="tab:cyan", label="E_z  pockets only")
    if len(on):
        a.plot(on.phi, on.Ex, "o-", color="tab:red", label="E_x  + channels")
        a.plot(on.phi, on.Ez, "^--", color="tab:orange", label="E_z  + channels")
    a.axvline(0.05, color="0.5", ls=":"); a.text(0.051, a.get_ylim()[1]*0.95,
              "rule of fives", fontsize=8, color="0.4")
    a.set_xlabel("brine volume fraction  phi_b"); a.set_ylabel("modulus (GPa)")
    a.set_title("(a) E_x, E_z vs brine fraction"); a.grid(alpha=0.3); a.legend(fontsize=8)

    a = ax[1]
    a.plot(off.phi, off.aniso, "o-", color="tab:blue", label="pockets only")
    if len(on):
        a.plot(on.phi, on.aniso, "o-", color="tab:red", label="+ vertical channels")
    a.axhline(1.0, color="0.6"); a.axvline(0.05, color="0.5", ls=":")
    a.set_xlabel("brine volume fraction  phi_b"); a.set_ylabel("E_z / E_x")
    a.set_title("(b) anisotropy onset at percolation"); a.grid(alpha=0.3); a.legend(fontsize=8)
    fig.tight_layout(); fig.savefig("study_percolation.png", dpi=160)
    print("wrote study_percolation.png")

# ===========================================================================
def do_monotonic():
    df = load("results_mono.csv")
    if df is None: return
    d = df.run_id.str.extract(r"z(\d+)")[0].astype(float).values / 100
    o = np.argsort(d); d = d[o]
    Eeff = (df.E_eff.values if "E_eff" in df else df.E_x.values)[o] / GPa
    Ez = df.E_z.values[o] / GPa
    Em = marchenko(d)
    zc = np.linspace(0.001, 1, 200); Emc = marchenko(1 - zc); dc = 1 - zc

    cmp = None
    if os.path.exists("results_marchenko.csv"):       # the old C-shape calibration
        c = load("results_marchenko.csv")
        dc2 = c.run_id.str.extract(r"z(\d+)")[0].astype(float).values / 100
        oc = np.argsort(dc2)
        cmp = (dc2[oc], (c.E_eff.values if "E_eff" in c else c.E_x.values)[oc] / GPa)

    rms = 100*np.sqrt(np.mean(((Eeff - Em)/Em)**2))
    fig, ax = plt.subplots(1, 2, figsize=(13, 6))
    a = ax[0]
    a.plot(Emc, dc, "-", color="tab:red", lw=2, label="Marchenko 2024 fit")
    a.plot(Em, d, "s", color="tab:red", ms=5, label="Marchenko @ slices")
    a.plot(Eeff, d, "o-", color="tab:green", label="SPAX monotonic-salinity E_eff")
    a.plot(Ez, d, "^--", color="tab:olive", lw=1, label="SPAX E_z")
    if cmp is not None:
        a.plot(cmp[1], cmp[0], "o:", color="0.6", lw=1, label="SPAX C-shape (prev.)")
    a.set_ylim(1, 0); a.set_xlim(0, 6)
    a.set_xlabel("Young's modulus (GPa)"); a.set_ylabel("depth z/H (0=top,1=base)")
    a.set_title("(a) Monotonic-salinity column vs Marchenko"); a.grid(alpha=0.3)
    a.legend(fontsize=8, loc="lower right")
    a = ax[1]
    a.plot(Eeff/Em, d, "o-", color="tab:green"); a.axvline(1, color="0.6")
    a.set_ylim(1, 0); a.set_xlim(0.6, 1.6)
    a.set_xlabel("E_SPAX / E_Marchenko"); a.set_ylabel("depth z/H")
    a.set_title(f"(b) match ratio (RMS {rms:.0f}%)"); a.grid(alpha=0.3)
    for i in range(len(d)):
        a.annotate(f"{Eeff[i]/Em[i]:.2f}", (Eeff[i]/Em[i], d[i]), fontsize=7,
                   xytext=(3, 3), textcoords="offset points")
    fig.tight_layout(); fig.savefig("study_monotonic.png", dpi=160)
    print(f"wrote study_monotonic.png  (RMS {rms:.0f}%)")

# ===========================================================================
def do_seasonal():
    df = load("results_seas.csv")
    if df is None: return
    scen = df.run_id.str.extract(r"SEAS_(w\d+)_")[0]
    d = df.run_id.str.extract(r"z(\d+)")[0].astype(float).values / 100
    Ex = df.E_x.values / GPa
    dd = pd.DataFrame({"scen": scen.values, "d": d, "Ex": Ex})
    colors = {"w20": "tab:blue", "w12": "tab:orange", "w06": "tab:red"}
    labels = {"w20": "T_top = -20 degC (mid-winter)",
              "w12": "T_top = -12 degC",
              "w06": "T_top = -6 degC (spring)"}
    fig, ax = plt.subplots(1, 2, figsize=(13, 6))
    a = ax[0]
    for s in ("w20", "w12", "w06"):
        g = dd[dd.scen == s].sort_values("d")
        if len(g):
            a.plot(g.Ex, g.d, "o-", color=colors[s], label=labels[s])
    a.set_ylim(1, 0); a.set_xlabel("Young's modulus E_x (GPa)")
    a.set_ylabel("depth z/H (0=top,1=base)")
    a.set_title("(a) column stiffness profile vs surface temperature")
    a.grid(alpha=0.3); a.legend(fontsize=8)

    # soft-layer thickness: fraction of column below a stiffness threshold
    a = ax[1]
    thr = 7.0
    bars = []
    for s in ("w20", "w12", "w06"):
        g = dd[dd.scen == s].sort_values("d")
        if len(g):
            frac = (g.Ex < thr).mean()
            bars.append((labels[s], frac))
    if bars:
        a.barh([b[0] for b in bars], [b[1] for b in bars],
               color=[colors[s] for s in ("w20", "w12", "w06")][:len(bars)])
        for i, b in enumerate(bars):
            a.text(b[1]+0.01, i, f"{b[1]*100:.0f}%", va="center", fontsize=9)
    a.set_xlim(0, 1); a.set_xlabel(f"fraction of column with E_x < {thr:.0f} GPa")
    a.set_title("(b) warm soft-layer thickness grows with warming")
    a.grid(alpha=0.3, axis="x")
    fig.tight_layout(); fig.savefig("study_seasonal.png", dpi=160)
    print("wrote study_seasonal.png")

# ===========================================================================
# Second battery
# ---------------------------------------------------------------------------
def do_channel():
    df = load("results_channel.csv")
    if df is None: return
    rch = df.run_id.str.extract(r"r(\d+)f")[0].astype(float) / 1000
    frc = df.run_id.str.extract(r"f(\d+)")[0].astype(float) / 100
    Ex, Ez = df.E_x.values / GPa, df.E_z.values / GPa
    d = pd.DataFrame({"rch": rch, "frc": frc, "Ex": Ex, "aniso": Ez / Ex})
    A = d.pivot_table(index="frc", columns="rch", values="aniso")
    E = d.pivot_table(index="frc", columns="rch", values="Ex")
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.4))
    for a, M, ttl, cb in ((ax[0], A, "(a) vertical anisotropy E_z/E_x", "E_z/E_x"),
                          (ax[1], E, "(b) horizontal modulus E_x (GPa)", "E_x (GPa)")):
        im = a.imshow(M.values, origin="lower", aspect="auto", cmap="viridis",
                      extent=[M.columns.min()-.003, M.columns.max()+.003,
                              M.index.min()-.1, M.index.max()+.1])
        a.set_xticks(M.columns); a.set_yticks(M.index)
        a.set_xlabel("channel radius (model units)")
        a.set_ylabel("fraction of brine in channels")
        a.set_title(ttl)
        for iy, gy in enumerate(M.index):
            for ix, sx in enumerate(M.columns):
                v = M.values[iy, ix]
                if np.isfinite(v):
                    a.text(sx, gy, f"{v:.3f}" if M is A else f"{v:.2f}",
                           ha="center", va="center", fontsize=9,
                           color="white" if v < np.nanmean(M.values) else "black")
        fig.colorbar(im, ax=a, fraction=0.046, label=cb)
    fig.tight_layout(); fig.savefig("study_channel.png", dpi=160)
    print("wrote study_channel.png")

def do_fymy():
    df = load("results_fymy.csv")
    if df is None: return
    typ = df.run_id.str.extract(r"FYMY_(fy|my)_")[0]
    d = df.run_id.str.extract(r"z(\d+)")[0].astype(float).values / 100
    Ex = df.E_x.values / GPa
    dd = pd.DataFrame({"typ": typ.values, "d": d, "Ex": Ex})
    fig, a = plt.subplots(figsize=(7, 6))
    for t, c, lab in (("fy", "tab:blue", "first-year (C-shape, brine-rich)"),
                      ("my", "tab:red", "multi-year (desalinated, gassy)")):
        g = dd[dd.typ == t].sort_values("d")
        if len(g): a.plot(g.Ex, g.d, "o-", color=c, label=lab)
    a.set_ylim(1, 0); a.set_xlabel("Young's modulus E_x (GPa)")
    a.set_ylabel("depth z/H (0=top,1=base)")
    a.set_title("First-year vs multi-year ice column stiffness")
    a.grid(alpha=0.3); a.legend(fontsize=9)
    fig.tight_layout(); fig.savefig("study_fymy.png", dpi=160)
    print("wrote study_fymy.png")

def do_mono2():
    df = load("results_mono2.csv")
    if df is None: return
    d = df.run_id.str.extract(r"z(\d+)")[0].astype(float).values / 100
    o = np.argsort(d); d = d[o]
    Eeff = (df.E_eff.values if "E_eff" in df else df.E_x.values)[o] / GPa
    Em = marchenko(d)
    zc = np.linspace(0.001, 1, 200); Emc = marchenko(1 - zc); dc = 1 - zc
    prev = None
    p = load("results_mono.csv")
    if p is not None:
        dp = p.run_id.str.extract(r"z(\d+)")[0].astype(float).values / 100
        op = np.argsort(dp)
        prev = (dp[op], (p.E_eff.values if "E_eff" in p else p.E_x.values)[op] / GPa)
    rms = 100*np.sqrt(np.mean(((Eeff - Em)/Em)**2))
    fig, ax = plt.subplots(1, 2, figsize=(13, 6))
    a = ax[0]
    a.plot(Emc, dc, "-", color="tab:red", lw=2, label="Marchenko 2024 fit")
    a.plot(Eeff, d, "o-", color="tab:green", label="SPAX mono2 (steeper S)")
    if prev is not None:
        a.plot(prev[1], prev[0], "o:", color="0.6", label="SPAX mono1 (gentle S)")
    a.set_ylim(1, 0); a.set_xlim(0, 6)
    a.set_xlabel("Young's modulus (GPa)"); a.set_ylabel("depth z/H")
    a.set_title("(a) steeper-salinity column vs Marchenko"); a.grid(alpha=0.3)
    a.legend(fontsize=8, loc="lower right")
    a = ax[1]
    a.plot(Eeff/Em, d, "o-", color="tab:green"); a.axvline(1, color="0.6")
    a.set_ylim(1, 0); a.set_xlim(0.6, 1.6)
    a.set_xlabel("E_SPAX / E_Marchenko"); a.set_ylabel("depth z/H")
    a.set_title(f"(b) match ratio (RMS {rms:.0f}%)"); a.grid(alpha=0.3)
    for i in range(len(d)):
        a.annotate(f"{Eeff[i]/Em[i]:.2f}", (Eeff[i]/Em[i], d[i]), fontsize=7,
                   xytext=(3, 3), textcoords="offset points")
    fig.tight_layout(); fig.savefig("study_mono2.png", dpi=160)
    print(f"wrote study_mono2.png  (RMS {rms:.0f}%)")

def _load_tensor(run_id):
    import os
    path = f"tensors/elasticity_tensor_{run_id}.csv"
    if not os.path.exists(path): return None
    rows = []
    with open(path) as f:
        for ln in f:
            p = ln.split(",")
            if p[0] in ("11", "22", "33", "12", "13", "23") and len(p) >= 7:
                rows.append([float(x) for x in p[1:7]])
            if len(rows) == 6: break
    return np.array(rows) if len(rows) == 6 else None

def do_basetensor():
    import glob, os
    ids = sorted(os.path.basename(p)[len("elasticity_tensor_"):-4]
                 for p in glob.glob("tensors/elasticity_tensor_BTEN_*.csv"))
    if not ids:
        print("[skip] no BTEN tensor csvs"); return
    rec = []
    for rid in ids:
        C = _load_tensor(rid)
        if C is None: continue
        S = np.linalg.inv(0.5*(C+C.T))
        rec.append(dict(d=int(rid.split("z")[1])/100,
                        Ex=1/S[0,0]/GPa, Ey=1/S[1,1]/GPa, Ez=1/S[2,2]/GPa,
                        Gyz=1/S[5,5]/GPa, Gxz=1/S[4,4]/GPa, Gxy=1/S[3,3]/GPa))
    if not rec: return
    r = pd.DataFrame(rec).sort_values("d")
    fig, ax = plt.subplots(1, 2, figsize=(13, 6))
    a = ax[0]
    a.plot(r.Ez, r.d, "o-", label="E_z (vertical)")
    a.plot(0.5*(r.Ex+r.Ey), r.d, "s-", label="E_xy (horizontal)")
    a.set_ylim(1, 0); a.set_xlabel("Young's modulus (GPa)"); a.set_ylabel("depth z/H")
    a.set_title("(a) E_z vs E_xy"); a.grid(alpha=0.3); a.legend(fontsize=9)
    a = ax[1]
    a.plot(r.Ez/(0.5*(r.Ex+r.Ey)), r.d, "o-", color="tab:purple", label="E_z/E_xy")
    a.plot(0.5*(r.Gxz+r.Gyz)/r.Gxy, r.d, "s-", color="tab:green", label="G_axial/G_xy")
    a.axvline(1, color="0.6"); a.set_ylim(1, 0)
    a.set_xlabel("anisotropy ratio"); a.set_ylabel("depth z/H")
    a.set_title("(b) transverse-isotropy ratios"); a.grid(alpha=0.3); a.legend(fontsize=9)
    fig.tight_layout(); fig.savefig("study_basetensor.png", dpi=160)
    print("wrote study_basetensor.png")

# ===========================================================================
# Third battery
# ---------------------------------------------------------------------------
def do_orient():
    df = load("results_orient.csv")
    if df is None: return
    df = df.set_index("run_id")
    order = ["ORI_rand", "ORI_X50", "ORI_X90", "ORI_Z50", "ORI_Z90",
             "ORI_granular", "ORI_columnar"]
    order = [r for r in order if r in df.index]
    lab = {"ORI_rand": "random\npockets", "ORI_X50": "X-aligned\nconc .5",
           "ORI_X90": "X-aligned\nconc .9", "ORI_Z50": "Z-aligned\nconc .5",
           "ORI_Z90": "Z-aligned\nconc .9", "ORI_granular": "granular\n(no chan)",
           "ORI_columnar": "columnar\n(Z chan)"}
    Ex = df.loc[order, "E_x"].astype(float).values / GPa
    Ez = df.loc[order, "E_z"].astype(float).values / GPa
    aniso = Ez / Ex
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.4))
    x = np.arange(len(order))
    ax[0].bar(x-0.2, Ex, 0.4, label="E_x", color="tab:blue")
    ax[0].bar(x+0.2, Ez, 0.4, label="E_z", color="tab:orange")
    ax[0].set_xticks(x); ax[0].set_xticklabels([lab[r] for r in order], fontsize=8)
    ax[0].set_ylabel("modulus (GPa)"); ax[0].set_title("(a) E_x vs E_z by texture")
    ax[0].legend(); ax[0].grid(alpha=0.3, axis="y")
    ax[1].bar(x, aniso, 0.6, color="tab:purple")
    ax[1].axhline(1, color="0.5")
    ax[1].set_xticks(x); ax[1].set_xticklabels([lab[r] for r in order], fontsize=8)
    ax[1].set_ylabel("E_z / E_x"); ax[1].set_title("(b) vertical anisotropy")
    ax[1].set_ylim(0.97, max(1.08, aniso.max()*1.02)); ax[1].grid(alpha=0.3, axis="y")
    for i, v in enumerate(aniso):
        ax[1].annotate(f"{v:.3f}", (i, v), ha="center", va="bottom", fontsize=8)
    fig.tight_layout(); fig.savefig("study_orient.png", dpi=160)
    print("wrote study_orient.png")

def do_gas():
    df = load("results_gas.csv")
    if df is None: return
    g = df.run_id.str.extract(r"v(\d+)")[0].astype(float).values / 100
    Ex = df.E_x.values / GPa
    o = np.argsort(g); g, Ex = g[o], Ex[o]
    fig, a = plt.subplots(figsize=(7.5, 5.6))
    a.plot(g*100, Ex, "o-", color="tab:blue")
    if len(Ex) and Ex[0] > 0:
        a.plot(g*100, Ex[0]*(1-1.9*g), "--", color="0.6",
               label="dilute-void guide  E0(1-1.9 phi)")
        a.legend(fontsize=8)
    a.set_xlabel("air-void fraction (%)"); a.set_ylabel("Young's modulus E_x (GPa)")
    a.set_title("Gas/porosity sweep (fixed brine phi_b=0.02)")
    a.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig("study_gas.png", dpi=160)
    print("wrote study_gas.png")

def do_seeds():
    df = load("results_seeds.csv")
    if df is None: return
    cfg = df.run_id.str.extract(r"SEED_(\w+?)_s")[0]
    df = df.assign(cfg=cfg, Ex=df.E_x.astype(float)/GPa, Ez=df.E_z.astype(float)/GPa)
    df["aniso"] = df.Ez / df.Ex
    order = [c for c in ("top", "mid", "base") if c in df.cfg.values]
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.4))
    print("\n  config   E_x mean+/-std (GPa)   CoV%    E_z/E_x mean+/-std")
    means, stds, an_m, an_s = [], [], [], []
    for c in order:
        g = df[df.cfg == c]
        # population s.d. (ddof=0), the convention declared in the manuscript;
        # pandas defaults to the sample s.d., which is ~12% wider at 5 packings
        m, s = g.Ex.mean(), g.Ex.std(ddof=0)
        am, as_ = g.aniso.mean(), g.aniso.std(ddof=0)
        means.append(m); stds.append(s); an_m.append(am); an_s.append(as_)
        print(f"  {c:6s}  {m:6.2f} +/- {s:.3f}     {100*s/m:4.1f}   {am:.3f} +/- {as_:.3f}")
    x = np.arange(len(order))
    ax[0].bar(x, means, yerr=stds, capsize=6, color="tab:blue")
    ax[0].set_xticks(x); ax[0].set_xticklabels(order)
    ax[0].set_ylabel("E_x (GPa)"); ax[0].set_title("(a) E_x: mean +/- 1 std (5 packings)")
    ax[0].grid(alpha=0.3, axis="y")
    ax[1].errorbar(x, an_m, yerr=an_s, fmt="o", capsize=6, ms=9, color="tab:purple")
    ax[1].axhline(1, color="0.5"); ax[1].set_xticks(x); ax[1].set_xticklabels(order)
    ax[1].set_ylabel("E_z/E_x"); ax[1].set_title("(b) anisotropy: mean +/- 1 std")
    lo = min(0.99, min(np.array(an_m)-np.array(an_s)))
    hi = max(1.05, max(np.array(an_m)+np.array(an_s)))
    ax[1].set_ylim(lo-0.005, hi+0.005); ax[1].grid(alpha=0.3, axis="y")
    for i, (m, s) in enumerate(zip(an_m, an_s)):
        ax[1].annotate(f"{m:.3f}", (i, m), xytext=(8, 0),
                       textcoords="offset points", va="center", fontsize=8)
    fig.tight_layout(); fig.savefig("study_seeds.png", dpi=160)
    print("wrote study_seeds.png")

def do_brine():
    df = load("results_brine.csv")
    if df is None: return
    G0, K0 = 440029.33528897085, 2.2e9
    Gmul = {"G0p1": 0.1, "G1": 1, "G10": 10, "G100": 100, "G1000": 1000}
    Kmul = {"K0p1": 0.1, "G1": 1, "K10": 10}
    df = df.assign(micro=df.run_id.str.extract(r"BRINE_(\w+?)_")[0],
                   tag=df.run_id.str.extract(r"BRINE_\w+?_(\w+)")[0],
                   Ex=df.E_x.astype(float)/GPa, Ez=df.E_z.astype(float)/GPa)
    df["aniso"] = df.Ez / df.Ex
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.6))
    for micro, c in (("iso", "tab:blue"), ("chan", "tab:red")):
        g = df[(df.micro == micro) & (df.tag.isin(Gmul))].copy()
        g["G"] = g.tag.map(Gmul) * G0
        g = g.sort_values("G")
        ax[0].plot(g.G, g.Ex, "o-", color=c, label=f"{micro} pockets")
        ax[1].plot(g.G, g.aniso, "o-", color=c, label=f"{micro}")
    ax[0].axvline(G0, color="0.6", ls=":"); ax[0].set_xscale("log")
    ax[0].annotate("physical\nbrine G", (G0, ax[0].get_ylim()[0]), fontsize=8, color="0.4")
    ax[0].set_xlabel("brine shear modulus G (Pa)"); ax[0].set_ylabel("E_x (GPa)")
    ax[0].set_title("(a) E sensitivity to brine stiffness"); ax[0].grid(alpha=0.3, which="both")
    ax[0].legend(fontsize=9)
    ax[1].axvline(G0, color="0.6", ls=":"); ax[1].set_xscale("log"); ax[1].axhline(1, color="0.7")
    ax[1].set_xlabel("brine shear modulus G (Pa)"); ax[1].set_ylabel("E_z / E_x")
    ax[1].set_title("(b) anisotropy sensitivity"); ax[1].grid(alpha=0.3, which="both")
    ax[1].legend(fontsize=9)
    fig.tight_layout(); fig.savefig("study_brine.png", dpi=160)
    print("wrote study_brine.png")
    # K sensitivity printout
    print("\n  K-sweep (E_x GPa, fixed G):")
    for micro in ("iso", "chan"):
        k = df[(df.micro == micro) & (df.tag.isin(Kmul))].copy()
        k["K"] = k.tag.map(Kmul) * K0
        k = k.sort_values("K")
        print("   " + micro + ": " + "  ".join(f"K={kk:.1e}:{ee:.2f}" for kk, ee in zip(k.K, k.Ex)))

def do_scf():
    df = load("results_scf.csv")
    if df is None: return
    order = ["GAS_v00", "BRINE_iso", "MORF_elong", "GAS_v10", "BRINE_chan", "SEAS_base"]
    lab = {"GAS_v00": "control\n(2% brine)", "BRINE_iso": "brine\npockets 5%",
           "MORF_elong": "elongated\npockets", "GAS_v10": "gas voids\n10%",
           "BRINE_chan": "brine\n+channels", "SEAS_base": "warm base\n(high VoF)"}
    df = df.set_index("run_id")
    order = [o for o in order if o in df.index]
    x = np.arange(len(order))
    p50 = df.loc[order, "SCF_p50"].astype(float).values
    p90 = df.loc[order, "SCF_p90"].astype(float).values
    p99 = df.loc[order, "SCF_p99"].astype(float).values
    mx = df.loc[order, "SCF_max"].astype(float).values
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.4))
    w = 0.26
    ax[0].bar(x-w, p50, w, label="P50", color="tab:green")
    ax[0].bar(x,   p90, w, label="P90", color="tab:orange")
    ax[0].bar(x+w, p99, w, label="P99", color="tab:red")
    ax[0].axhline(2.0, color="0.5", ls=":")
    ax[0].set_xticks(x); ax[0].set_xticklabels([lab[o] for o in order], fontsize=8)
    ax[0].set_ylabel("SCF = $\\sigma_1^{max}/\\bar\\sigma_{11}$")
    ax[0].set_title("(a) matrix stress-concentration distribution")
    ax[0].legend(); ax[0].grid(alpha=0.3, axis="y")
    fg2 = df.loc[order, "volfrac_SCF_gt2"].astype(float).values * 100
    ax[1].bar(x, fg2, color="tab:purple")
    ax[1].set_xticks(x); ax[1].set_xticklabels([lab[o] for o in order], fontsize=8)
    ax[1].set_ylabel("matrix volume with SCF > 2  (\\%)")
    ax[1].set_title("(b) extent of stress localization")
    ax[1].grid(alpha=0.3, axis="y")
    for i, v in enumerate(fg2):
        ax[1].annotate(f"{v:.1f}%", (i, v), ha="center", va="bottom", fontsize=8)
    fig.tight_layout(); fig.savefig("study_scf.png", dpi=160)
    print("wrote study_scf.png")

if __name__ == "__main__":
    do_scf()
    do_morphology()
    do_percolation()
    do_monotonic()
    do_seasonal()
    do_channel()
    do_fymy()
    do_mono2()
    do_basetensor()
    do_orient()
    do_gas()
    do_seeds()
    do_brine()
    print("done.")
