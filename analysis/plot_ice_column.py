"""Plots for the sea-ice column RVE study (parametric_sea_ice_column.csv ->
results_column.csv). Produces depth profiles and cross-plots."""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

GPa = 1e9
df = pd.read_csv("results_column.csv")
df = df.sort_values("run_id").reset_index(drop=True)
z = np.array([int(r[5:]) for r in df.run_id]) / 100.0   # depth fraction 0..1

# --- design drivers (winter FY profile used to build the CSV) -----------------
T = -20.0 + (-1.8 - (-20.0)) * z                         # linear T(z), degC
S = np.array([7.0,5.5,4.8,4.5,4.3,4.3,4.5,5.0,6.0,8.0])  # C-shaped salinity ppt
phi_b = S * (-49.185 / T + 0.532) / 1000.0               # Frankenstein-Garner
phi_total = df.VoF_sphere.values                         # brine pockets + gas
PHI_C = 0.05                                             # rule-of-fives threshold
z_perc = np.interp(PHI_C, phi_b, z)                      # percolation depth

Ex, Ey, Ez = df.E_x/GPa, df.E_y/GPa, df.E_z/GPa
Eeff = df.E_eff/GPa
Exy = 0.5*(Ex+Ey)
Gxy, Gxz, Gyz = df.G_xy/GPa, df.G_xz/GPa, df.G_yz/GPa
aniso = df.E_anisotropy.values
ezxy = df.E_z_over_xy.values
Emat = df.E_matrix.values/GPa

def style_depth(ax):
    ax.set_ylim(1.0, 0.0)                 # depth increases downward
    ax.axhspan(z_perc, 1.0, color="0.92", zorder=0)
    ax.axhline(z_perc, color="tab:red", ls="--", lw=1)
    ax.set_ylabel("depth  z/H   (0 = top / air, 1 = base / ocean)")
    ax.grid(alpha=0.3)

# ============================ FIGURE 1: profiles ==============================
fig, axes = plt.subplots(1, 4, figsize=(16, 6.5), sharey=True)

# (1) drivers: T and brine fraction
ax = axes[0]; style_depth(ax)
ax.plot(T, z, "o-", color="tab:blue", label="temperature")
ax.set_xlabel("temperature  (°C)", color="tab:blue")
ax.tick_params(axis="x", labelcolor="tab:blue")
ax2 = ax.twiny()
ax2.plot(phi_b*100, z, "s-", color="tab:green", label="brine vol. fraction")
ax2.set_xlabel("brine volume fraction  φ_b (%)", color="tab:green")
ax2.tick_params(axis="x", labelcolor="tab:green")
ax.set_title("(a) Drivers: temperature & brine")

# (2) Young's moduli E(z)
ax = axes[1]; style_depth(ax)
ax.plot(Ez, z, "o-", color="tab:red",   label="E_z (vertical)")
ax.plot(Exy, z, "s-", color="tab:blue", label="E_xy (horizontal)")
ax.plot(Emat, z, ":", color="0.4",      label="E_matrix (pure ice)")
ax.set_xlabel("Young's modulus  (GPa)")
ax.set_title("(b) Stiffness profile  E(z)")
ax.legend(loc="lower left", fontsize=8)

# (3) shear moduli
ax = axes[2]; style_depth(ax)
ax.plot(Gxz, z, "o-", label="G_xz")
ax.plot(Gyz, z, "^-", label="G_yz")
ax.plot(Gxy, z, "s-", label="G_xy (in-plane)")
ax.set_xlabel("shear modulus  (GPa)")
ax.set_title("(c) Shear profile  G(z)")
ax.legend(loc="lower left", fontsize=8)

# (4) anisotropy
ax = axes[3]; style_depth(ax)
ax.plot(ezxy, z, "o-", color="tab:purple", label="E_z / E_xy")
ax.plot(aniso, z, "s-", color="tab:orange", label="E_max / E_min")
ax.axvline(1.0, color="0.6", lw=1)
ax.set_xlabel("anisotropy ratio")
ax.set_title("(d) Anisotropy profile")
ax.legend(loc="lower left", fontsize=8)
ax.text(0.98, z_perc-0.02, "percolation\n(rule of fives)", color="tab:red",
        fontsize=8, ha="right", va="bottom")

fig.suptitle("Sea-ice column: depth-resolved RVE homogenization "
             "(first-year, winter T-profile)", fontsize=13)
fig.tight_layout(rect=[0,0,1,0.96])
fig.savefig("ice_column_profiles.png", dpi=160)
print("wrote ice_column_profiles.png")

# ============================ FIGURE 2: cross-plots ===========================
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

# (A) normalized stiffness vs brine fraction, with Voigt/Reuss-style context
ax = axes[0]
sc = ax.scatter(phi_b*100, Eeff/Emat, c=T, cmap="coolwarm_r", s=90,
                edgecolor="k", zorder=3)
for i in range(len(z)):
    ax.annotate(f"z{int(z[i]*100):02d}", (phi_b[i]*100, (Eeff/Emat)[i]),
                fontsize=7, xytext=(4,4), textcoords="offset points")
# dilute estimate E/E0 ~ 1 - k*phi (k~2 for soft/void inclusions), guide only
xx = np.linspace(0, phi_b.max()*100, 50)
ax.plot(xx, 1 - 0.019*xx, "--", color="0.5",
        label="dilute guide  E/E₀ ≈ 1 − 1.9·φ_b")
cb = fig.colorbar(sc, ax=ax); cb.set_label("temperature (°C)")
ax.set_xlabel("brine volume fraction  φ_b (%)")
ax.set_ylabel("normalized stiffness  E_eff / E_matrix")
ax.set_title("(a) Stiffness knockdown vs brine content")
ax.legend(fontsize=8); ax.grid(alpha=0.3)

# (B) anisotropy onset vs brine fraction
ax = axes[1]
chan = df.generate_channels.str.strip().str.lower().eq("yes").values
ax.scatter(phi_b[~chan]*100, ezxy[~chan], s=90, color="tab:blue",
           edgecolor="k", label="isolated pockets", zorder=3)
ax.scatter(phi_b[chan]*100, ezxy[chan], s=120, marker="D", color="tab:red",
           edgecolor="k", label="percolating channels", zorder=3)
ax.axhline(1.0, color="0.6", lw=1)
ax.axvline(PHI_C*100, color="tab:red", ls="--", lw=1)
ax.text(PHI_C*100+0.2, 1.001, "φ_c ≈ 5%\n(rule of fives)", color="tab:red", fontsize=8)
for i in range(len(z)):
    ax.annotate(f"z{int(z[i]*100):02d}", (phi_b[i]*100, ezxy[i]),
                fontsize=7, xytext=(4,4), textcoords="offset points")
ax.set_xlabel("brine volume fraction  φ_b (%)")
ax.set_ylabel("transverse-isotropy index  E_z / E_xy")
ax.set_title("(b) Vertical stiffening appears only after percolation")
ax.legend(fontsize=8); ax.grid(alpha=0.3)

fig.suptitle("Sea-ice column: stiffness & anisotropy vs brine content", fontsize=13)
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig("ice_column_crossplots.png", dpi=160)
print("wrote ice_column_crossplots.png")

# summary table to stdout
out = pd.DataFrame({"z/H":z, "T_C":np.round(T,1), "S_ppt":S,
                    "phi_b_%":np.round(phi_b*100,2),
                    "E_eff_GPa":np.round(Eeff,3), "E_z/E_xy":np.round(ezxy,4),
                    "channels":df.generate_channels.values})
print(out.to_string(index=False))
print(f"percolation depth z/H ~ {z_perc:.2f}")
