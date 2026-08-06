import numpy as np, pandas as pd, matplotlib.pyplot as plt


GPa=1e9
df = pd.read_csv("results_column.csv").sort_values("run_id").reset_index(drop=True)
d = np.array([int(r[5:]) for r in df.run_id])/100.0     # depth from top, 0..1
Eeff = df.E_eff.values/GPa
Ez, Exy = df.E_z.values/GPa, 0.5*(df.E_x.values+df.E_y.values)/GPa

# Marchenko (2024) Eq. (17), p.11 -- a Kerr & Palmer (1972) form fitted to
# E(z) computed from measured brine content via his empirical formula (5),
# not to the vibrating-beam moduli themselves.  d is depth from the top.
E0, alpha, n = 4.4, 0.38, 0.6
E_march = E0*(1.0 - (1.0 - alpha)*d**n)
M = 1.0/alpha                        # surface/base ratio, for the caption

k = np.mean(Eeff/E_march)            # mean overprediction factor
fig, ax = plt.subplots(1,2, figsize=(13,6))

# (a) depth profiles overlaid
a=ax[0]
a.plot(Eeff, d, "o-", color="tab:blue", label="SPAX RVE  E_eff (sim)")
a.plot(Ez,   d, "^--", color="tab:cyan", lw=1, label="SPAX  E_z (vertical)")
a.plot(E_march, d, "s-", color="tab:red", label="Marchenko 2024 fit")
a.plot(Eeff/k, d, ":", color="tab:blue", lw=1.5,
       label=f"SPAX shape, rescaled /{k:.1f}")
a.set_ylim(1,0); a.set_xlim(0,10)
a.set_xlabel("Young's modulus  (GPa)")
a.set_ylabel("depth  z/H  (0=top, 1=bottom)")
a.set_title("(a) Simulated vs Marchenko (2024)")
a.grid(alpha=0.3); a.legend(fontsize=8, loc="lower right")

# (b) ratio sim/experiment vs depth
a=ax[1]
a.plot(Eeff/E_march, d, "o-", color="tab:purple")
a.axvline(1, color="0.6"); a.axvline(2, color="0.6", ls=":")
a.set_ylim(1,0); a.set_xlim(0,3.5)
a.set_xlabel("E_sim / E_Marchenko")
a.set_ylabel("depth  z/H")
a.set_title("(b) Stiffness overprediction factor")
a.grid(alpha=0.3)
for i in range(len(d)):
    a.annotate(f"{Eeff[i]/E_march[i]:.1f}×",(Eeff[i]/E_march[i],d[i]),
               fontsize=7,xytext=(3,3),textcoords="offset points")

fig.suptitle("Sea-ice column stiffness: SPAX homogenization vs Marchenko (2024) field fit",
             fontsize=13)
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig("ice_column_marchenko.png", dpi=160)
print("wrote ice_column_marchenko.png")

print("\n depth  E_sim   E_March  ratio")
for i in range(len(d)):
    print(f" {d[i]:.2f}  {Eeff[i]:6.2f}  {E_march[i]:6.2f}  {Eeff[i]/E_march[i]:.2f}")
print(f"\n sim range : {Eeff.min():.2f}-{Eeff.max():.2f} GPa  (top/bot ratio {Eeff[0]/Eeff[-1]:.2f})")
print(f" March range: {E_march.min():.2f}-{E_march.max():.2f} GPa  (M={M})")
print(f" mean overprediction factor k = {k:.2f}")
