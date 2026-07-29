import numpy as np, pandas as pd, matplotlib.pyplot as plt
GPa=1e9
old = pd.read_csv("results_column.csv").sort_values("run_id").reset_index(drop=True)
new = pd.read_csv("results_marchenko.csv").sort_values("run_id").reset_index(drop=True)
d = np.array([int(r[5:]) for r in new.run_id])/100.0
Eold, Enew = old.E_eff.values/GPa, new.E_eff.values/GPa

Ebot, M, n = 1.67, 2.63, 0.5
zeta = 1.0 - d
E_march = Ebot*((M-1)*zeta**n + 1.0)
zz=np.linspace(0.001,1,200); E_march_c = Ebot*((M-1)*zz**n+1.0); d_c=1-zz

rms = np.sqrt(np.mean((Enew-E_march)**2))
rms_pct = 100*np.sqrt(np.mean(((Enew-E_march)/E_march)**2))

fig, ax = plt.subplots(1,2, figsize=(13,6))
a=ax[0]
a.plot(E_march_c, d_c, "-", color="tab:red", lw=2, label="Marchenko 2024 fit")
a.plot(E_march, d, "s", color="tab:red", ms=5, label="Marchenko @ slice depths")
a.plot(Enew, d, "o-", color="tab:green", label="SPAX recalibrated (E_mat≈4.6 GPa)")
a.plot(Eold, d, "o--", color="0.7", lw=1, label="SPAX original (E_mat≈9.4 GPa)")
a.set_ylim(1,0); a.set_xlim(0,10)
a.set_xlabel("Young's modulus  (GPa)"); a.set_ylabel("depth  z/H  (0=top, 1=bottom)")
a.set_title("(a) Recalibrated SPAX vs Marchenko (2024)")
a.grid(alpha=0.3); a.legend(fontsize=8, loc="lower right")

a=ax[1]
a.plot(Enew/E_march, d, "o-", color="tab:green")
a.axvline(1, color="0.6")
a.set_ylim(1,0); a.set_xlim(0.8,1.5)
a.set_xlabel("E_SPAX(recal) / E_Marchenko"); a.set_ylabel("depth  z/H")
a.set_title(f"(b) Match ratio   (RMS error {rms_pct:.0f}%)")
a.grid(alpha=0.3)
for i in range(len(d)):
    a.annotate(f"{Enew[i]/E_march[i]:.2f}",(Enew[i]/E_march[i],d[i]),
               fontsize=7,xytext=(3,3),textcoords="offset points")
fig.suptitle("Marchenko-matched scenario: vibrating-beam-effective matrix modulus",fontsize=13)
fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig("ice_column_marchenko_match.png",dpi=160)
print("wrote ice_column_marchenko_match.png")
print("\n depth  E_recal  E_March  ratio")
for i in range(len(d)):
    print(f" {d[i]:.2f}  {Enew[i]:6.2f}  {E_march[i]:6.2f}  {Enew[i]/E_march[i]:.2f}")
print(f"\n RMS error vs Marchenko: {rms:.2f} GPa ({rms_pct:.0f}%)")
print(f" ends: top {Enew[0]:.2f} vs {E_march[0]:.2f} | bot {Enew[-1]:.2f} vs {E_march[-1]:.2f}")
