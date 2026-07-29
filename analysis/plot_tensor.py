import numpy as np, matplotlib.pyplot as plt

def load_C(run_id):
    rows = []
    with open(f"tensors/elasticity_tensor_{run_id}.csv") as f:
        for ln in f:
            if ln.startswith("#") or ln.strip()=="" or ln.startswith(","): continue
            p = ln.split(",")
            if p[0] in ("11","22","33","12","13","23") and len(p)>=7:
                rows.append([float(x) for x in p[1:7]])
            if len(rows)==6: break
    return np.array(rows)

lab = ["11","22","33","12","13","23"]
for rid in ["ICE_z05","ICE_z95"]:
    C = load_C(rid)/1e9
    Cs = 0.5*(C+C.T)                     # symmetrize (noise removal)
    S = np.linalg.inv(Cs)               # compliance (GPa^-1)
    Ex,Ey,Ez = 1/S[0,0],1/S[1,1],1/S[2,2]
    Gyz,Gxz,Gxy = 1/S[5,5],1/S[4,4],1/S[3,3]  # 23,13,12
    nuxy,nuxz,nuyz = -S[1,0]/S[0,0], -S[2,0]/S[0,0], -S[2,1]/S[1,1]
    # max normal-shear coupling as fraction of mean normal diagonal
    coup = np.abs(Cs[:3,3:]).max()/np.mean(np.diag(Cs)[:3])
    print(f"\n=== {rid} ===")
    print(f" E_x,E_y,E_z   = {Ex:.3f} {Ey:.3f} {Ez:.3f} GPa   (E_z/E_xy={Ez/(0.5*(Ex+Ey)):.3f})")
    print(f" G_yz,G_xz,G_xy= {Gyz:.3f} {Gxz:.3f} {Gxy:.3f} GPa")
    print(f" nu_xy,nu_xz,nu_yz = {nuxy:.3f} {nuxz:.3f} {nuyz:.3f}")
    print(f" normal-shear coupling (max/diag) = {coup*100:.2f}%  -> orthotropic")

# heatmap comparison
fig, axes = plt.subplots(1,2, figsize=(12,5.2))
for ax, rid, ttl in zip(axes, ["ICE_z05","ICE_z95"],
        ["z05  cold top  (T=-19°C, φ_b=2%)","z95  warm base  (T=-3°C, φ_b=15%)"]):
    C = load_C(rid)/1e9
    Cs = 0.5*(C+C.T)
    im = ax.imshow(Cs, cmap="viridis", vmin=0, vmax=13)
    ax.set_xticks(range(6)); ax.set_yticks(range(6))
    ax.set_xticklabels(lab); ax.set_yticklabels(lab)
    ax.set_title(ttl, fontsize=11)
    for i in range(6):
        for j in range(6):
            v=Cs[i,j]
            ax.text(j,i,f"{v:.2f}",ha="center",va="center",
                    color="white" if v<7 else "black",fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.046, label="C_ij (GPa)")
fig.suptitle("Effective stiffness tensor C_ij (Voigt): isotropic top → orthotropic base",
             fontsize=13)
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig("ice_column_tensor.png", dpi=160)
print("\nwrote ice_column_tensor.png")
