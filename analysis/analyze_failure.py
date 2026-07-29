"""Offline failure-onset analysis (pandas + matplotlib, NO Abaqus).

Reads results_failure.csv (one row per column slice ICE_zNN, produced by
failure_extract.py) and turns the strength-independent stress-concentration
percentiles into a depth-resolved first-failure prediction:

  tensile (max-principal):  sigma_fail(z) = sigma_t / SCF_p99(z)
  Mohr-Coulomb:             sigma_fail(z) = 2 c cos(phi) / MCnorm_p99(z)

The slice with the SMALLEST sigma_fail cracks first -> that depth is the
failure-onset depth. The *ranking* of depths is independent of the chosen
strengths (they only rescale the curve); sigma_t, c just set the absolute MPa.

Strengths (override via CLI): sigma_t [Pa]  c [Pa].  Defaults: 1.0 MPa, 0.6 MPa.
phi is read from the MC_phi_deg column. Uses P99 (mesh-robust), not the peak.

Usage:  python3 analyze_failure.py [results_failure.csv] [sigma_t_Pa] [c_Pa]
"""
import sys, math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fn = sys.argv[1] if len(sys.argv) > 1 else "results_failure.csv"
SIGMA_T = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0e6   # ice tensile strength
COH = float(sys.argv[3]) if len(sys.argv) > 3 else 0.6e6       # cohesion

df = pd.read_csv(fn)
# depth fraction from run_id ICE_zNN  (z05 -> 0.05 ... z95 -> 0.95; 0=surface,1=base)
df["depth"] = df["run_id"].str.extract(r"z(\d+)").astype(float) / 100.0
df = df.sort_values("depth").reset_index(drop=True)
phi = float(df["MC_phi_deg"].iloc[0]) if "MC_phi_deg" in df else 30.0

df["sigfail_tensile_MPa"] = (SIGMA_T / df["SCF_p99"]) / 1e6
df["sigfail_MC_MPa"] = (2.0 * COH * math.cos(math.radians(phi)) / df["MCnorm_p99"]) / 1e6

it = df["sigfail_tensile_MPa"].idxmin()
imc = df["sigfail_MC_MPa"].idxmin()
print("=== failure-onset (lower sigma_fail = fails first) ===")
print(df[["run_id", "depth", "SCF_p99", "MCnorm_p99",
          "sigfail_tensile_MPa", "sigfail_MC_MPa"]].to_string(index=False,
          float_format=lambda x: "%.3f" % x))
print("\nTensile first-failure : %s (depth %.0f%%)  sigma_fail=%.3f MPa"
      % (df.run_id[it], 100*df.depth[it], df.sigfail_tensile_MPa[it]))
print("Mohr-Coulomb first-failure (phi=%.0f): %s (depth %.0f%%)  sigma_fail=%.3f MPa"
      % (phi, df.run_id[imc], 100*df.depth[imc], df.sigfail_MC_MPa[imc]))

fig, (a0, a1) = plt.subplots(1, 2, figsize=(11, 4.6))
a0.plot(df.SCF_p99, df.depth, "o-", label="tensile SCF P99")
a0.plot(df.MCnorm_p99, df.depth, "s--", label="Mohr-Coulomb P99 (phi=%.0f)" % phi)
a0.invert_yaxis(); a0.set_ylabel("depth fraction (0=surface, 1=base)")
a0.set_xlabel("stress-concentration P99 (per unit macro stress)")
a0.set_title("(a) local stress concentration vs depth"); a0.grid(alpha=0.3); a0.legend()

a1.plot(df.sigfail_tensile_MPa, df.depth, "o-",
        label="tensile (sigma_t=%.1f MPa)" % (SIGMA_T/1e6))
a1.plot(df.sigfail_MC_MPa, df.depth, "s--",
        label="Mohr-Coulomb (c=%.1f MPa)" % (COH/1e6))
a1.scatter([df.sigfail_tensile_MPa[it]], [df.depth[it]], s=120, facecolors="none",
           edgecolors="C0", linewidths=2, zorder=5)
a1.scatter([df.sigfail_MC_MPa[imc]], [df.depth[imc]], s=120, facecolors="none",
           edgecolors="C1", linewidths=2, zorder=5)
a1.invert_yaxis(); a1.set_xlabel("macro stress to first failure [MPa]")
a1.set_title("(b) failure-onset stress vs depth (circles = first to fail)")
a1.grid(alpha=0.3); a1.legend()
fig.tight_layout()
fig.savefig("study_failure.png", dpi=160)
print("\nwrote study_failure.png")
