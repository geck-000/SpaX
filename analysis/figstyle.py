"""Shared figure styling for the depth-resolved study figures in main_rev.tex.

Okabe-Ito colourblind-safe palette, enlarged fonts, and the depth-on-vertical
convention requested by the reviewers (Reviewer 2: "put z/H on the y-axis in all
related figures; check font sizes"). Imported by make_rev_figs.py and the
analyze_*.py study figures so every depth profile shares one coordinate frame,
palette, and type size. z/H runs 0 (cold surface) at the top to 1 (warm base) at
the bottom, matching the ice-sheet schematic of Fig.2.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Okabe-Ito CVD-safe palette (matching viz/render_rve.py) ---------------
BLUE, ORANGE, GREEN = "#0072B2", "#E69F00", "#009E73"
VERM, PURPLE, SKY, BLACK = "#D55E00", "#CC79A7", "#56B4E9", "#222222"

RC = {
    "font.size": 14, "axes.titlesize": 15, "axes.labelsize": 15,
    "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 12.5,
    "axes.linewidth": 1.0, "lines.linewidth": 2.2, "lines.markersize": 7,
    "figure.dpi": 120, "savefig.bbox": "tight", "axes.grid": True,
    "grid.alpha": 0.30, "grid.linewidth": 0.7,
}


def apply():
    """Install the shared rcParams (call once at import time in each script)."""
    plt.rcParams.update(RC)


def depth_axis(ax, label=True):
    """Put z/H on the vertical axis: 0 (cold surface) at top, 1 (warm base) at
    bottom, with the warm percolated bottom band lightly shaded."""
    ax.set_ylim(1.0, 0.0)                 # 0 (surface) at top, 1 (base) at bottom
    if label:
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
