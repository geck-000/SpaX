r"""Shared figure styling for the depth-resolved study figures.

Figures are authored at the size they are printed. The paper's text block is
415.13pt = 5.74in wide, so a figure drawn at ``FIG_W`` and included at
``\linewidth`` is reproduced 1:1 and its type is the size set here. Authoring
wider and letting LaTeX scale down is what made the earlier figures illegible:
a 12.6in figure at \linewidth is shrunk to 46%, taking a 14pt label to 6.4pt.

Type is Times to match the body font (the class loads times/mathptmx, which
ships as NimbusRomNo9L), with STIX for maths, which is Times-metric compatible.
Colours are the Okabe-Ito colourblind-safe palette used by the TikZ figures.

z/H runs 0 (cold surface) at the top to 1 (warm base) at the bottom, matching
the ice-sheet schematic.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Okabe-Ito CVD-safe palette (matching viz/render_rve.py) ---------------
BLUE, ORANGE, GREEN = "#0072B2", "#E69F00", "#009E73"
VERM, PURPLE, SKY, BLACK = "#D55E00", "#CC79A7", "#56B4E9", "#222222"
GREY = "#6E6E6E"

# Text-block width of the paper, in inches (415.13pt / 72.27).
FIG_W = 5.74

SERIF = ["Nimbus Roman", "TeX Gyre Termes", "STIXGeneral", "DejaVu Serif"]

RC = {
    "font.family": "serif", "font.serif": SERIF,
    "mathtext.fontset": "stix",
    "font.size": 8.5, "axes.titlesize": 9.0, "axes.labelsize": 9.0,
    "xtick.labelsize": 8.0, "ytick.labelsize": 8.0, "legend.fontsize": 7.5,
    "axes.linewidth": 0.7, "lines.linewidth": 1.5, "lines.markersize": 4.0,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 3.0, "ytick.major.size": 3.0,
    "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,
    "axes.grid": True, "grid.alpha": 0.5, "grid.linewidth": 0.4,
    "grid.color": "#C8C8C8",
    "legend.frameon": True, "legend.framealpha": 0.92,
    "legend.edgecolor": "#BBBBBB", "legend.borderpad": 0.4,
    "legend.handlelength": 1.8, "legend.handletextpad": 0.6,
    "legend.labelspacing": 0.35,
    "figure.dpi": 120, "savefig.dpi": 400,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "axes.axisbelow": True,
}


def apply():
    """Install the shared rcParams (call once at import time in each script)."""
    plt.rcParams.update(RC)


def size(aspect, frac=1.0):
    """Figure size for a figure printed at ``frac`` of the text width.

    ``aspect`` is height/width. Returns the (w, h) inches to author at, so the
    saved figure needs no scaling in the document.
    """
    w = FIG_W * frac
    return (w, w * aspect)


def clean(ax):
    """Drop the top and right spines and lighten what is left."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#888888")
    ax.tick_params(colors="#444444", labelcolor="black")


def panel(ax, letter, x=0.022, y=0.975, ha="left", va="top", **kw):
    """Panel letter, placed inside the axes in the paper's (a)/(b) style."""
    ax.text(x, y, "(%s)" % letter, transform=ax.transAxes, fontsize=9.0,
            fontweight="bold", va=va, ha=ha, zorder=6,
            bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.4), **kw)


def depth_axis(ax, label=True, band=True):
    """Put z/H on the vertical axis: 0 (cold surface) at top, 1 (warm base) at
    bottom, with the warm percolated bottom band lightly shaded."""
    ax.set_ylim(1.0, 0.0)
    if label:
        ax.set_ylabel(r"normalised depth $z/H$")
    if band:
        ax.axhspan(0.80, 1.0, color=SKY, alpha=0.13, lw=0, zorder=0)


def orient_labels(ax, fontsize=7.0):
    """Cold-surface / warm-bottom orientation, as ticks outside the right edge."""
    tw = ax.twinx()
    tw.set_ylim(1.0, 0.0)
    tw.set_yticks([0.0, 1.0])
    tw.set_yticklabels(["surface (cold)", "base (warm)"], fontsize=fontsize,
                       color="0.35")
    tw.tick_params(axis="y", length=0)
    tw.grid(False)
    for spine in tw.spines.values():
        spine.set_visible(False)
    return tw


def logticks(ax, axis, ticks, fmt="%g"):
    """Plain decimal ticks on a log axis, with the minor labels suppressed.

    Matplotlib's default log formatter writes every decade and subdivision as
    ``2 x 10^-2``; at 8pt in a 2in panel those overlap into an unreadable band.
    """
    from matplotlib.ticker import FixedLocator, NullFormatter, FuncFormatter
    a = ax.xaxis if axis == "x" else ax.yaxis
    a.set_major_locator(FixedLocator(ticks))
    a.set_major_formatter(FuncFormatter(lambda v, _p: fmt % v))
    a.set_minor_locator(FixedLocator([]))
    a.set_minor_formatter(NullFormatter())
