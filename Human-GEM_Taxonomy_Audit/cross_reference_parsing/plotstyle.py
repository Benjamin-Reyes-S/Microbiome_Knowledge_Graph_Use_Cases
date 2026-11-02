"""Figure style of the cross-reference parsing comparison (same tokens as ../plotstyle.py, own colours for the three outcomes)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import FIGURES

SURFACE, INK, INK2, RULE = "#fcfcfb", "#0b0b0b", "#52514e", "#d9d8d3"
BLUE, ORANGE, GRAY = "#2a78d6", "#eb6834", "#9a9a95"        # human / other taxa only / cannot be mapped
plt.rcParams.update({"figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "font.size": 10,
                     "axes.edgecolor": RULE, "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK})


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.xaxis.grid(True, color=RULE, linewidth=0.5)
    ax.set_axisbelow(True)


def finish(fig, name):
    fig.tight_layout()
    fig.savefig(FIGURES / name, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("figure:", FIGURES / name)
