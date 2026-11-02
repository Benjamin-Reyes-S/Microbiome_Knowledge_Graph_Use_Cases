"""Shared figure style (dataviz skill): one hue per match status, thin 2 px-gapped segments, text in ink tokens, no grid clutter.
Light surface only (PNG for the walkthrough)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import FIGURES

SURFACE, INK, INK2, RULE = "#fcfcfb", "#0b0b0b", "#52514e", "#d9d8d3"
COLOR = {"matched_human": "#2a78d6", "mapped_not_human": "#eb6834", "unmappable": "#9a9a95"}   # blue, orange, neutral gray
NICE = {"matched_human": "matched to human", "mapped_not_human": "in KG, no human evidence", "unmappable": "no id in KG"}
plt.rcParams.update({"figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
                     "font.size": 10, "axes.edgecolor": RULE, "text.color": INK, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK})


def stacked_h(ax, names, values, series, totals=None, labels=True, colors=None, nice=None):
    """Horizontal stacked bars, share of each row's total (0-100 %), 2 px surface gap between segments."""
    colors, nice = colors or COLOR, nice or NICE
    import numpy as np
    v = np.asarray(values, float)
    tot = v.sum(axis=1, keepdims=True)
    pct = 100 * v / np.where(tot == 0, 1, tot)
    left = np.zeros(len(names))
    for k, s in enumerate(series):
        ax.barh(range(len(names)), pct[:, k], left=left, color=colors[s], height=0.62, edgecolor=SURFACE, linewidth=2, label=nice[s])
        if labels:
            for i in range(len(names)):
                if pct[i, k] >= 7:
                    ax.text(left[i] + pct[i, k] / 2, i, f"{int(v[i, k]):,}", ha="center", va="center", color="white", fontsize=9)
        left += pct[:, k]
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([f"{n}  (n={int(t):,})" if totals is not None else n for n, t in zip(names, totals if totals is not None else names)])
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of entities")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.xaxis.grid(True, color=RULE, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=9, loc="upper left", ncol=len(series), bbox_to_anchor=(0, -0.18 if len(names) < 8 else -0.08))


def finish(fig, name):
    fig.tight_layout()
    fig.savefig(FIGURES / name, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("figure:", FIGURES / name)
