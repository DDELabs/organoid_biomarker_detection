"""Static figures (matplotlib), validated reference palette, thin marks, recessive axes."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from . import survival as S

BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e4e3df"


def _axes(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def km_plot(t, e, responder, title, path):
    fig, ax = plt.subplots(figsize=(4.6, 3.6), dpi=150)
    _axes(ax)
    for mask, name, col in ((responder, "Predicted responder", BLUE), (~responder, "Predicted non-responder", ORANGE)):
        x, y = S.km_curve(t[mask], e[mask])
        ax.step(x, y, where="post", color=col, linewidth=2, label=f"{name} (n={mask.sum()})")
        ax.text(x[-1], y[-1], f"  {y[-1]:.2f}", color=INK, fontsize=7, va="center")
    ax.axvline(60, color=MUTED, linestyle="--", linewidth=0.8)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Overall survival (months)", color=INK, fontsize=8)
    ax.set_ylabel("Survival probability", color=INK, fontsize=8)
    ax.set_title(title, color=INK, fontsize=9, loc="left")
    ax.legend(frameon=False, fontsize=7, loc="lower left")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def forest_plot(rows, title, path):
    """rows: list of (label, HR, low, high)."""
    fig, ax = plt.subplots(figsize=(5.2, 0.45 * len(rows) + 1.2), dpi=150)
    _axes(ax)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", color=GRID, linewidth=0.6)
    for i, (label, hr, lo, hi) in enumerate(rows):
        y = len(rows) - i
        if np.isfinite(hr):
            ax.plot([max(lo, 0.05), min(hi, 20)], [y, y], color=BLUE, linewidth=2)
            ax.plot(hr, y, "o", color=BLUE, markersize=7, markeredgecolor="white", markeredgewidth=2)
            ax.text(21, y, f"{hr:.2f} ({lo:.2f}-{hi:.2f})", fontsize=7, va="center", color=INK)
    ax.set_yticks(range(len(rows), 0, -1))
    ax.set_yticklabels([r[0] for r in rows], fontsize=7, color=INK)
    ax.axvline(1, color=MUTED, linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xlim(0.05, 20)
    ax.set_xlabel("Adjusted hazard ratio per SD of predicted IC50 (log scale)", color=INK, fontsize=8)
    ax.set_title(title, color=INK, fontsize=9, loc="left")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def stability_plot(table, title, path, top=15):
    t = table.sort_values("selection_freq", ascending=False).head(top)[::-1]
    fig, ax = plt.subplots(figsize=(6.4, 0.32 * len(t) + 1.0), dpi=150)
    _axes(ax)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", color=GRID, linewidth=0.6)
    cols = [BLUE if c < 0 else ORANGE for c in t["mean_coef"]]
    ax.barh(range(len(t)), t["selection_freq"], color=cols, height=0.6)
    ax.set_yticks(range(len(t)))
    ax.set_yticklabels([p.replace("REACTOME_", "").replace("_", " ").lower()[:55] for p in t.index], fontsize=7, color=INK)
    ax.axvline(0.5, color=MUTED, linestyle="--", linewidth=0.8)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Selection frequency (blue: sensitising, orange: resistance)", color=INK, fontsize=8)
    ax.set_title(title, color=INK, fontsize=9, loc="left")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
