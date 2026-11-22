"""Static figures: volcano plot of genes and a bar chart of the top pathways."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

UP, DOWN, NEUTRAL = "#e34948", "#2a78d6", "#c9c8c2"
INK, MUTED = "#0b0b0b", "#52514e"


def _style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.grid(axis="x", color="#e6e5e0", linewidth=0.8)
    ax.set_axisbelow(True)


def volcano(de: pd.DataFrame, path: str | Path, title: str = "", fdr: float = 0.05,
            min_abs_lfc: float = 1.0, label_top: int = 8) -> Path:
    x = de["log2fc"]
    y = -np.log10(de["pvalue"].clip(lower=1e-300))
    sig = (de["fdr"] < fdr) & (x.abs() >= min_abs_lfc)
    up, down = sig & (x > 0), sig & (x < 0)

    fig, ax = plt.subplots(figsize=(7, 5.5), dpi=150)
    _style(ax)
    ax.scatter(x[~sig], y[~sig], s=6, c=NEUTRAL, linewidths=0, rasterized=True)
    ax.scatter(x[up], y[up], s=9, c=UP, linewidths=0, label=f"up ({int(up.sum())})")
    ax.scatter(x[down], y[down], s=9, c=DOWN, linewidths=0, label=f"down ({int(down.sum())})")
    for v in (-min_abs_lfc, min_abs_lfc):
        ax.axvline(v, color=MUTED, linewidth=0.6, linestyle=":")
    for g in de[sig].sort_values("pvalue").index[:label_top]:
        ax.annotate(g, (x[g], y[g]), xytext=(3, 2), textcoords="offset points", fontsize=7, color=INK)
    ax.set_xlabel("log2 fold change", color=INK)
    ax.set_ylabel("-log10 p-value", color=INK)
    ax.set_title(title, loc="left", fontsize=11, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2,
              title=f"FDR < {fdr}, |log2FC| ≥ {min_abs_lfc}", title_fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return Path(path)


def top_pathways(gsea: pd.DataFrame, path: str | Path, title: str = "", n: int = 20) -> Path:
    """Top pathways by GSEA FDR. An FDR of 0 means below the permutation floor."""
    top = gsea.dropna(subset=["nes"]).nsmallest(n, ["fdr"]).sort_values("nes")
    fig, ax = plt.subplots(figsize=(8, 0.32 * len(top) + 1.2), dpi=150)
    _style(ax)
    colors = [UP if v > 0 else DOWN for v in top["nes"]]
    ax.barh(range(len(top)), top["nes"], color=colors, height=0.7)
    labels = [p if len(p) <= 60 else p[:57] + "..." for p in top["pathway"]]
    ax.set_yticks(range(len(top)), labels, fontsize=8, color=INK)
    for i, (nes, q) in enumerate(zip(top["nes"], top["fdr"])):
        text = f"FDR {q:.2g}" if q > 0 else "FDR < 0.001"
        ax.annotate(text, (nes, i), xytext=(4 if nes > 0 else -4, 0), textcoords="offset points",
                    ha="left" if nes > 0 else "right", va="center", fontsize=7, color=MUTED)
    ax.axvline(0, color=MUTED, linewidth=0.8)
    lim = max(abs(top["nes"]).max() * 1.35, 1)
    ax.set_xlim(-lim, lim)
    ax.set_xlabel("normalized enrichment score (red = up in case, blue = down)", color=INK, fontsize=9)
    fig.suptitle(title, x=0.01, ha="left", fontsize=11, color=INK)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return Path(path)
