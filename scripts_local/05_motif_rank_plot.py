# Created: 2026-04-23 20:00
# Updated: 2026-04-23 20:00
"""
Rank-sorted motif enrichment plot (Kathiriya 2026 Fig 5E/F style).

Input:
    outputs/dar_rnaCM_0423/motif_HOM_vs_WT_up.tsv
    outputs/dar_rnaCM_0423/motif_HOM_vs_WT_down.tsv

Output:
    outputs/dar_rnaCM_0423/figures/fig5E_motif_rank_HOM_up.pdf
    outputs/dar_rnaCM_0423/figures/fig5F_motif_rank_HOM_down.pdf
    outputs/dar_rnaCM_0423/figures/fig5EF_combined_motif_rank_HOM.pdf
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(__file__).resolve().parents[1]
IN_DIR = ROOT / "outputs" / "dar_rnaCM_0423"
OUT_DIR = IN_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TOP_N = 30
DPI = 300

# ArchR "comet" palette (blue → purple → orange → yellow)
COMET_COLORS = ["#1B1B3A", "#3B3363", "#753D80", "#C04F62", "#F68B43", "#FDD84E"]
COMET = LinearSegmentedColormap.from_list("comet", COMET_COLORS, N=256)


def _clean_tf(name: str) -> str:
    """GATA2_388 → GATA2, ENSG00000234254_595 → ENSG00000234254."""
    return re.sub(r"_\d+$", "", name)


def _load(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df.sort_values("mlog10Padj", ascending=False).reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)
    df["label"] = df["TF"].map(_clean_tf)
    return df


def _single_panel(df: pd.DataFrame, title: str, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.0, 7.5))
    vmax = df["mlog10Padj"].max()
    sc = ax.scatter(
        df["rank"],
        df["mlog10Padj"],
        c=df["mlog10Padj"],
        cmap=COMET,
        s=12,
        vmin=0,
        vmax=vmax,
        linewidths=0,
    )
    top = df.head(TOP_N)
    texts = [
        ax.text(r, y, lbl, fontsize=8)
        for r, y, lbl in zip(top["rank"], top["mlog10Padj"], top["label"])
    ]
    adjust_text(
        texts,
        ax=ax,
        arrowprops=dict(arrowstyle="-", color="grey", lw=0.4),
        expand=(1.2, 1.6),
        force_text=(0.4, 0.8),
    )
    ax.set_xlabel("Rank sorted TFs enriched")
    ax.set_ylabel("-log10(P-adj) motif enrichment")
    ax.set_title(title)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("-log10(P-adj)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".png"), dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out_path.name}  (+ .png)  (top TF: {top['label'].iloc[0]})")


def _combined_panel(df_up: pd.DataFrame, df_down: pd.DataFrame, out_path: Path) -> None:
    up = df_up.copy()
    down = df_down.copy()
    up["signed_rank"] = up["rank"]
    down["signed_rank"] = -down["rank"]
    up["direction"] = "HOM_up"
    down["direction"] = "HOM_down"
    combined = pd.concat([up, down], ignore_index=True)
    combined.loc[combined["mlog10Padj"] == 0, "mlog10Padj"] = 0.001

    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    vmax = combined["mlog10Padj"].max()
    sc = ax.scatter(
        combined["signed_rank"],
        combined["mlog10Padj"],
        c=combined["mlog10Padj"],
        cmap=COMET,
        s=10,
        vmin=0,
        vmax=vmax,
        linewidths=0,
    )
    top_up = up.head(TOP_N)
    top_down = down.head(TOP_N)
    texts = []
    for r, y, lbl in zip(top_up["signed_rank"], top_up["mlog10Padj"], top_up["label"]):
        texts.append(ax.text(r, y, lbl, fontsize=7))
    for r, y, lbl in zip(top_down["signed_rank"], top_down["mlog10Padj"], top_down["label"]):
        texts.append(ax.text(r, y, lbl, fontsize=7))
    adjust_text(
        texts,
        ax=ax,
        arrowprops=dict(arrowstyle="-", color="grey", lw=0.3),
        expand=(1.1, 1.5),
        force_text=(0.3, 0.6),
    )
    ax.axvline(0, linestyle="--", color="grey", lw=0.8)
    ax.set_xlabel("Rank sorted TFs enriched  (DOWN ← | → UP)")
    ax.set_ylabel("-log10(P-adj) motif enrichment")
    ax.set_title("HOM vs WT — motif enrichment (combined)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("-log10(P-adj)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".png"), dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out_path.name}  (+ .png)")


def main() -> None:
    up = _load(IN_DIR / "motif_HOM_vs_WT_up.tsv")
    down = _load(IN_DIR / "motif_HOM_vs_WT_down.tsv")
    print(f"Loaded {len(up)} up-motif TFs, {len(down)} down-motif TFs")
    print(f"Top UP:   {', '.join(up['label'].head(10).tolist())}")
    print(f"Top DOWN: {', '.join(down['label'].head(10).tolist())}")

    _single_panel(
        up,
        "HOM vs WT — up-peak motif enrichment (Fig 5E style)",
        OUT_DIR / "fig5E_motif_rank_HOM_up.pdf",
    )
    _single_panel(
        down,
        "HOM vs WT — down-peak motif enrichment (Fig 5F style)",
        OUT_DIR / "fig5F_motif_rank_HOM_down.pdf",
    )
    _combined_panel(up, down, OUT_DIR / "fig5EF_combined_motif_rank_HOM.pdf")
    print("Done.")


if __name__ == "__main__":
    main()
