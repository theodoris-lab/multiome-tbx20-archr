# Created: 2026-04-24 17:25
# Updated: 2026-04-24 17:25
"""
Plot Het chromVAR projection results.

Inputs (from outputs/chromvar_proj_rnaCM_0424/):
  - chromvar_cell_projection.tsv   (per-cell projection score along WT->HOM axis)
  - chromvar_projection_summary.tsv (group means by genotype x sample)
  - chromvar_axis_loadings.tsv      (per-motif contribution to axis)

Outputs (-> outputs/chromvar_proj_rnaCM_0424/plots/):
  - projection_density.{png,pdf}        density overlay WT/HET/HOM
  - projection_violin_by_sample.{png,pdf} violin per sample (batch check)
  - axis_loadings_top.{png,pdf}         top 20 WT-side + top 20 HOM-side bar
  - projection_summary_table.md         report-ready markdown table
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INDIR = ROOT / "outputs/chromvar_proj_rnaCM_0424"
OUTDIR = INDIR / "plots"
OUTDIR.mkdir(parents=True, exist_ok=True)

GENO_ORDER = ["WT", "HET", "HOM"]
GENO_COLOR = {"WT": "#1f77b4", "HET": "#ff7f0e", "HOM": "#d62728"}


def save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(OUTDIR / f"{stem}.{ext}", dpi=200, bbox_inches="tight")
    print(f"  saved: {stem}.png/.pdf")


# ---------------- 1. Density overlay ----------------
proj = pd.read_csv(INDIR / "chromvar_cell_projection.tsv", sep="\t")
proj = proj[proj["genotype"].isin(GENO_ORDER)].copy()
print(f"loaded per-cell projections: n={len(proj)}")
print(proj.groupby("genotype")["projection"].describe()[["count", "mean", "50%", "std"]])

fig, ax = plt.subplots(figsize=(6.4, 4.0))
xmin, xmax = proj["projection"].quantile([0.001, 0.999])
xs = np.linspace(xmin, xmax, 400)
for g in GENO_ORDER:
    vals = proj.loc[proj["genotype"] == g, "projection"].values
    if len(vals) < 5:
        continue
    # gaussian kde via numpy histogram smoothing
    from scipy.stats import gaussian_kde

    kde = gaussian_kde(vals, bw_method=0.25)
    ax.fill_between(xs, kde(xs), alpha=0.25, color=GENO_COLOR[g])
    ax.plot(xs, kde(xs), color=GENO_COLOR[g], lw=1.8,
            label=f"{g} (n={len(vals)})")
ax.axvline(0, color="grey", ls="--", lw=0.8)
ax.axvline(1, color="grey", ls="--", lw=0.8)
ax.text(0, ax.get_ylim()[1] * 0.97, " WT centroid", fontsize=8, color="grey",
        ha="left", va="top")
ax.text(1, ax.get_ylim()[1] * 0.97, " HOM centroid", fontsize=8, color="grey",
        ha="left", va="top")
ax.set_xlabel("Projection score on WT->HOM motif-deviation axis\n(0 = WT mean, 1 = HOM mean)")
ax.set_ylabel("Density")
ax.set_title("Per-cell chromVAR projection (rnaCM)")
ax.legend(frameon=False, loc="upper right")
save(fig, "projection_density")
plt.close(fig)

# ---------------- 2. Violin per sample ----------------
samples = sorted(proj["sample"].unique())
fig, ax = plt.subplots(figsize=(max(5, 1.4 * len(samples) * 3), 4.0))

positions, ticks, labels = [], [], []
for i, s in enumerate(samples):
    for j, g in enumerate(GENO_ORDER):
        vals = proj.loc[(proj["sample"] == s) & (proj["genotype"] == g),
                        "projection"].values
        if len(vals) < 5:
            continue
        pos = i * 4 + j
        positions.append(pos)
        parts = ax.violinplot(vals, positions=[pos], widths=0.85,
                              showmeans=False, showmedians=True,
                              showextrema=False)
        for body in parts["bodies"]:
            body.set_facecolor(GENO_COLOR[g])
            body.set_alpha(0.7)
            body.set_edgecolor("black")
            body.set_linewidth(0.5)
        if "cmedians" in parts:
            parts["cmedians"].set_color("black")
            parts["cmedians"].set_linewidth(1.0)
    ticks.append(i * 4 + 1)
    labels.append(s)

ax.axhline(0, color="grey", ls="--", lw=0.8)
ax.axhline(1, color="grey", ls="--", lw=0.8)
ax.set_xticks(ticks)
ax.set_xticklabels(labels)
ax.set_ylabel("Projection score")
ax.set_title("Per-cell chromVAR projection by sample")

handles = [plt.Rectangle((0, 0), 1, 1, color=GENO_COLOR[g], alpha=0.7,
                         label=g) for g in GENO_ORDER]
ax.legend(handles=handles, frameon=False, loc="upper left")
save(fig, "projection_violin_by_sample")
plt.close(fig)

# ---------------- 3. Axis loadings top 20 each side ----------------
load = pd.read_csv(INDIR / "chromvar_axis_loadings.tsv", sep="\t")
load["TF_clean"] = load["TF"].str.replace(r"_\d+$", "", regex=True)
top_neg = load.sort_values("delta").head(20).iloc[::-1]   # WT-side (negative delta)
top_pos = load.sort_values("delta", ascending=False).head(20).iloc[::-1]  # HOM-side

fig, axes = plt.subplots(1, 2, figsize=(10.5, 6.5), sharex=True)
for ax, df, title, color in [
    (axes[0], top_neg, "Top WT-side TFs (delta < 0)", GENO_COLOR["WT"]),
    (axes[1], top_pos, "Top HOM-side TFs (delta > 0)", GENO_COLOR["HOM"]),
]:
    ax.barh(df["TF_clean"], df["delta"], color=color, alpha=0.8,
            edgecolor="black", linewidth=0.4)
    ax.axvline(0, color="black", lw=0.5)
    ax.set_xlabel("HOM mean z  -  WT mean z")
    ax.set_title(title)
    ax.tick_params(axis="y", labelsize=8)

fig.suptitle("chromVAR axis loadings (rnaCM WT->HOM)", fontsize=11)
fig.tight_layout()
save(fig, "axis_loadings_top")
plt.close(fig)

# ---------------- 4. Markdown summary table ----------------
summ = pd.read_csv(INDIR / "chromvar_projection_summary.tsv", sep="\t")
summ = summ.sort_values(["sample", "genotype"])
md_lines = ["| sample | genotype | n_cells | proj_mean | proj_median | proj_sd |",
            "|--------|----------|---------|-----------|-------------|---------|"]
for _, r in summ.iterrows():
    md_lines.append(
        f"| {r['sample']} | {r['genotype']} | {int(r['n_cells'])} | "
        f"{r['proj_mean']:.3f} | {r['proj_median']:.3f} | {r['proj_sd']:.3f} |"
    )
overall = proj.groupby("genotype")["projection"].agg(["count", "mean", "median", "std"])
md_lines.append("")
md_lines.append("Overall (pooled across samples):")
md_lines.append("")
md_lines.append("| genotype | n_cells | proj_mean | proj_median | proj_sd |")
md_lines.append("|----------|---------|-----------|-------------|---------|")
for g in GENO_ORDER:
    if g not in overall.index:
        continue
    r = overall.loc[g]
    md_lines.append(
        f"| {g} | {int(r['count'])} | {r['mean']:.3f} | "
        f"{r['median']:.3f} | {r['std']:.3f} |"
    )
(OUTDIR / "projection_summary_table.md").write_text("\n".join(md_lines) + "\n")
print("  saved: projection_summary_table.md")

print(f"\n[done] all outputs in {OUTDIR}")
