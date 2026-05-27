# Created: 2026-04-24 18:20
# Updated: 2026-04-24 18:20
"""
Plot per-TF Fisher enrichment in HOM-DOWN DAR peaks
(TALE family + WT-side positive controls).

Inputs (from outputs/meis_validation_rnaCM_0424/):
  - tale_perTF_enrichment.tsv
  - wt_side_perTF_enrichment.tsv

Outputs (-> outputs/meis_validation_rnaCM_0424/plots/):
  - perTF_OR_forest.{png,pdf}
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
INDIR = ROOT / "outputs/meis_validation_rnaCM_0424"
OUTDIR = INDIR / "plots"
OUTDIR.mkdir(parents=True, exist_ok=True)


def load(path):
    df = pd.read_csv(path, sep="\t")
    # 95% CI for OR via Wald on log scale
    a = df["n_target_hit"].astype(float)
    b = (df["n_target"] - df["n_target_hit"]).astype(float)
    c = df["n_bg_hit"].astype(float)
    d = (df["n_bg"] - df["n_bg_hit"]).astype(float)
    se_log = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    log_or = np.log(df["OR"].replace(0, np.nan))
    df["ci_lo"] = np.exp(log_or - 1.96 * se_log)
    df["ci_hi"] = np.exp(log_or + 1.96 * se_log)
    return df


tale = load(INDIR / "tale_perTF_enrichment.tsv").sort_values("OR", ascending=True)
wt = load(INDIR / "wt_side_perTF_enrichment.tsv").sort_values("OR", ascending=True)

fig, axes = plt.subplots(1, 2, figsize=(11, 7), sharex=True)

for ax, df, title, family_color in [
    (axes[0], tale, "TALE family (vs HOM-DOWN DAR enrichment)", "#377eb8"),
    (axes[1], wt, "WT-side controls (chromVAR top)", "#984ea3"),
]:
    y = np.arange(len(df))
    sig = (df["fdr"] < 0.05).values
    colors = np.where(sig, family_color, "lightgrey")
    ax.errorbar(
        df["OR"], y,
        xerr=[df["OR"] - df["ci_lo"], df["ci_hi"] - df["OR"]],
        fmt="o", ecolor="grey", elinewidth=0.8, capsize=2,
        markerfacecolor="none", markeredgecolor="none",
    )
    ax.scatter(df["OR"], y, c=colors, s=55, edgecolor="black", linewidth=0.6, zorder=3)
    ax.set_yticks(y)
    labels = []
    for tf, fdr in zip(df["TF_clean"], df["fdr"]):
        marker = " *" if fdr < 0.05 else ""
        labels.append(f"{tf}{marker}")
    ax.set_yticklabels(labels, fontsize=9)
    ax.axvline(1, color="black", lw=0.5, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("Odds ratio  (DOWN peak motif occurrence vs NS background)")
    ax.set_title(title, fontsize=11)
    ax.tick_params(axis="x", labelsize=9)

# annotate MEIS1 specifically on TALE panel
tale_idx = list(tale["TF_clean"]).index("MEIS1")
axes[0].annotate(
    "MEIS1 OR=0.97, NS\n(not enriched in HOM-DOWN)",
    xy=(tale.iloc[tale_idx]["OR"], tale_idx),
    xytext=(0.55, tale_idx - 4),
    fontsize=8, color="darkred",
    arrowprops=dict(arrowstyle="->", color="darkred", lw=0.8),
)

fig.suptitle(
    "Per-TF Fisher enrichment in HOM-DOWN DAR peaks (n=1,310) vs NS background (n=255,581)\n"
    "* = FDR < 0.05",
    fontsize=11,
)
fig.tight_layout(rect=[0, 0, 1, 0.94])
for ext in ("png", "pdf"):
    fig.savefig(OUTDIR / f"perTF_OR_forest.{ext}", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"saved: {OUTDIR}/perTF_OR_forest.png/.pdf")
