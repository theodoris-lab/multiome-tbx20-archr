# Created: 2026-04-21 00:15
# Quick local plot: L2FC × position scatter at chr7:34.8-35.8Mb (TBX20 ±500kb)
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

OUT = Path("/Users/bkim/vscode/V2_multiome_2026-04-14/outputs/dar_0420")
TBX20_START = 35_242_030
TBX20_END = 35_293_758
WIN_START = 34_700_000
WIN_END = 35_950_000

up = pd.read_csv(OUT / "DAR_HOM_vs_WT_up.tsv", sep="\t")
dn = pd.read_csv(OUT / "DAR_HOM_vs_WT_down.tsv", sep="\t")
up["dir"] = "UP"
dn["dir"] = "DOWN"
all_dar = pd.concat([up, dn], ignore_index=True)

locus = all_dar[(all_dar["seqnames"] == "chr7")
                & (all_dar["start"] >= WIN_START)
                & (all_dar["end"] <= WIN_END)].copy()
locus["mid"] = (locus["start"] + locus["end"]) / 2
locus["-log10FDR"] = -pd.Series(locus["FDR"]).apply(lambda x: pd.np.log10(x) if hasattr(pd, "np") else __import__("math").log10(x))

fig, ax = plt.subplots(figsize=(11, 4.5))

# TBX20 gene body
ax.axvspan(TBX20_START, TBX20_END, color="#ffd58a", alpha=0.5, zorder=0, label="TBX20 gene body")

# 0 line
ax.axhline(0, color="grey", lw=0.5)

# points
for _, r in locus.iterrows():
    c = "#d62728" if r["dir"] == "UP" else "#1f77b4"
    ax.scatter(r["mid"], r["Log2FC"], s=80, c=c, edgecolor="black", lw=0.5, zorder=3)
    ax.annotate(f"FDR={r['FDR']:.1e}",
                xy=(r["mid"], r["Log2FC"]),
                xytext=(5, 5), textcoords="offset points",
                fontsize=8)

# background: all peaks in window (grey dots, even non-significant if we had them)
# Note: we only have filtered DAR here, so just show what's there

ax.set_xlim(WIN_START, WIN_END)
ax.set_xlabel("chr7 position (bp)")
ax.set_ylabel("Log2FC (HOM vs WT)")
ax.set_title("TBX20 locus DARs (±500 kb) — HOM vs WT")
ax.ticklabel_format(style="plain", axis="x")
ax.grid(alpha=0.3)

up_patch = mpatches.Patch(color="#d62728", label="UP in HOM")
dn_patch = mpatches.Patch(color="#1f77b4", label="DOWN in HOM")
gene_patch = mpatches.Patch(color="#ffd58a", alpha=0.5, label="TBX20 gene body")
ax.legend(handles=[up_patch, dn_patch, gene_patch], loc="upper left")

plt.tight_layout()
out_png = OUT / "TBX20_locus_DAR_scatter.png"
plt.savefig(out_png, dpi=150)
print(f"saved: {out_png}")
print(locus[["seqnames", "start", "end", "Log2FC", "FDR", "dir"]].to_string(index=False))
