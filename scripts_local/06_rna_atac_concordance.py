# Created: 2026-04-24 09:30
# Updated: 2026-04-24 09:55
"""
RNA-ATAC concordance via proximity-based peak2gene linking.

Input:
  - outputs/dar_rnaCM_0423/DAR_HOM_vs_WT_{up,down}.tsv  (peak coords + Log2FC/FDR)
  - outputs/deg_4res_compare_0422/res0.1_c0_7/DEG_D{10,15,30}_HOM_vs_WT.csv (mature CM)
  - downloads/refFlat_hg38.txt.gz  (UCSC refFlat for TSS coords)

Method:
  For each DAR peak, assign gene via: (1) gene body containment takes
  priority, (2) else nearest TSS within ±100 kb. If multiple gene bodies
  contain the peak, pick the one whose TSS is closest. Collapse refFlat
  transcripts per gene to canonical TSS (most 5' per strand).
  Aggregate DEG across D10/D15/D30 via mean log2FC.
  Merge peak × gene × RNA log2FC. Plot concordance scatter.

Output:
  outputs/concordance_rnaCM_0424/
    peak2gene_nearest.tsv        — per-peak gene link
    concordance_scatter.{png,pdf} — ATAC l2fc (x) vs RNA l2fc (y)
    meis1_locus_peaks.tsv         — MEIS1 locus peaks + linked genes
"""
from pathlib import Path
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from adjustText import adjust_text

ROOT = Path("/Users/bkim/vscode/V2_multiome_2026-04-14")
DAR_DIR = ROOT / "outputs/dar_rnaCM_0423"
DEG_DIR = ROOT / "outputs/deg_4res_compare_0422/res0.1_c0_7"
REFFLAT = ROOT / "downloads/refFlat_hg38.txt.gz"
LVNC_CSV = ROOT / "downloads/LVNC_all_genes.csv"
OUT = ROOT / "outputs/concordance_rnaCM_0424"
OUT.mkdir(parents=True, exist_ok=True)

WINDOW = 100_000  # ±100 kb around peak midpoint
DEG_TIMEPOINTS = ["D10", "D15", "D30"]  # mature CM window

# ---------- 1. DAR peaks ----------
up = pd.read_csv(DAR_DIR / "DAR_HOM_vs_WT_up.tsv", sep="\t")
dn = pd.read_csv(DAR_DIR / "DAR_HOM_vs_WT_down.tsv", sep="\t")
up["dir"] = "UP"
dn["dir"] = "DOWN"
dar = pd.concat([up, dn], ignore_index=True)
dar["mid"] = ((dar["start"] + dar["end"]) // 2).astype(int)
dar["peak_id"] = dar["seqnames"] + ":" + dar["start"].astype(str) + "-" + dar["end"].astype(str)
print(f"DAR peaks: {len(dar)} (UP {len(up)}, DOWN {len(dn)})")

# ---------- 2. refFlat → per-gene canonical TSS + body ----------
refflat_cols = ["gene", "refseq", "chrom", "strand",
                "txStart", "txEnd", "cdsStart", "cdsEnd",
                "exonCount", "exonStarts", "exonEnds"]
with gzip.open(REFFLAT, "rt") as fh:
    rf = pd.read_csv(fh, sep="\t", names=refflat_cols,
                     usecols=["gene", "chrom", "strand", "txStart", "txEnd"])
rf = rf[rf["chrom"].str.match(r"^chr[0-9XYM]+$")]
# canonical per gene: 5'-most TSS + outermost body
def canonical_row(g):
    strand = g["strand"].iloc[0]
    return pd.Series({
        "chrom": g["chrom"].iloc[0],
        "strand": strand,
        "tss": g["txStart"].min() if strand == "+" else g["txEnd"].max(),
        "body_start": int(g["txStart"].min()),
        "body_end": int(g["txEnd"].max()),
    })
canon = rf.groupby("gene", group_keys=False).apply(canonical_row).reset_index()
print(f"refFlat genes: {canon['gene'].nunique()}")

# ---------- 3. peak2gene: body containment > nearest TSS ----------
by_chrom_body = {c: g.sort_values("body_start").reset_index(drop=True)
                 for c, g in canon.groupby("chrom")}
by_chrom_tss = {c: g.sort_values("tss").reset_index(drop=True)
                for c, g in canon.groupby("chrom")}

def assign_gene(peak_chr, peak_start, peak_end, peak_mid, window=WINDOW):
    if peak_chr not in by_chrom_body:
        return None, None, None, None
    sub = by_chrom_body[peak_chr]
    # gene body containment (peak midpoint inside gene body)
    inside = sub[(sub["body_start"] <= peak_mid) & (sub["body_end"] >= peak_mid)]
    if not inside.empty:
        # pick gene whose TSS is closest to peak_mid
        inside = inside.assign(_d=(inside["tss"] - peak_mid).abs())
        best = inside.sort_values("_d").iloc[0]
        return best["gene"], int(best["tss"]), int(peak_mid - best["tss"]), "body"
    # else nearest TSS within window
    tss_sub = by_chrom_tss[peak_chr]
    idx = tss_sub["tss"].searchsorted(peak_mid)
    cands = []
    if idx > 0:
        cands.append(tss_sub.iloc[idx - 1])
    if idx < len(tss_sub):
        cands.append(tss_sub.iloc[idx])
    best = None
    best_d = None
    for r in cands:
        d = peak_mid - r["tss"]
        if abs(d) <= window and (best_d is None or abs(d) < abs(best_d)):
            best = r
            best_d = d
    if best is None:
        return None, None, None, None
    return best["gene"], int(best["tss"]), int(best_d), "tss"

assigned = dar.apply(
    lambda r: pd.Series(assign_gene(r["seqnames"], r["start"], r["end"], r["mid"])), axis=1)
assigned.columns = ["gene", "gene_tss", "peak_to_tss", "link_type"]
dar = pd.concat([dar, assigned], axis=1)
linked = dar.dropna(subset=["gene"]).copy()
print(f"peaks assigned (body+tss): {len(linked)}/{len(dar)}")
print(f"  by gene body: {(linked['link_type']=='body').sum()}")
print(f"  by nearest TSS: {(linked['link_type']=='tss').sum()}")

# ---------- 4. DEG aggregation ----------
deg_frames = []
for tp in DEG_TIMEPOINTS:
    f = DEG_DIR / f"DEG_{tp}_HOM_vs_WT.csv"
    t = pd.read_csv(f)[["names", "logfoldchanges", "pvals_adj"]]
    t["tp"] = tp
    deg_frames.append(t)
deg_all = pd.concat(deg_frames)
deg_mean = deg_all.groupby("names", as_index=False).agg(
    rna_l2fc_mean=("logfoldchanges", "mean"),
    rna_l2fc_max_abs=("logfoldchanges", lambda x: x.loc[x.abs().idxmax()]),
    rna_fdr_min=("pvals_adj", "min"),
    n_tp_sig=("pvals_adj", lambda x: (x < 0.05).sum()),
)

# ---------- 5. merge peak × gene × RNA ----------
merged = linked.merge(deg_mean, left_on="gene", right_on="names", how="left")
merged["rna_l2fc_mean"] = merged["rna_l2fc_mean"].fillna(0)
merged["has_rna"] = merged["names"].notna()
merged.to_csv(OUT / "peak2gene_nearest.tsv", sep="\t", index=False)
print(f"linked peaks w/ RNA expression: {merged['has_rna'].sum()}/{len(merged)}")

# concordance categories
def cat(row):
    a = row["Log2FC"]
    r = row["rna_l2fc_mean"]
    if not row["has_rna"]:
        return "no_RNA"
    if abs(r) < 0.25:
        return "RNA_flat"
    same = (a > 0) == (r > 0)
    return "concordant" if same else "discordant"
merged["concordance"] = merged.apply(cat, axis=1)
print("\nConcordance breakdown:")
print(merged["concordance"].value_counts())

# ---------- 6. MEIS1 locus dump ----------
meis1_peaks = merged[merged["gene"] == "MEIS1"].copy()
meis1_peaks.to_csv(OUT / "meis1_locus_peaks.tsv", sep="\t", index=False)
print(f"\nMEIS1 locus peaks: {len(meis1_peaks)}")
if len(meis1_peaks):
    print(meis1_peaks[["seqnames", "start", "end", "Log2FC", "FDR", "dir",
                       "peak_to_tss", "rna_l2fc_mean", "rna_fdr_min"]].to_string(index=False))

# ---------- 7. concordance scatter ----------
# load LVNC tiers
lvnc = pd.read_csv(LVNC_CSV)
lvnc["gene"] = lvnc["gene"].replace({"NKX2.5": "NKX2-5"})
lvnc_map = dict(zip(lvnc["gene"], lvnc["category"]))
merged["lvnc_tier"] = merged["gene"].map(lvnc_map)

fig, ax = plt.subplots(figsize=(9, 8))
plot_df = merged[merged["has_rna"]].copy()
# thin bg (|rna|<0.25, |atac|<0.5)
bg = plot_df[(plot_df["rna_l2fc_mean"].abs() < 0.25) & (plot_df["Log2FC"].abs() < 0.5)]
fg = plot_df.drop(bg.index)

ax.scatter(bg["Log2FC"], bg["rna_l2fc_mean"], s=6, c="#d0d0d0", alpha=0.5, linewidths=0, zorder=1)

for tier, col, sz in [(None, "#555", 12), ("Limited", "#bdbdbd", 20),
                      ("Moderate", "#ff9f1c", 40), ("Definitive", "#d62728", 60)]:
    if tier is None:
        sub = fg[fg["lvnc_tier"].isna()]
    else:
        sub = fg[fg["lvnc_tier"] == tier]
    if sub.empty:
        continue
    ax.scatter(sub["Log2FC"], sub["rna_l2fc_mean"], s=sz, c=col,
               edgecolor="black" if tier in ("Definitive", "Moderate") else "none",
               linewidths=0.4, alpha=0.85,
               zorder=4 if tier == "Definitive" else 3 if tier == "Moderate" else 2,
               label=f"LVNC {tier}" if tier else "other fg")

# reference axes
ax.axhline(0, color="grey", lw=0.5)
ax.axvline(0, color="grey", lw=0.5)
lim_a = max(plot_df["Log2FC"].abs().max(), 3) + 0.5
lim_r = max(plot_df["rna_l2fc_mean"].abs().max(), 3) + 0.5
ax.set_xlim(-lim_a, lim_a)
ax.set_ylim(-lim_r, lim_r)

# shade concordant quadrants
ax.axhspan(0, lim_r, xmin=0.5, xmax=1, alpha=0.05, color="green", zorder=0)
ax.axhspan(-lim_r, 0, xmin=0, xmax=0.5, alpha=0.05, color="green", zorder=0)

# MEIS1 highlight
meis1_row = plot_df[plot_df["gene"] == "MEIS1"]
if not meis1_row.empty:
    for _, r in meis1_row.iterrows():
        ax.scatter(r["Log2FC"], r["rna_l2fc_mean"], s=180,
                   facecolors="none", edgecolors="blue", linewidths=1.5, zorder=6)

# labels: Definitive + Moderate + MEIS1 (always)
texts = []
to_label = plot_df[(plot_df["lvnc_tier"].isin(["Definitive", "Moderate"])) |
                    (plot_df["gene"].isin(["MEIS1", "MEIS2", "MEIS3", "TGIF1", "PKNOX1", "GATA4", "NKX2-5"]))]
# dedupe gene × keep strongest peak
to_label = to_label.sort_values(
    "Log2FC", key=lambda s: s.abs(), ascending=False
).drop_duplicates("gene")
for _, r in to_label.iterrows():
    col = "#d62728" if r["lvnc_tier"] == "Definitive" else \
          "#ff9f1c" if r["lvnc_tier"] == "Moderate" else "black"
    texts.append(ax.text(r["Log2FC"], r["rna_l2fc_mean"], r["gene"],
                         fontsize=8, color=col, fontweight="bold"))
if texts:
    adjust_text(texts, ax=ax,
                only_move={"text": "xy"},
                arrowprops=dict(arrowstyle="-", color="black", lw=0.3, alpha=0.5),
                expand=(1.05, 1.15), force_text=(0.3, 0.4), time_lim=1.0)

ax.set_xlabel("ATAC Log2FC  (HOM vs WT, rnaCM peak)")
ax.set_ylabel("RNA Log2FC  (HOM vs WT, mean of D10/D15/D30)")
ax.set_title("RNA-ATAC concordance — proximity peak2gene (±100 kb)")
ax.grid(alpha=0.25)

legend = [
    mpatches.Patch(color="#d62728", label="LVNC Definitive"),
    mpatches.Patch(color="#ff9f1c", label="LVNC Moderate"),
    mpatches.Patch(color="#bdbdbd", label="LVNC Limited"),
    mpatches.Patch(color="#d0d0d0", label="bg (|RNA|<.25 & |ATAC|<.5)"),
    mpatches.Patch(edgecolor="blue", facecolor="none", label="MEIS1 peaks"),
]
ax.legend(handles=legend, loc="lower right", fontsize=8)

fig.tight_layout()
for ext in ["png", "pdf"]:
    fp = OUT / f"concordance_scatter.{ext}"
    fig.savefig(fp, dpi=200 if ext == "png" else None, bbox_inches="tight")
    print(f"saved: {fp}")

# ---------- 8. summary print ----------
print("\n=== Concordance summary ===")
print("Total linked peaks :", len(plot_df))
print("Concordant         :", (merged["concordance"] == "concordant").sum())
print("Discordant         :", (merged["concordance"] == "discordant").sum())
print("RNA-flat           :", (merged["concordance"] == "RNA_flat").sum())
print("\nTop discordant peaks (|ATAC|>1, |RNA|<0.25):")
disc = plot_df[(plot_df["Log2FC"].abs() > 1) &
               (plot_df["rna_l2fc_mean"].abs() < 0.25)]
print(disc[["peak_id", "gene", "Log2FC", "rna_l2fc_mean", "dir"]]
      .sort_values("Log2FC", key=lambda s: s.abs(), ascending=False)
      .head(15).to_string(index=False))
