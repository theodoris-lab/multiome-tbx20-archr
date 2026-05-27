# Created: 2026-04-23 17:50
# Updated: 2026-04-23 17:50
# ==============================================================================
# Phase 8 of DAR_analysis_plan_rnaCM_0423 — A/B comparison between:
#   A) Original ATAC-CM DAR  (outputs/dar_0420/,         14,823 cells)
#   B) rnaCM-gated DAR      (outputs/dar_rnaCM_0423/,   16,789 cells)
#
# Deliverables:
#   1. CM cell-set overlap (Venn)
#   2. DAR count table (HOM/HET up/down both runs)
#   3. DAR peak coord overlap (Jaccard + venn counts) for HOM_up, HOM_down
#   4. Top-N motif Spearman rank-corr
#   5. TBX20 locus replication check
#
# Writes:
#   outputs/dar_rnaCM_0423/ab_comparison_summary.md
#   outputs/dar_rnaCM_0423/ab_comparison_cm_venn.png
#   outputs/dar_rnaCM_0423/ab_comparison_dar_overlap.png
# ==============================================================================

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

ROOT = Path("/Users/bkim/vscode/V2_multiome_2026-04-14")
A_DIR = ROOT / "outputs/dar_0420"              # original ATAC-CM
B_DIR = ROOT / "outputs/dar_rnaCM_0423"        # rnaCM-gated
OUT = B_DIR
SUMMARY = OUT / "ab_comparison_summary.md"


def log(s=""):
    print(s)
    lines.append(s)


lines: list[str] = []
log(f"# A/B comparison: ATAC-CM (original) vs rnaCM-gated DAR")
log(f"Generated: 2026-04-23\n")
log(f"- A (original) : `{A_DIR.relative_to(ROOT)}`")
log(f"- B (rnaCM)    : `{B_DIR.relative_to(ROOT)}`\n")


# ------------------------------------------------------------------ 1) CM cell overlap
log("## 1. CM cell-set overlap\n")
a_cells = set((B_DIR / "atacCM_cellNames.txt").read_text().splitlines())
b_cells = set((B_DIR / "rnaCM_cellNames.txt").read_text().splitlines())
inter = a_cells & b_cells
only_a = a_cells - b_cells
only_b = b_cells - a_cells

log(f"- ATAC-CM only: {len(only_a):,}")
log(f"- rnaCM only : {len(only_b):,}")
log(f"- both       : {len(inter):,}")
log(f"- Jaccard    : {len(inter)/len(a_cells | b_cells):.3f}")
log(f"- A recall in B : {len(inter)/len(a_cells):.1%}  "
    f"(fraction of original ATAC-CM retained by RNA gate)")
log(f"- B precision vs A : {len(inter)/len(b_cells):.1%}\n")

fig, ax = plt.subplots(figsize=(5, 4))
venn2([a_cells, b_cells], set_labels=("ATAC-CM", "rnaCM"), ax=ax)
ax.set_title("CM cell overlap")
fig.tight_layout()
fig.savefig(OUT / "ab_comparison_cm_venn.png", dpi=150)
plt.close(fig)


# ------------------------------------------------------------------ 2) DAR counts
log("## 2. DAR counts\n")


def count_rows(p: Path) -> int:
    if not p.exists():
        return 0
    return max(sum(1 for _ in p.open()) - 1, 0)


count_rows_table = pd.DataFrame(
    {
        "A_ATAC-CM": [
            count_rows(A_DIR / "DAR_HOM_vs_WT_up.tsv"),
            count_rows(A_DIR / "DAR_HOM_vs_WT_down.tsv"),
            count_rows(A_DIR / "DAR_HET_vs_WT_up.tsv"),
            count_rows(A_DIR / "DAR_HET_vs_WT_down.tsv"),
        ],
        "B_rnaCM": [
            count_rows(B_DIR / "DAR_HOM_vs_WT_up.tsv"),
            count_rows(B_DIR / "DAR_HOM_vs_WT_down.tsv"),
            count_rows(B_DIR / "DAR_HET_vs_WT_up.tsv"),
            count_rows(B_DIR / "DAR_HET_vs_WT_down.tsv"),
        ],
    },
    index=["HOM_up", "HOM_down", "HET_up", "HET_down"],
)
count_rows_table["B/A"] = (count_rows_table["B_rnaCM"] /
                            count_rows_table["A_ATAC-CM"].replace(0, np.nan)).round(3)
log(count_rows_table.to_markdown())
log("")


# ------------------------------------------------------------------ 3) DAR peak coord overlap
log("## 3. DAR peak coord overlap (HOM)\n")


def load_dar(p: Path) -> pd.DataFrame:
    if not p.exists():
        return pd.DataFrame(columns=["seqnames", "start", "end"])
    return pd.read_csv(p, sep="\t")[["seqnames", "start", "end"]].copy()


def overlap_counts(a: pd.DataFrame, b: pd.DataFrame):
    """Return (only_a, only_b, both_a, both_b) — peaks with any overlap."""
    if a.empty or b.empty:
        return len(a), len(b), 0, 0
    a = a.sort_values(["seqnames", "start"]).reset_index(drop=True)
    b = b.sort_values(["seqnames", "start"]).reset_index(drop=True)
    matched_a = np.zeros(len(a), dtype=bool)
    matched_b = np.zeros(len(b), dtype=bool)
    for chrom in a["seqnames"].unique():
        ai = np.where(a["seqnames"].values == chrom)[0]
        bi = np.where(b["seqnames"].values == chrom)[0]
        if len(bi) == 0:
            continue
        b_starts = b["start"].values[bi]
        b_ends = b["end"].values[bi]
        order = np.argsort(b_starts)
        b_starts_s = b_starts[order]
        b_ends_s = b_ends[order]
        bi_s = bi[order]
        for idx in ai:
            s, e = a["start"].iloc[idx], a["end"].iloc[idx]
            # candidates where b_start < e; take all whose b_end > s
            lo = np.searchsorted(b_starts_s, e, side="left")
            cand = np.arange(lo)
            cand = cand[b_ends_s[cand] > s]
            if cand.size:
                matched_a[idx] = True
                matched_b[bi_s[cand]] = True
    return (
        int(np.sum(~matched_a)),
        int(np.sum(~matched_b)),
        int(np.sum(matched_a)),
        int(np.sum(matched_b)),
    )


peak_summary = []
for label in ["HOM_vs_WT_up", "HOM_vs_WT_down"]:
    a = load_dar(A_DIR / f"DAR_{label}.tsv")
    b = load_dar(B_DIR / f"DAR_{label}.tsv")
    only_a, only_b, both_a, both_b = overlap_counts(a, b)
    total = only_a + only_b + both_a  # union (both_a and both_b count overlap from each side)
    jaccard = both_a / (only_a + only_b + both_a) if total > 0 else 0.0
    recall_a = both_a / len(a) if len(a) else 0.0
    recall_b = both_b / len(b) if len(b) else 0.0
    peak_summary.append(
        dict(
            contrast=label,
            A=len(a),
            B=len(b),
            A_only=only_a,
            B_only=only_b,
            overlap_from_A=both_a,
            overlap_from_B=both_b,
            jaccard=round(jaccard, 3),
            recall_A_in_B=round(recall_a, 3),
            recall_B_in_A=round(recall_b, 3),
        )
    )

peak_df = pd.DataFrame(peak_summary).set_index("contrast")
log(peak_df.to_markdown())
log("")
log("_interpretation: `recall_A_in_B` = what fraction of the original DAR peaks "
    "are recovered by the rnaCM run (peak coord ∩); Jaccard = shared / union._\n")

# Venn-ish bar plot
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for ax, row in zip(axes, peak_summary):
    vals = [row["A_only"], row["B_only"], row["overlap_from_A"]]
    ax.bar(["A only", "B only", "overlap"], vals,
           color=["#4c72b0", "#dd8452", "#55a868"])
    for i, v in enumerate(vals):
        ax.text(i, v, str(v), ha="center", va="bottom", fontsize=10)
    ax.set_title(row["contrast"])
    ax.set_ylabel("# peaks")
fig.tight_layout()
fig.savefig(OUT / "ab_comparison_dar_overlap.png", dpi=150)
plt.close(fig)


# ------------------------------------------------------------------ 4) Motif Spearman
log("## 4. Motif enrichment overlap (HOM)\n")


def load_motif(p: Path) -> pd.DataFrame:
    if not p.exists():
        return pd.DataFrame()
    return pd.read_csv(p, sep="\t")


for direction in ["up", "down"]:
    a = load_motif(A_DIR / f"motif_HOM_vs_WT_{direction}.tsv")
    b = load_motif(B_DIR / f"motif_HOM_vs_WT_{direction}.tsv")
    if a.empty or b.empty:
        log(f"- HOM {direction}: missing motif table, skipped.")
        continue
    merged = a.merge(b, on="TF", suffixes=("_A", "_B"))
    if len(merged) < 10:
        log(f"- HOM {direction}: only {len(merged)} shared TFs, skipping corr.")
        continue

    rho_p, pval_p = spearmanr(merged["mlog10Padj_A"], merged["mlog10Padj_B"])
    rho_e, _ = spearmanr(merged["Enrichment_A"], merged["Enrichment_B"])

    top_a = set(a.sort_values("mlog10Padj", ascending=False).head(20)["TF"])
    top_b = set(b.sort_values("mlog10Padj", ascending=False).head(20)["TF"])
    overlap20 = top_a & top_b
    top_a_50 = set(a.sort_values("mlog10Padj", ascending=False).head(50)["TF"])
    top_b_50 = set(b.sort_values("mlog10Padj", ascending=False).head(50)["TF"])

    log(f"### HOM_vs_WT {direction}")
    log(f"- Shared TFs (in both tables): **{len(merged)}**")
    log(f"- Spearman ρ (mlog10Padj): **{rho_p:.3f}** (p={pval_p:.2e})")
    log(f"- Spearman ρ (Enrichment): **{rho_e:.3f}**")
    log(f"- Top-20 overlap: **{len(overlap20)}/20**")
    log(f"- Top-50 overlap: **{len(top_a_50 & top_b_50)}/50**")
    log(f"- Top-20 ATAC-CM: {sorted(top_a)}")
    log(f"- Top-20 rnaCM  : {sorted(top_b)}")
    log(f"- Top-20 shared : {sorted(overlap20)}\n")


# ------------------------------------------------------------------ 5) TBX20 locus
log("## 5. TBX20 locus replication\n")
# TBX20 gene body ~chr7:35,240-35,295 kb; regulatory window ~chr7:35,000-36,000 kb
# From original report: HOM-up peaks flagged at chr7:35.42Mb and chr7:35.89Mb
TBX20_WINDOW = ("chr7", 35_000_000, 36_000_000)


def in_window(df: pd.DataFrame, chrom, lo, hi):
    if df.empty:
        return df
    return df[(df["seqnames"] == chrom) &
              (df["start"] >= lo) & (df["end"] <= hi)].copy()


for label in ["HOM_vs_WT_up", "HOM_vs_WT_down"]:
    a_all = load_dar(A_DIR / f"DAR_{label}.tsv")
    b_all = load_dar(B_DIR / f"DAR_{label}.tsv")
    a_tbx = in_window(a_all, *TBX20_WINDOW)
    b_tbx = in_window(b_all, *TBX20_WINDOW)
    log(f"### {label}  (chr7:35.0-36.0 Mb)")
    log(f"- A (ATAC-CM) peaks in window: {len(a_tbx)}")
    if len(a_tbx):
        log("  ```")
        for _, r in a_tbx.iterrows():
            log(f"  {r['seqnames']}:{r['start']:,}-{r['end']:,}")
        log("  ```")
    log(f"- B (rnaCM)  peaks in window: {len(b_tbx)}")
    if len(b_tbx):
        log("  ```")
        for _, r in b_tbx.iterrows():
            log(f"  {r['seqnames']}:{r['start']:,}-{r['end']:,}")
        log("  ```")
    log("")


# ------------------------------------------------------------------ write
SUMMARY.write_text("\n".join(lines) + "\n")
print(f"\nWrote {SUMMARY}")
