# Created: 2026-04-17 13:57
"""
Phase 0: Extract valid barcodes + metadata for ArchR (Wynton/Hive).

Input:  data/processed_wtD0/metadata_wtD0.tsv.gz   (43,014 doublet-removed cells)
Output: scripts_hive/phase0_out/
    - valid_barcodes_G1.txt, G2.txt, G3.txt   (ATAC barcodes, one per line)
    - cell_metadata_for_archr.csv             (archr_cell_id = "G#<atac_barcode>")

Why ATAC barcodes (not GEX): ArchR's createArrowFiles() reads fragment files,
which are keyed by ATAC cell barcodes. In 10x Multiome, the GEX and ATAC
barcodes differ per-cell, so validBarcodes must be the ATAC side.

Run:
    python scripts_hive/00_extract_barcodes.py
"""
from pathlib import Path
import pandas as pd

PROJ = Path("/Users/bkim/vscode/V2_multiome_2026-04-14")
META = PROJ / "data/processed_wtD0/metadata_wtD0.tsv.gz"
OUT = PROJ / "scripts_hive/phase0_out"
OUT.mkdir(parents=True, exist_ok=True)

KEEP_COLS = [
    "arc_library", "multiseq_group", "multiseq_sample",
    "gex_barcode", "atac_barcode",
    "seurat_clusters", "wsnn_res.1",
]

df = pd.read_csv(META, sep="\t")
print(f"Loaded metadata: {df.shape[0]:,} cells × {df.shape[1]} columns")

# Sanity checks
assert df.shape[0] == 43014, f"Expected 43,014 cells, got {df.shape[0]}"
assert df["arc_library"].isin(["G1", "G2", "G3"]).all()
assert df["atac_barcode"].str.endswith("-1").all(), "ATAC barcodes must end in -1"
assert df["multiseq_group"].notna().all(), "multiseq_group has NaNs"

# archr_cell_id matches ArchR's cellNames: "G1#ACGTACGT-1"
df["archr_cell_id"] = df["arc_library"].astype(str) + "#" + df["atac_barcode"].astype(str)
assert df["archr_cell_id"].is_unique, "archr_cell_id collision"

# Per-sample valid barcodes (raw ATAC barcode, with -1 suffix)
for sample in ["G1", "G2", "G3"]:
    bcs = df.loc[df["arc_library"] == sample, "atac_barcode"].tolist()
    out_path = OUT / f"valid_barcodes_{sample}.txt"
    out_path.write_text("\n".join(bcs) + "\n")
    print(f"  {sample}: {len(bcs):>6,} barcodes → {out_path.name}")

# Metadata CSV for ArchR annotation
out_meta = df[["archr_cell_id"] + KEEP_COLS].copy()
out_meta_path = OUT / "cell_metadata_for_archr.csv"
out_meta.to_csv(out_meta_path, index=False)
print(f"  metadata: {len(out_meta):,} rows → {out_meta_path.name}")

# Summary tables
print("\n=== arc_library × multiseq_group ===")
print(pd.crosstab(df["arc_library"], df["multiseq_group"], margins=True))

print("\n=== arc_library × multiseq_sample ===")
print(pd.crosstab(df["arc_library"], df["multiseq_sample"]))

print("\nDone.")
