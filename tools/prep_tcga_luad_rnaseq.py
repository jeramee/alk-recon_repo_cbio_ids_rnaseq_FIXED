import numpy as np
import pandas as pd
from pathlib import Path

IN_PATH  = Path(r"data\3\Human__TCGA_LUAD__UNC__RNAseq__GA_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct.txt")
OUT_EXPR = Path(r"data\real\rnaseq_expression_tcga_luad.tsv")
OUT_MAP  = Path(r"data\real\sample_map.tcga_luad.tsv")

df = pd.read_csv(IN_PATH, sep="\t")

# First column is the gene identifier in these Firehose/LinkedOmics-style matrices
gene_col = df.columns[0]
df = df.rename(columns={gene_col: "gene"})

# Keep a small subset of samples for iteration speed
sample_cols = list(df.columns[1:101])  # first 100 samples
df = df[["gene", *sample_cols]]

# Convert log2(RSEM) -> pseudo-count-ish non-negative numbers
vals = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
vals = (np.power(2.0, vals) - 1.0)
vals = np.maximum(vals, 0.0)
df[sample_cols] = vals

OUT_EXPR.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(OUT_EXPR, sep="\t", index=False)

# Build sample_map
rows = []
for sid in sample_cols:
    sid = str(sid).strip()
    rows.append({
        "sample_id": sid,
        "case_id": sid[:12],          # TCGA patient barcode
        "study_id": "TCGA-LUAD",
        "timepoint_id": "T0",
        "source": "LINKEDOMICS"
    })

pd.DataFrame(rows).to_csv(OUT_MAP, sep="\t", index=False)

print("Wrote:", OUT_EXPR)
print("Wrote:", OUT_MAP)
