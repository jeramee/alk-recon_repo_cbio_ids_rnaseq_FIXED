"""Convert a LinkedOmics expression matrix into ALK-RECON rnaseq inputs.

LinkedOmics TCGA files (e.g., ...RSEM_log2.cct.txt / .gz) are already log2-scaled.
We don't need to "unlog" them; we just need consistent IDs and a clean TSV.

Usage (PowerShell):
  python tools\prep_linkedomics_expression.py \
    --in  data\3\Human__TCGA_LUAD__...RSEM_log2.cct.gz \
    --out-counts data\real\rnaseq_expression_tcga_luad.full.tsv \
    --out-meta   data\real\sample_map.tcga_luad.full.tsv \
    --study-id TCGA-LUAD \
    --timepoint-id T0 \
    --source LINKEDOMICS \
    --id-mode patient12

id-mode:
  - patient12: normalize sample columns to TCGA-XX-YYYY (first 12 chars)
  - full: keep the column names after '.'->'-' replacement
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd


def _read_any(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            return pd.read_csv(f, sep="\t")
    return pd.read_csv(path, sep="\t")


def _norm_id(s: str) -> str:
    return str(s).strip().replace(".", "-")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out-counts", required=True)
    ap.add_argument("--out-meta", required=True)
    ap.add_argument("--study-id", default="TCGA-LUAD")
    ap.add_argument("--timepoint-id", default="T0")
    ap.add_argument("--source", default="LINKEDOMICS")
    ap.add_argument("--id-mode", choices=["patient12", "full"], default="patient12")
    args = ap.parse_args()

    inp = Path(args.inp)
    out_counts = Path(args.out_counts)
    out_meta = Path(args.out_meta)
    out_counts.parent.mkdir(parents=True, exist_ok=True)
    out_meta.parent.mkdir(parents=True, exist_ok=True)

    df = _read_any(inp)
    if df.shape[1] < 3:
        raise SystemExit(f"Expression matrix looks too small: {inp} -> shape {df.shape}")

    # First column is gene symbol in LinkedOmics matrices
    gene_col = df.columns[0]
    df[gene_col] = df[gene_col].astype(str)
    df = df.rename(columns={gene_col: "gene"})

    # Normalize sample IDs from headers
    sample_cols = list(df.columns[1:])
    norm_cols = []
    for c in sample_cols:
        c2 = _norm_id(c)
        if args.id_mode == "patient12":
            c2 = c2[:12]
        norm_cols.append(c2)

    df.columns = ["gene"] + norm_cols

    # De-duplicate sample columns (common when collapsing patient12)
    # Keep mean across duplicates.
    # pandas groupby(axis=1) is deprecated; do transpose-groupby-transpose
    numeric = df.drop(columns=["gene"]).apply(pd.to_numeric, errors="coerce")
    numeric = numeric.T.groupby(level=0).mean().T
    out_df = pd.concat([df[["gene"]], numeric], axis=1)

    out_df.to_csv(out_counts, sep="\t", index=False)

    # Build sample map from columns
    rows = []
    for sid in list(numeric.columns):
        sid = str(sid)
        rows.append(
            {
                "sample_id": sid,
                "case_id": sid,
                "study_id": args.study_id,
                "timepoint_id": args.timepoint_id,
                "source": args.source,
            }
        )
    meta = pd.DataFrame(rows)
    meta.to_csv(out_meta, sep="\t", index=False)

    print("Wrote:", out_counts)
    print("Wrote:", out_meta)
    print("Counts shape:", out_df.shape, "samples:", len(numeric.columns))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
