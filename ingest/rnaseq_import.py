"""RNA-seq ingest (MVP).

This module is intentionally minimal and auditable.

What it does
- Load a count matrix (genes × samples) OR (samples × genes) and normalize orientation.
- Load a sample metadata table.
- Compute a simple 'persister signature score' per sample (z-scored across samples).
- Build an index mapping a CaseSnapshot-style key -> expression_summary dict.

What it does NOT do
- Differential expression (DESeq2), batch correction, or clinical interpretation.

Scoring
- For a set of marker genes, compute per-gene z-scores across samples.
- The sample score is the mean z-score across signature genes.

Why z-score?
- It produces roughly standardized values (~-2..+2) where a default threshold like 1.0 is meaningful.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Default signature file (one gene per line; '#' comments allowed)
DEFAULT_SIGNATURE_PATH = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "persistence_signatures"
    / "example_persister_genes.txt"
)


def _load_signature_genes(path: Optional[str] = None) -> List[str]:
    p = Path(path) if path else DEFAULT_SIGNATURE_PATH
    if not p.exists():
        return []
    genes: List[str] = []
    for raw in p.read_text(encoding="utf-8").splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        genes.append(s)
    # de-dup while preserving order
    seen = set()
    out: List[str] = []
    for g in genes:
        if g not in seen:
            out.append(g)
            seen.add(g)
    return out


def _ensure_counts_orientation(df: pd.DataFrame) -> pd.DataFrame:
    """Return counts as a DataFrame with columns: gene + sample_id columns.

    Accepts:
    - genes × samples (first column is gene id)
    - samples × genes (first column is sample_id)
    """
    if df.shape[1] < 2:
        raise ValueError("Counts matrix must have at least 2 columns.")

    first_col = df.columns[0]
    # Heuristic: if first column values look like sample IDs, transpose.
    head_vals = df[first_col].astype(str).head(5).tolist()
    looks_like_sample_ids = any(
        v.startswith("TCGA-") or v.startswith("SAMPLE") or v.startswith("SRR") for v in head_vals
    )

    if looks_like_sample_ids:
        # samples × genes -> transpose to genes × samples
        df_t = df.set_index(first_col).T.reset_index()
        df_t = df_t.rename(columns={"index": "gene"})
        return df_t

    # genes × samples -> normalize gene column name
    if first_col != "gene":
        df = df.rename(columns={first_col: "gene"})
    return df


def load_rnaseq_counts(path: str, *, delimiter: str = "\t") -> pd.DataFrame:
    df = pd.read_csv(path, sep=delimiter, dtype=str, keep_default_na=False, comment="#")
    df = _ensure_counts_orientation(df)

    # Convert sample columns to numeric.
    for c in df.columns:
        if c == "gene":
            continue
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
    return df


def load_rnaseq_metadata(path: str, *, delimiter: str = "\t") -> pd.DataFrame:
    df = pd.read_csv(path, sep=delimiter, dtype=str, keep_default_na=False, comment="#")
    if "sample_id" not in df.columns:
        raise ValueError("RNA-seq metadata must include a 'sample_id' column.")
    return df


def compute_persister_signature_scores(
    counts_df: pd.DataFrame,
    *,
    signature_genes: Optional[List[str]] = None,
) -> Dict[str, float]:
    """Compute a persister signature score per sample.

    Returns {sample_id: score}. If scoring can't be computed, returns {}.
    """
    signature_genes = signature_genes or _load_signature_genes()
    if not signature_genes:
        return {}

    if "gene" not in counts_df.columns:
        raise ValueError("Counts dataframe must include a 'gene' column.")

    # counts_df: columns = gene + samples
    sub = counts_df[counts_df["gene"].astype(str).isin(signature_genes)].copy()
    if sub.empty:
        return {}

    sub = sub.set_index("gene")
    # gene x sample numeric matrix
    x = sub.values.astype(float)
    # per-gene stats across samples (rows)
    mu = x.mean(axis=1, keepdims=True)
    sigma = x.std(axis=1, keepdims=True)
    sigma[sigma == 0] = 1.0
    z = (x - mu) / sigma  # gene x sample

    # score per sample: mean z across genes
    scores = z.mean(axis=0)
    sample_ids = list(sub.columns)
    return {sid: float(scores[i]) for i, sid in enumerate(sample_ids)}


def build_expression_index(
    counts_path: str,
    meta_path: str,
    *,
    delimiter: str = "\t",
    signature_path: Optional[str] = None,
) -> Dict[Tuple[Optional[str], Optional[str], str, Optional[str]], Dict]:
    """Build an expression index keyed like (study_id, case_id, sample_id, timepoint_id).

    Also includes a fallback key (None, None, sample_id, None) so callers can match
    by sample_id only when higher-level IDs are missing.
    """
    counts = load_rnaseq_counts(counts_path, delimiter=delimiter)
    meta = load_rnaseq_metadata(meta_path, delimiter=delimiter)

    signature_genes = _load_signature_genes(signature_path)
    score_map = compute_persister_signature_scores(counts, signature_genes=signature_genes)

    expr_index: Dict[Tuple[Optional[str], Optional[str], str, Optional[str]], Dict] = {}

    for _, row in meta.iterrows():
        sample_id = str(row.get("sample_id")).strip()
        if not sample_id:
            continue
        study_id = row.get("study_id") if "study_id" in meta.columns else None
        case_id = row.get("case_id") if "case_id" in meta.columns else None
        timepoint_id = row.get("timepoint_id") if "timepoint_id" in meta.columns else None

        summary = {
            "has_expression_data": True,
            "persister_signature_score": score_map.get(sample_id),
            "signature_genes_present": signature_genes,
        }

        key_full = (study_id, case_id, sample_id, timepoint_id)
        expr_index[key_full] = summary

        # Smart fallback matching: sample_id-only key
        key_sample_only = (None, None, sample_id, None)
        expr_index.setdefault(key_sample_only, summary)

    return expr_index
