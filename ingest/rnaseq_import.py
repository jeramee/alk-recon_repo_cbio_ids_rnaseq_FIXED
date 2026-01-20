from __future__ import annotations

"""RNA-seq counts + metadata import (research/education only).

This module is intentionally conservative:
- It does NOT attempt to replicate DESeq2/edgeR.
- It provides a simple, auditable pathway to compute *signature scores*
  (e.g., a generic 'persister_score') from a counts matrix.

The output is designed to attach to :class:`schema.case_snapshot.ExpressionSummary`.

Input conventions
-----------------
Counts matrix (CSV/TSV):
  - Preferred: rows = genes, columns = sample IDs.
  - Alternate: rows = sample IDs, columns = genes (auto-detected and transposed).
  - One column must represent gene names (default auto-detect: first column or 'gene').

Metadata (CSV/TSV):
  - Must include a sample identifier column.
  - May include: case_id, study_id, timepoint_id.

This module is not medical advice and is not intended for clinical use.
"""

from dataclasses import asdict
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import numpy as np
import pandas as pd

from schema.case_snapshot import ExpressionSummary
import math
...
top = []
for gene, val in vec.items():
    try:
        f = float(val)
    except:
        continue
    if math.isfinite(f):
        top.append((gene, f))

top = sorted(top, key=lambda x: x[1], reverse=True)[:top_n]

def _read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".tab"}:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def _auto_pick_col(df: pd.DataFrame, candidates: Iterable[str]) -> Optional[str]:
    cols = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols:
            return cols[cand.lower()]
    return None


def _load_signature_genes(path: str | Path | None) -> list[str]:
    """Return a newline-delimited gene list. If path is None, use the bundled example list."""
    if path is None:
        bundled = Path(__file__).resolve().parents[1] / "data" / "persistence_signatures" / "example_persister_genes.txt"
        path = bundled
    p = Path(path)
    genes: list[str] = []
    for line in p.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        genes.append(s)
    # unique, keep order
    seen = set()
    out = []
    for g in genes:
        if g not in seen:
            out.append(g)
            seen.add(g)
    return out


def _ensure_counts_orientation(counts: pd.DataFrame, meta_sample_ids: set[str]) -> Tuple[pd.DataFrame, str]:
    """Ensure counts columns are sample IDs; transpose if needed.

    Returns (counts_df, orientation_note).
    """
    # If any metadata sample IDs match columns, we assume columns are samples.
    if meta_sample_ids and len(set(counts.columns).intersection(meta_sample_ids)) >= max(1, min(3, len(meta_sample_ids))):
        return counts, "genes_by_samples"

    # If sample IDs match the index, assume rows are samples and transpose.
    if meta_sample_ids and len(set(map(str, counts.index)).intersection(meta_sample_ids)) >= max(1, min(3, len(meta_sample_ids))):
        return counts.T, "samples_by_genes_transposed"

    # Fallback: keep as-is.
    return counts, "unknown_assumed_genes_by_samples"


def compute_signature_scores(
    counts_gene_by_sample: pd.DataFrame,
    signature_genes: list[str],
    top_k_markers: int = 10,
) -> Tuple[pd.Series, pd.DataFrame, float]:
    """Compute a simple signature score per sample.

    Pipeline:
      1) CPM normalize
      2) log2(CPM+1)
      3) z-score each gene across samples
      4) per sample, score = mean(z) across signature genes present

    Returns:
      scores: pd.Series indexed by sample_id
      marker_z: pd.DataFrame with z-scores for signature genes (rows=gene, cols=sample)
      coverage: fraction of signature genes present in matrix
    """
    # Coerce numeric, drop junk, disallow negatives
    c = counts_gene_by_sample.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    c = c.clip(lower=0.0)

    # CPM normalize
    lib = c.sum(axis=0).replace(0, np.nan)
    cpm = c.div(lib, axis=1) * 1e6

    # log2(CPM+1) — force a DataFrame for both runtime + Pylance typing
    log2cpm = pd.DataFrame(
        np.log2(cpm.to_numpy(dtype=float) + 1.0),
        index=cpm.index,
        columns=cpm.columns,
    )


    # Standardize per gene across samples (pandas skips NaNs by default)
    mu = log2cpm.mean(axis=1)
    sd = log2cpm.std(axis=1)

    # Avoid divide-by-zero and all-NaN genes
    sd = sd.replace(0, 1e-6).fillna(1.0)

    z = log2cpm.sub(mu, axis=0).div(sd, axis=0)

    genes_present = [g for g in signature_genes if g in z.index]
    coverage = (len(genes_present) / max(1, len(signature_genes)))
    if not genes_present:
        scores = pd.Series({sid: float("nan") for sid in z.columns})
        marker_z = pd.DataFrame(index=[], columns=z.columns)
        return scores, marker_z, coverage

    marker_z = z.loc[genes_present]
    scores = marker_z.mean(axis=0)
    return scores, marker_z, coverage


def build_expression_index(
        counts_path: str | Path,
        metadata_path: str | Path,
        gene_col: str | None = None,
        sample_id_col: str | None = None,
        case_id_col: str | None = None,
        study_id_col: str | None = None,
        timepoint_id_col: str | None = None,
        signature_path: str | Path | None = None,
        top_k_markers: int = 10,   # <--- ADD THIS
    ) -> Dict[Tuple[Optional[str], Optional[str], str, Optional[str]], ExpressionSummary]:
    """Load RNA-seq counts + metadata and return an index to attach to CaseSnapshots.

    Key:
      (study_id, case_id, sample_id, timepoint_id)

    Notes:
      - study_id/case_id/timepoint_id may be None in public data; sample_id is required.
      - You can attach by exact key match, or by falling back to (None, case_id, sample_id, None).
    """
    counts_raw = _read_table(counts_path)
    meta = _read_table(metadata_path)

    # Auto-detect ID columns in metadata
    sample_id_col = sample_id_col or _auto_pick_col(meta, ["sample_id", "sample", "sampleid", "sampleId"])
    if not sample_id_col:
        raise ValueError("RNA-seq metadata must include a sample_id column (or provide --rnaseq-sample-id-col).")

    case_id_col = case_id_col or _auto_pick_col(meta, ["case_id", "case", "patient_id", "patient", "patientid", "patientId"])
    study_id_col = study_id_col or _auto_pick_col(meta, ["study_id", "study", "studyId"])
    timepoint_id_col = timepoint_id_col or _auto_pick_col(meta, ["timepoint_id", "timepoint", "time", "visit"])

    meta[sample_id_col] = meta[sample_id_col].astype(str)
    meta_sample_ids = set(meta[sample_id_col].tolist())

    # Auto-detect gene column in counts
    gene_col = gene_col or _auto_pick_col(counts_raw, ["gene", "gene_id", "gene_name", "symbol", "Gene"])
    if gene_col is None:
        gene_col = counts_raw.columns[0]

    counts = counts_raw.copy()
    counts[gene_col] = counts[gene_col].astype(str)
    counts = counts.set_index(gene_col)

    # Ensure columns are samples (transpose if needed)
    counts, orientation = _ensure_counts_orientation(counts, meta_sample_ids)

    # Keep only samples present in metadata
    shared_samples = [s for s in counts.columns if str(s) in meta_sample_ids]
    if not shared_samples:
        raise ValueError(
            "No overlapping sample IDs between counts matrix and metadata. "
            "Check that your counts columns (or rows) match metadata sample IDs."
        )
    counts = counts[shared_samples]

    signature_genes = _load_signature_genes(signature_path)
    scores, marker_z, coverage = compute_signature_scores(counts, signature_genes)

    # Build index
    index: Dict[Tuple[Optional[str], Optional[str], str, Optional[str]], ExpressionSummary] = {}
    for _, row in meta.iterrows():
        sid = str(row[sample_id_col])
        if sid not in scores.index:
            continue
        csid = str(row[case_id_col]) if case_id_col and not pd.isna(row[case_id_col]) else None
        stid = str(row[study_id_col]) if study_id_col and not pd.isna(row[study_id_col]) else None
        tpid = str(row[timepoint_id_col]) if timepoint_id_col and not pd.isna(row[timepoint_id_col]) else None

        s = float(scores.loc[sid]) if np.isfinite(scores.loc[sid]) else float("nan")
        # Top markers: signature genes with highest z in this sample
        top: list[dict] = []
        if sid in marker_z.columns and len(marker_z.index) > 0:
            vec = marker_z[sid].sort_values(ascending=False).head(10)
            top = [{"gene": g, "z": float(v)} for g, v in vec.items()]

        expr = ExpressionSummary(
            platform="rnaseq_counts",
            contrast=None,
            signature_scores={
                "persister_score": s,
                "signature_coverage": float(coverage),
            },
            top_markers=top,
        )
        index[(stid, csid, sid, tpid)] = expr

        # Add a few "alias" keys to make attachment robust when a downstream input
        # (e.g., a variant table) lacks study_id or timepoint_id.
        #
        # NOTE: This is an intentional trade-off for MVP ergonomics. If you ingest
        # multiple studies that share case/sample IDs, you should include study_id
        # in ALL inputs and avoid relying on these aliases.
        if stid is not None:
            index[("unspecified", csid, sid, tpid)] = expr
            index[(None, csid, sid, tpid)] = expr
        if tpid is not None:
            index[(stid, csid, sid, "")] = expr
            index[("unspecified", csid, sid, "")] = expr

        # Another common convention is using the literal string "unspecified"
        # when an upstream file lacks a timepoint. Provide that alias as well.
        index[(stid, csid, sid, "unspecified")] = expr
        index[("unspecified", csid, sid, "unspecified")] = expr

    return index
