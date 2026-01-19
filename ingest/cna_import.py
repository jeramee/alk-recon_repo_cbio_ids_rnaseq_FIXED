from __future__ import annotations

"""Copy-number (CNA) import from cBioPortal-style matrices.

Supports:
  - data_cna.txt (GISTIC thresholded, typically -2..2)
  - data_linear_cna.txt (continuous copy-ratio like log2)

Design goals
------------
* Conservative: we do not attempt full GISTIC interpretation or segmentation.
* Auditable: outputs are simple "events" attached to CaseSnapshot.genomic.
* Practical: works with the same sample_map TSV you already use for RNA-seq.

Output
------
Returns an index:
  (study_id, case_id, sample_id, timepoint_id) -> list[dict]

Each dict has:
  gene, kind ("amp"|"del"|"gain"|"loss"), value (float), source

This module is not medical advice and is not intended for clinical use.
"""

from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import numpy as np
import pandas as pd


def _read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".tab", ".txt"}:
        return pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    return pd.read_csv(path, comment="#", low_memory=False)


def _auto_pick_col(df: pd.DataFrame, candidates: Iterable[str]) -> Optional[str]:
    cols = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols:
            return cols[cand.lower()]
    return None


def _normalize_id(s: str) -> str:
    # LinkedOmics often uses dots; cBioPortal uses hyphens.
    return str(s).strip().replace(".", "-")


def _to_patient_barcode(sample_id: str) -> str:
    # TCGA patient barcode is first 12 chars, e.g. TCGA-44-6146
    sid = _normalize_id(sample_id)
    return sid[:12]


def _load_sample_map(metadata_path: str | Path,
                     sample_id_col: str | None,
                     case_id_col: str | None,
                     study_id_col: str | None,
                     timepoint_id_col: str | None) -> pd.DataFrame:
    meta = _read_table(metadata_path)

    sample_id_col = sample_id_col or _auto_pick_col(meta, ["sample_id", "sample", "sampleid", "sampleId"])
    if not sample_id_col:
        raise ValueError("CNA metadata must include a sample_id column (or provide --cna-sample-id-col).")

    case_id_col = case_id_col or _auto_pick_col(meta, ["case_id", "case", "patient_id", "patient"])
    study_id_col = study_id_col or _auto_pick_col(meta, ["study_id", "study"])
    timepoint_id_col = timepoint_id_col or _auto_pick_col(meta, ["timepoint_id", "timepoint", "time", "visit"])

    # Normalize IDs
    meta[sample_id_col] = meta[sample_id_col].astype(str).map(_normalize_id)
    if case_id_col:
        meta[case_id_col] = meta[case_id_col].astype(str).map(_normalize_id)
    if study_id_col:
        meta[study_id_col] = meta[study_id_col].astype(str)
    if timepoint_id_col:
        meta[timepoint_id_col] = meta[timepoint_id_col].astype(str)

    meta.attrs["sample_id_col"] = sample_id_col
    meta.attrs["case_id_col"] = case_id_col
    meta.attrs["study_id_col"] = study_id_col
    meta.attrs["timepoint_id_col"] = timepoint_id_col
    return meta


def _iter_cna_events_thresholded(
    df: pd.DataFrame,
    genes_of_interest: set[str] | None,
    source: str,
) -> Dict[str, list[dict]]:
    """Return sample_id -> events for thresholded CNA."""
    gene_col = _auto_pick_col(df, ["Hugo_Symbol", "gene", "Gene", "symbol"])
    if gene_col is None:
        gene_col = df.columns[0]

    df = df.copy()
    df[gene_col] = df[gene_col].astype(str)
    df = df.set_index(gene_col)

    # Drop the Entrez column if present
    if "Entrez_Gene_Id" in df.columns:
        df = df.drop(columns=["Entrez_Gene_Id"])

    # Numeric coercion; keep NaN (means missing)
    mat = df.apply(pd.to_numeric, errors="coerce")

    out: Dict[str, list[dict]] = {}
    for gene, row in mat.iterrows():
        if genes_of_interest and gene not in genes_of_interest:
            continue
        for sid, v in row.items():
            if pd.isna(v) or float(v) == 0.0:
                continue
            vv = float(v)
            kind = "amp" if vv >= 2 else "gain" if vv > 0 else "del" if vv <= -2 else "loss"
            out.setdefault(_normalize_id(sid), []).append(
                {"gene": str(gene), "kind": kind, "value": vv, "source": source}
            )
    return out


def _iter_cna_events_linear(
    df: pd.DataFrame,
    genes_of_interest: set[str] | None,
    source: str,
    amp_threshold: float = 0.7,
    del_threshold: float = -0.7,
) -> Dict[str, list[dict]]:
    """Return sample_id -> events for continuous CNA (log2 ratios, etc.)."""
    gene_col = _auto_pick_col(df, ["Hugo_Symbol", "gene", "Gene", "symbol"])
    if gene_col is None:
        gene_col = df.columns[0]

    df = df.copy()
    df[gene_col] = df[gene_col].astype(str)
    df = df.set_index(gene_col)
    if "Entrez_Gene_Id" in df.columns:
        df = df.drop(columns=["Entrez_Gene_Id"])

    mat = df.apply(pd.to_numeric, errors="coerce")

    out: Dict[str, list[dict]] = {}
    for gene, row in mat.iterrows():
        if genes_of_interest and gene not in genes_of_interest:
            continue
        for sid, v in row.items():
            if pd.isna(v) or float(v) == 0.0:
                continue
            vv = float(v)
            if vv >= amp_threshold:
                kind = "amp"
            elif vv <= del_threshold:
                kind = "del"
            else:
                continue
            out.setdefault(_normalize_id(sid), []).append(
                {"gene": str(gene), "kind": kind, "value": vv, "source": source}
            )
    return out


def build_cna_index(
    metadata_path: str | Path,
    thresholded_path: str | Path | None = None,
    linear_path: str | Path | None = None,
    sample_id_col: str | None = None,
    case_id_col: str | None = None,
    study_id_col: str | None = None,
    timepoint_id_col: str | None = None,
    genes_of_interest: set[str] | None = None,
) -> Dict[Tuple[Optional[str], Optional[str], str, Optional[str]], list[dict]]:
    """Build an index that can be attached to CaseSnapshots.

    Key: (study_id, case_id, sample_id, timepoint_id)

    We add a couple alias keys ("unspecified" and None) to match pipeline fallbacks.
    """
    if thresholded_path is None and linear_path is None:
        return {}

    meta = _load_sample_map(metadata_path, sample_id_col, case_id_col, study_id_col, timepoint_id_col)
    sid_col = meta.attrs["sample_id_col"]
    cid_col = meta.attrs["case_id_col"]
    stid_col = meta.attrs["study_id_col"]
    tpid_col = meta.attrs["timepoint_id_col"]

    by_sample: Dict[str, list[dict]] = {}
    if thresholded_path is not None:
        df = _read_table(thresholded_path)
        m = _iter_cna_events_thresholded(df, genes_of_interest, source=Path(thresholded_path).name)
        for sid, evs in m.items():
            by_sample.setdefault(sid, []).extend(evs)

    if linear_path is not None:
        df = _read_table(linear_path)
        m = _iter_cna_events_linear(df, genes_of_interest, source=Path(linear_path).name)
        for sid, evs in m.items():
            by_sample.setdefault(sid, []).extend(evs)

    index: Dict[Tuple[Optional[str], Optional[str], str, Optional[str]], list[dict]] = {}
    for _, row in meta.iterrows():
        sid = _normalize_id(row[sid_col])
        evs = by_sample.get(sid)
        if not evs:
            continue

        case_id = _normalize_id(row[cid_col]) if cid_col and not pd.isna(row[cid_col]) else None
        study_id = str(row[stid_col]) if stid_col and not pd.isna(row[stid_col]) else None
        time_id = str(row[tpid_col]) if tpid_col and not pd.isna(row[tpid_col]) else None

        key = (study_id, case_id, sid, time_id)
        index[key] = list(evs)

        # Common fallbacks/aliases
        index[("unspecified" if study_id else None, case_id, sid, time_id)] = list(evs)
        index[(study_id, case_id, _to_patient_barcode(sid), time_id)] = list(evs)
        index[("unspecified" if study_id else None, case_id, _to_patient_barcode(sid), time_id)] = list(evs)
        index[(study_id, _to_patient_barcode(case_id or ""), _to_patient_barcode(sid), time_id)] = list(evs)
        index[("unspecified" if study_id else None, _to_patient_barcode(case_id or ""), _to_patient_barcode(sid), time_id)] = list(evs)

    return index
