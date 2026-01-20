# ingest/variant_table_import.py

from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from schema.case_snapshot import (
    CaseSnapshot,
    ClinicalContext,
    EvidenceItem,
    EvidenceLevel,
    ExpressionSummary,
    GenomicAlterations,
    ProteinVariant,
    Provenance,
    SampleLevel,
    VariantEffect,
)

# --- Column aliases (many NGS exports use different headers)
DEFAULT_COLUMN_ALIASES: dict[str, List[str]] = {
    # Optional cross-reference fields (helpful if your source is cBioPortal)
    "study_id": [
        "study_id",
        "study",
        "cbio_study",
        "cbio_study_id",
        "cbioportal_study_id",
        "cbioportal:studyid",
    ],
    "case_id": ["case_id", "patient", "patient_id", "subject", "id"],
    "sample_id": ["sample_id", "sample", "specimen", "tumor_sample", "sample_name"],
    "cbio_sample_id": [
        "cbio_sample_id",
        "cbioportal_sample_id",
        "cbioportal:sampleid",
        "sampleid",
    ],
    "cbio_patient_id": [
        "cbio_patient_id",
        "cbioportal_patient_id",
        "cbioportal:patientid",
        "patientid",
    ],
    "timepoint_id": ["timepoint", "timepoint_id", "visit", "collection", "draw"],
    "gene": ["gene", "hugo_symbol", "symbol"],
    "protein_change": [
        "protein_change",
        "hgvsp",
        "protein",
        "aa_change",
        "amino_acid_change",
        "hgvs_p",
    ],
    "cDNA_change": ["cdna_change", "hgvsc", "cdna", "hgvs_c"],
    "variant_class": ["variant_class", "variant_type", "type", "classification"],
    "effect": ["effect", "consequence"],
    "vaf": ["vaf", "allele_fraction", "tumor_vaf", "af"],
    "copy_number": ["copy_number", "cn", "cna", "copy_number_change"],
    "fusion": ["fusion", "fusion_name", "rearrangement"],
    "persister_score": [
        "persister_score",
        "persistence_score",
        "tolerance_score",
        "persister_signature",
        "drug_tolerant_score",
        "dtp_score",
    ],
}

AA_SUB_RE = re.compile(r"(?i)^p?\.?(?P<ref>[A-Z])(?P<pos>\d+)(?P<alt>[A-Z])$")


def _now_utc() -> str:
    return datetime.utcnow().replace(microsecond=0).isoformat() + "Z"


def _find_column(df: pd.DataFrame, canonical: str, explicit: Optional[str]) -> Optional[str]:
    if explicit and explicit in df.columns:
        return explicit
    aliases = DEFAULT_COLUMN_ALIASES.get(canonical, [])
    lower_map = {c.lower(): c for c in df.columns}
    for a in aliases:
        if a.lower() in lower_map:
            return lower_map[a.lower()]
    return None


def _coerce_vaf(x: Any) -> str:
    """Normalize VAF field.
    Keeps dtype=str, but ensures percent-like inputs become fractions.
    Example: "22" -> "0.22"
    """
    if x is None:
        return ""
    s = str(x).strip()
    if s == "":
        return ""
    try:
        v = float(s)
    except Exception:
        return s
    if v > 1.0:
        v = v / 100.0
    return str(v)


def load_variant_table(
    path: str | Path,
    delimiter: str | None = None,
    case_col: str | None = None,
    sample_col: str | None = None,
    timepoint_col: str | None = None,
    min_vaf: float | None = None,
) -> pd.DataFrame:
    """Load a CSV/TSV ctDNA/NGS variant table.

    Permissive loader: attempts to locate key columns via aliases and
    normalizes them to canonical names.

    Returns a DataFrame with canonical columns present when possible.
    """
    path = Path(path)
    if delimiter is None:
        delimiter = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","

    df = pd.read_csv(path, sep=delimiter, dtype=str, keep_default_na=False)

    # Normalize key columns
    c_study = _find_column(df, "study_id", None)
    c_case = _find_column(df, "case_id", case_col)
    c_sample = _find_column(df, "sample_id", sample_col)
    c_time = _find_column(df, "timepoint_id", timepoint_col)
    c_gene = _find_column(df, "gene", None)
    c_persist = _find_column(df, "persister_score", None)

    if c_gene is None:
        raise ValueError("Could not find a gene column. Add one or pass a mapping.")

    # Build missing ids if necessary
    if c_case is None:
        df["case_id"] = "CASE_" + (df.index.astype(str))
        c_case = "case_id"
    if c_sample is None:
        df["sample_id"] = df[c_case].astype(str)
        c_sample = "sample_id"
    if c_time is None:
        df["timepoint_id"] = "unspecified"
        c_time = "timepoint_id"
    if c_study is None:
        df["study_id"] = "unspecified"
        c_study = "study_id"

    rename_map: Dict[str, str] = {
        c_study: "study_id",
        c_case: "case_id",
        c_sample: "sample_id",
        c_time: "timepoint_id",
        c_gene: "gene",
    }
    if c_persist is not None:
        rename_map[c_persist] = "persister_score"

    df = df.rename(columns=rename_map)

    # --- VAF normalization (22 -> 0.22) ---
    c_vaf_pre = _find_column(df, "vaf", None)
    if c_vaf_pre and c_vaf_pre != "vaf":
        df = df.rename(columns={c_vaf_pre: "vaf"})
    if "vaf" in df.columns:
        df["vaf"] = df["vaf"].map(_coerce_vaf)

    # vaf filter (optional)
    if min_vaf is not None and "vaf" in df.columns:
        def _to_float(x: str) -> float:
            try:
                return float(x)
            except Exception:
                return float("nan")

        v = df["vaf"].map(_to_float)
        df = df.loc[(v.isna()) | (v >= float(min_vaf))].copy()

    return df


def _parse_protein_change(raw: str) -> Tuple[Optional[str], Optional[int], Optional[str], Optional[str]]:
    if not raw:
        return None, None, None, None
    raw = raw.strip()
    m = AA_SUB_RE.match(raw)
    if not m:
        # allow strings like "p.G1202R" or unknown formats; keep as-is
        return (f"p.{raw}" if not raw.lower().startswith("p.") else raw), None, None, None
    ref = m.group("ref").upper()
    pos = int(m.group("pos"))
    alt = m.group("alt").upper()
    return f"p.{ref}{pos}{alt}", pos, ref, alt


def _effect_from_row(row: pd.Series) -> VariantEffect:
    # attempt to infer effect from variant_class/effect text
    txt = " ".join([str(row.get(c, "")) for c in row.index]).lower()
    if "fusion" in txt or "rearr" in txt:
        return VariantEffect.FUSION
    if "amp" in txt or "amplification" in txt:
        return VariantEffect.AMPLIFICATION
    if "del" in txt and "deletion" in txt:
        return VariantEffect.DELETION
    if "frameshift" in txt:
        return VariantEffect.FRAMESHIFT
    if "nonsense" in txt or "stop_gained" in txt or "stop" in txt:
        return VariantEffect.NONSENSE
    if "splice" in txt:
        return VariantEffect.SPLICE
    if "inframe" in txt:
        return VariantEffect.INFRAME_INDEL
    if "missense" in txt or "nonsyn" in txt:
        return VariantEffect.MISSENSE
    return VariantEffect.UNKNOWN


def dataframe_to_case_snapshots(
    df: pd.DataFrame,
    source_id: str = "variant_table",
) -> list[CaseSnapshot]:
    """Group rows into CaseSnapshot objects.

    One snapshot per (case_id, sample_id, timepoint_id).
    """
    fetched_at = _now_utc()
    prov = Provenance(source="local_file", source_id=source_id, fetched_at_utc=fetched_at)

    c_prot = _find_column(df, "protein_change", None)
    c_fusion = _find_column(df, "fusion", None)
    c_cn = _find_column(df, "copy_number", None)
    c_persist = _find_column(df, "persister_score", None)

    snaps: List[CaseSnapshot] = []

    group_cols = (
        ["study_id", "case_id", "sample_id", "timepoint_id"]
        if "study_id" in df.columns
        else ["case_id", "sample_id", "timepoint_id"]
    )

    for keys, sub in df.groupby(group_cols, dropna=False):
        if len(group_cols) == 4:
            study_id, case_id, sample_id, timepoint_id = keys
        else:
            study_id = None
            case_id, sample_id, timepoint_id = keys

        genomic = GenomicAlterations()
        evidence: List[EvidenceItem] = []

        # Collect alterations
        for _, row in sub.reset_index(drop=True).iterrows():
            gene = str(row.get("gene", "")).strip()
            if not gene:
                continue

            raw_prot = str(row.get(c_prot, "")).strip() if c_prot else ""
            raw_fusion = str(row.get(c_fusion, "")).strip() if c_fusion else ""
            raw_cn = str(row.get(c_cn, "")).strip() if c_cn else ""

            eff = _effect_from_row(row)

            eid = f"E{len(evidence) + 1}"
            ev_val: Dict[str, Any] = {"gene": gene}

            if raw_prot:
                norm, pos, ref, alt = _parse_protein_change(raw_prot)
                ev_val.update({"protein_change": raw_prot, "normalized": norm})
            if raw_fusion:
                ev_val.update({"fusion": raw_fusion})
            if raw_cn:
                ev_val.update({"copy_number": raw_cn})

            evidence.append(
                EvidenceItem(
                    id=eid,
                    level=EvidenceLevel.KNOWN,
                    kind="variant_row",
                    label=f"Variant table row: {gene}",
                    value=ev_val,
                    provenance=prov,
                )
            )

            # Populate genomic structured fields for ALK / bypass
            if gene.upper() == "ALK":
                if eff == VariantEffect.FUSION and (raw_fusion or raw_prot):
                    genomic.alk_fusion = raw_fusion or raw_prot
                elif raw_prot:
                    norm, pos, ref, alt = _parse_protein_change(raw_prot)
                    pv = ProteinVariant(
                        gene="ALK",
                        raw=raw_prot,
                        normalized=norm,
                        effect=eff if eff != VariantEffect.UNKNOWN else VariantEffect.MISSENSE,
                        aa_pos=pos,
                        ref_aa=ref,
                        alt_aa=alt,
                    )
                    genomic.alk_variants.append(pv)
            else:
                if raw_cn:
                    genomic.copy_number_events.append({"gene": gene, "copy_number": raw_cn})
                if raw_prot or raw_fusion or raw_cn:
                    genomic.bypass_events.append(
                        {
                            "gene": gene,
                            "protein_change": raw_prot,
                            "fusion": raw_fusion,
                            "copy_number": raw_cn,
                            "effect": eff.value,
                        }
                    )

        # Optional persistence score at this snapshot level
        expr = None
        if c_persist:
            vals = sub[c_persist].astype(str).tolist()
            pval = None
            for v in vals:
                v = v.strip()
                if not v or v.lower() == "nan":
                    continue
                try:
                    pval = float(v)
                    break
                except Exception:
                    continue

            if pval is not None:
                expr = ExpressionSummary(
                    platform="unknown",
                    contrast=None,
                    signature_scores={"persister_score": pval},
                    top_markers=[],
                )
                evidence.append(
                    EvidenceItem(
                        id=f"E{len(evidence) + 1}",
                        level=EvidenceLevel.KNOWN,
                        kind="expression_signature",
                        label="Persister/tolerance signature score provided",
                        value={"persister_score": pval},
                        provenance=prov,
                    )
                )

        snap = CaseSnapshot(
            schema_version="0.1.0",
            level=SampleLevel.TIMEPOINT,
            case_id=str(case_id),
            study_id=str(study_id) if study_id is not None else None,
            patient_id=str(case_id),
            sample_id=str(sample_id),
            timepoint_id=str(timepoint_id),
            collected_at_utc=None,
            genomic=genomic,
            expression=expr,
            clinical=ClinicalContext(),
            evidence=evidence,
            mechanism_calls=[],
            routing=None,
            tags=["ctDNA/NGS"],
            provenance=prov,
        )
        snaps.append(snap)

    return snaps
