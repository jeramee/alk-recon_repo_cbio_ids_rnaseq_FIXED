"""ALK-RECON — CaseSnapshot schema (research / education only)

This file defines the canonical "unit of analysis" for the ALK-RECON pipeline:
a single CaseSnapshot representing one sample (optionally with timepoint) plus
structured evidence and mechanism calls.

Design goals:
- Consistent, machine-readable fields (for rules + ML + reporting)
- Explicit uncertainty + provenance (what we know vs infer)
- JSON Schema export for validation + interoperability

No medical advice is provided or implied by this schema.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from datetime import datetime
from enum import Enum
from typing import Any, Dict, List, Literal, Optional, Union
import json


# ---------------------------------------------------------------------
# Enums / controlled vocab
# ---------------------------------------------------------------------

class EvidenceLevel(str, Enum):
    """How solid is this evidence? (NOT clinical guidance; just data quality)"""
    KNOWN = "known"          # Directly observed in the input data
    INFERRED = "inferred"    # Derived/engineered from known data (e.g., category flag)
    SPECULATIVE = "speculative"  # Hypothesis / annotation, not observed


class VariantEffect(str, Enum):
    """Best-effort classification of variant effect (can be 'unknown')."""
    MISSENSE = "missense"
    NONSENSE = "nonsense"
    FRAMESHIFT = "frameshift"
    SPLICE = "splice"
    INFRAME_INDEL = "inframe_indel"
    FUSION = "fusion"
    AMPLIFICATION = "amplification"
    DELETION = "deletion"
    UNKNOWN = "unknown"


class MechanismType(str, Enum):
    """Top-level resistance mechanism hypotheses."""
    ON_TARGET_ALK = "on_target_alk"
    BYPASS = "bypass"
    PERSISTENCE = "persistence"


class StrategyBucket(str, Enum):
    """Strategy bucket routing output (non-prescriptive, research framing)."""
    A = "A_nextgen_atp_site_inhibitor"
    B = "B_allosteric_conformation_lock"
    C = "C_targeted_degradation"
    D = "D_persistence_epigenetic_adjunct"
    E = "E_sequencing_monitoring_logic"


class SampleLevel(str, Enum):
    """Whether the snapshot is at patient / sample / timepoint level."""
    PATIENT = "patient"
    SAMPLE = "sample"
    TIMEPOINT = "timepoint"


# ---------------------------------------------------------------------
# Core evidence objects
# ---------------------------------------------------------------------

@dataclass
class Provenance:
    """Where did this datum come from? Keep it boring and explicit."""
    source: str                       # e.g., "cbioportal", "local_file", "manual_note"
    source_id: Optional[str] = None    # e.g., studyId, file path, accession
    fetched_at_utc: Optional[str] = None  # ISO timestamp
    raw_ref: Optional[Dict[str, Any]] = None  # optional: endpoint params, row id, etc.


@dataclass
class EvidenceItem:
    """
    A single piece of evidence that supports (or contradicts) a mechanism call.
    This is what makes the system auditable and LLM-safe (narrate what exists).
    """
    id: str
    level: EvidenceLevel
    kind: str  # e.g., "mutation", "copy_number", "expression_signature", "clinical", "literature"
    label: str  # human-readable tag, e.g. "ALK G1202R detected", "MET amplification"
    value: Optional[Union[str, float, int, bool, Dict[str, Any]]] = None
    units: Optional[str] = None
    direction: Optional[Literal["up", "down", "present", "absent", "unknown"]] = None
    confidence: Optional[float] = None  # 0..1, if your pipeline assigns it
    provenance: Optional[Provenance] = None
    notes: Optional[str] = None


# ---------------------------------------------------------------------
# Genomic / expression / clinical substructures
# ---------------------------------------------------------------------

@dataclass
class ProteinVariant:
    """
    A protein-level variant representation (keep original + normalized).
    Example: raw="G1202R" -> normalized="p.G1202R"
    """
    gene: str
    raw: str
    normalized: Optional[str] = None
    effect: VariantEffect = VariantEffect.UNKNOWN
    aa_pos: Optional[int] = None
    ref_aa: Optional[str] = None
    alt_aa: Optional[str] = None
    @property
    def protein_change(self) -> Optional[str]:
        return self.normalized or self.raw



@dataclass
class GenomicAlterations:
    """Structured alterations relevant to ALK resistance mapping."""
    alk_fusion: Optional[str] = None  # e.g., "EML4-ALK" if known
    alk_variants: List[ProteinVariant] = field(default_factory=list)
    copy_number_events: List[Dict[str, Any]] = field(default_factory=list)  # amps/dels
    bypass_events: List[Dict[str, Any]] = field(default_factory=list)       # MET amp, EGFR, KRAS, etc.

    # Convenience engineered flags (derived; store as INFERRED evidence too)
    flags: Dict[str, bool] = field(default_factory=dict)


@dataclass
class ExpressionSummary:
    """
    Optional. Used for persistence/tolerance hypotheses.
    Store computed signature scores, DESeq2 contrasts, etc.
    """
    platform: Optional[str] = None  # e.g., "rnaseq", "microarray"
    contrast: Optional[str] = None  # e.g., "treated_vs_control"
    signature_scores: Dict[str, float] = field(default_factory=dict)  # e.g., {"persister_score": 1.7}
    top_markers: List[Dict[str, Any]] = field(default_factory=list)   # e.g., [{"gene":"...", "log2FC":...}]


@dataclass
class ClinicalContext:
    """
    Minimal, non-exhaustive clinical context fields.
    This schema supports public datasets where treatment history is often missing.
    """
    diagnosis: Optional[str] = None
    cancer_type: Optional[str] = None
    stage: Optional[str] = None
    age: Optional[float] = None
    sex: Optional[str] = None

    # Treatment history is often absent in cBioPortal; keep optional.
    prior_alk_inhibitors: List[str] = field(default_factory=list)  # free text, normalized later
    response_pattern: Optional[str] = None  # e.g., "early_progression", "late_progression", "mixed"
    toxicity_constraints: List[str] = field(default_factory=list)  # non-prescriptive program constraints


# ---------------------------------------------------------------------
# Outputs: mechanism calls + strategy routing
# ---------------------------------------------------------------------

@dataclass
class MechanismCall:
    """A ranked hypothesis about the resistance mechanism for this snapshot."""
    mechanism: MechanismType
    score: float  # 0..1 relative confidence within your engine (not clinical)
    rationale: str  # short human-readable reason
    supporting_evidence_ids: List[str] = field(default_factory=list)
    contradicting_evidence_ids: List[str] = field(default_factory=list)

    @property
    def type(self):
        """Compatibility alias for tests: `.type` mirrors `.mechanism`."""
        return self.mechanism


@dataclass
class StrategyRouting:
    """Non-prescriptive routing: which strategy buckets are worth exploring next, given evidence."""
    ranked_buckets: List[StrategyBucket] = field(default_factory=list)
    what_to_test_next: List[str] = field(default_factory=list)  # research validation steps
    what_to_avoid: List[str] = field(default_factory=list)      # guardrails to avoid confusion


# ---------------------------------------------------------------------
# The canonical unit: CaseSnapshot
# ---------------------------------------------------------------------

@dataclass
class CaseSnapshot:
    """Canonical data object representing a single patient/sample/timepoint snapshot."""
    schema_version: str

    # Identity
    level: SampleLevel
    case_id: str                   # stable id for the patient/case
    study_id: Optional[str] = None # e.g., cBioPortal studyId
    patient_id: Optional[str] = None
    sample_id: Optional[str] = None
    timepoint_id: Optional[str] = None  # e.g., "baseline", "progression_1"
    collected_at_utc: Optional[str] = None

    # Main data
    genomic: GenomicAlterations = field(default_factory=GenomicAlterations)
    expression: Optional[ExpressionSummary] = None
    clinical: Optional[ClinicalContext] = None

    # Evidence ledger (auditable)
    evidence: List[EvidenceItem] = field(default_factory=list)

    # Outputs (filled by your mechanism engine)
    mechanism_calls: List[MechanismCall] = field(default_factory=list)
    routing: Optional[StrategyRouting] = None

    # Meta
    tags: List[str] = field(default_factory=list)  # e.g., ["ALK+", "ctDNA"]
    notes: Optional[str] = None
    provenance: Optional[Provenance] = None

    def to_dict(self) -> Dict[str, Any]:
        """Dataclass -> JSON-serializable dict (deep)."""
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, sort_keys=False, default=str)


# ---------------------------------------------------------------------
# JSON Schema (self-contained; no external deps required)
# ---------------------------------------------------------------------

JSON_SCHEMA: Dict[str, Any] = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "https://alk-recon.local/schema/case_snapshot.schema.json",
    "title": "CaseSnapshot",
    "type": "object",
    "required": ["schema_version", "level", "case_id", "genomic", "evidence", "mechanism_calls"],
    "properties": {
        "schema_version": {"type": "string"},
        "level": {"type": "string", "enum": [e.value for e in SampleLevel]},
        "case_id": {"type": "string"},
        "study_id": {"type": ["string", "null"]},
        "patient_id": {"type": ["string", "null"]},
        "sample_id": {"type": ["string", "null"]},
        "timepoint_id": {"type": ["string", "null"]},
        "collected_at_utc": {"type": ["string", "null"], "description": "ISO 8601 UTC timestamp"},

        "genomic": {
            "type": "object",
            "required": ["alk_variants", "copy_number_events", "bypass_events", "flags"],
            "properties": {
                "alk_fusion": {"type": ["string", "null"]},
                "alk_variants": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "required": ["gene", "raw", "effect"],
                        "properties": {
                            "gene": {"type": "string"},
                            "raw": {"type": "string"},
                            "normalized": {"type": ["string", "null"]},
                            "effect": {"type": "string", "enum": [e.value for e in VariantEffect]},
                            "aa_pos": {"type": ["integer", "null"]},
                            "ref_aa": {"type": ["string", "null"]},
                            "alt_aa": {"type": ["string", "null"]},
                        },
                        "additionalProperties": False,
                    },
                },
                "copy_number_events": {"type": "array", "items": {"type": "object"}},
                "bypass_events": {"type": "array", "items": {"type": "object"}},
                "flags": {"type": "object", "additionalProperties": {"type": "boolean"}},
            },
            "additionalProperties": False,
        },

        "expression": {
            "type": ["object", "null"],
            "properties": {
                "platform": {"type": ["string", "null"]},
                "contrast": {"type": ["string", "null"]},
                "signature_scores": {"type": "object", "additionalProperties": {"type": "number"}},
                "top_markers": {"type": "array", "items": {"type": "object"}},
            },
            "additionalProperties": False,
        },

        "clinical": {
            "type": ["object", "null"],
            "properties": {
                "diagnosis": {"type": ["string", "null"]},
                "cancer_type": {"type": ["string", "null"]},
                "stage": {"type": ["string", "null"]},
                "age": {"type": ["number", "null"]},
                "sex": {"type": ["string", "null"]},
                "prior_alk_inhibitors": {"type": "array", "items": {"type": "string"}},
                "response_pattern": {"type": ["string", "null"]},
                "toxicity_constraints": {"type": "array", "items": {"type": "string"}},
            },
            "additionalProperties": False,
        },

        "evidence": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["id", "level", "kind", "label"],
                "properties": {
                    "id": {"type": "string"},
                    "level": {"type": "string", "enum": [e.value for e in EvidenceLevel]},
                    "kind": {"type": "string"},
                    "label": {"type": "string"},
                    "value": {"type": ["string", "number", "integer", "boolean", "object", "array", "null"]},
                    "units": {"type": ["string", "null"]},
                    "direction": {"type": ["string", "null"], "enum": ["up", "down", "present", "absent", "unknown", None]},
                    "confidence": {"type": ["number", "null"], "minimum": 0.0, "maximum": 1.0},
                    "provenance": {
                        "type": ["object", "null"],
                        "properties": {
                            "source": {"type": "string"},
                            "source_id": {"type": ["string", "null"]},
                            "fetched_at_utc": {"type": ["string", "null"]},
                            "raw_ref": {"type": ["object", "null"]},
                        },
                        "additionalProperties": False,
                    },
                    "notes": {"type": ["string", "null"]},
                },
                "additionalProperties": False,
            },
        },

        "mechanism_calls": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["mechanism", "score", "rationale", "supporting_evidence_ids", "contradicting_evidence_ids"],
                "properties": {
                    "mechanism": {"type": "string", "enum": [e.value for e in MechanismType]},
                    "score": {"type": "number", "minimum": 0.0, "maximum": 1.0},
                    "rationale": {"type": "string"},
                    "supporting_evidence_ids": {"type": "array", "items": {"type": "string"}},
                    "contradicting_evidence_ids": {"type": "array", "items": {"type": "string"}},
                },
                "additionalProperties": False,
            },
        },

        "routing": {
            "type": ["object", "null"],
            "properties": {
                "ranked_buckets": {"type": "array", "items": {"type": "string", "enum": [e.value for e in StrategyBucket]}},
                "what_to_test_next": {"type": "array", "items": {"type": "string"}},
                "what_to_avoid": {"type": "array", "items": {"type": "string"}},
            },
            "additionalProperties": False,
        },

        "tags": {"type": "array", "items": {"type": "string"}},
        "notes": {"type": ["string", "null"]},
        "provenance": {
            "type": ["object", "null"],
            "properties": {
                "source": {"type": "string"},
                "source_id": {"type": ["string", "null"]},
                "fetched_at_utc": {"type": ["string", "null"]},
                "raw_ref": {"type": ["object", "null"]},
            },
            "additionalProperties": False,
        },
    },
    "additionalProperties": False,
}


def get_json_schema() -> Dict[str, Any]:
    """Return the canonical JSON Schema dict for CaseSnapshot."""
    return JSON_SCHEMA


def save_json_schema(path: str) -> None:
    """Write schema to disk (e.g., schema/case_snapshot.schema.json)."""
    with open(path, "w", encoding="utf-8") as f:
        json.dump(JSON_SCHEMA, f, indent=2, ensure_ascii=False)


# ---------------------------------------------------------------------
# Minimal example instance (handy for tests / docs)
# ---------------------------------------------------------------------

def example_case_snapshot():
    now = datetime.utcnow().replace(microsecond=0).isoformat() + "Z"
    ev1 = EvidenceItem(
        id="E1",
        level=EvidenceLevel.KNOWN,
        kind="mutation",
        label="ALK variant detected",
        value={"gene": "ALK", "raw": "G1202R", "normalized": "p.G1202R"},
        provenance=Provenance(source="local_file", fetched_at_utc=now),
    )
    genomic = GenomicAlterations(
        alk_fusion="EML4-ALK",
        alk_variants=[ProteinVariant(
            gene="ALK",
            raw="G1202R",
            normalized="p.G1202R",
            effect=VariantEffect.MISSENSE,
            aa_pos=1202,
            ref_aa="G",
            alt_aa="R",
        )],
        flags={"has_G1202R": True, "has_L1196M": False, "has_compound_alk_mutations": False},
    )
    call = MechanismCall(
        mechanism=MechanismType.ON_TARGET_ALK,
        score=0.85,
        rationale="Solvent-front ALK mutation pattern consistent with on-target resistance.",
        supporting_evidence_ids=["E1"],
        contradicting_evidence_ids=[],
    )
    routing = StrategyRouting(
        ranked_buckets=[StrategyBucket.A_ATP_MUTATION_AWARE, StrategyBucket.C_DEGRADATION],
        what_to_test_next=[
            "Confirm mutant panel sensitivity in vitro (WT + solvent-front + gatekeeper + compound).",
            "Check bypass events (MET/EGFR/MAPK) in the same NGS panel if available.",
        ],
        what_to_avoid=[
            "Do not assume persistence dominates without expression/persistence evidence.",
        ],
    )
    return CaseSnapshot(
        schema_version="0.1.0",
        level=SampleLevel.SAMPLE,
        case_id="CASE_0001",
        study_id="example_study",
        patient_id="PAT_0001",
        sample_id="SAMP_0001",
        collected_at_utc=now,
        genomic=genomic,
        expression=None,
        clinical=None,
        evidence=[ev1],
        mechanism_calls=[call],
        routing=routing,
        tags=["ALK+", "example"],
        provenance=Provenance(source="manual_note", source_id="example", fetched_at_utc=now),
    )

# --- Backwards-compatible alias expected by older tests ---
MechanismCallType = MechanismType
