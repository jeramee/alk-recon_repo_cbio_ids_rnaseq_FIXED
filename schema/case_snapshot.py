from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional


@dataclass
class MechanismCall:
    mechanism: str  # on_target_alk | bypass | persistence
    score: float
    rationale: List[str] = field(default_factory=list)


@dataclass
class CaseSnapshot:
    """Unit of analysis: one case/sample/(optional)timepoint.

    This object is intentionally small and serialization-friendly.
    """
    study_id: Optional[str] = None
    case_id: Optional[str] = None
    sample_id: Optional[str] = None
    timepoint_id: Optional[str] = None

    variants: List[Dict[str, Any]] = field(default_factory=list)
    expression_summary: Dict[str, Any] = field(default_factory=dict)

    flags: Dict[str, Any] = field(default_factory=dict)
    evidence_ledger: List[Dict[str, Any]] = field(default_factory=list)
    mechanism_calls: List[MechanismCall] = field(default_factory=list)

    def _next_evidence_id(self) -> str:
        return f"E{len(self.evidence_ledger) + 1:03d}"

    def add_evidence(
        self,
        *,
        layer: str,
        feature_name: str,
        feature_value: Any,
        source_type: str,
        source_ref: str,
        gene: Optional[str] = None,
        note: Optional[str] = None,
    ) -> str:
        eid = self._next_evidence_id()
        self.evidence_ledger.append(
            {
                "evidence_id": eid,
                "study_id": self.study_id,
                "case_id": self.case_id,
                "sample_id": self.sample_id,
                "timepoint_id": self.timepoint_id,
                "layer": layer,
                "feature_name": feature_name,
                "feature_value": feature_value,
                "source_type": source_type,
                "source_ref": source_ref,
                "gene": gene,
                "note": note or "",
            }
        )
        return eid

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        # asdict will expand MechanismCall dataclasses into dicts automatically.
        return d
