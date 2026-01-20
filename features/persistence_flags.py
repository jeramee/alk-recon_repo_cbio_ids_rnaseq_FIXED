"""Persistence/tolerance flag computation.

The ingest layer can attach a persister_score into:
  cs.expression.signature_scores["persister_score"]

Tests rely on:
- has_persister_score
- has_persister_score_high

We also keep older keys for back-compat.
"""

from __future__ import annotations

from schema.case_snapshot import CaseSnapshot


def _extract_persister_score(cs: CaseSnapshot):
    if not cs.expression or not cs.expression.signature_scores:
        return None
    # be permissive about casing
    for k, v in cs.expression.signature_scores.items():
        if str(k).lower() == "persister_score":
            try:
                return float(v)
            except Exception:
                return None
    return None


def apply_persistence_flags(cs: CaseSnapshot, high_threshold: float = 2.0) -> CaseSnapshot:
    if cs.genomic is None:
        return cs

    flags = dict(cs.genomic.flags or {})

    score = _extract_persister_score(cs)
    has_score = score is not None
    is_high = bool(has_score and score >= float(high_threshold))

    flags["has_persister_score"] = has_score
    flags["has_persister_score_high"] = is_high

    # older keys (keep them if downstream uses them)
    flags["persister_signature_score"] = score if has_score else None
    flags["persister_signature_score_high"] = is_high

    cs.genomic.flags = flags
    return cs
