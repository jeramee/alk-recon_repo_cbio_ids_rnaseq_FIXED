"""Persistence / tolerance flags.

This module stays deliberately simple for the MVP:
- It looks for an expression_summary on the CaseSnapshot.
- It derives boolean flags used by the rule engine.

Terminology:
- persister_signature_score: numeric score (higher == more 'persister-like')
- has_expression_data: we have any RNA/expression attached
- has_persister_score_high: score >= threshold (default 1.0)
"""

from __future__ import annotations

from typing import Any, Dict

from schema.case_snapshot import CaseSnapshot


def apply_persistence_flags(cs: CaseSnapshot, *, threshold: float = 1.0) -> None:
    """Attach persistence-related flags to cs.flags in-place.

    Backwards-compatible: we also set older flag names if they exist in downstream code.
    """
    flags: Dict[str, Any] = cs.flags

    expr = cs.expression_summary or {}
    # If expression_summary is absent/empty, mark as no expression evidence.
    if not isinstance(expr, dict) or len(expr) == 0:
        flags.setdefault("has_expression_data", False)
        flags.setdefault("has_persister_score_high", False)
        # legacy names
        flags.setdefault("has_persistence_evidence", False)
        flags.setdefault("persister_signature_score_high", False)
        return

    flags["has_expression_data"] = True
    flags["has_persistence_evidence"] = True  # legacy-friendly umbrella

    score = expr.get("persister_signature_score")
    try:
        score_f = float(score) if score is not None else None
    except Exception:
        score_f = None

    if score_f is None:
        flags["has_persister_score_high"] = False
        flags["persister_signature_score_high"] = False  # legacy
        return

    flags["persister_signature_score"] = score_f  # convenience mirror
    flags["has_persister_score_high"] = score_f >= float(threshold)
    flags["persister_signature_score_high"] = flags["has_persister_score_high"]  # legacy
