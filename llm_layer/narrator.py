from __future__ import annotations

from typing import Any, Dict

from schema.case_snapshot import CaseSnapshot


SAFE_NARRATOR_SYSTEM = (
    "You are a scientific assistant. You ONLY narrate what is present in the provided JSON. "
    "You do not give medical advice and you do not recommend treatment for any individual. "
    "You can propose research validation steps at a high level, but never dosing, prescribing, or patient guidance. "
    "If information is missing, say so explicitly."
)


def build_safe_narrator_prompt(case: CaseSnapshot) -> Dict[str, str]:
    """Return a {system, user} message pair for an LLM.

    This is a *safe* layer: structured-in -> structured-out.
    """
    payload = case.to_json(indent=2)
    user = (
        "Given the following CaseSnapshot JSON, write a concise mechanism summary and evidence list. "
        "Rules: (1) cite evidence item IDs; (2) do not invent data; (3) no medical advice; "
        "(4) separate known vs inferred; (5) keep it under 250 words.\n\n"
        f"CASE_SNAPSHOT_JSON:\n{payload}"
    )
    return {"system": SAFE_NARRATOR_SYSTEM, "user": user}


def deterministic_narration(case: CaseSnapshot) -> str:
    """LLM-free narration (always available)."""
    top = case.mechanism_calls[0] if case.mechanism_calls else None
    lines = []
    lines.append(f"Case: {case.case_id}")
    if top:
        lines.append(f"Top mechanism: {top.mechanism.value} (score={top.score:.2f})")
        lines.append(f"Rationale: {top.rationale}")
        if top.supporting_evidence_ids:
            lines.append("Supporting evidence IDs: " + ", ".join(top.supporting_evidence_ids))
    else:
        lines.append("No mechanism calls present (pipeline not run?).")
    return "\n".join(lines)
