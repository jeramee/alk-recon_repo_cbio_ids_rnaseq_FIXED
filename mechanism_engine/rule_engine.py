from __future__ import annotations

from dataclasses import replace

from schema.case_snapshot import (
    CaseSnapshot,
    MechanismCall,
    MechanismType,
    StrategyBucket,
    StrategyRouting,
)

on_target_score = 0.0
bypass_score = 0.0
persistence_score = 0.0

if flags.get("has_any_alk_variant") or flags.get("has_any_alk_mutation"):
    on_target_score = 0.5

if flags.get("has_MET_amp") or flags.get("has_any_bypass_event"):
    bypass_score = 0.6

if flags.get("has_persister_score_high"):
    persistence_score = 0.6


def _score_on_target(cs):
    """
    Minimal scorer. Returns a plain dict (no new types).
    """
    flags = {}
    g = getattr(cs, "genomic", None)
    if g is not None:
        flags = getattr(g, "flags", {}) or {}

    score = 0.0
    reasons = []

    if flags.get("has_alk_fusion"):
        score += 3.0; reasons.append("ALK fusion")
    if flags.get("has_alk_mutation"):
        score += 1.0; reasons.append("ALK mutation")
    if flags.get("has_gatekeeper_mutation"):
        score += 2.0; reasons.append("ALK gatekeeper")
    if flags.get("has_solvent_front_mutation"):
        score += 3.0; reasons.append("ALK solvent-front")

    return {
        "bucket": "on_target",
        "score": float(score),
        "reasons": reasons,
    }


def _score_bypass(cs):
    flags = {}
    g = getattr(cs, "genomic", None)
    if g is not None:
        flags = getattr(g, "flags", {}) or {}

    score = 0.0
    reasons = []

    # Use whatever flags you already set in your feature layer
    if flags.get("has_met_amp_or_high") or flags.get("has_met_alt"):
        score += 2.0; reasons.append("MET alteration")
    if flags.get("has_egfr_alt"):
        score += 1.5; reasons.append("EGFR alteration")
    if flags.get("has_kras_alt"):
        score += 1.5; reasons.append("KRAS alteration")
    if flags.get("has_erbb2_alt"):
        score += 1.0; reasons.append("ERBB2 alteration")
    if flags.get("has_ret_alt"):
        score += 1.0; reasons.append("RET alteration")

    return {
        "bucket": "bypass",
        "score": float(score),
        "reasons": reasons,
    }

def _score_persistence(cs):
    """
    Uses expression.signature_scores['persister_score'] if present.
    Returns a dict.
    """
    expr = getattr(cs, "expression", None)
    sig = getattr(expr, "signature_scores", {}) if expr is not None else {}
    sig = sig or {}

    pers = sig.get("persister_score", None)

    score = 0.0
    reasons = []

    if pers is None:
        reasons.append("no expression signature")
        score = 0.0
    else:
        # very conservative thresholds; adjust later
        try:
            p = float(pers)
        except Exception:
            p = None

        if p is None:
            reasons.append("persister_score not numeric")
            score = 0.0
        elif p >= 0.5:
            score = 2.0; reasons.append(f"high persister_score={p:.3f}")
        elif p >= 0.0:
            score = 1.0; reasons.append(f"mid persister_score={p:.3f}")
        else:
            score = 0.0; reasons.append(f"low persister_score={p:.3f}")

    return {
        "bucket": "persistence",
        "score": float(score),
        "reasons": reasons,
    }


def score_mechanisms_and_route(cs: CaseSnapshot) -> CaseSnapshot:
    """Score mechanisms (rule-based) and compute strategy routing.

    This is intentionally interpretable and auditable for early MVP.
    You can later replace scoring with an ML model, but keep the same
    inputs/outputs.
    """

    flags = cs.genomic.flags or {}

    # --- Mechanism scoring (0..1, heuristic) ---
    on_target = 0.0
    bypass = 0.0
    persistence = 0.0

    # On-target ALK evidence
    if flags.get("has_any_alk_mutation"):
        on_target += 0.40
    if flags.get("has_G1202R"):
        on_target += 0.35
    if flags.get("has_L1196M"):
        on_target += 0.25
    if flags.get("has_compound_alk_mutations"):
        on_target += 0.35

    # Bypass evidence
    if flags.get("has_any_bypass_event"):
        bypass += 0.35
    if flags.get("has_MET_event"):
        bypass += 0.30
    if flags.get("has_EGFR_event"):
        bypass += 0.25
    if flags.get("has_MAPK_event"):
        bypass += 0.20

    # Persistence evidence
    if flags.get("has_expression_data"):
        persistence += 0.15
    if flags.get("has_persister_score_high"):
        persistence += 0.60

    # Normalize into 0..1 by simple cap
    on_target = min(1.0, on_target)
    bypass = min(1.0, bypass)
    persistence = min(1.0, persistence)

    # Build calls (ranked)
    calls = []
    calls.append(
        MechanismCall(
            mechanism=MechanismType.ON_TARGET_ALK,
            score=on_target,
            rationale=_rationale_on_target(cs),
            supporting_evidence_ids=_supporting_evidence_ids(cs, "on_target_alk"),
            contradicting_evidence_ids=[],
        )
    )
    calls.append(
        MechanismCall(
            mechanism=MechanismType.BYPASS,
            score=bypass,
            rationale=_rationale_bypass(cs),
            supporting_evidence_ids=_supporting_evidence_ids(cs, "bypass"),
            contradicting_evidence_ids=[],
        )
    )
    calls.append(
        MechanismCall(
            mechanism=MechanismType.PERSISTENCE,
            score=persistence,
            rationale=_rationale_persistence(cs),
            supporting_evidence_ids=_supporting_evidence_ids(cs, "persistence"),
            contradicting_evidence_ids=[],
        )
    )

    calls_sorted = sorted(calls, key=lambda c: c.score, reverse=True)
    cs.mechanism_calls = calls_sorted

    cs.routing = _route(cs)
    return cs


def _route(cs: CaseSnapshot) -> StrategyRouting:
    flags = cs.genomic.flags or {}
    top = cs.mechanism_calls[0] if cs.mechanism_calls else None

    ranked: list[StrategyBucket] = []
    what_to_test: list[str] = []
    what_to_avoid: list[str] = []

    # Some stable guardrails (keeps reports honest)
    what_to_avoid.append(
        "Do not treat this as a clinical recommendation; this is a research summary of evidence." 
    )

    # If we have on-target mutations, prioritize mutation-aware + degradation
    if flags.get("has_any_alk_mutation"):
        ranked.append(StrategyBucket.A)
        ranked.append(StrategyBucket.C)
        what_to_test.append("Confirm effect across a WT + mutant panel (include solvent-front, gatekeeper, compound if possible).")

        if flags.get("has_compound_alk_mutations"):
            what_to_test.append("If compound mutations exist, explicitly test compound constructs; single-mutant results can mislead.")

    # If bypass features, add bypass-aware logic
    if flags.get("has_any_bypass_event"):
        ranked.append(StrategyBucket.E)
        what_to_test.append("Validate bypass signal (e.g., MET/EGFR/MAPK) with orthogonal evidence when available.")

    # If persistence evidence, include persistence adjunct bucket
    if flags.get("has_persister_score_high"):
        ranked.append(StrategyBucket.D)
        what_to_test.append("Run a durability assay (days–weeks) to quantify a persister fraction and rebound kinetics.")

    # If weak evidence overall, emphasize sequencing/monitoring logic
    if top is not None and top.score < 0.35:
        ranked = [StrategyBucket.E]
        what_to_test = [
            "Acquire higher-yield evidence: ctDNA/tissue NGS for ALK mutations and key bypass events; optional expression signatures if available."
        ]

    # Deduplicate while preserving order
    ranked_unique: list[StrategyBucket] = []
    for b in ranked:
        if b not in ranked_unique:
            ranked_unique.append(b)

    return StrategyRouting(
        ranked_buckets=ranked_unique,
        what_to_test_next=what_to_test,
        what_to_avoid=what_to_avoid,
    )

def _rationale_on_target(cs: CaseSnapshot) -> str:
    f = cs.genomic.flags or {}
    if not f.get("has_any_alk_mutation"):
        return "No ALK kinase-domain variant was detected in the provided variant table."
    bits = []
    if f.get("has_G1202R"):
        bits.append("solvent-front pattern (G1202R)")
    if f.get("has_L1196M"):
        bits.append("gatekeeper pattern (L1196M)")
    if f.get("has_compound_alk_mutations"):
        bits.append("multiple ALK variants (possible compound context)")
    if not bits:
        bits.append("ALK variant(s) detected")
    return "On-target ALK resistance is supported by: " + ", ".join(bits) + "."


def _rationale_bypass(cs: CaseSnapshot) -> str:
    f = cs.genomic.flags or {}
    if not f.get("has_any_bypass_event"):
        return "No obvious bypass-related alterations were detected in the provided variant table."
    bits = []
    if f.get("has_MET_event"):
        bits.append("MET-related alteration")
    if f.get("has_EGFR_event"):
        bits.append("EGFR/ERBB-family alteration")
    if f.get("has_MAPK_event"):
        bits.append("MAPK-pathway alteration")
    if not bits:
        bits.append("bypass-related alteration")
    return "Bypass signaling is supported by: " + ", ".join(bits) + "."


def _rationale_persistence(cs: CaseSnapshot) -> str:
    f = cs.genomic.flags or {}
    if not f.get("has_expression_data"):
        return "No expression-level persistence evidence was provided (expression block is missing)."
    if f.get("has_persister_score_high"):
        return "Expression signature score suggests a possible drug-tolerant persister program (threshold exceeded)."
    return "Expression data was provided, but no persistence signature crossed the current threshold."


def _supporting_evidence_ids(cs: CaseSnapshot, mechanism: str) -> list[str]:
    """Heuristic mapping from flags to evidence items.

    For the MVP we keep this simple: any evidence label containing a key term
    is considered supporting. Later you can make this explicit and exact.
    """
    mech_terms = {
        "on_target_alk": ["ALK", "G1202R", "L1196M"],
        "bypass": ["MET", "EGFR", "ERBB", "KRAS", "NRAS", "BRAF", "MAP2K", "MAPK"],
        "persistence": ["persister", "signature", "expression"],
    }
    terms = mech_terms.get(mechanism, [])
    out: list[str] = []
    for ev in cs.evidence:
        text = f"{ev.label} {ev.value}".upper()
        if any(t.upper() in text for t in terms):
            out.append(ev.id)
    return out

# --- Backwards-compatible wrappers for older tests / docs ---

# --- Backwards-compatible wrappers for older tests / docs ---

def compute_mechanism_calls(cs: CaseSnapshot):
    """
    Compatibility wrapper for older tests.
    Ensures cs.mechanism_calls is populated and returns it.
    """
    score_mechanisms_and_route(cs)
    return cs.mechanism_calls or []


def route_strategy(cs: CaseSnapshot, calls=None):
    """Compatibility wrapper.

    Some tests call route_strategy(cs, calls).
    Internally we route based on cs.mechanism_calls.
    """
    if calls is not None:
        cs.mechanism_calls = list(calls)

    # If no calls provided, compute them.
    if not getattr(cs, "mechanism_calls", None):
        cs.mechanism_calls = compute_mechanism_calls(cs)

    # Ensure ranked best-to-worst.
    cs.mechanism_calls = sorted(cs.mechanism_calls, key=lambda c: c.score, reverse=True)

    cs.routing = _route(cs)
    return cs.routing


