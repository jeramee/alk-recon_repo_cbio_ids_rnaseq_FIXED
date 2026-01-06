from __future__ import annotations

from typing import List

from schema.case_snapshot import CaseSnapshot, MechanismCall


MECHANISMS = ("on_target_alk", "bypass", "persistence")


def compute_mechanism_calls(cs: CaseSnapshot) -> List[MechanismCall]:
    f = cs.flags

    on_target = 0.0
    if f.get("has_any_alk_mutation"):
        on_target += 1.0
    if f.get("has_G1202R"):
        on_target += 3.0
    if f.get("has_L1196M"):
        on_target += 2.0
    if f.get("has_compound_alk_mutations"):
        on_target += 2.0

    bypass = 0.0
    if f.get("has_met_event"):
        bypass += 2.0
    if f.get("has_egfr_event"):
        bypass += 1.5
    if f.get("has_mapk_event"):
        bypass += 1.0
    if f.get("has_bypass_evidence"):
        bypass += 0.5

    persistence = 0.0
    if f.get("has_persister_score_high"):
        persistence += 2.0
    # weaker: if expression exists and score is moderately positive
    score = f.get("persister_signature_score") or cs.expression_summary.get("persister_signature_score")
    try:
        score_f = float(score) if score is not None else None
    except Exception:
        score_f = None
    if f.get("has_expression_data") and score_f is not None and score_f >= 0.5:
        persistence += 1.0

    calls = [
        MechanismCall("on_target_alk", on_target, rationale=[]),
        MechanismCall("bypass", bypass, rationale=[]),
        MechanismCall("persistence", persistence, rationale=[]),
    ]

    # Add short rationale lines referencing available evidence IDs.
    ev = cs.evidence_ledger
    # Map feature_name -> first evidence_id
    first_eid = {}
    for e in ev:
        fn = e.get("feature_name")
        if fn and fn not in first_eid:
            first_eid[fn] = e.get("evidence_id")

    def add_rationale(mech: str, lines: List[str]):
        for c in calls:
            if c.mechanism == mech:
                c.rationale.extend(lines)

    if f.get("has_G1202R"):
        add_rationale("on_target_alk", [f"ALK solvent-front mutation (EID {first_eid.get('has_G1202R','?')})."])
    if f.get("has_L1196M"):
        add_rationale("on_target_alk", [f"ALK gatekeeper mutation (EID {first_eid.get('has_L1196M','?')})."])
    if f.get("has_compound_alk_mutations"):
        add_rationale("on_target_alk", [f"Compound ALK changes (EID {first_eid.get('has_compound_alk_mutations','?')})."])

    if f.get("has_met_event"):
        add_rationale("bypass", [f"MET bypass evidence (EID {first_eid.get('has_met_event','?')})."])
    if f.get("has_egfr_event"):
        add_rationale("bypass", [f"EGFR bypass evidence (EID {first_eid.get('has_egfr_event','?')})."])
    if f.get("has_mapk_event"):
        add_rationale("bypass", [f"MAPK-pathway event (EID {first_eid.get('has_mapk_event','?')})."])

    if f.get("has_persister_score_high"):
        add_rationale("persistence", ["Persister-like signature score above threshold (RNA-derived)." ])

    # Sort high-to-low score
    calls.sort(key=lambda c: c.score, reverse=True)
    cs.mechanism_calls = calls
    return calls
