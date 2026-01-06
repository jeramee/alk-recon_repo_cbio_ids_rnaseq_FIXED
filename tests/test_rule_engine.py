from schema.case_snapshot import CaseSnapshot
from features.apply_features import apply_all_features
from mechanism_engine.rule_engine import compute_mechanism_calls


def top_mech(cs):
    calls = compute_mechanism_calls(cs)
    return calls[0].mechanism


def test_rule_engine_on_target_wins_for_g1202r():
    cs = CaseSnapshot(case_id="C1", sample_id="S1", variants=[
        {"gene": "ALK", "protein_change": "G1202R", "variant_type": "SNV"},
    ])
    apply_all_features(cs)
    assert top_mech(cs) == "on_target_alk"


def test_rule_engine_bypass_wins_for_met_amp():
    cs = CaseSnapshot(case_id="C2", sample_id="S2", variants=[
        {"gene": "MET", "variant_type": "AMP", "copy_number": "8"},
    ])
    apply_all_features(cs)
    assert top_mech(cs) == "bypass"


def test_rule_engine_persistence_wins_for_high_signature():
    cs = CaseSnapshot(case_id="C3", sample_id="S3", variants=[])
    # Provide expression summary directly (simulate RNA attach)
    cs.expression_summary = {"has_expression_data": True, "persister_signature_score": 1.2}
    apply_all_features(cs, persister_threshold=1.0)
    assert top_mech(cs) == "persistence"
