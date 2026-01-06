from schema.case_snapshot import CaseSnapshot
from features.apply_features import apply_all_features


def test_on_target_flags_detect_g1202r():
    cs = CaseSnapshot(case_id="C1", sample_id="S1", timepoint_id="t0", variants=[
        {"gene": "ALK", "protein_change": "G1202R", "variant_type": "SNV"},
    ])
    apply_all_features(cs)
    assert cs.flags["has_G1202R"] is True
    assert cs.flags["has_L1196M"] is False


def test_bypass_flags_detect_met_amp():
    cs = CaseSnapshot(case_id="C2", sample_id="S2", variants=[
        {"gene": "MET", "variant_type": "AMP", "copy_number": "8"},
    ])
    apply_all_features(cs)
    assert cs.flags["has_met_event"] is True
    assert cs.flags["has_bypass_evidence"] is True
