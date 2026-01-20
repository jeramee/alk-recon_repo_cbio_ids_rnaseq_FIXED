from ingest.variant_table_import import load_variant_table, dataframe_to_case_snapshots
from features.apply_features import apply_features_to_snapshot
from mechanism_engine.rule_engine import compute_mechanism_calls, route_strategy
from schema.case_snapshot import MechanismCallType, StrategyBucket


def _snapshot(tmp_path, tsv_text: str):
    p = tmp_path / "v.tsv"
    p.write_text(tsv_text)
    df = load_variant_table(str(p))
    s = dataframe_to_case_snapshots(df)[0]
    return apply_features_to_snapshot(s)


def test_rule_engine_on_target_alk_wins_for_solvent_front(tmp_path):
    s = _snapshot(
        tmp_path,
        "case_id\tsample_id\ttimepoint_id\tgene\tprotein_change\tvariant_type\tvaf\n"
        "CASE_1\tS1\tT1\tALK\tG1202R\tSNV\t0.20\n",
    )
    calls = compute_mechanism_calls(s)
    top = max(calls, key=lambda c: c.score)

    assert top.type == MechanismCallType.ON_TARGET_ALK
    buckets = route_strategy(s, calls)
    assert StrategyBucket.A in buckets  # mutation-aware ATP-site
    assert StrategyBucket.C in buckets  # degrader is often relevant


def test_rule_engine_bypass_wins_for_met_amp_without_alk(tmp_path):
    s = _snapshot(
        tmp_path,
        "case_id\tsample_id\ttimepoint_id\tgene\tvariant_type\tcopy_number\n"
        "CASE_2\tS2\tT1\tMET\tAMP\t8\n",
    )
    calls = compute_mechanism_calls(s)
    top = max(calls, key=lambda c: c.score)

    assert top.type == MechanismCallType.BYPASS
    buckets = route_strategy(s, calls)
    assert StrategyBucket.E in buckets  # sequencing/monitoring logic is typical here


def test_rule_engine_persistence_wins_for_high_persister_score(tmp_path):
    s = _snapshot(
        tmp_path,
        "case_id\tsample_id\ttimepoint_id\tgene\tvariant_type\tpersister_score\n"
        "CASE_3\tS3\tT1\tEGFR\tSNV\t2.1\n",
    )
    calls = compute_mechanism_calls(s)
    top = max(calls, key=lambda c: c.score)

    assert top.type == MechanismCallType.PERSISTENCE
    buckets = route_strategy(s, calls)
    assert StrategyBucket.D in buckets  # persistence suppression adjunct logic
