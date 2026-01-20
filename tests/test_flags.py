from ingest.variant_table_import import load_variant_table, dataframe_to_case_snapshots
from features.apply_features import apply_features_to_snapshot


def _build_snapshots_from_tsv_text(tmp_path, text: str):
    p = tmp_path / "variants.tsv"
    p.write_text(text)
    df = load_variant_table(str(p))
    return dataframe_to_case_snapshots(df)


def test_flags_detect_solvent_front_gatekeeper_and_compound(tmp_path):
    snapshots = _build_snapshots_from_tsv_text(
        tmp_path,
        "case_id\tsample_id\ttimepoint_id\tgene\tprotein_change\tvariant_type\tvaf\n"
        "CASE_C\tS1\tT1\tALK\tG1202R\tSNV\t0.20\n"
        "CASE_C\tS1\tT1\tALK\tL1196M\tSNV\t0.10\n",
    )
    s = apply_features_to_snapshot(snapshots[0])

    flags = s.genomic.flags
    assert flags.get("has_G1202R") is True
    assert flags.get("has_L1196M") is True
    assert flags.get("has_compound_alk_mutations") is True
    assert flags.get("has_any_alk_variant") is True


def test_flags_detect_bypass_driver_from_amp_event(tmp_path):
    snapshots = _build_snapshots_from_tsv_text(
        tmp_path,
        "case_id\tsample_id\ttimepoint_id\tgene\tprotein_change\tvariant_type\tcopy_number\n"
        "CASE_B\tS2\tT1\tMET\t\tAMP\t8\n",
    )
    s = apply_features_to_snapshot(snapshots[0])

    flags = s.genomic.flags
    assert flags.get("has_MET_amp") is True
    assert flags.get("has_bypass_driver") is True


def test_flags_detect_persistence_from_persister_score_column(tmp_path):
    # The ingest layer supports persister_score and stores it in expression.signature_scores
    snapshots = _build_snapshots_from_tsv_text(
        tmp_path,
        "case_id\tsample_id\ttimepoint_id\tgene\tprotein_change\tvariant_type\tpersister_score\n"
        "CASE_P\tS3\tT1\tEGFR\t\tSNV\t2.5\n",
    )
    s = apply_features_to_snapshot(snapshots[0])

    flags = s.genomic.flags
    assert flags.get("has_persister_score") is True
    assert flags.get("has_persister_score_high") is True
