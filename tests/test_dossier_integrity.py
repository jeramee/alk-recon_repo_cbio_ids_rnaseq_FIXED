import json
import re

from ingest.variant_table_import import load_variant_table, dataframe_to_case_snapshots
from features.apply_features import apply_features_to_snapshot
from mechanism_engine.rule_engine import compute_mechanism_calls, route_strategy
from reports.dossier import build_json_dossier, build_markdown_dossier


def test_dossier_has_evidence_ids_and_machine_readable_structure(tmp_path):
    p = tmp_path / "variants.tsv"
    p.write_text(
        "study_id\tcase_id\tsample_id\ttimepoint_id\tgene\tprotein_change\tvariant_type\tvaf\tcopy_number\tpersister_score\n"
        "luad_tcga_pan_can_atlas_2018\tCASE_X\tSAMPLE_X\tT1\tALK\tG1202R\tSNV\t0.2\t\t\n"
        "luad_tcga_pan_can_atlas_2018\tCASE_X\tSAMPLE_X\tT1\tMET\t\tAMP\t\t8\t\n"
    )

    df = load_variant_table(str(p))
    s = dataframe_to_case_snapshots(df)[0]
    s = apply_features_to_snapshot(s)

    calls = compute_mechanism_calls(s)
    routing = route_strategy(s, calls)

    j = build_json_dossier(s, calls, routing)
    assert j["case_id"] == "CASE_X"
    assert "mechanism_calls" in j and len(j["mechanism_calls"]) >= 3
    assert "evidence_ledger" in j and len(j["evidence_ledger"]) > 0

    # Evidence IDs should exist and look like "E1", "E2", ...
    eids = [e.get("evidence_id") for e in j["evidence_ledger"]]
    assert any(re.fullmatch(r"E\d+", str(x)) for x in eids)

    md = build_markdown_dossier(s, calls, routing)
    assert "Mechanism" in md  # headline section presence
    assert re.search(r"\bE\d+\b", md)  # evidence IDs appear in markdown

    # JSON should be serializable (this catches sneaky non-serializable objects)
    json.dumps(j)

