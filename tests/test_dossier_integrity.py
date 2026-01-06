from pathlib import Path
import json

from schema.case_snapshot import CaseSnapshot
from features.apply_features import apply_all_features
from mechanism_engine.rule_engine import compute_mechanism_calls
from reports.dossier import write_dossier_md, write_dossier_json


def test_dossier_outputs_are_written_and_auditable(tmp_path: Path):
    cs = CaseSnapshot(study_id="STUDY1", case_id="CASE_0001", sample_id="S1", timepoint_id="baseline", variants=[
        {"gene": "ALK", "protein_change": "G1202R", "variant_type": "SNV"},
        {"gene": "MET", "variant_type": "AMP", "copy_number": "8"},
    ])
    apply_all_features(cs)
    compute_mechanism_calls(cs)

    md_path = tmp_path / "dossiers" / "CASE_0001.md"
    js_path = tmp_path / "dossiers" / "CASE_0001.json"
    write_dossier_md(cs, md_path)
    write_dossier_json(cs, js_path)

    assert md_path.exists()
    assert js_path.exists()

    md = md_path.read_text(encoding="utf-8")
    assert "Mechanism calls" in md
    assert "| evidence_id |" in md

    d = json.loads(js_path.read_text(encoding="utf-8"))
    assert "evidence_ledger" in d
    assert len(d["evidence_ledger"]) >= 1
    assert d["evidence_ledger"][0]["evidence_id"].startswith("E")
