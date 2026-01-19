from __future__ import annotations

import json
from pathlib import Path

from alk_recon.pipeline import run_pipeline


def test_smoke(tmp_path: Path):
    input_path = Path(__file__).parent / "data" / "variants_sample.tsv"
    outdir = tmp_path / "out"
    run_pipeline(str(input_path), str(outdir))

    # Expect one dossier per case
    cases = {"CASE_G1202R", "CASE_L1196M", "CASE_COMPOUND", "CASE_BYPASS", "CASE_PERSIST"}
    for case in cases:
        js = outdir / case / f"{case}.case_snapshot.json"
        md = outdir / case / f"{case}.dossier.md"
        assert js.exists(), js
        assert md.exists(), md

        obj = json.loads(js.read_text(encoding="utf-8"))
        assert obj["case_id"] == case
        assert len(obj["mechanism_calls"]) >= 1

    # Quick sanity on two cases
    g1202r = json.loads((outdir / "CASE_G1202R" / "CASE_G1202R.case_snapshot.json").read_text(encoding="utf-8"))
    flags = g1202r["genomic"]["flags"]
    assert flags.get("has_G1202R") is True

    persist = json.loads((outdir / "CASE_PERSIST" / "CASE_PERSIST.case_snapshot.json").read_text(encoding="utf-8"))
    pflags = persist["genomic"]["flags"]
    assert pflags.get("persister_signature_score_high") is True
