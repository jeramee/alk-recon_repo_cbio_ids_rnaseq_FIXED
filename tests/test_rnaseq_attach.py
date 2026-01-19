from __future__ import annotations

import json
from pathlib import Path

from alk_recon.pipeline import run_pipeline


def test_rnaseq_counts_attach_drives_persistence_flags(tmp_path: Path) -> None:
    """Smoke test: RNA-seq counts+meta can drive persistence flags when variant table lacks persister_score."""

    outdir = tmp_path / "out"
    run_pipeline(
        input_path=Path(__file__).parent / "data" / "variants_sample_noexpr.tsv",
        out_dir=outdir,
        rnaseq_counts_path=Path(__file__).parent / "data" / "rnaseq_counts_test.tsv",
        rnaseq_metadata_path=Path(__file__).parent / "data" / "rnaseq_meta_test.tsv",
        rnaseq_sample_id_col="sample_id",
        rnaseq_case_id_col="case_id",
        rnaseq_study_id_col="study_id",
        rnaseq_timepoint_id_col="timepoint_id",
    )

    cs_path = outdir / "CASE_PERSIST" / "CASE_PERSIST.case_snapshot.json"
    assert cs_path.exists()

    cs = json.loads(cs_path.read_text(encoding="utf-8"))

    # Expression summary is attached
    assert cs.get("expression") is not None
    assert "signature_scores" in cs["expression"]
    assert "persister_score" in cs["expression"]["signature_scores"]

    # And a derived persistence flag is present (thresholding depends on cohort size / scaling)
    assert cs["genomic"]["flags"].get("persister_signature_score_high") in {True, False}
