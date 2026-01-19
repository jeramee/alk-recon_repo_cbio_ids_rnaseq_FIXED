from __future__ import annotations

import json, math
from pathlib import Path

from alk_recon.pipeline import run_pipeline


def test_rnaseq_top_markers_no_nan(tmp_path: Path) -> None:
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
    expr = cs.get("expression")
    assert expr is not None, "Expected expression to be attached"

    # persister_score should be finite
    score = (expr.get("signature_scores") or {}).get("persister_score")
    assert score is not None
    assert isinstance(score, (int, float))
    assert math.isfinite(float(score))

    # top_markers z should never be NaN/Inf
    for m in (expr.get("top_markers") or []):
        z = m.get("z")
        assert z is not None
        assert math.isfinite(float(z)), f"Non-finite z for {m.get('gene')}: {z}"
