from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd
import pytest


def _norm_tcga_id(s: str) -> str:
    """Normalize TCGA-like ids across sources (dots vs hyphens, patient vs sample)."""
    s = str(s).strip().replace(".", "-")
    return s[:12]  # patient barcode


@pytest.mark.skipif(
    not Path("data/real/variants_cbio_luad_tcga.tsv").exists(),
    reason="TCGA smoke test requires local real data files (data/real/variants_cbio_luad_tcga.tsv)",
)
def test_tcga_smoke_filtered(tmp_path: Path):
    """Local regression test: ensure TCGA LUAD data attaches expression + CNA without NaNs.

    This test is meant for *your workstation* (not CI). It auto-skips if the files
    are missing.
    """
    from alk_recon.pipeline import run_pipeline

    variants_path = Path("data/real/variants_cbio_luad_tcga.tsv")
    rnaseq_counts = Path("data/real/rnaseq_expression_tcga_luad.tsv")
    rnaseq_meta = Path("data/real/sample_map.tcga_luad.tsv")

    if not rnaseq_counts.exists() or not rnaseq_meta.exists():
        pytest.skip("TCGA smoke test requires data/real/rnaseq_expression_tcga_luad.tsv and sample_map.tcga_luad.tsv")

    # Bundle CNA paths (typical extraction layout)
    cna_thresh = None
    for candidate in [
        Path("data/1/luad_tcga_bundle/luad_tcga/data_cna.txt"),
        Path("data/1/luad_tcga_bundle/luad_tcga_pub/data_cna.txt"),
    ]:
        if candidate.exists():
            cna_thresh = candidate
            break
    cna_lin = None
    for candidate in [
        Path("data/1/luad_tcga_bundle/luad_tcga/data_linear_cna.txt"),
        Path("data/1/luad_tcga_bundle/luad_tcga_pub/data_linear_cna.txt"),
    ]:
        if candidate.exists():
            cna_lin = candidate
            break

    # Filter variants down to cases we have RNA for (keeps the test fast and makes attachment expectations strict)
    v = pd.read_csv(variants_path, sep="\t", dtype=str).fillna("")
    m = pd.read_csv(rnaseq_meta, sep="\t", dtype=str).fillna("")

    m["case_norm"] = m["case_id"].map(_norm_tcga_id)
    keep_cases = set(m["case_norm"].tolist())

    v["case_norm"] = v["case_id"].map(_norm_tcga_id)
    v_f = v[v["case_norm"].isin(keep_cases)].drop(columns=["case_norm"], errors="ignore")
    assert len(v_f) > 0, "No variant rows overlap the rnaseq sample map; check your IDs"

    filtered_variants = tmp_path / "variants_tcga_smoke_filtered.tsv"
    v_f.to_csv(filtered_variants, sep="\t", index=False)

    outdir = tmp_path / "out_tcga_smoke_filtered"
    run_pipeline(
        input_variants_path=filtered_variants,
        outdir=outdir,
        rnaseq_counts_path=rnaseq_counts,
        rnaseq_metadata_path=rnaseq_meta,
        rnaseq_study_id_col="study_id",
        rnaseq_timepoint_id_col="timepoint_id",
        cna_thresholded_path=str(cna_thresh) if cna_thresh else None,
        cna_linear_path=str(cna_lin) if cna_lin else None,
        cna_metadata_path=str(rnaseq_meta),
        cna_study_id_col="study_id",
        cna_timepoint_id_col="timepoint_id",
    )

    snaps = list(outdir.glob("**/*.case_snapshot.json"))
    assert len(snaps) > 0

    # All snapshots in this filtered run should have expression attached and be finite
    saw_cna = False
    for p in snaps:
        d = json.loads(p.read_text(encoding="utf-8"))
        e = d.get("expression")
        assert e is not None, f"Missing expression in {p}"

        score = (e.get("signature_scores") or {}).get("persister_score")
        assert score is not None
        assert isinstance(score, (int, float))
        assert math.isfinite(float(score))

        for mkr in e.get("top_markers") or []:
            z = mkr.get("z")
            assert z is not None
            assert math.isfinite(float(z)), f"Non-finite z-score in top_markers for {p}: {mkr}"

        g = d.get("genomic") or {}
        cna_events = g.get("copy_number_events") or []
        if cna_events:
            saw_cna = True

    # If CNA files exist locally, ensure we attached at least one CNA event.
    if cna_thresh or cna_lin:
        assert saw_cna, "CNA input provided but no CNA events were attached"
