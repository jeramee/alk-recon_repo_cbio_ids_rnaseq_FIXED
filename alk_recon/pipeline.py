from __future__ import annotations

"""End-to-end pipeline wrapper.

The repo is designed to be auditable and testable. The unit tests in this repo
expect a small, stable public API:

- run_pipeline(...)
  * accepts either `input_path` or the alias `input_variants_path`
  * writes dossiers to: <outdir>/dossiers/<case_id>.dossier.md and ...json
  * writes an index file: <outdir>/index.json with a top-level `cases` key

This is research/education tooling – not medical advice.
"""

from pathlib import Path
from typing import Optional

from ingest.variant_table_import import load_variant_table, dataframe_to_case_snapshots
from ingest.rnaseq_import import build_expression_index
from ingest.cna_import import build_cna_index
from features.apply_features import apply_features_to_snapshot
from mechanism_engine.rule_engine import compute_mechanism_calls, route_strategy
from reports.dossier import write_dossier_bundle, write_dossier_index


def _as_path(p: str | Path) -> Path:
    return p if isinstance(p, Path) else Path(str(p))


def run_pipeline(
    input_path: str | Path | None = None,
    outdir: str | Path | None = None,
    *args: object,
    # Back-compat aliases used in some tests / older callers
    input_variants_path: str | Path | None = None,
    out_dir: str | Path | None = None,
    rnaseq_counts_path: str | Path | None = None,
    rnaseq_metadata_path: str | Path | None = None,
    rnaseq_study_id_col: str | None = None,
    rnaseq_timepoint_id_col: str | None = None,
    cna_thresholded_path: str | Path | None = None,
    cna_linear_path: str | Path | None = None,
    cna_metadata_path: str | Path | None = None,
    cna_study_id_col: str | None = None,
    cna_timepoint_id_col: str | None = None,
    **_ignored: object,
):
    """Run the full pipeline and return index records.

    The extra kwargs are accepted intentionally to keep the API stable for local
    experiments and smoke tests.
    """

    # --- Positional back-compat ---
    # Some tests call: run_pipeline(str(input_path), str(outdir))
    if args:
        if input_path is None and len(args) >= 1:
            input_path = args[0]  # type: ignore[assignment]
        if outdir is None and len(args) >= 2:
            outdir = args[1]  # type: ignore[assignment]

    # --- Keyword aliases ---
    if outdir is None and out_dir is not None:
        outdir = out_dir
    if input_path is None and input_variants_path is not None:
        input_path = input_variants_path
    if input_path is None:
        raise ValueError("run_pipeline requires input_path (or input_variants_path)")

    input_path = _as_path(input_path)
    outdir = _as_path(outdir or Path("out"))
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_variant_table(input_path)
    snapshots = dataframe_to_case_snapshots(df)

    # Optional: attach RNA-seq expression
    if rnaseq_counts_path and rnaseq_metadata_path:
        expr_index = build_expression_index(
            counts_path=rnaseq_counts_path,
            metadata_path=rnaseq_metadata_path,
            study_id_col=rnaseq_study_id_col,
            timepoint_id_col=rnaseq_timepoint_id_col,
        )
        for cs in snapshots:
            key = (cs.study_id, cs.case_id, cs.sample_id, cs.timepoint_id)
            if key in expr_index:
                cs.expression = expr_index[key]

    # Optional: attach CNA events
    if (cna_thresholded_path or cna_linear_path) and cna_metadata_path:
        cna_index = build_cna_index(
            metadata_path=cna_metadata_path,
            thresholded_path=cna_thresholded_path,
            linear_path=cna_linear_path,
            study_id_col=cna_study_id_col,
            timepoint_id_col=cna_timepoint_id_col,
        )
        for cs in snapshots:
            key = (cs.study_id, cs.case_id, cs.sample_id, cs.timepoint_id)
            evs = cna_index.get(key)
            if evs:
                if cs.genomic is not None:
                    cs.genomic.copy_number_events = list(evs)

    dossiers_dir = outdir / "dossiers"
    dossiers_dir.mkdir(parents=True, exist_ok=True)

    records = []
    for cs in snapshots:
        cs = apply_features_to_snapshot(cs)
        calls = compute_mechanism_calls(cs)
        route_strategy(cs, calls)
        records.append(write_dossier_bundle(cs, dossiers_dir))

    write_dossier_index(records, outdir / "index.json")
    return records
