from __future__ import annotations

from pathlib import Path
from typing import Optional

from ingest.variant_table_import import load_variant_table, table_to_snapshots
from ingest.rnaseq_import import build_expression_index
from features.apply_features import apply_all_features
from mechanism_engine.rule_engine import compute_mechanism_calls
from reports.dossier import write_dossier_json, write_dossier_md


def run_pipeline(
    *,
    input_path: str,
    outdir: str,
    delimiter: Optional[str] = None,
    case_col: str = "case_id",
    sample_col: str = "sample_id",
    timepoint_col: str = "timepoint_id",
    study_col: str = "study_id",
    min_vaf: Optional[float] = None,
    rnaseq_counts_path: Optional[str] = None,
    rnaseq_meta_path: Optional[str] = None,
    persister_threshold: float = 1.0,
) -> None:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    df = load_variant_table(input_path, delimiter=delimiter)

    # Optional VAF filter (expects 0-1 or 0-100; keeps rows that parse)
    if min_vaf is not None and "vaf" in df.columns:
        def _ok(v):
            try:
                f = float(v)
                if f > 1.0:
                    f = f / 100.0
                return f >= float(min_vaf)
            except Exception:
                return True
        df = df[df["vaf"].apply(_ok)]

    cases = table_to_snapshots(
        df,
        case_col=case_col,
        sample_col=sample_col,
        timepoint_col=timepoint_col,
        study_col=study_col,
    )

    expr_index = None
    if rnaseq_counts_path and rnaseq_meta_path:
        expr_index = build_expression_index(
            rnaseq_counts_path,
            rnaseq_meta_path,
            delimiter=delimiter or "\t",
        )

    index_rows = []
    for cs in cases:
        # Attach RNA-seq summary if present (explicit key + smart fallback)
        if expr_index is not None and cs.sample_id:
            key_full = (cs.study_id, cs.case_id, cs.sample_id, cs.timepoint_id)
            expr = expr_index.get(key_full)
            if expr is None:
                expr = expr_index.get((None, None, cs.sample_id, None))
            if expr is not None:
                cs.expression_summary = expr

        apply_all_features(cs, persister_threshold=persister_threshold)
        compute_mechanism_calls(cs)

        case_id = cs.case_id or cs.sample_id or "UNKNOWN"
        md_path = out / "dossiers" / f"{case_id}.md"
        json_path = out / "dossiers" / f"{case_id}.json"
        write_dossier_md(cs, md_path)
        write_dossier_json(cs, json_path)

        index_rows.append(
            {
                "case_id": cs.case_id,
                "sample_id": cs.sample_id,
                "timepoint_id": cs.timepoint_id,
                "md_path": str(md_path.relative_to(out)),
                "json_path": str(json_path.relative_to(out)),
                "top_mechanism": cs.mechanism_calls[0].mechanism if cs.mechanism_calls else None,
            }
        )

    (out / "index.json").write_text(__import__("json").dumps(index_rows, indent=2), encoding="utf-8")
