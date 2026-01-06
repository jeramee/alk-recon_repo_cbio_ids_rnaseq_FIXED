from __future__ import annotations

import argparse

from alk_recon.pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="alk-recon", description="ALK-RECON: auditable ALK resistance mechanism dossiers (research/education)")
    p.add_argument("--input", required=True, help="Variant table CSV/TSV path.")
    p.add_argument("--outdir", required=True, help="Output directory.")
    p.add_argument("--delimiter", default=None, help="Delimiter override (default: auto-detect).")
    p.add_argument("--case-col", default="case_id", help="Name of case_id column.")
    p.add_argument("--sample-col", default="sample_id", help="Name of sample_id column.")
    p.add_argument("--timepoint-col", default="timepoint_id", help="Name of timepoint_id column.")
    p.add_argument("--study-col", default="study_id", help="Name of study_id column.")
    p.add_argument("--min-vaf", type=float, default=None, help="Optional minimum VAF filter (0-1).")

    # RNA-seq front door
    p.add_argument("--rnaseq-counts", default=None, help="Optional RNA-seq counts TSV.")
    p.add_argument("--rnaseq-meta", default=None, help="Optional RNA-seq metadata TSV (must include sample_id).")
    p.add_argument("--persister-threshold", type=float, default=1.0, help="Threshold for persister signature score (default: 1.0).")

    return p


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)

    run_pipeline(
        input_path=args.input,
        outdir=args.outdir,
        delimiter=args.delimiter,
        case_col=args.case_col,
        sample_col=args.sample_col,
        timepoint_col=args.timepoint_col,
        study_col=args.study_col,
        min_vaf=args.min_vaf,
        rnaseq_counts_path=args.rnaseq_counts,
        rnaseq_meta_path=args.rnaseq_meta,
        persister_threshold=args.persister_threshold,
    )


if __name__ == "__main__":
    main()
