# ALK-RECON

**ALK-RECON** is a research/education pipeline starter kit that turns ALK-related ctDNA/NGS findings into an **auditable mechanism report**.

This project is **not medical advice**.

## What you get
- A canonical, JSON-serializable `CaseSnapshot` data object (`schema/case_snapshot.py`).
- A flexible CSV/TSV importer for ctDNA/NGS variant tables (`ingest/variant_table_import.py`).
- Feature engineering for ALK resistance flags + bypass flags (`features/`).
- A transparent rule-based mechanism scorer and strategy-bucket router (`mechanism_engine/`).
- Markdown + JSON dossier output (`reports/`).
- A safe LLM "narrator" prompt template that only narrates existing evidence (`llm_layer/`).

## Quick start

### 1) Install

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

### 2) Run on an example table

```bash
python -m alk_recon.cli \
  --input examples/variant_table_template.tsv \
  --outdir out
```

Outputs:
- `out/CASE_<id>/dossier.json`
- `out/CASE_<id>/dossier.md`

## Input format
See `examples/variant_table_template.tsv` (and remove the leading `#` comment lines before using it).

## Optional RNA-seq module (counts matrix + metadata)

If you have RNA-seq **raw counts** (genes × samples) and a **metadata table** that maps those samples
to `case_id` / `sample_id` (and optionally `study_id` / `timepoint_id`), ALK-RECON can attach an
`ExpressionSummary` to each matching `CaseSnapshot` and compute a simple, auditable signature score
intended to support **persistence/tolerance** hypotheses.

This is **not DESeq2** and not meant to be publishable differential expression on its own — it’s a
minimal scoring pathway you can swap for your lab’s preferred method.

### Expected files

1) **Counts matrix** (`--rnaseq-counts`):
   - TSV/CSV
   - One gene identifier column (default: first column, or `gene`/`gene_id`/`symbol`)
   - Remaining columns are sample IDs

2) **Metadata** (`--rnaseq-meta`):
   - TSV/CSV
   - Must include a sample id column (default autodetect: `sample_id`, `SampleID`, etc.)
   - Recommended: `case_id` (and optionally `study_id`, `timepoint_id`)

### Run with RNA-seq inputs

```bash
python -m alk_recon.cli \
  --input examples/variant_table_template.tsv \
  --rnaseq-counts examples/rnaseq_counts_example.tsv \
  --rnaseq-meta examples/rnaseq_meta_example.tsv \
  --outdir out
```

Defaults:
- A small **example** persistence gene list is bundled at `data/persistence_signatures/example_persister_genes.txt`.
- Override with `--signature-genes path/to/your_genes.txt`.

### cBioPortal traceability (study + sample IDs)

If your variants came from cBioPortal, add the optional `study_id` column (the cBioPortal `studyId`).

Practical mapping:
- `study_id`  = cBioPortal studyId (e.g., `luad_tcga_pan_can_atlas_2018`)
- `sample_id` = cBioPortal sampleId
- `case_id`   = a stable per-patient id (often the cBioPortal `patientId`)

When present, `study_id` is carried into each `CaseSnapshot` and appears in JSON/Markdown dossiers so you can backtrack to the originating cohort and sample.

## Roadmap
- Add cBioPortal ingestion as a first-class input.
- Improve expression module: plug-in DESeq2/edgeR outputs + richer persistence assays.
- Add a small set of curated "golden" cases for regression tests.

