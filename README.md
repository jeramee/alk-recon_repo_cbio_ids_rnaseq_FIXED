# ALK-RECON (cbio_ids + rnaseq)

**ALK-RECON** is a research/education pipeline starter kit that turns ALK-related ctDNA/NGS findings into an **auditable resistance-mechanism dossier**.

This repo flavor supports:

- **cBioPortal traceability** via `study_id` + `sample_id` columns (“cbio_ids”)
- **Optional RNA-seq**: count matrix + metadata (“rnaseq”)

> **Not medical advice.** This project is for research/education only.

---

## Project title

**ALK Lockpick** — *Stop resistance from changing the locks.*

---

## What you get

- A canonical, JSON-serializable **`CaseSnapshot`** data object (`schema/case_snapshot.py`)
- A flexible CSV/TSV importer for ctDNA/NGS variant tables (`ingest/variant_table_import.py`)
- Feature engineering for ALK resistance flags + bypass flags (`features/`)
- A transparent rule-based mechanism scorer + strategy-bucket router (`mechanism_engine/`)
- Markdown + JSON dossier output (`reports/`)
- A safe LLM “narrator” layer that **only** narrates existing evidence (`llm_layer/`)
- Unit tests + tiny fixtures (synthetic; no real patient data required) (`tests/`)

---

## Quick start (run the pipeline)

### 1) Install (editable)

> Run these commands **from the repo root** (the folder containing `pyproject.toml`, `Makefile`, and `alk_recon/`).

**Windows (PowerShell)**

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip uninstall -y alk-recon
python -m pip install -e .
```

**macOS / Linux**

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip uninstall -y alk-recon
python -m pip install -e .
```

### 2) Run (variants only)

```bash
python -m alk_recon.cli \
  --input examples/variant_table_template.tsv \
  --outdir out
```

### 3) Run (variants + RNA-seq)

```bash
python -m alk_recon.cli \
  --input examples/variant_table_template.tsv \
  --rnaseq-counts examples/rnaseq_counts_example.tsv \
  --rnaseq-meta examples/rnaseq_meta_example.tsv \
  --outdir out_rnaseq
```

### 4) Outputs

You should see dossiers written under your output directory, for example:

- `out/dossiers/<case_id>.md`
- `out/dossiers/<case_id>.json`
- `out/index.json`

---

## Input formats

### A) Variant table (CSV/TSV) — required

Minimum required columns:

- `case_id`
- `gene`

Strongly recommended:

- `sample_id` (helps link variants to RNA if you have it)
- `timepoint_id` (baseline vs progression is huge for resistance logic)
- `variant_type` (SNV / INDEL / AMP / DEL / FUSION)
- `protein_change` (e.g., `G1202R`, `L1196M`)
- `vaf` (0–1 or 0–100; the pipeline can normalize if configured)
- `copy_number` (optional)
- `notes` (optional)

Example (TSV):

```tsv
study_id	case_id	sample_id	timepoint_id	gene	protein_change	variant_type	vaf
luad_tcga_pan_can_atlas_2018	CASE_0001	TCGA-XX-YYYY-01	baseline	ALK	G1202R	SNV	0.22
luad_tcga_pan_can_atlas_2018	CASE_0001	TCGA-XX-YYYY-01	baseline	MET		AMP	
```

> The importer is designed to ignore leading `#` comment lines in templates.

---

### B) RNA-seq count matrix (optional)

Expected shape:

- **genes × samples**
- First column: gene identifier (symbol / gene_id)
- Remaining columns: sample IDs (must match the metadata `sample_id` column)

Example:

```tsv
gene	TCGA-XX-YYYY-01	TCGA-XX-YYYY-02
ALK	10	8
MET	200	350
EGFR	50	45
```

### C) RNA-seq metadata table (optional)

Required:

- `sample_id`

Recommended:

- `case_id`
- `study_id`
- `timepoint_id`
- `condition`
- `batch`

---

## cBioPortal traceability (study + sample IDs)

If your variants came from cBioPortal, add:

- `study_id` = cBioPortal **studyId** (e.g., `luad_tcga_pan_can_atlas_2018`)
- `sample_id` = cBioPortal **sampleId**
- `case_id` = your stable per-patient/per-case id (often maps to cBioPortal patientId)

These IDs are carried into the dossier so you can backtrack to the originating cohort/sample later.

---

## Optional RNA-seq module (what it does)

If you provide RNA-seq counts + metadata, ALK-RECON can attach an **ExpressionSummary** to each matching `CaseSnapshot` and compute a simple, auditable signature score intended to support **persistence/tolerance** hypotheses.

- This is **not DESeq2** and not meant to be publishable differential expression by itself.
- It’s a minimal scoring path you can later swap for your lab’s preferred method.

---

## Running tests (pytest)

### 1) Unzip + open a terminal in the repo root

You want to be in the folder that contains `pyproject.toml`, `Makefile`, and `alk_recon/`.

### 2) Create + activate a virtual environment

**Windows (PowerShell)**

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
```

**macOS / Linux**

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```

### 3) Install the project (editable)

This makes sure `alk-recon` runs from *this folder*.

```bash
python -m pip uninstall -y alk-recon
python -m pip install -e .
```

### 4) Install pytest (dev dependency)

```bash
python -m pip install pytest
```

> `pyproject.toml` does not include pytest by default, so install it explicitly.

### 5) Run tests

Fast way:

```bash
python -m pytest -q
```

Run one test file:

```bash
python -m pytest -q tests/test_cli.py
```

Run tests matching a keyword:

```bash
python -m pytest -q -k rnaseq
```

### 6) Using the Makefile (optional but nice)

Your `Makefile` provides quick commands:

```bash
make dev     # installs project + pytest
make test    # runs pytest -q
make clean   # wipes pytest caches + common output dirs
```

**Windows note:** `make` works if you’re using **Git Bash**, **WSL**, or you have `make` installed.  
If not, just use the `python -m pytest -q` commands above.

### 7) If your tests folder came as a separate zip (tests.zip)

Unzip it into the repo root so you end up with:

```text
<repo>/
  tests/
    test_*.py
```

Then re-run:

```bash
python -m pytest -q
```

### 8) Sanity checks (am I running the right repo flavor?)

These help when multiple ALK-RECON flavors were installed previously:

```bash
python -c "import alk_recon; print(alk_recon.__file__)"
```

Then:

- Windows:

  ```powershell
  where alk-recon
  ```

- macOS/Linux:

  ```bash
  which alk-recon
  ```

If `alk_recon.__file__` points inside your current repo path, you’re good.

---

## Roadmap (practical)

- Tighten input validators + clearer error messages
- Expand golden test cases (compound, bypass-dominant, persistence-dominant)
- Add a simple report viewer (static HTML or Streamlit)
- Optional: cBioPortal fetch mode (study → sample lists → molecular profiles)
- Optional: plug-in DESeq2/edgeR outputs for a richer expression/persistence module

---

## License

MIT
