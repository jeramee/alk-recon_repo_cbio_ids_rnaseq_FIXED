# Quick Start — ALK-RECON

← Back to [root README](../README.md)  
→ Next: [Synthetic test data](README_test_data.md)

Run the pipeline end-to-end with minimal setup.

```bash
python -m alk_recon.cli run   --variants data/test/variants_synthetic_alk.tsv   --out out
```

This ingests variants, builds snapshots, and emits dossiers.
