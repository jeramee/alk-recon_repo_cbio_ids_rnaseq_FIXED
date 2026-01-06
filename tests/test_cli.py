from pathlib import Path
import subprocess
import sys


def test_cli_end_to_end_with_rnaseq(tmp_path: Path):
    fixtures = Path(__file__).parent / "fixtures"
    variants = fixtures / "variants_with_comments.tsv"
    counts = fixtures / "rnaseq_counts.tsv"
    meta = fixtures / "rnaseq_meta.tsv"

    outdir = tmp_path / "out"
    cmd = [
        sys.executable, "-m", "alk_recon.cli",
        "--input", str(variants),
        "--rnaseq-counts", str(counts),
        "--rnaseq-meta", str(meta),
        "--outdir", str(outdir),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr

    assert (outdir / "index.json").exists()
    assert (outdir / "dossiers" / "CASE_0001.md").exists()
