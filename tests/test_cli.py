import json
import subprocess
import sys


def test_cli_writes_index_and_dossiers(tmp_path):
    variants = tmp_path / "variants.tsv"
    variants.write_text(
        "study_id\tcase_id\tsample_id\ttimepoint_id\tgene\tprotein_change\tvariant_type\tvaf\n"
        "luad_tcga_pan_can_atlas_2018\tCASE_CLI\tS_CLI\tbaseline\tALK\tL1196M\tSNV\t0.12\n"
    )

    outdir = tmp_path / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        "-m",
        "alk_recon.cli",
        "--input",
        str(variants),
        "--out",
        str(outdir),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, f"CLI failed:\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"

    index_path = outdir / "index.json"
    assert index_path.exists()

    index = json.loads(index_path.read_text())
    assert "cases" in index
    assert len(index["cases"]) == 1

    # Expect at least one dossier file written
    dossier_dir = outdir / "dossiers"
    assert dossier_dir.exists()
    md_files = list(dossier_dir.glob("*.md"))
    json_files = list(dossier_dir.glob("*.json"))
    assert md_files and json_files
