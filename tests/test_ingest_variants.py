import pandas as pd

from ingest.variant_table_import import load_variant_table, dataframe_to_case_snapshots


def test_load_variant_table_accepts_cbioportal_ids(tmp_path):
    # Synthetic TSV includes study_id + sample_id to match the cbio_ids flavor.
    tsv = tmp_path / "variants_cbio.tsv"
    tsv.write_text(
        "study_id\tcase_id\tsample_id\ttimepoint_id\tgene\tprotein_change\tvariant_type\tvaf\tcopy_number\tnotes\n"
        "luad_tcga_pan_can_atlas_2018\tCASE_0001\tTCGA-XX-YYYY-01\tbaseline\tALK\tG1202R\tSNV\t22\t\tctDNA\n"
        "luad_tcga_pan_can_atlas_2018\tCASE_0001\tTCGA-XX-YYYY-01\tbaseline\tMET\t\tAMP\t\t8\tbypass\n"
    )

    df = load_variant_table(str(tsv))

    # Basic shape/columns
    assert len(df) == 2
    for col in ["study_id", "case_id", "sample_id", "timepoint_id", "gene", "variant_type"]:
        assert col in df.columns

    # VAF normalization: 22 -> 0.22 (the loader should normalize 0–100 to 0–1)
    alk_row = df[df["gene"] == "ALK"].iloc[0]
    assert 0.0 <= float(alk_row["vaf"]) <= 1.0
    assert abs(float(alk_row["vaf"]) - 0.22) < 1e-9


def test_dataframe_to_case_snapshots_builds_expected_identity(tmp_path):
    tsv = tmp_path / "variants_one_case.tsv"
    tsv.write_text(
        "study_id\tcase_id\tsample_id\ttimepoint_id\tgene\tprotein_change\tvariant_type\tvaf\n"
        "luad_tcga_pan_can_atlas_2018\tCASE_0002\tSAMPLE_A\tprogression_1\tALK\tL1196M\tSNV\t0.15\n"
    )

    df = load_variant_table(str(tsv))
    snapshots = dataframe_to_case_snapshots(df)

    assert len(snapshots) == 1
    s = snapshots[0]

    assert s.case_id == "CASE_0002"
    assert s.sample_id == "SAMPLE_A"
    assert s.timepoint_id == "progression_1"
    assert s.study_id == "luad_tcga_pan_can_atlas_2018"

    # ALK variant should be captured and normalized
    assert s.genomic is not None
    assert len(s.genomic.alk_variants) == 1
    pv = s.genomic.alk_variants[0]
    # We accept either raw "L1196M" or normalized "p.L1196M" depending on your current normalizer.
    assert "L1196M" in (pv.protein_change or "")
