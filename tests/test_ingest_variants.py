from pathlib import Path

from ingest.variant_table_import import load_variant_table, table_to_snapshots


def test_variant_table_ignores_comment_lines():
    p = Path(__file__).parent / "fixtures" / "variants_with_comments.tsv"
    df = load_variant_table(str(p))
    assert "case_id" in df.columns
    assert len(df) == 3

    snaps = table_to_snapshots(df)
    # CASE_0001 + baseline, CASE_0002 + baseline
    assert len(snaps) == 2
    s1 = [s for s in snaps if s.case_id == "CASE_0001"][0]
    assert len(s1.variants) == 2
