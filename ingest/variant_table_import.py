from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import pandas as pd

from schema.case_snapshot import CaseSnapshot


def _guess_delimiter(path: str) -> str:
    # Look at the first non-comment line.
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        return "\t" if "\t" in s else "," if "," in s else "\t"
    return "\t"


def load_variant_table(path: str, *, delimiter: Optional[str] = None) -> pd.DataFrame:
    delim = delimiter or _guess_delimiter(path)
    df = pd.read_csv(path, sep=delim, dtype=str, keep_default_na=False, comment="#")
    # Normalize column names (strip whitespace)
    df.columns = [c.strip() for c in df.columns]
    return df


def table_to_snapshots(
    df: pd.DataFrame,
    *,
    case_col: str = "case_id",
    sample_col: str = "sample_id",
    timepoint_col: str = "timepoint_id",
    study_col: str = "study_id",
) -> List[CaseSnapshot]:
    if case_col not in df.columns:
        raise ValueError(f"Variant table missing required column: {case_col}")
    if "gene" not in df.columns:
        raise ValueError("Variant table missing required column: gene")

    # Fill missing optional columns to avoid KeyError during row -> dict
    for c in [study_col, sample_col, timepoint_col, "protein_change", "variant_type", "vaf", "copy_number", "notes"]:
        if c and c not in df.columns:
            df[c] = ""

    group_cols = [study_col, case_col, sample_col, timepoint_col]
    # Use a stable grouping even if study/sample/timepoint are empty.
    grouped = df.groupby(group_cols, dropna=False, sort=False)

    snapshots: List[CaseSnapshot] = []
    for (study_id, case_id, sample_id, timepoint_id), g in grouped:
        # Empty strings -> None for identity fields
        sid = study_id if str(study_id).strip() else None
        cid = case_id if str(case_id).strip() else None
        samp = sample_id if str(sample_id).strip() else None
        tp = timepoint_id if str(timepoint_id).strip() else None

        cs = CaseSnapshot(study_id=sid, case_id=cid, sample_id=samp, timepoint_id=tp)

        for i, row in g.reset_index(drop=True).iterrows():
            v = {k: (row[k] if k in row else "") for k in df.columns}
            v["_row_index"] = int(i)
            cs.variants.append(v)

        snapshots.append(cs)
    return snapshots
