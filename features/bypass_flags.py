from __future__ import annotations

from typing import Any, Dict, List

from schema.case_snapshot import CaseSnapshot

MAPK_GENES = {"KRAS", "NRAS", "BRAF", "MAP2K1", "MAP2K2"}


def apply_bypass_flags(cs: CaseSnapshot) -> None:
    flags: Dict[str, Any] = cs.flags

    genes = [str(r.get("gene", "")).upper() for r in cs.variants]
    vtypes = [str(r.get("variant_type", "")).upper() for r in cs.variants]
    pcs = [str(r.get("protein_change", "")) for r in cs.variants]

    def has_amp(gene: str) -> bool:
        for r in cs.variants:
            if str(r.get("gene", "")).upper() != gene:
                continue
            if str(r.get("variant_type", "")).upper() == "AMP":
                return True
            cn = r.get("copy_number", "")
            try:
                if cn and float(cn) >= 6:
                    return True
            except Exception:
                pass
        return False

    flags["has_met_event"] = has_amp("MET") or ("MET" in genes and "AMP" in vtypes)
    flags["has_egfr_event"] = has_amp("EGFR") or ("EGFR" in genes and "AMP" in vtypes)
    flags["has_mapk_event"] = any(g in MAPK_GENES for g in genes)

    flags["has_bypass_evidence"] = any([flags["has_met_event"], flags["has_egfr_event"], flags["has_mapk_event"]])

    if flags["has_met_event"]:
        cs.add_evidence(
            layer="KNOWN",
            feature_name="has_met_event",
            feature_value=True,
            source_type="variants",
            source_ref="variant_table",
            gene="MET",
            note="MET event detected (AMP or high copy_number).",
        )
    if flags["has_egfr_event"]:
        cs.add_evidence(
            layer="KNOWN",
            feature_name="has_egfr_event",
            feature_value=True,
            source_type="variants",
            source_ref="variant_table",
            gene="EGFR",
            note="EGFR event detected (AMP or high copy_number).",
        )
    if flags["has_mapk_event"]:
        cs.add_evidence(
            layer="INFERRED",
            feature_name="has_mapk_event",
            feature_value=True,
            source_type="variants",
            source_ref="variant_table",
            gene="MAPK",
            note="MAPK pathway gene event present (e.g., KRAS/BRAF/MAP2K1...).",
        )
