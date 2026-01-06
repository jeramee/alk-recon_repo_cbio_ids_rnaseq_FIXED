from __future__ import annotations

from typing import Any, Dict, List

from schema.case_snapshot import CaseSnapshot


def _norm_pc(x: str) -> str:
    s = (x or "").strip()
    if s.startswith("p."):
        s = s[2:]
    return s.upper()


def apply_on_target_flags(cs: CaseSnapshot) -> None:
    flags: Dict[str, Any] = cs.flags

    alk_rows: List[Dict[str, Any]] = [r for r in cs.variants if str(r.get("gene", "")).upper() == "ALK"]
    pcs = [_norm_pc(r.get("protein_change", "")) for r in alk_rows if r.get("protein_change", "")]
    pcs = [p for p in pcs if p]

    flags["alk_mutations"] = sorted(set(pcs))
    flags["has_any_alk_mutation"] = len(flags["alk_mutations"]) > 0

    flags["has_G1202R"] = any("G1202R" in p for p in pcs)
    flags["has_L1196M"] = any("L1196M" in p for p in pcs)

    # Compound = 2+ distinct ALK protein changes
    flags["has_compound_alk_mutations"] = len(set(pcs)) >= 2

    # Evidence ledger
    if flags["has_G1202R"]:
        cs.add_evidence(
            layer="KNOWN",
            feature_name="has_G1202R",
            feature_value=True,
            source_type="variants",
            source_ref="variant_table",
            gene="ALK",
            note="Solvent-front mutation detected in ALK.",
        )
    if flags["has_L1196M"]:
        cs.add_evidence(
            layer="KNOWN",
            feature_name="has_L1196M",
            feature_value=True,
            source_type="variants",
            source_ref="variant_table",
            gene="ALK",
            note="Gatekeeper mutation detected in ALK.",
        )
    if flags["has_compound_alk_mutations"]:
        cs.add_evidence(
            layer="INFERRED",
            feature_name="has_compound_alk_mutations",
            feature_value=True,
            source_type="variants",
            source_ref="variant_table",
            gene="ALK",
            note=f"Multiple distinct ALK protein changes observed: {', '.join(sorted(set(pcs)))}",
        )
