from __future__ import annotations

import re

from schema.case_snapshot import CaseSnapshot, EvidenceItem, EvidenceLevel


# Common ALK resistance labels seen in literature. This list is not exhaustive.
SOLVENT_FRONT = {"G1202R"}
GATEKEEPER = {"L1196M"}


def _extract_protein_changes(cs: CaseSnapshot) -> list[str]:
    out: list[str] = []
    for v in cs.genomic.alk_variants:
        raw = (v.normalized or v.raw or "").strip()
        # normalize "p.G1202R" -> "G1202R"
        m = re.match(r"^p\.?([A-Za-z])(\d+)([A-Za-z\*])$", raw)
        if m:
            out.append(f"{m.group(1).upper()}{int(m.group(2))}{m.group(3).upper()}")
            continue
        # already looks like "G1202R"
        m2 = re.match(r"^([A-Za-z])(\d+)([A-Za-z\*])$", raw)
        if m2:
            out.append(f"{m2.group(1).upper()}{int(m2.group(2))}{m2.group(3).upper()}")
    return out


def apply_alk_flags(cs: CaseSnapshot) -> None:
    """Populate cs.genomic.flags with ALK-specific derived signals."""
    flags = cs.genomic.flags

    changes = _extract_protein_changes(cs)
    uniq = sorted(set(changes))

    flags["has_any_alk_mutation"] = bool(uniq)
    # Back-compat: tests expect this name.
    flags["has_any_alk_variant"] = bool(cs.genomic.alk_variants) if cs.genomic else False
    flags["has_G1202R"] = any(c in SOLVENT_FRONT for c in uniq)
    flags["has_L1196M"] = any(c in GATEKEEPER for c in uniq)
    flags["has_compound_alk_mutations"] = len(uniq) >= 2

    # add auditable inferred evidence items for any flags we compute
    # keep IDs stable-ish within snapshot by using a deterministic key
    for k in ["has_any_alk_mutation", "has_G1202R", "has_L1196M", "has_compound_alk_mutations"]:
        cs.evidence.append(
            EvidenceItem(
                id=f"F_{k}",
                level=EvidenceLevel.INFERRED,
                kind="feature_flag",
                label=f"Derived flag {k}",
                value=bool(flags.get(k, False)),
            )
        )

    if uniq:
        cs.evidence.append(
            EvidenceItem(
                id="F_alk_variant_list",
                level=EvidenceLevel.INFERRED,
                kind="feature",
                label="Derived list of ALK protein changes",
                value=uniq,
            )
        )
