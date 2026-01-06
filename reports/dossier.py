from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict

from schema.case_snapshot import CaseSnapshot


def build_dossier(cs: CaseSnapshot) -> Dict[str, Any]:
    return {
        "identifiers": {
            "study_id": cs.study_id,
            "case_id": cs.case_id,
            "sample_id": cs.sample_id,
            "timepoint_id": cs.timepoint_id,
        },
        "mechanism_calls": [c.__dict__ for c in cs.mechanism_calls],
        "flags": cs.flags,
        "evidence_ledger": cs.evidence_ledger,
        "routing_buckets": _routing_buckets(cs),
    }


def _routing_buckets(cs: CaseSnapshot):
    f = cs.flags
    buckets = []
    # A: mutation-aware ATP-site
    if f.get("has_any_alk_mutation") or f.get("has_G1202R") or f.get("has_L1196M"):
        buckets.append({"bucket": "A", "label": "Mutation-aware ATP-site inhibitor", "why": "On-target ALK evidence present."})
    # B: allosteric clamp thesis
    if f.get("has_G1202R") or f.get("has_compound_alk_mutations"):
        buckets.append({"bucket": "B", "label": "Allosteric / conformation-lock", "why": "Hard mutations may benefit from non-ATP-site strategies."})
    # C: degradation
    if f.get("has_compound_alk_mutations") or f.get("has_G1202R"):
        buckets.append({"bucket": "C", "label": "Targeted degradation", "why": "Degrader strategies can bypass some binding-site liabilities."})
    # D: persistence
    if f.get("has_persister_score_high") or f.get("has_expression_data"):
        buckets.append({"bucket": "D", "label": "Persistence suppression adjuncts", "why": "Expression suggests tolerance/persistence is plausible."})
    # E: sequencing/monitoring (always)
    buckets.append({"bucket": "E", "label": "Sequencing / monitoring logic", "why": "Track evolving mechanisms over time."})
    return buckets


def write_dossier_json(cs: CaseSnapshot, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(build_dossier(cs), f, indent=2)


def write_dossier_md(cs: CaseSnapshot, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    d = build_dossier(cs)

    lines = []
    ids = d["identifiers"]
    lines.append(f"# ALK-RECON Dossier: {ids.get('case_id') or 'UNKNOWN'}")
    lines.append("")
    lines.append("## Identifiers")
    lines.append(f"- study_id: {ids.get('study_id')}")
    lines.append(f"- case_id: {ids.get('case_id')}")
    lines.append(f"- sample_id: {ids.get('sample_id')}")
    lines.append(f"- timepoint_id: {ids.get('timepoint_id')}")
    lines.append("")
    lines.append("## Mechanism calls (ranked)")
    for c in d["mechanism_calls"]:
        lines.append(f"- **{c['mechanism']}** — score={c['score']}")
        for r in c.get("rationale", []):
            lines.append(f"  - {r}")
    lines.append("")
    lines.append("## Top evidence")
    for e in d["evidence_ledger"][:10]:
        lines.append(f"- {e['evidence_id']}: {e['feature_name']} = {e['feature_value']} ({e['layer']})")
    lines.append("")
    lines.append("## Evidence ledger")
    lines.append("| evidence_id | layer | feature | value | source_type | gene | note |")
    lines.append("|---|---|---|---|---|---|---|")
    for e in d["evidence_ledger"]:
        lines.append(f"| {e['evidence_id']} | {e['layer']} | {e['feature_name']} | {e['feature_value']} | {e['source_type']} | {e.get('gene','')} | {e.get('note','')} |")
    lines.append("")
    lines.append("## Routing buckets")
    for b in d["routing_buckets"]:
        lines.append(f"- **{b['bucket']}** {b['label']}: {b['why']}")
    lines.append("")
    out_path.write_text("\n".join(lines), encoding="utf-8")
