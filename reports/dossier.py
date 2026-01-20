from __future__ import annotations

from pathlib import Path
from typing import List
import json

from schema.case_snapshot import CaseSnapshot


def _md_escape(s: str) -> str:
    return s.replace("\n", " ").strip()


def render_markdown_dossier(cs: CaseSnapshot) -> str:
    """One-case, one-page-ish dossier (Markdown)."""

    lines: List[str] = []
    lines.append(f"# ALK-RECON Dossier — {cs.case_id}")
    lines.append("")
    lines.append("> Research/education use only. Not medical advice.")
    lines.append("")

    # Identity
    lines.append("## Identity")
    lines.append("")
    lines.append(f"- Level: `{cs.level}`")
    if cs.study_id:
        lines.append(f"- Study: `{cs.study_id}`")
    if cs.patient_id:
        lines.append(f"- Patient: `{cs.patient_id}`")
    if cs.sample_id:
        lines.append(f"- Sample: `{cs.sample_id}`")
    if cs.timepoint_id:
        lines.append(f"- Timepoint: `{cs.timepoint_id}`")
    if cs.collected_at_utc:
        lines.append(f"- Collected: `{cs.collected_at_utc}`")

    # Genomics summary
    lines.append("")
    lines.append("## Genomics summary")
    lines.append("")
    if cs.genomic.alk_fusion:
        lines.append(f"- ALK fusion: **{_md_escape(cs.genomic.alk_fusion)}**")
    if cs.genomic.alk_variants:
        vtxt = ", ".join([v.normalized or v.raw for v in cs.genomic.alk_variants])
        lines.append(f"- ALK variants: **{_md_escape(vtxt)}**")
    else:
        lines.append("- ALK variants: *(none recorded in this snapshot)*")

    if cs.genomic.bypass_events:
        lines.append("- Bypass events:")
        for ev in cs.genomic.bypass_events[:12]:
            gene = ev.get("gene") or "?"
            et = ev.get("event_type") or ev.get("type") or "event"
            val = ev.get("value")
            suffix = f" ({val})" if val is not None else ""
            lines.append(f"  - {gene}: {et}{suffix}")

    # Optional expression summary (typically used for persistence/tolerance hypotheses)
    if cs.expression is not None:
        lines.append("")
        lines.append("## Expression summary (RNA-seq / signatures)")
        lines.append("")
        if cs.expression.platform:
            lines.append(f"- Platform: `{_md_escape(cs.expression.platform)}`")
        if cs.expression.contrast:
            lines.append(f"- Contrast: `{_md_escape(cs.expression.contrast)}`")
        if cs.expression.signature_scores:
            lines.append("- Signature scores:")
            for k, v in sorted(cs.expression.signature_scores.items()):
                try:
                    lines.append(f"  - `{k}`: {float(v):.3f}")
                except Exception:
                    lines.append(f"  - `{k}`: {_md_escape(str(v))}")
        else:
            lines.append("- Signature scores: *(none)*")

        if cs.expression.top_markers:
            lines.append("- Top signature markers (within-sample):")
            for m in cs.expression.top_markers[:10]:
                g = m.get("gene") or "?"
                zv = m.get("z")
                if zv is not None:
                    try:
                        lines.append(f"  - {g}: z={float(zv):.2f}")
                    except Exception:
                        lines.append(f"  - {g}: z={_md_escape(str(zv))}")
                else:
                    lines.append(f"  - {g}")

    # Flags
    lines.append("")
    lines.append("## Engineered flags")
    lines.append("")
    if cs.genomic.flags:
        for k in sorted(cs.genomic.flags.keys()):
            lines.append(f"- `{k}` = {cs.genomic.flags[k]}")
    else:
        lines.append("- *(no flags set)*")

    # Mechanism calls
    lines.append("")
    lines.append("## Mechanism calls (rule engine)")
    lines.append("")
    if cs.mechanism_calls:
        for mc in sorted(cs.mechanism_calls, key=lambda x: x.score, reverse=True):
            lines.append(f"- **{mc.mechanism}** — score {mc.score:.2f} — {_md_escape(mc.rationale)}")
    else:
        lines.append("- *(no mechanism calls produced)*")

    # Routing
    lines.append("")
    lines.append("## Strategy routing (non-prescriptive)")
    lines.append("")
    if cs.routing:
        if cs.routing.ranked_buckets:
            lines.append("Ranked buckets:")
            for b in cs.routing.ranked_buckets:
                lines.append(f"- `{b}`")
        if cs.routing.what_to_test_next:
            lines.append("")
            lines.append("What to test next (research):")
            for w in cs.routing.what_to_test_next:
                lines.append(f"- {_md_escape(w)}")
        if cs.routing.what_to_avoid:
            lines.append("")
            lines.append("Guardrails:")
            for w in cs.routing.what_to_avoid:
                lines.append(f"- {_md_escape(w)}")
    else:
        lines.append("- *(no routing computed)*")

    # Evidence ledger (small)
    lines.append("")
    lines.append("## Evidence ledger (top items)")
    lines.append("")
    if cs.evidence:
        for ev in cs.evidence[:25]:
            lines.append(f"- [{ev.id}] ({ev.level}/{ev.kind}) {_md_escape(ev.label)}")
    else:
        lines.append("- *(no evidence items)*")

    return "\n".join(lines) + "\n"


def write_dossier_files(cs: CaseSnapshot, outdir: Path) -> None:
    """Write both JSON and Markdown outputs for a CaseSnapshot."""
    outdir.mkdir(parents=True, exist_ok=True)

    # JSON
    json_path = outdir / f"{cs.case_id}.case_snapshot.json"
    json_path.write_text(json.dumps(cs.to_dict(), indent=2, ensure_ascii=False), encoding="utf-8")

    # Markdown
    md_path = outdir / f"{cs.case_id}.dossier.md"
    md_path.write_text(render_markdown_dossier(cs), encoding="utf-8")


def write_dossier_bundle(cs: CaseSnapshot, out_root: Path) -> dict:
    """Write dossier outputs for a single snapshot.

    IMPORTANT: tests expect a *flat* layout under out_root:
      out_root/<CASE>.dossier.md
      out_root/<CASE>.case_snapshot.json

    (No per-case subdirectory.)
    """
    out_root.mkdir(parents=True, exist_ok=True)

    md_path = out_root / f"{cs.case_id}.dossier.md"
    json_path = out_root / f"{cs.case_id}.case_snapshot.json"

    # render + write
    md_path.write_text(render_markdown_dossier(cs), encoding="utf-8")
    json_path.write_text(json.dumps(cs.to_dict(), indent=2, ensure_ascii=False), encoding="utf-8")

    record = {
        "case_id": cs.case_id,
        "study_id": cs.study_id,
        "sample_id": cs.sample_id,
        "timepoint_id": cs.timepoint_id,
        "case_dir": str(out_root.resolve()),
        "dossier_md": str(md_path.resolve()),
        "snapshot_json": str(json_path.resolve()),
    }
    return record


# ---- Back-compat aliases + index writer ----

def build_json_dossier(cs, calls=None, routing=None):
    return build_dossier_json(cs)

def build_markdown_dossier(cs, calls=None, routing=None):
    return build_dossier_markdown(cs)

def build_dossier_markdown(cs: CaseSnapshot) -> str:
    """Back-compat alias used by older code/tests."""
    return render_markdown_dossier(cs)

def build_dossier_json(cs: CaseSnapshot) -> dict:
    """Machine-readable dossier payload (JSON-serializable)."""
    # Keep structure stable for tests; the snapshot already contains the evidence ledger.
    return cs.to_dict()

def write_dossier_index(records: list[dict], index_path: Path) -> Path:
    """Write an index.json that wraps records under a `cases` key."""
    index_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"cases": records}
    index_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    return index_path
