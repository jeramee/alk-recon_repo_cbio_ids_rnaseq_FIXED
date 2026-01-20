"""Bypass-pathway flag computation.

Tests rely on these keys:
- has_MET_amp: True for MET amplification / AMP CNA events
- has_bypass_driver: True when any bypass driver is present

We keep the logic permissive because different ingest sources encode CNA/AMP
in different fields.
"""

from __future__ import annotations

from schema.case_snapshot import CaseSnapshot

# A lightweight set of common bypass drivers. This is not meant to be exhaustive.
BYPASS_DRIVER_GENES = {
    "MET",
    "EGFR",
    "ERBB2",
    "HER2",
    "KRAS",
    "NRAS",
    "BRAF",
    "MAP2K1",
    "PIK3CA",
    "KIT",
    "RET",
    "ROS1",
    "NTRK1",
    "NTRK2",
    "NTRK3",
}


def _is_amp_event(ev: dict) -> bool:
    """Return True if a bypass/cna event looks like an amplification."""
    effect = str(
        ev.get("effect")
        or ev.get("event_type")
        or ev.get("type")
        or ev.get("variant_type")
        or ""
    ).upper()
    cn_raw = ev.get("copy_number") or ev.get("cn") or ""
    cn = str(cn_raw).strip()

    if "AMP" in effect or "AMPL" in effect:
        return True

    if cn and cn.lower() not in {"nan", "none", "null"}:
        # If a copy-number is present, treat it as a likely CNA event.
        # Prefer a numeric threshold when possible.
        try:
            return float(cn) >= 6.0
        except Exception:
            return True

    return False


def apply_bypass_flags(cs: CaseSnapshot) -> CaseSnapshot:
    """Update cs.genomic.flags with bypass-driver booleans."""
    if cs.genomic is None:
        return cs

    flags = dict(cs.genomic.flags or {})

    bypass_events = list(cs.genomic.bypass_events or [])
    cna_events = list(getattr(cs.genomic, "copy_number_events", []) or [])

    has_met_amp = False
    for ev in bypass_events:
        gene = str(ev.get("gene") or "").upper()
        if gene == "MET" and _is_amp_event(ev):
            has_met_amp = True
            break

    if not has_met_amp:
        for ev in cna_events:
            gene = str(ev.get("gene") or "").upper()
            if gene == "MET" and _is_amp_event(ev):
                has_met_amp = True
                break

    flags["has_MET_amp"] = has_met_amp

    # Compatibility aliases used by the rule engine
    flags["has_met_amp_or_high"] = has_met_amp
    flags["has_met_alt"] = has_met_amp

    has_bypass_driver = False
    for ev in bypass_events:
        gene = str(ev.get("gene") or "").upper()
        if gene in BYPASS_DRIVER_GENES:
            has_bypass_driver = True
            break

    if not has_bypass_driver:
        for ev in cna_events:
            gene = str(ev.get("gene") or "").upper()
            if gene in BYPASS_DRIVER_GENES:
                has_bypass_driver = True
                break

    flags["has_bypass_driver"] = bool(has_bypass_driver or has_met_amp)
    # Back-compat naming used elsewhere
    flags.setdefault("has_bypass_event", bool(bypass_events or cna_events))

    # Existing semantic flag (kept for backward compatibility)
    flags.setdefault("has_bypass_event", bool(bypass_events or cna_events))

    cs.genomic.flags = flags
    return cs
