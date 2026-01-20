from __future__ import annotations

from schema.case_snapshot import CaseSnapshot

from .alk_flags import apply_alk_flags
from .bypass_flags import apply_bypass_flags
from .persistence_flags import apply_persistence_flags


def apply_all_features(cs: CaseSnapshot) -> CaseSnapshot:
    """Mutates CaseSnapshot in-place and returns it for convenience."""
    apply_alk_flags(cs)
    apply_bypass_flags(cs)
    apply_persistence_flags(cs)
    return cs

# --- Backwards-compatible alias for older tests / docs ---
def apply_features_to_snapshot(cs):
    """Compatibility wrapper. Older tests import this name."""
    return apply_all_features(cs)


