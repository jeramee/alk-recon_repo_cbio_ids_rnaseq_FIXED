from __future__ import annotations

from schema.case_snapshot import CaseSnapshot

from features.alk_flags import apply_on_target_flags
from features.bypass_flags import apply_bypass_flags
from features.persistence_flags import apply_persistence_flags


def apply_all_features(cs: CaseSnapshot, *, persister_threshold: float = 1.0) -> None:
    apply_on_target_flags(cs)
    apply_bypass_flags(cs)
    apply_persistence_flags(cs, threshold=persister_threshold)
