# ALK-RECON Dossier — CASE_PERSIST

> Research/education use only. Not medical advice.

## Identity

- Level: `SampleLevel.TIMEPOINT`
- Study: `unspecified`
- Patient: `CASE_PERSIST`
- Sample: `S1`
- Timepoint: `unspecified`

## Genomics summary

- ALK fusion: **EML4-ALK**
- ALK variants: *(none recorded in this snapshot)*

## Expression summary (RNA-seq / signatures)

- Platform: `rnaseq_counts`
- Signature scores:
  - `persister_score`: 0.589
  - `signature_coverage`: 0.800
- Top signature markers (within-sample):
  - HSPA1A: z=0.71
  - SLC2A1: z=0.71
  - JUN: z=0.71
  - ALDH1A1: z=0.71
  - NGFR: z=0.71
  - FOS: z=0.71
  - HSPB1: z=0.71
  - AXL: z=0.71
  - ZEB1: z=0.71
  - VIM: z=0.71

## Engineered flags

- `has_G1202R` = False
- `has_L1196M` = False
- `has_MET_amp` = False
- `has_any_alk_mutation` = False
- `has_any_alk_variant` = False
- `has_bypass_driver` = False
- `has_bypass_event` = False
- `has_compound_alk_mutations` = False
- `has_met_alt` = False
- `has_met_amp_or_high` = False
- `has_persister_score` = True
- `has_persister_score_high` = False
- `persister_signature_score` = 0.5892556509887897
- `persister_signature_score_high` = False

## Mechanism calls (rule engine)

- **MechanismType.ON_TARGET_ALK** — score 0.00 — No ALK kinase-domain variant was detected in the provided variant table.
- **MechanismType.BYPASS** — score 0.00 — No obvious bypass-related alterations were detected in the provided variant table.
- **MechanismType.PERSISTENCE** — score 0.00 — No expression-level persistence evidence was provided (expression block is missing).

## Strategy routing (non-prescriptive)

Ranked buckets:
- `StrategyBucket.E_SEQUENCING_LOGIC`

What to test next (research):
- Acquire higher-yield evidence: ctDNA/tissue NGS for ALK mutations and key bypass events; optional expression signatures if available.

Guardrails:
- Do not treat this as a clinical recommendation; this is a research summary of evidence.

## Evidence ledger (top items)

- [E1] (EvidenceLevel.KNOWN/variant_row) Variant table row: ALK
- [F_has_any_alk_mutation] (EvidenceLevel.INFERRED/feature_flag) Derived flag has_any_alk_mutation
- [F_has_G1202R] (EvidenceLevel.INFERRED/feature_flag) Derived flag has_G1202R
- [F_has_L1196M] (EvidenceLevel.INFERRED/feature_flag) Derived flag has_L1196M
- [F_has_compound_alk_mutations] (EvidenceLevel.INFERRED/feature_flag) Derived flag has_compound_alk_mutations
