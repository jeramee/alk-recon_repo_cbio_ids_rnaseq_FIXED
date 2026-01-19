import pandas as pd
from pathlib import Path

BUNDLE = Path(r"data\1\luad_tcga_bundle\luad_tcga")
OUT    = Path(r"data\real\variants_cbio_luad_tcga.tsv")
OUT.parent.mkdir(parents=True, exist_ok=True)

def pick(df, candidates):
    low = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    return None

rows = []

# ---- Mutations (MAF-like) ----
mut_file = None
for name in ["data_mutations_extended.txt", "data_mutations.txt", "data_mutations_extended.maf", "data_mutations.maf"]:
    p = BUNDLE / name
    if p.exists():
        mut_file = p
        break

if mut_file:
    m = pd.read_csv(mut_file, sep="\t", comment="#", low_memory=False)

    gene_c   = pick(m, ["Hugo_Symbol", "gene", "Gene"])
    samp_c   = pick(m, ["Tumor_Sample_Barcode", "sample_id", "Sample_ID"])
    pat_c    = pick(m, ["Patient_ID", "case_id", "patient"])
    prot_c   = pick(m, ["HGVSp_Short", "Protein_Change", "Amino_Acid_Change"])
    vc_c     = pick(m, ["Variant_Classification"])

    for _, r in m.iterrows():
        gene = str(r.get(gene_c, "")).strip()
        sid  = str(r.get(samp_c, "")).strip()
        if not gene or not sid:
            continue

        case = str(r.get(pat_c, "")).strip() if pat_c else ""
        if not case:
            case = sid[:12]  # TCGA patient barcode fallback

        prot = str(r.get(prot_c, "")).strip() if prot_c else ""
        notes = str(r.get(vc_c, "")).strip() if vc_c else ""

        # Keep it focused for now: ALK + common bypass genes
        if gene not in {"ALK","MET","EGFR","KRAS","BRAF","ERBB2","ROS1","RET"}:
            continue

        rows.append({
            "study_id": "TCGA-LUAD",
            "case_id": case,
            "sample_id": sid,
            "timepoint_id": "T0",
            "gene": gene,
            "variant_type": "SNV",
            "protein_change": prot,
            "fusion_partner": "",
            "copy_number": "",
            "vaf": "",
            "notes": notes
        })

# ---- Fusions ----
fus_file = BUNDLE / "data_fusions.txt"
if fus_file.exists():
    f = pd.read_csv(fus_file, sep="\t", comment="#", low_memory=False)

    samp_c = pick(f, ["Tumor_Sample_Barcode", "Sample_ID", "sample_id"])
    pat_c  = pick(f, ["Patient_ID", "case_id"])
    fus_c  = pick(f, ["Fusion", "fusion", "Event_Info", "event_info"])

    # If there isn't a single "Fusion" string column, we try common gene1/gene2 patterns
    g1 = pick(f, ["Site1_Hugo_Symbol", "Gene_1", "gene1", "Hugo_Symbol"])
    g2 = pick(f, ["Site2_Hugo_Symbol", "Gene_2", "gene2", "Partner"])

    for _, r in f.iterrows():
        sid = str(r.get(samp_c, "")).strip() if samp_c else ""
        if not sid:
            continue
        case = str(r.get(pat_c, "")).strip() if pat_c else ""
        if not case:
            case = sid[:12]

        fusion_str = str(r.get(fus_c, "")).strip() if fus_c else ""
        a = b = ""

        if fusion_str and ("-" in fusion_str):
            parts = fusion_str.replace("--","-").split("-")
            if len(parts) >= 2:
                a, b = parts[0].strip(), parts[1].strip()
        elif g1 and g2:
            a, b = str(r.get(g1,"")).strip(), str(r.get(g2,"")).strip()

        # Normalize to ALK-centric row
        if "ALK" not in {a,b}:
            continue

        partner = b if a == "ALK" else a
        rows.append({
            "study_id": "TCGA-LUAD",
            "case_id": case,
            "sample_id": sid,
            "timepoint_id": "T0",
            "gene": "ALK",
            "variant_type": "FUSION",
            "protein_change": "",
            "fusion_partner": partner,
            "copy_number": "",
            "vaf": "",
            "notes": fusion_str or f"{a}-{b}"
        })

df = pd.DataFrame(rows).drop_duplicates()
df.to_csv(OUT, sep="\t", index=False)
print("Wrote:", OUT, "rows:", len(df))
