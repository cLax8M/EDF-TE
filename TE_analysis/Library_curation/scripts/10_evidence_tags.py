#!/usr/bin/env python3
"""
Add evidence tags to each sequence in the quality-filtered library:
  _PE  = Protein Evidence (domain hit confirms classification)
  _DA  = Disagreement (domain hit contradicts classification)
  _NE  = No Evidence (no domain hits)

Classification is read from the FASTA header (e.g. #LTR/Gypsy, #DNA/TIR, #Unknown)
Domain-to-class mapping follows MycoMobilome Table 1.
"""
import re, sys
from Bio import SeqIO

WORK = "/data/carlos/EDF-TE/TE_analysis/Library_curation"
IN_FA   = f"{WORK}/09_quality/consensus_quality_filtered.fa"
EVID_F  = f"{WORK}/07_decision/domain_evidence.tsv"
OUT_FA  = f"{WORK}/10_tagged/consensus_tagged.fa"

import os; os.makedirs(f"{WORK}/10_tagged", exist_ok=True)

# Domain → high-level TE class mapping (REPET profile names contain these substrings)
CLASS_DOMAINS = {
    "DNA":        {"Tase", "DDE", "HTH", "ATP", "INT", "AP"},
    "RC":         {"HEL", "EN", "RPA"},
    "LTR":        {"RT", "INT", "RH", "GAG", "AP", "VirusRelated", "LTRrelated",
                   "Caulimovirus", "ClassIrelated", "ENV"},
    "LINE":       {"RT", "EN", "RH", "GAG", "ClassIrelated", "LINErelated"},
    "PLE":        {"RT", "EN", "ClassIrelated"},
    "Retroposon": {"RT", "INT", "RH", "GAG", "AP", "VirusRelated", "LTRrelated",
                   "Caulimovirus", "ClassIrelated", "ENV", "EN", "LINErelated"},
}
# RepeatPeps hit name prefixes that map to high-level class
RPEP_CLASS = {
    "DNA": "DNA",
    "RC":  "RC",
    "LTR": "LTR",
    "LINE": "LINE",
}

def get_seq_class(header):
    """Extract high-level class from RepeatMasker-style classification in header."""
    m = re.search(r'#(\w+)', header)
    if not m:
        return "Unknown"
    top = m.group(1).upper()
    for cls in CLASS_DOMAINS:
        if top.startswith(cls.upper()):
            return cls
    return "Unknown"

def domain_supports_class(domain_str, seq_class):
    """Check if any domain string suggests the expected class."""
    if not domain_str:
        return None   # no evidence
    expected = CLASS_DOMAINS.get(seq_class, set())
    for d in domain_str.split(";"):
        if any(e.upper() in d.upper() for e in expected):
            return True
    return False

def rpep_supports_class(rpep_str, seq_class):
    if not rpep_str:
        return None
    expected_prefix = seq_class.upper()
    for hit in rpep_str.split(";"):
        for prefix, cls in RPEP_CLASS.items():
            if cls.upper() == expected_prefix and prefix.upper() in hit.upper():
                return True
    return False

# Load evidence table
evidence = {}
with open(EVID_F) as f:
    next(f)  # header
    for line in f:
        parts = line.rstrip("\n").split("\t")
        seq_id = parts[0]
        evidence[seq_id] = {
            "hmmer": parts[1],
            "rpeps": parts[2],
            "decision": parts[4]
        }

with open(OUT_FA, "w") as out:
    for rec in SeqIO.parse(IN_FA, "fasta"):
        seq_id = rec.id
        seq_class = get_seq_class(rec.description)
        ev = evidence.get(seq_id, {"hmmer": "", "rpeps": "", "decision": "keep_no_evidence"})

        hmmer_support = domain_supports_class(ev["hmmer"], seq_class)
        rpep_support  = rpep_supports_class(ev["rpeps"], seq_class)

        has_any_domain = bool(ev["hmmer"] or ev["rpeps"])

        if not has_any_domain:
            tag = "_NE"
        elif hmmer_support or rpep_support:
            tag = "_PE"
        else:
            tag = "_DA"

        new_id = rec.id + tag
        out.write(f">{new_id} {rec.description.split(' ', 1)[1] if ' ' in rec.description else ''}\n")
        out.write(str(rec.seq) + "\n")

print(f"Tagging complete → {OUT_FA}")
