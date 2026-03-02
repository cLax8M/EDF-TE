#!/usr/bin/env python3
"""
Decide host vs TE for each consensus sequence.
Outputs: keep_ids.txt, reject_ids.txt, domain_evidence.tsv
Logic follows MycoMobilome:
- keep if: hits RepeatPeps OR hits HMMER (fseq_evalue<=0.001, bitscore>=50)
- reject if: hits RefSeq AND (no TE evidence OR TE evidence doesn't overlap RefSeq hit)
- Uncharacterized/hypothetical RefSeq hits are NOT treated as host
"""
import sys, os
from collections import defaultdict

WORK = "/data/carlos/EDF-TE/TE_analysis/Library_curation"
REFSEQ_F   = f"{WORK}/04_host_filter/refseq_hits.tsv"
HMMER_F    = f"{WORK}/05_te_domains/te_domains.tbl"
RPEPS_F    = f"{WORK}/06_repeatpeps/repeatpeps_hits.tsv"
OUT_DIR    = f"{WORK}/07_decision"

os.makedirs(OUT_DIR, exist_ok=True)

KEEP_KEYWORDS = {"uncharacterized", "hypothetical", "low quality", "predicted protein"}

# ── Load RefSeq hits (strip _1.._6 frame suffix from query ID) ───────────
host = {}   # query_base -> list of (qseqid, stitle, evalue)
with open(REFSEQ_F) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        qid, stitle, evalue = parts[0], parts[6] if len(parts) > 6 else "", float(parts[4])
        qbase = "_".join(qid.split("_")[:-1])  # remove _1.._6 frame suffix
        stitle_lower = stitle.lower()
        if any(kw in stitle_lower for kw in KEEP_KEYWORDS):
            continue   # treat as potential TE, ignore
        host.setdefault(qbase, []).append((qid, stitle, evalue))

# ── Load HMMER hits (filter: fseq_evalue <= 0.001, bitscore >= 50) ───────
hmmer = defaultdict(list)  # query_base -> list of domain names
with open(HMMER_F) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 10:
            continue
        target_name = parts[0]   # HMM profile (domain name)
        query_name  = parts[2]   # query sequence
        fseq_evalue = float(parts[4])
        fseq_bitscore = float(parts[5])
        if fseq_evalue <= 0.001 and fseq_bitscore >= 50:
            qbase = "_".join(query_name.split("_")[:-1])
            hmmer[qbase].append(target_name)

# ── Load RepeatPeps hits ──────────────────────────────────────────────────
repeatpeps = defaultdict(list)   # query_base -> list of (sseqid, qstart, qend)
with open(RPEPS_F) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        qid, sseqid = parts[0], parts[1]
        qstart, qend = int(parts[4]), int(parts[5])
        qbase = "_".join(qid.split("_")[:-1])
        repeatpeps[qbase].append((sseqid, qstart, qend))

# ── Collect all query bases ───────────────────────────────────────────────
all_queries = set(host) | set(hmmer) | set(repeatpeps)

keep = set()
reject = set()
evidence = {}  # qbase -> dict with evidence info

for qbase in all_queries:
    has_host    = qbase in host
    has_hmmer   = qbase in hmmer
    has_rpeps   = qbase in repeatpeps
    has_te      = has_hmmer or has_rpeps

    if has_te:
        keep.add(qbase)
        evidence[qbase] = {
            "hmmer_domains": ";".join(hmmer.get(qbase, [])),
            "repeatpeps_hits": ";".join(h[0] for h in repeatpeps.get(qbase, [])),
            "refseq_conflict": "yes" if has_host else "no",
            "decision": "keep_te_evidence"
        }
    elif has_host:
        reject.add(qbase)
        evidence[qbase] = {
            "hmmer_domains": "",
            "repeatpeps_hits": "",
            "refseq_conflict": "yes",
            "decision": "reject_host"
        }
    # sequences with no hits at all → keep (no evidence against)
    else:
        keep.add(qbase)
        evidence[qbase] = {
            "hmmer_domains": "",
            "repeatpeps_hits": "",
            "refseq_conflict": "no",
            "decision": "keep_no_evidence"
        }

with open(f"{OUT_DIR}/keep_ids.txt", "w") as o:
    for k in sorted(keep):
        o.write(k + "\n")

with open(f"{OUT_DIR}/reject_ids.txt", "w") as o:
    for r in sorted(reject):
        o.write(r + "\n")

with open(f"{OUT_DIR}/domain_evidence.tsv", "w") as o:
    o.write("seq_id\thmmer_domains\trepeatpeps_hits\trefseq_conflict\tdecision\n")
    for qbase in sorted(evidence):
        ev = evidence[qbase]
        o.write(f"{qbase}\t{ev['hmmer_domains']}\t{ev['repeatpeps_hits']}\t{ev['refseq_conflict']}\t{ev['decision']}\n")

print(f"Keep: {len(keep)}  |  Reject: {len(reject)}")
print(f"  - keep with TE evidence: {sum(1 for e in evidence.values() if e['decision']=='keep_te_evidence')}")
print(f"  - keep no evidence:      {sum(1 for e in evidence.values() if e['decision']=='keep_no_evidence')}")
print(f"  - reject (host):         {sum(1 for e in evidence.values() if e['decision']=='reject_host')}")
