#!/usr/bin/env python3
"""
03c_align_and_codon.py

Replacement for HTT pipeline steps 3+4.

Root cause of original failure: the OrthoFinder MSA has trimmed/gapped sequences
that don't match the full CDS length, causing PAL2NAL inconsistency errors.

Fix: translate CDS -> protein (guarantees protein/CDS consistency), align 2
full proteins with MAFFT, then run PAL2NAL on aligned protein + original CDS.

Usage:
    python3 03c_align_and_codon.py [--threads N]
"""

import argparse
import json
import logging
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

RESULTS    = Path("/legserv/Temp/Carlos/EDF-TE/results/HTT/03_ks_divergence")
BUSCO_RUNS = Path("/legserv/Temp/Carlos/EDF-TE/results/HTT/01_phylo_tree/busco_runs")
ALIGN_DIR  = RESULTS / "alignments"
CODON_DIR  = RESULTS / "codon_alignments"
PAIRS_FILE = RESULTS / "node_analyses" / "node_gene_pairs.json"
TMPDIR     = Path("/legserv/Temp/Carlos/tmp_R/03c_tmp")

MAFFT    = "/home/carlos/micromamba/envs/esptreeEDFTE/bin/mafft"
PAL2NAL  = "/home/carlos/micromamba/envs/HTT_pal2nal/bin/pal2nal.pl"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("/data/carlos/EDF-TE/logs/HTT/03c_align_codon.log"),
    ],
)
log = logging.getLogger()

# Standard genetic code
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate(cds: str) -> str:
    """Translate CDS nucleotides to protein, stopping at stop codon."""
    cds = cds.upper().replace('-', '').replace('\n', '')
    protein = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i + 3]
        if len(codon) < 3:
            break
        aa = GENETIC_CODE.get(codon, 'X')
        if aa == '*':
            break
        protein.append(aa)
    return ''.join(protein)


def read_fasta_first(path: Path) -> str:
    """Return sequence of the first record in a FASTA file (no gaps)."""
    seq_lines = []
    in_seq = False
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if in_seq:
                    break
                in_seq = True
            elif in_seq:
                seq_lines.append(line.replace('-', ''))
    return ''.join(seq_lines)


def process_pair(pair: dict, sample_busco_dirs: dict, tmpdir: Path) -> str:
    """
    Process one node-gene pair: MAFFT alignment + PAL2NAL codon alignment.
    Returns 'ok' or an error message.
    """
    node_name = pair["node_name"]
    gene      = pair["gene"]
    reps      = pair["representatives"]  # [[name, msa_seq], ...]

    protein_seqs: dict[str, str] = {}
    cds_seqs:     dict[str, str] = {}

    for rep_name, _msa_seq in reps:
        if rep_name not in sample_busco_dirs:
            return f"no_busco:{rep_name}"
        nuc_file = sample_busco_dirs[rep_name] / f"{gene}.fna"
        if not nuc_file.exists():
            return f"no_nuc:{rep_name}:{gene}"
        cds  = read_fasta_first(nuc_file)
        prot = translate(cds)
        if len(prot) < 10:
            return f"short_prot:{rep_name}:{gene}:{len(prot)}"
        protein_seqs[rep_name] = prot
        cds_seqs[rep_name]     = cds

    if len(protein_seqs) < 2:
        return "too_few_reps"

    # Write protein FASTA for MAFFT
    prot_file  = tmpdir / f"{node_name}_{gene}_prot.faa"
    align_file = ALIGN_DIR / f"{node_name}_{gene}_aligned.faa"
    nuc_tmp    = tmpdir / f"{node_name}_{gene}_nuc.fna"
    codon_file = CODON_DIR / f"{node_name}_{gene}_codon_aligned.fasta"

    try:
        with open(prot_file, 'w') as fh:
            for name, seq in protein_seqs.items():
                fh.write(f">{name}\n{seq}\n")

        # MAFFT alignment (2 sequences — very fast)
        mafft_result = subprocess.run(
            [MAFFT, '--auto', '--quiet', str(prot_file)],
            capture_output=True, text=True, check=False
        )
        prot_file.unlink(missing_ok=True)

        if mafft_result.returncode != 0 or not mafft_result.stdout:
            return f"mafft_fail:{mafft_result.returncode}"

        align_file.write_text(mafft_result.stdout)

        # Write nucleotide sequences (no gaps, no stop codon)
        with open(nuc_tmp, 'w') as fh:
            for name, cds in cds_seqs.items():
                # Remove stop codon if present (last 3 nt)
                if cds[-3:].upper() in ('TAA', 'TAG', 'TGA'):
                    cds = cds[:-3]
                fh.write(f">{name}\n{cds}\n")

        # PAL2NAL codon alignment
        pal2nal_result = subprocess.run(
            [PAL2NAL, str(align_file), str(nuc_tmp), '-output', 'fasta', '-nogap'],
            capture_output=True, text=True, check=False
        )
        nuc_tmp.unlink(missing_ok=True)

        if pal2nal_result.returncode != 0:
            return f"pal2nal_fail:{pal2nal_result.returncode}"

        stdout = pal2nal_result.stdout
        if not stdout or 'ERROR' in stdout or not stdout.strip().startswith('>'):
            return f"pal2nal_empty_or_error"

        codon_file.write_text(stdout)
        if codon_file.stat().st_size == 0:
            return "codon_empty"

        return "ok"

    except Exception as exc:
        for p in [prot_file, nuc_tmp]:
            p.unlink(missing_ok=True)
        return f"exception:{exc}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    ALIGN_DIR.mkdir(parents=True, exist_ok=True)
    CODON_DIR.mkdir(parents=True, exist_ok=True)
    TMPDIR.mkdir(parents=True, exist_ok=True)

    # Verify tools
    for tool, path in [("mafft", MAFFT), ("pal2nal.pl", PAL2NAL)]:
        if not Path(path).exists():
            log.error(f"{tool} not found at {path}")
            sys.exit(1)
    log.info(f"MAFFT:    {MAFFT}")
    log.info(f"PAL2NAL:  {PAL2NAL}")

    # Build species → busco CDS dir mapping
    sample_busco_dirs: dict[str, Path] = {}
    for busco_dir in BUSCO_RUNS.glob("*_dothideomycetes_odb10"):
        sp = busco_dir.name.replace("_dothideomycetes_odb10", "")
        seq_dir = (busco_dir / "run_dothideomycetes_odb10"
                             / "busco_sequences"
                             / "single_copy_busco_sequences")
        if seq_dir.exists():
            sample_busco_dirs[sp] = seq_dir
    log.info(f"Species with CDS: {len(sample_busco_dirs)}")

    # Load pairs
    with open(PAIRS_FILE) as fh:
        pairs = json.load(fh)
    log.info(f"Total node-gene pairs: {len(pairs)}")

    success = 0
    failed  = 0
    error_counts: dict[str, int] = {}

    log.info(f"Processing with {args.threads} thread(s)…")

    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = {
            pool.submit(process_pair, pair, sample_busco_dirs, TMPDIR): i
            for i, pair in enumerate(pairs)
        }
        for i, fut in enumerate(as_completed(futures)):
            result = fut.result()
            if result == "ok":
                success += 1
            else:
                failed += 1
                key = result.split(':')[0]
                error_counts[key] = error_counts.get(key, 0) + 1
            if (i + 1) % 2000 == 0:
                log.info(f"  [{i+1}/{len(pairs)}] success={success} failed={failed}")

    log.info(f"Done. Success={success} | Failed={failed}")
    if error_counts:
        for k, v in sorted(error_counts.items(), key=lambda x: -x[1]):
            log.info(f"  error [{k}]: {v}")
    log.info(f"Codon alignments: {CODON_DIR}")


if __name__ == "__main__":
    main()
