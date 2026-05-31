#!/usr/bin/env python3
"""
Step 1: Extract sex locus windows from selected genomes for synteny figure.
Selected species:
  Mucci3   - scaffold_1   (JGI, sexP)
  mu402    - scaffold_1   (Extra, sexP Mucor)
  Phybl4   - scaffold_10  (JGI, sexM Phycomyces)
  PhyblU21 - scaffold_17  (Extra, sexP Phycomyces)
  Rarr_CCUG - JBSZQK010000002.1  (Extra, Rhizopus arrhizus)
  RhiMicr1  - CDSBDI010000001.1  (Extra, R. microsporus oligosporus)
"""
import subprocess
from pathlib import Path
import gzip

EXTRA   = Path("/legserv/Temp/Carlos/EDF-TE/data/Extra")
JGI     = Path("/legserv/Temp/Carlos/EDF-TE/data/genomes")
OUTDIR  = Path("/data/carlos/EDF-TE/TE_analysis/sex_locus/synteny_figure")
OUTDIR.mkdir(exist_ok=True)

BEDTOOLS = "/opt/easybuild/software/BEDTools/2.31.1-GCC-13.2.0/bin/bedtools"

# Genome definitions:  label -> (fasta, [(scaffold, sex_center, window_half)])
WINDOWS = 50000  # ±50kb around sex gene center

# (scaffold, sex_gene_center, strand_inversion)
LOCI = {
    "Mucci3":   (JGI/"MU_Mucci3.fasta",    "scaffold_1",                 2213094, "+"),
    "mu402":    (EXTRA/"mu402GenomeJGI.fa", "scaffold_1",                 2225952, "+"),
    "Phybl4":   (JGI/"MU_Phybl4.fasta",    "scaffold_10",                 726371, "+"),
    "PhyblU21": (EXTRA/"PhyblU21_2_AssemblyScaffolds.fasta.gz",
                                             "scaffold_17",                332700, "+"),
    "Rarr_CCUG":(EXTRA/"RarrCCUG61147.fna", "gb|JBSZQK010000002.1|",    2552119, "+"),
    "RhiMicr1": (EXTRA/"RhiMicr1.fna",      "emb|CDSBDI010000001.1|",   2979312, "+"),
}

def extract_fasta_region(fasta, scaffold, start, end, out_fasta):
    """Extract region from genome FASTA using bedtools. Handles gzipped inputs."""
    if str(fasta).endswith(".gz"):
        # decompress first
        tmp = OUTDIR / f"tmp_{Path(str(fasta).replace('.gz','')).name}"
        if not tmp.exists():
            with gzip.open(fasta,"rb") as fi, open(tmp,"wb") as fo:
                fo.write(fi.read())
        fasta = tmp

    # Create genome index if needed
    fai = Path(str(fasta) + ".fai")
    if not fai.exists():
        subprocess.run(["samtools", "faidx", str(fasta)], check=True,
                       capture_output=True)

    # Use samtools faidx region extraction
    region = f"{scaffold}:{start+1}-{end}"  # 1-based for samtools
    result = subprocess.run(
        ["samtools", "faidx", str(fasta), region],
        capture_output=True, text=True
    )
    if result.returncode != 0 or not result.stdout.strip():
        print(f"  WARNING: samtools faidx failed for {scaffold}:{start}-{end}, trying bedtools")
        # Fallback: bedtools getfasta
        bed_str = f"{scaffold}\t{start}\t{end}\n"
        bed_tmp = OUTDIR / "tmp.bed"
        bed_tmp.write_text(bed_str)
        result = subprocess.run(
            [BEDTOOLS, "getfasta", "-fi", str(fasta), "-bed", str(bed_tmp)],
            capture_output=True, text=True
        )
    seq = result.stdout
    # Rename header
    lines = seq.split("\n")
    label = out_fasta.stem
    lines[0] = f">{label}_{scaffold}:{start}-{end}"
    out_fasta.write_text("\n".join(lines))
    return len(seq) > 100

# Check samtools
samtools = subprocess.run(["which","samtools"], capture_output=True, text=True).stdout.strip()
if not samtools:
    samtools = "/opt/easybuild/software/SAMtools/1.20-GCC-13.2.0/bin/samtools"
    import os; os.environ["PATH"] += f":{Path(samtools).parent}"
print(f"samtools: {samtools}")

for label, (fasta, scaffold, center, strand) in LOCI.items():
    start = max(0, center - WINDOWS)
    end   = center + WINDOWS
    out   = OUTDIR / f"{label}_locus.fasta"
    print(f"Extracting {label}: {scaffold}:{start}-{end}")
    ok = extract_fasta_region(fasta, scaffold, start, end, out)
    print(f"  -> {out} {'OK' if ok else 'FAILED'}")

print("\nDone extracting windows.")
