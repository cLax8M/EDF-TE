#!/usr/bin/env Rscript
# run_dnashape.R
# Compute DNA shape parameters for one genome FASTA using DNAshapeR.
# Called by run_nuclpred_array.sh per species.
#
# Usage:
#   Rscript run_dnashape.R <genome.fasta> <output_dir>
#
# Output dir will contain: genome.MGW, genome.Roll, genome.ProT,
#                          genome.HelT, genome.EP
# (5 shape files required by NuclPred generate_matrix.py)
#
# Requires: DNAshapeR (Bioconductor)
#   R -e "BiocManager::install('DNAshapeR')"

suppressPackageStartupMessages(library(DNAshapeR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: run_dnashape.R <genome.fasta> <output_dir>")
}

fasta_path <- args[1]
out_dir    <- args[2]

if (!file.exists(fasta_path)) stop(paste("FASTA not found:", fasta_path))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Copy FASTA to output dir so shape files land there
fasta_basename <- basename(fasta_path)
dest_fasta <- file.path(out_dir, fasta_basename)
if (!file.exists(dest_fasta)) file.symlink(fasta_path, dest_fasta)

message("[DNAshapeR] Computing shape parameters for: ", fasta_basename)
message("[DNAshapeR] Output dir: ", out_dir)

pred <- getShape(dest_fasta)

message("[OK] Shape files written:")
for (nm in names(pred)) {
  message("  ", file.path(out_dir, paste0(fasta_basename, ".", nm)))
}
