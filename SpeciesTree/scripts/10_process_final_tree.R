#!/usr/bin/env Rscript
# 10_process_final_tree.R
#
# PURPOSE: Post-process the IQ-TREE output treefile for use in downstream analyses.
#
# WHAT THIS DOES:
#   1. Reads the raw IQ-TREE treefile (*.treefile from 09_final_tree/)
#   2. Cleans tip labels: removes ".aa" suffixes (e.g. MU_Mucmuc1.aa -> MU_Mucmuc1)
#      Exception: CH_Batde5 has no .aa suffix — handled explicitly (no-op)
#   3. Saves the cleaned non-ultrametric tree (for analyses that accept branch lengths,
#      e.g. treeWAS, some phytools functions)
#   4. Ultrametricizes via chronos() for analyses requiring equal root-to-tip distances:
#      - B02_pagel_lambda.R (Pagel's λ)
#      - D02_pgls_regression.R (PGLS with corBrownian)
#      - A01_ancestral_reconstruction.R (phytools::fastAnc)
#   5. Saves the ultrametric tree in newick format
#
# USAGE:
#   Rscript 10_process_final_tree.R <treefile> <output_dir>
#   Example:
#   Rscript SpeciesTree/scripts/10_process_final_tree.R \
#       SpeciesTree/09_final_tree/EDF-TE.treefile \
#       SpeciesTree/10_processed_tree/
#
# OUTPUTS (in output_dir/):
#   - final_tree_cleaned.nwk        — raw IQ-TREE tree, tip labels cleaned
#   - final_tree_ultrametric.nwk    — ultrametric tree (chronos, λ=1)
#   - tree_processing_summary.txt   — log: n tips, label check, ultrametric check
#
# DEPENDENCIES: ape (CRAN)

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript 10_process_final_tree.R <treefile> <output_dir>\n")
  quit(status = 1)
}

treefile   <- args[1]
output_dir <- args[2]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_lines <- c()
log <- function(msg) { cat(msg, "\n"); log_lines <<- c(log_lines, msg) }

# ── 1. Read tree ──────────────────────────────────────────────────────────────
log(paste("Reading tree:", treefile))
tree <- read.tree(treefile)
log(paste("Tips in raw tree:", length(tree$tip.label)))

# ── 2. Clean tip labels ───────────────────────────────────────────────────────
# Remove ".aa" suffix; CH_Batde5 has no suffix so sub() is a no-op for it
original_labels  <- tree$tip.label
tree$tip.label   <- sub("\\.aa$", "", tree$tip.label)
cleaned_labels   <- tree$tip.label

n_cleaned <- sum(original_labels != cleaned_labels)
log(paste("Labels cleaned (removed .aa suffix):", n_cleaned))
log(paste("Labels unchanged (no .aa suffix, e.g. CH_Batde5):", sum(original_labels == cleaned_labels)))

# Verify no duplicates after cleaning
if (anyDuplicated(tree$tip.label)) {
  stop("ERROR: Duplicate tip labels after cleaning — check treefile.")
}

# ── 3. Save cleaned non-ultrametric tree ─────────────────────────────────────
cleaned_nwk <- file.path(output_dir, "final_tree_cleaned.nwk")
write.tree(tree, file = cleaned_nwk)
log(paste("Saved cleaned tree:", cleaned_nwk))

# ── 4. Ultrametricize via chronos() ──────────────────────────────────────────
# chronos() uses penalized likelihood; lambda=1 is a good default for divergence time estimation
# This does NOT require calibration points — it just makes the tree clocklike
log("Ultrametricizing with chronos(lambda=1)...")
tree_ultra <- chronos(tree, lambda = 1, quiet = TRUE)

# Verify ultrametric
is_ultra <- is.ultrametric(tree_ultra, tol = 1e-6)
log(paste("Is ultrametric:", is_ultra))

if (!is_ultra) {
  warning("Tree is not perfectly ultrametric after chronos() — check branch lengths")
}

# ── 5. Save ultrametric tree ──────────────────────────────────────────────────
ultra_nwk <- file.path(output_dir, "final_tree_ultrametric.nwk")
write.tree(tree_ultra, file = ultra_nwk)
log(paste("Saved ultrametric tree:", ultra_nwk))

# ── 6. Summary ────────────────────────────────────────────────────────────────
log_lines <- c(log_lines,
  "",
  "── Downstream usage ──",
  "  B02_pagel_lambda.R       → use final_tree_ultrametric.nwk",
  "  D02_pgls_regression.R    → use final_tree_ultrametric.nwk",
  "  A01_ancestral_recon.R    → use final_tree_ultrametric.nwk",
  "  TreeWAS (01_run_treewas.R) → use final_tree_cleaned.nwk (branch lengths preserved)"
)

summary_file <- file.path(output_dir, "tree_processing_summary.txt")
writeLines(log_lines, summary_file)
log(paste("Summary written:", summary_file))

log("Done.")
