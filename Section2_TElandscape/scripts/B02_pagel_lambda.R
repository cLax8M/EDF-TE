#!/usr/bin/env Rscript
# B02_pagel_lambda.R
# Compute Pagel's lambda (phylogenetic signal) for Shannon H' of TE diversity.
#
# Usage:
#   Rscript B02_pagel_lambda.R <species_shannon_tsv> <tree_newick> <output_tsv>
#
# Requires: ape, phytools (available in treewasEDFTE conda env)

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: B02_pagel_lambda.R <shannon_tsv> <tree_newick> <output_tsv>")
}

shannon_file <- args[1]
tree_file    <- args[2]
output_file  <- args[3]

# ── load data ─────────────────────────────────────────────────
df   <- read.table(shannon_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
tree <- read.tree(tree_file)

# ── prune tree to species with Shannon data ───────────────────
common <- intersect(df$species, tree$tip.label)
missing_tree <- setdiff(df$species, tree$tip.label)
missing_data <- setdiff(tree$tip.label, df$species)

if (length(missing_tree) > 0)
  message("[WARN] Species in data but not in tree (", length(missing_tree), "): ",
          paste(head(missing_tree, 5), collapse = ", "))
if (length(missing_data) > 0)
  message("[INFO] Tree tips pruned (no Shannon data): ", length(missing_data))

tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, common))
df_sub      <- df[df$species %in% common, ]
rownames(df_sub) <- df_sub$species

# Named vectors for phytools
h_vec <- setNames(df_sub$H_prime, df_sub$species)
j_vec <- setNames(as.numeric(df_sub$J_evenness), df_sub$species)
j_vec <- j_vec[!is.na(j_vec)]

results <- list()

for (trait_name in c("H_prime", "J_evenness")) {
  vec <- if (trait_name == "H_prime") h_vec else j_vec
  # Align to tree
  vec <- vec[intersect(names(vec), tree_pruned$tip.label)]
  tr  <- drop.tip(tree_pruned, setdiff(tree_pruned$tip.label, names(vec)))

  message("[INFO] Running phylosig (", trait_name, ") on ", length(vec), " species")
  ps <- phylosig(tr, vec, method = "lambda", test = TRUE)

  results[[trait_name]] <- data.frame(
    trait     = trait_name,
    lambda    = round(ps$lambda, 4),
    logL      = round(ps$logL, 4),
    logL0     = round(ps$logL0, 4),
    LR        = round(2 * (ps$logL - ps$logL0), 4),
    p_value   = signif(ps$P, 4),
    n_species = length(vec)
  )
}

out_df <- do.call(rbind, results)
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write.table(out_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

message("[OK] Pagel's lambda results: ", output_file)
print(out_df)
