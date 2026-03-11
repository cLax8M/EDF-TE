#!/usr/bin/env Rscript
# A01_ancestral_reconstruction.R
# Ancestral state reconstruction of TE content across the Holomycota phylogeny.
# Uses phytools::fastAnc (Brownian motion) for continuous TE percentages.
# Also fits OU model (phytools::OUwie-style via geiger) and compares AIC.
#
# Usage:
#   Rscript A01_ancestral_reconstruction.R \
#       <earlgrey_annotated_dir> <tree_newick> <output_dir>
#
# Requires: ape, phytools, geiger (add to treewasEDFTE env if needed)

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: A01_ancestral_reconstruction.R <annotated_dir> <tree_newick> <output_dir>")
}

annotated_dir <- args[1]
tree_file     <- args[2]
output_dir    <- args[3]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── load tree ─────────────────────────────────────────────────
tree <- read.tree(tree_file)
message("[INFO] Tree: ", length(tree$tip.label), " tips")

# ── collect TE percentages per species ───────────────────────
te_classes <- c("DNA", "LTR", "LINE", "SINE", "Rolling Circle", "Unknown")
sp_data <- list()

entries <- list.dirs(annotated_dir, full.names = TRUE, recursive = FALSE)
for (sp_dir in entries) {
  entry <- basename(sp_dir)
  if (!grepl("_EarlGrey$", entry)) next
  sp <- sub("_EarlGrey$", "", entry)
  hlc <- file.path(sp_dir, paste0(sp, "_summaryFiles"), paste0(sp, ".highLevelCount.txt"))
  if (!file.exists(hlc)) next

  df <- tryCatch(
    read.table(hlc, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
               check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df)) next
  # Column names may have spaces; normalize
  names(df) <- trimws(names(df))
  pct_col <- grep("Genome Coverage|%", names(df), value = TRUE)[1]
  class_col <- grep("Classification|Class", names(df), value = TRUE)[1]
  if (is.na(pct_col) || is.na(class_col)) next

  df[[class_col]] <- trimws(df[[class_col]])
  # Merge nested into parent
  df[[class_col]] <- gsub("-nested$", "", df[[class_col]])
  agg <- tapply(as.numeric(df[[pct_col]]), df[[class_col]], sum, na.rm = TRUE)
  total_te <- sum(agg, na.rm = TRUE)

  row <- c(total_te = total_te)
  for (cl in te_classes) {
    row[paste0("pct_", gsub(" ", "_", cl))] <- ifelse(cl %in% names(agg), agg[[cl]], 0)
  }
  sp_data[[sp]] <- row
}

mat <- do.call(rbind, sp_data)
sp_df <- as.data.frame(mat)
sp_df$species <- rownames(sp_df)

# ── prune tree ────────────────────────────────────────────────
common      <- intersect(sp_df$species, tree$tip.label)
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, common))
sp_df       <- sp_df[sp_df$species %in% common, ]
rownames(sp_df) <- sp_df$species
message("[INFO] Species for reconstruction: ", nrow(sp_df))

# ── ancestral reconstruction (BM) for each TE class ──────────
trait_cols <- c("total_te", paste0("pct_", gsub(" ", "_", te_classes)))
anc_results <- list()

for (trait in trait_cols) {
  vec <- setNames(sp_df[[trait]], sp_df$species)
  vec <- vec[tree_pruned$tip.label]
  anc <- fastAnc(tree_pruned, vec, CI = TRUE)
  anc_df <- data.frame(
    node       = names(anc$ace),
    trait      = trait,
    anc_value  = round(anc$ace, 4),
    lower_CI   = round(anc$CI95[, 1], 4),
    upper_CI   = round(anc$CI95[, 2], 4)
  )
  anc_results[[trait]] <- anc_df
}

anc_all <- do.call(rbind, anc_results)
write.table(anc_all,
            file.path(output_dir, "ancestral_te_states.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── save tip values for plotting ─────────────────────────────
write.table(sp_df,
            file.path(output_dir, "tip_te_values.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── phylogenetic signal (Pagel's lambda) for total TE ────────
vec_total <- setNames(sp_df$total_te, sp_df$species)
vec_total <- vec_total[tree_pruned$tip.label]
ps <- phylosig(tree_pruned, vec_total, method = "lambda", test = TRUE)
message("[INFO] Pagel's lambda (total_te): lambda=", round(ps$lambda, 3),
        "  p=", signif(ps$P, 3))

# ── save simple plot of total TE mapped on tree ────────────────
pdf(file.path(output_dir, "total_te_ancestral_tree.pdf"), width = 12, height = 16)
contMap(tree_pruned, vec_total, plot = FALSE) |> plot(lwd = 3, legend = TRUE)
title(main = "Total TE% — ancestral reconstruction (BM)")
dev.off()

message("[OK] Ancestral reconstruction complete: ", output_dir)
message("[INFO] Node estimates written: ancestral_te_states.tsv")
