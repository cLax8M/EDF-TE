#!/usr/bin/env Rscript
# 01_run_treewas.R
# ----------------
# Run treeWAS for one phenotype × one genotype matrix.
#
# Usage (called by 01_submit_treewas.sh):
#   Rscript 01_run_treewas.R \
#       --phenotype log10_Genome_Size \
#       --matrix    OG \
#       --datadir   /path/to/TreeWAS/results \
#       --outdir    /path/to/TreeWAS/results
#
# Outputs (<outdir>/<phenotype>_<matrix>/):
#   treeWAS_results.rds        full treeWAS S4 object
#   terminal_hits.tsv
#   simultaneous_hits.tsv
#   subsequent_hits.tsv
#   summary.tsv                union of hits, annotated with n_tests_significant

suppressPackageStartupMessages({
  library(optparse)
  library(treeWAS)
  library(ape)
  library(data.table)
})

# ── argument parsing ─────────────────────────────────────────────────────────
option_list <- list(
  make_option("--phenotype",  type = "character", help = "Phenotype column name"),
  make_option("--matrix",     type = "character", help = "OG or IPR"),
  make_option("--datadir",    type = "character", help = "Directory with prepared inputs"),
  make_option("--outdir",     type = "character", help = "Output base directory"),
  make_option("--n_sim",      type = "integer",   default = 10000L,
              help = "Number of null simulations [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$phenotype), !is.null(opt$matrix),
          !is.null(opt$datadir),   !is.null(opt$outdir))

matrix_label <- toupper(opt$matrix)
stopifnot(matrix_label %in% c("OG", "IPR"))

run_label <- paste0(opt$phenotype, "_", matrix_label)
cat(sprintf("[treeWAS] phenotype=%s  matrix=%s\n", opt$phenotype, matrix_label))

# ── output directory ─────────────────────────────────────────────────────────
run_outdir <- file.path(opt$outdir, run_label)
dir.create(run_outdir, showWarnings = FALSE, recursive = TRUE)

# ── load phenotype vector ────────────────────────────────────────────────────
pheno_path <- file.path(opt$datadir, "phenotypes.tsv")
pheno_df   <- as.data.frame(fread(pheno_path, header = TRUE))
rownames(pheno_df) <- pheno_df[[1]]
pheno_df   <- pheno_df[, -1, drop = FALSE]

if (!(opt$phenotype %in% colnames(pheno_df))) {
  stop(sprintf("Phenotype '%s' not found. Available: %s",
               opt$phenotype, paste(colnames(pheno_df), collapse = ", ")))
}

# treeWAS requires a named numeric vector with no NAs
phen_raw <- pheno_df[, opt$phenotype]
names(phen_raw) <- rownames(pheno_df)
phen <- phen_raw[!is.na(phen_raw)]
cat(sprintf("    Phenotype vector: %d species (dropped %d NAs)\n",
            length(phen), sum(is.na(phen_raw))))

# ── load genotype matrix ──────────────────────────────────────────────────────
x_file <- file.path(opt$datadir, paste0("X_", matrix_label, ".tsv"))
x_raw  <- as.data.frame(fread(x_file, header = TRUE))
rownames(x_raw) <- x_raw[[1]]
x_raw  <- x_raw[, -1, drop = FALSE]
cat(sprintf("    Genotype matrix: %d × %d\n", nrow(x_raw), ncol(x_raw)))

# ── load pruned tree ──────────────────────────────────────────────────────────
tree_path <- file.path(opt$datadir, "tree_pruned.nwk")
tree <- read.tree(tree_path)
cat(sprintf("    Tree: %d tips\n", length(tree$tip.label)))

# ── align species ─────────────────────────────────────────────────────────────
# intersect species present in all three: phen, X, tree
sp_common <- Reduce(intersect, list(names(phen), rownames(x_raw), tree$tip.label))
cat(sprintf("    Common species after alignment: %d\n", length(sp_common)))

if (length(sp_common) < 20) {
  stop("Too few common species (<20). Check that species names match across files.")
}

phen  <- phen[sp_common]
x_mat <- as.matrix(x_raw[sp_common, ])
tree  <- keep.tip(tree, sp_common)

# treeWAS expects snps as species × features (confirmed format)
# ensure numeric
storage.mode(x_mat) <- "integer"

# ── run treeWAS ──────────────────────────────────────────────────────────────
cat(sprintf("    Running treeWAS (n.snps.sim=%d) …\n", opt$n_sim))
set.seed(42)
res <- treeWAS(
  snps             = x_mat,
  phen             = phen,
  tree             = tree,
  n.snps.sim       = opt$n_sim,
  test             = c("terminal", "simultaneous", "subsequent"),
  p.value          = 0.01,
  p.value.correct  = "bonf",
  p.value.by       = "count",
  plot.tree        = FALSE,
  plot.manhattan   = FALSE,
  plot.null.dist   = FALSE,
  plot.dist        = FALSE
)

# ── save full object ──────────────────────────────────────────────────────────
rds_path <- file.path(run_outdir, "treeWAS_results.rds")
saveRDS(res, rds_path)
cat(sprintf("    Saved: %s\n", rds_path))

# ── extract hits per test ─────────────────────────────────────────────────────
extract_hits <- function(res_obj, test_name) {
  r <- res_obj[[test_name]]
  if (is.null(r)) return(data.frame())

  sig_df <- r$sig.snps
  # when no hits, treeWAS sets sig.snps = "No significant SNPs found." (character)
  if (is.null(sig_df) || is.character(sig_df) || !is.data.frame(sig_df) || nrow(sig_df) == 0) {
    return(data.frame())
  }

  data.frame(
    feature  = rownames(sig_df),
    p_value  = sig_df$p.value,
    score    = sig_df$score,
    test     = test_name,
    stringsAsFactors = FALSE
  )
}

tests <- c("terminal", "simultaneous", "subsequent")
hits_list <- lapply(tests, function(tt) extract_hits(res, tt))
names(hits_list) <- tests

for (tt in tests) {
  out_path <- file.path(run_outdir, paste0(tt, "_hits.tsv"))
  write.table(hits_list[[tt]], out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("    %s hits: %d  → %s\n", tt, nrow(hits_list[[tt]]), out_path))
}

# ── summary: robust hits (significant in ≥2 tests) ───────────────────────────
all_hits <- do.call(rbind, hits_list)
if (nrow(all_hits) > 0) {
  test_counts <- aggregate(test ~ feature, data = all_hits, FUN = length)
  colnames(test_counts)[2] <- "n_tests_significant"
  test_names  <- aggregate(test ~ feature, data = all_hits,
                            FUN = function(x) paste(sort(x), collapse = ","))
  colnames(test_names)[2] <- "tests_significant"
  pvals_agg   <- aggregate(p_value ~ feature, data = all_hits, FUN = min, na.rm = TRUE)
  colnames(pvals_agg)[2] <- "min_p_value"
  scores_agg  <- aggregate(score ~ feature, data = all_hits, FUN = mean, na.rm = TRUE)
  colnames(scores_agg)[2] <- "mean_score"

  summary_df <- Reduce(function(a, b) merge(a, b, by = "feature"),
                       list(test_counts, test_names, pvals_agg, scores_agg))
  summary_df$robust <- summary_df$n_tests_significant >= 2
  summary_df$phenotype <- opt$phenotype
  summary_df$matrix    <- matrix_label

  summary_path <- file.path(run_outdir, "summary.tsv")
  write.table(summary_df, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("    Robust hits (≥2 tests): %d\n", sum(summary_df$robust)))
  cat(sprintf("    Summary: %s\n", summary_path))
} else {
  cat("    No significant hits found.\n")
  summary_df <- data.frame()
  write.table(summary_df, file.path(run_outdir, "summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

cat(sprintf("[treeWAS] DONE: %s\n", run_label))
