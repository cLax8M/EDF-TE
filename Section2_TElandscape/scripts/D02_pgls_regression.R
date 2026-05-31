#!/usr/bin/env Rscript
# D02_pgls_regression.R
# PGLS regression: TE class percentages ~ genomic features.
# Uses Brownian motion covariance (ape::corBrownian + nlme::gls).
# Bonferroni correction across all tests.
#
# Usage:
#   Rscript D02_pgls_regression.R <merged_tsv> <tree_newick> <output_dir>
#
# Requires: ape, nlme, caper (available or add to treewasEDFTE env)

suppressPackageStartupMessages({
  library(ape)
  library(nlme)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: D02_pgls_regression.R <merged_tsv> <tree_newick> <output_dir>")
}

merged_file <- args[1]
tree_file   <- args[2]
output_dir  <- args[3]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── load ──────────────────────────────────────────────────────
df   <- read.table(merged_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
tree <- read.tree(tree_file)

# ── prune tree ────────────────────────────────────────────────
common      <- intersect(df$species, tree$tip.label)
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, common))
df          <- df[df$species %in% common, ]
rownames(df) <- df$species

# ── define response and predictor variables ───────────────────
te_responses <- c("pct_TE_total", "pct_DNA", "pct_LTR",
                  "pct_LINE", "pct_RC", "pct_PLE", "pct_Unknown")

genomic_predictors <- c("Genome_Size_bp", "Pct_Intronic", "Pct_Intergenic",
                        "Avg_Introns_Multiexon", "Total_Genes",
                        "Avg_Exons_Per_Gene")

# ── run PGLS for each response ~ predictor pair ───────────────
results <- list()
n_tests <- 0

for (response in te_responses) {
  if (!(response %in% names(df))) next
  y <- df[[response]]

  for (pred in genomic_predictors) {
    if (!(pred %in% names(df))) next
    x <- df[[pred]]

    # drop rows with NA in either variable
    keep <- !is.na(y) & !is.na(x) & is.finite(y) & is.finite(x)
    if (sum(keep) < 10) next

    df_sub  <- df[keep, ]
    tr_sub  <- drop.tip(tree_pruned, setdiff(tree_pruned$tip.label, rownames(df_sub)))
    cor_bm  <- corBrownian(phy = tr_sub)

    tryCatch({
      fit <- gls(
        as.formula(paste(response, "~", pred)),
        data       = df_sub,
        correlation = cor_bm,
        method     = "ML"
      )
      s  <- summary(fit)
      ct <- coef(s)
      # row for the predictor (row 2 — row 1 is intercept)
      if (nrow(ct) >= 2) {
        n_tests <- n_tests + 1
        results[[paste(response, pred, sep = "~")]] <- data.frame(
          response  = response,
          predictor = pred,
          estimate  = round(ct[2, "Value"], 6),
          se        = round(ct[2, "Std.Error"], 6),
          t_value   = round(ct[2, "t-value"], 4),
          p_raw     = signif(ct[2, "p-value"], 4),
          n         = sum(keep),
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) {
      message("[WARN] PGLS failed for ", response, " ~ ", pred, ": ", conditionMessage(e))
    })
  }
}

if (length(results) == 0) {
  stop("[ERROR] No PGLS tests completed. Check merged_tsv columns and tree tip labels.")
}

out_df <- do.call(rbind, results)

# Bonferroni correction
out_df$p_bonferroni <- pmin(out_df$p_raw * n_tests, 1.0)
out_df$significant  <- out_df$p_bonferroni < 0.05

out_df <- out_df[order(out_df$p_raw), ]

out_file <- file.path(output_dir, "pgls_results.tsv")
write.table(out_df, out_file, sep = "\t", row.names = FALSE, quote = FALSE)

message("[OK] PGLS results: ", out_file)
message("[INFO] Tests run: ", n_tests, " | Significant (Bonferroni p<0.05): ",
        sum(out_df$significant, na.rm = TRUE))
print(head(out_df[out_df$significant, ], 20))
