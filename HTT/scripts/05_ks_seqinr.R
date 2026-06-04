#!/usr/bin/env Rscript
# 05_ks_seqinr.R
# Calculate Ks divergence per gene pair and per node using seqinr::kaks()
# Usage: Rscript 05_ks_seqinr.R <results_dir> <min_genes>

suppressPackageStartupMessages({
    library(seqinr)
    library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("Usage: Rscript 05_ks_seqinr.R <results_dir> <min_genes>\n")
    quit(status = 1)
}

results_dir <- args[1]
min_genes   <- as.integer(args[2])

codon_dir <- file.path(results_dir, "codon_alignments")
ks_dir    <- file.path(results_dir, "ks_calculations")
dir.create(ks_dir, showWarnings = FALSE, recursive = TRUE)

codon_files <- list.files(codon_dir, pattern = "_codon_aligned\\.fasta$", full.names = TRUE)
cat("Found", length(codon_files), "codon alignments\n")

if (length(codon_files) == 0) {
    cat("ERROR: No codon alignments found in", codon_dir, "\n")
    quit(status = 1)
}

ks_results <- list()
n_ok    <- 0L
n_skip  <- 0L

for (i in seq_along(codon_files)) {
    codon_file <- codon_files[[i]]
    base_name  <- basename(codon_file)
    pair_name  <- sub("_codon_aligned\\.fasta$", "", base_name)

    if (i %% 500 == 0) {
        cat("Progress:", i, "/", length(codon_files),
            "- ok:", n_ok, "skip:", n_skip, "\n", flush = TRUE)
    }

    tryCatch({
        aln <- read.alignment(codon_file, format = "fasta")

        if (length(aln$seq) != 2) {
            n_skip <- n_skip + 1L
            return(NULL)
        }

        res <- kaks(aln)

        # seqinr >= 4.x returns list with $ks and $ka
        ks_val <- tryCatch(as.numeric(res$ks)[1], error = function(e) NA_real_)
        ka_val <- tryCatch(as.numeric(res$ka)[1], error = function(e) NA_real_)

        # Parse node and gene from pair name: <node1>_<node2>_<gene>
        tokens    <- strsplit(pair_name, "_", fixed = TRUE)[[1]]
        node_name <- paste(tokens[1:2], collapse = "_")
        gene_name <- paste(tokens[3:length(tokens)], collapse = "_")

        ks_results[[pair_name]] <- list(
            node  = node_name,
            gene  = gene_name,
            ks    = ks_val,
            ka    = ka_val,
            ka_ks = if (!is.na(ks_val) && is.finite(ks_val) && ks_val > 0)
                        ka_val / ks_val
                    else NA_real_
        )
        n_ok <- n_ok + 1L

    }, error = function(e) {
        n_skip <<- n_skip + 1L
    })
}

cat("Calculated Ks for", length(ks_results), "gene pairs\n")
cat("Skipped:", n_skip, "\n")

if (length(ks_results) == 0) {
    cat("ERROR: No Ks values calculated\n")
    quit(status = 1)
}

# ── Save raw per-pair Ks ────────────────────────────────────────────────────
write_json(ks_results, file.path(ks_dir, "ks_values.json"), pretty = TRUE)

ks_raw_df <- do.call(rbind, lapply(names(ks_results), function(pair) {
    data.frame(
        node = ks_results[[pair]]$node,
        gene = ks_results[[pair]]$gene,
        ks   = ks_results[[pair]]$ks,
        ka   = ks_results[[pair]]$ka,
        stringsAsFactors = FALSE
    )
}))
write.table(ks_raw_df,
            file.path(ks_dir, "ks_values_by_node.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── Node-level statistics ───────────────────────────────────────────────────
node_stats <- list()
for (result in ks_results) {
    node   <- result$node
    ks_val <- result$ks
    if (!is.na(ks_val) && is.finite(ks_val) && ks_val > 0) {
        node_stats[[node]] <- c(node_stats[[node]], ks_val)
    }
}

node_summary <- list()
for (node in names(node_stats)) {
    vals <- node_stats[[node]]
    if (length(vals) >= min_genes) {
        node_summary[[node]] <- list(
            n_genes         = length(vals),
            ks_mean         = mean(vals),
            ks_median       = median(vals),
            ks_0.5_quantile = as.numeric(quantile(vals, 0.005)),
            ks_1_quantile   = as.numeric(quantile(vals, 0.01)),
            ks_5_quantile   = as.numeric(quantile(vals, 0.05)),
            ks_min          = min(vals),
            ks_max          = max(vals)
        )
        cat("Node", node, "- Genes:", length(vals),
            "- Ks 0.5%:", round(quantile(vals, 0.005), 4), "\n")
    }
}

write_json(node_summary,
           file.path(ks_dir, "node_ks_summary.json"),
           pretty = TRUE)

if (length(node_summary) > 0) {
    summary_df <- do.call(rbind, lapply(names(node_summary), function(node) {
        s <- node_summary[[node]]
        data.frame(
            node            = node,
            n_genes         = s$n_genes,
            ks_min          = s$ks_min,
            ks_0.5_quantile = s$ks_0.5_quantile,
            ks_1_quantile   = s$ks_1_quantile,
            ks_5_quantile   = s$ks_5_quantile,
            ks_median       = s$ks_median,
            ks_mean         = s$ks_mean,
            ks_max          = s$ks_max,
            stringsAsFactors = FALSE
        )
    }))

    summary_df$pass_0.01  <- summary_df$ks_0.5_quantile > 0.01
    summary_df$pass_0.005 <- summary_df$ks_0.5_quantile > 0.005
    summary_df$pass_0.001 <- summary_df$ks_0.5_quantile > 0.001

    write.table(summary_df,
                file.path(ks_dir, "node_ks_summary.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)

    cat("\nNode summary (", nrow(summary_df), "nodes with >=", min_genes, "genes):\n")
    print(summary_df)
}

cat("\nKs calculation complete\n")
