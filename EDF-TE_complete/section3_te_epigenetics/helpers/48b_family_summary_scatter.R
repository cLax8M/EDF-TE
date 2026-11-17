#!/usr/bin/env Rscript
# =============================================================================
# 48b_family_summary_scatter.R
#
# Reproducible replacement for 48_family_summary_scatter_v2.pdf (Supp. Fig. 6).
#
# The v2 figure was generated in an interactive R session and its code was
# never saved, so the published panel could not be reproduced from the data
# on disk (the reported "n" did not match the TSV under any definition).
# This script rebuilds the analysis from scratch, with an explicit definition
# of every quantity, and writes BOTH the figure and the statistics table.
#
# WHAT IS COMPUTED
# ----------------
#   x = te_log2fc_meth  (log2 fold-change of 6mA/5mC over the TE copy)
#   y = te_log2fc_rna   (log2 fold-change of TE expression)
#
#   Primary statistic : Pearson correlation on ONE OBSERVATION PER TE COPY
#                       (log2FC averaged across contrasts), which avoids the
#                       pseudoreplication of counting the same TE once per
#                       contrast.
#   p_one             : one-tailed p-value (directional hypothesis r > 0)
#   95% CI            : Fisher z-transformation
#   n_copies          : number of independent TE copies tested   <-- use this
#   n_pairs           : number of TE-gene/contrast observations  (for reference)
#
# OUTPUTS (written to $OUTDIR)
#   ../../48b_family_stats.tsv          per-group, per-superfamily statistics
#   48b_family_summary_scatter.pdf      new figure (Supp. Fig. 6 replacement)
#   48b_family_summary_scatter.png      same, 300 dpi
#   48b_reproducibility.txt             sessionInfo() + input checksums
#
# Usage:
#   Rscript 48b_family_summary_scatter.R [input.tsv] [outdir]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
OUTDIR_DEF <- "/legserv/Temp/Carlos/EDF-TE/results/TE_dynamics/48_adjacent_genes"
INFILE <- if (length(args) >= 1) args[1] else file.path(OUTDIR_DEF, "48_adjacent_gene_expr.tsv")
OUTDIR <- if (length(args) >= 2) args[2] else OUTDIR_DEF

MIN_COPIES <- 10     # minimum independent TE copies to test a superfamily
SEED       <- 42

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Load and aggregate
# ---------------------------------------------------------------------------
dat <- fread(INFILE)
cat(sprintf("Input: %s\n  rows = %d\n", INFILE, nrow(dat)))

req <- c("species", "mark", "te_id", "superfam",
         "te_log2fc_meth", "te_log2fc_rna")
miss <- setdiff(req, names(dat))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

dat <- dat[!is.na(te_log2fc_meth) & !is.na(te_log2fc_rna)]
cat(sprintf("  rows with both methylation and TE-RNA = %d\n", nrow(dat)))

# one observation per TE copy: mean log2FC across contrasts
copies <- dat[, .(meth    = mean(te_log2fc_meth),
                  rna     = mean(te_log2fc_rna),
                  n_pairs = .N),
              by = .(species, mark, te_id, superfam)]
cat(sprintf("  independent TE copies = %d\n", nrow(copies)))

# display labels
GRP_LEVELS <- c("Phycomyces\n5mC", "Phycomyces\n6mA", "Rhizopus\n6mA")
copies[, group := factor(paste0(species, "\n", mark), levels = GRP_LEVELS)]

# ---------------------------------------------------------------------------
# 2. Statistics
# ---------------------------------------------------------------------------
fisher_ci <- function(r, n, conf = 0.95) {
  if (!is.finite(r) || n < 4) return(c(NA_real_, NA_real_))
  z  <- atanh(max(min(r, 0.999), -0.999))
  se <- 1 / sqrt(n - 3)
  zc <- qnorm(1 - (1 - conf) / 2)
  c(tanh(z - zc * se), tanh(z + zc * se))
}

stats <- copies[, {
  n <- .N
  if (n < MIN_COPIES) {
    .(n_copies = n, n_pairs = sum(n_pairs),
      r = NA_real_, p_two = NA_real_, p_one = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_)
  } else {
    ct <- suppressWarnings(cor.test(meth, rna, method = "pearson"))
    r  <- unname(ct$estimate)
    p2 <- ct$p.value
    p1 <- if (r > 0) p2 / 2 else 1 - p2 / 2
    ci <- fisher_ci(r, n)
    .(n_copies = n, n_pairs = sum(n_pairs),
      r = r, p_two = p2, p_one = p1,
      ci_lo = ci[1], ci_hi = ci[2])
  }
}, by = .(species, mark, superfam)]

stats[, group := factor(paste0(species, "\n", mark), levels = GRP_LEVELS)]
setcolorder(stats, c("group", "species", "mark", "superfam",
                     "n_copies", "n_pairs", "r", "p_two", "p_one",
                     "ci_lo", "ci_hi"))
stats <- stats[order(group, -r)]

stats_file <- file.path(OUTDIR, "48b_family_stats.tsv")
fwrite(stats, stats_file, sep = "\t", na = "NA")
cat(sprintf("\nStatistics table -> %s\n", stats_file))

# per-group summary
summ <- stats[!is.na(r), .(n_superfam  = .N,
                           n_positive  = sum(r > 0),
                           n_sig_one   = sum(p_one < 0.05),
                           n_sig_two   = sum(p_two < 0.05)),
              by = group]
print(summ)
summ_file <- file.path(OUTDIR, "48b_family_group_summary.tsv")
fwrite(summ, summ_file, sep = "\t")
cat(sprintf("Group summary     -> %s\n", summ_file))

# ---------------------------------------------------------------------------
# 3. Figure
# ---------------------------------------------------------------------------
pl <- stats[!is.na(r)]
pl[, superfam := factor(superfam, levels = unique(superfam[order(group, r)]))]
pl[, sig := fifelse(p_one < 0.05, "p < 0.05 (one-tailed)", "n.s.")]

# annotation: n superfamilies positive / significant, per group
ann <- pl[, .(label = sprintf("%d/%d r > 0 | %d significant (one-tailed)",
                              sum(r > 0), .N, sum(p_one < 0.05))),
          by = group]

forest <- ggplot(pl, aes(x = r, y = superfam)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, colour = sig),
                 height = 0, linewidth = 0.5) +
  geom_point(aes(size = n_copies, colour = sig), alpha = 0.9) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = -0.05, vjust = 1.4,
            size = 2.6, colour = "grey25") +
  facet_wrap(~ group, nrow = 1) +
  scale_colour_manual(values = c("p < 0.05 (one-tailed)" = "#c0392b",
                                 "n.s." = "#9aa0a6"), name = NULL) +
  scale_size_continuous(range = c(0.6, 4.5), name = "TE copies") +
  labs(x = "Pearson r (Δmethylation ↔ ΔTE expression), 95% CI",
       y = NULL,
       title = "TE methylation–expression correlation by superfamily") +
  theme_bw(base_size = 9) +
  theme(legend.position = "right",
        strip.background = element_rect(fill = "grey92"),
        panel.grid.minor = element_blank())

# representative scatter panels: prefer significant superfamilies, then largest |r|
rep_sf <- pl[order(group, p_one >= 0.05, -abs(r)), head(.SD, 2), by = group]
scat_data <- copies[rep_sf[, .(group, superfam)],
                    on = .(group, superfam), nomatch = 0]
lab <- rep_sf[, .(group, superfam,
                  label = sprintf("r = %+.2f\np = %s (one-tailed)\nn = %d copies",
                                  r,
                                  formatC(p_one, format = "g", digits = 2),
                                  n_copies))]

scatter <- ggplot(scat_data, aes(meth, rna)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_point(alpha = 0.45, size = 0.9, colour = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "firebrick", linewidth = 0.6) +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = -0.08, vjust = 1.25,
            size = 2.5, colour = "grey20") +
  facet_wrap(~ group + superfam, nrow = 2) +
  labs(x = "log₂FC (TE methylation)", y = "log₂FC (TE RNA)") +
  theme_bw(base_size = 9) +
  theme(strip.background = element_rect(fill = "grey92"),
        panel.grid.minor = element_blank())

fig <- forest / scatter +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(
    caption = paste0("One observation per TE copy (log2FC averaged across contrasts). ",
                     "Pearson correlation; 95% CI by Fisher z-transformation. ",
                     "Generated by 48b_family_summary_scatter.R"),
    tag_levels = "A")

pdf_file <- file.path(OUTDIR, "48b_family_summary_scatter.pdf")
png_file <- file.path(OUTDIR, "48b_family_summary_scatter.png")
ggsave(pdf_file, fig, width = 12, height = 8.5, device = cairo_pdf)
ggsave(png_file, fig, width = 12, height = 8.5, dpi = 300)
cat(sprintf("Figure            -> %s\n", pdf_file))
cat(sprintf("                  -> %s\n", png_file))

# ---------------------------------------------------------------------------
# 4. Reproducibility record
# ---------------------------------------------------------------------------
rec <- file.path(OUTDIR, "48b_reproducibility.txt")
con <- file(rec, "w")
writeLines(c(
  "48b_family_summary_scatter.R — run record",
  paste("Date      :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("Input     :", INFILE),
  paste("MD5 input :", unname(tools::md5sum(INFILE))),
  paste("Rows used :", nrow(dat)),
  paste("TE copies :", nrow(copies)),
  paste("Min copies per superfamily:", MIN_COPIES),
  paste("Seed      :", SEED),
  "", "--- sessionInfo() ---", ""
), con)
sink(con, append = TRUE); print(sessionInfo()); sink()
close(con)
cat(sprintf("Run record        -> %s\n", rec))
cat("\nDone.\n")
