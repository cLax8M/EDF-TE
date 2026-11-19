#!/usr/bin/env Rscript
# 37b_6mA_expression_combined.R
# Combined visualization: RNA expression in 6mA+ vs 6mA- TE copies
# for Phycomyces blakesleeanus (Phybl2) and Rhizopus microsporus (Rhimi59)
# Requires: te_phybl2_rna_6mA.tsv (existing) + te_rhimi59_rna_6mA.tsv (from script 37)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Cairo)
})

# ---- Paths ---------------------------------------------------------------
PHY_TSV  <- "/legserv/Temp/Carlos/EDF-TE/results/TE_dynamics/20_chiapello_extension/data/te_phybl2_rna_6mA.tsv"
RHI_TSV  <- "/legserv/Temp/Carlos/EDF-TE/results/TE_dynamics/37_6mA_expression/te_rhimi59_rna_6mA.tsv"
PHY_STAT <- "/legserv/Temp/Carlos/EDF-TE/results/TE_dynamics/20_chiapello_extension/data/F5_superfamily_stats.tsv"
RHI_STAT <- "/legserv/Temp/Carlos/EDF-TE/results/TE_dynamics/37_6mA_expression/rhimi59_superfam_stats.tsv"
OUTDIR   <- "/legserv/Temp/Carlos/EDF-TE/results/TE_dynamics/37_6mA_expression"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---- Broad TE category function (6 categories) --------------------------
# Handles both full RepBase paths (LINE/L1) AND fine names (L1, CR1, Academ)
classify_sf <- function(x) {
  ifelse(grepl("Gypsy|Copia|DIRS|BelPao|^LTR|LTR/", x, ignore.case = TRUE), "LTR",
  ifelse(grepl(paste0("^LINE|LINE/|^Retroposon",
                      "|^L1|^L2|^CR1|^RTE|^I-|^R1$|^R2",
                      "|^Rex|Jockey|^Tad1|^SLACS|^TRAS|^Hero"),
               x, ignore.case = TRUE), "LINE",
  ifelse(grepl("^SINE|SINE/|^Alu/|^Alu$|^B1$|^B2$|^B4$|^ID$|^MIR",
               x, ignore.case = TRUE), "SINE",
  ifelse(grepl(paste0("^DNA|DNA/|TcMar|CMC|hAT|Harbinger|Maverick|Kolobok|",
                      "MuLE|Transib|Polinton|Crypton|Zisupton|ISRm|Merlin|Academ"),
               x, ignore.case = TRUE), "DNA",
  ifelse(grepl("^RC|RC/|Helitron", x, ignore.case = TRUE), "Rolling Circle",
                                                            "Unclassified")))))
}

# ---- Load per-copy tables ------------------------------------------------
phy <- fread(PHY_TSV)
rhi <- fread(RHI_TSV)

# Reclassify using fine names (handles full paths like LINE/L1 and bare names like L1, CR1)
phy[, superfam := classify_sf(superfam_fine)]
rhi[, superfam := classify_sf(superfam_fine)]   # re-classify; fixes LINE/CR1/etc.
rhi[, species := "Rhizopus"]
phy[, species := "Phycomyces"]

# Use a common condition label
phy[, cond_label := sub("Micelium",     "Mycelium", condition)]
phy[, cond_label := sub("Sporangiophore", "Sporangiophore", cond_label)]
rhi[, cond_label := condition]

# ---- Pick superfamilies to display --------------------------------------
# Superfamilies with enough 6mA+ copies in at least one species/condition
MIN_WITH <- 20L

sf_keep <- function(dt) {
  k <- dt[, .(n_with = sum(has_6mA), n_tot = .N), by = .(superfam, cond_label)]
  unique(k[n_with >= MIN_WITH, superfam])
}
phy_sf <- sf_keep(phy)
rhi_sf <- sf_keep(rhi)

phy_plt <- phy[superfam %in% phy_sf & RNA_signal >= 0]
rhi_plt <- rhi[superfam %in% rhi_sf & RNA_signal >= 0]

SF_ORDER <- c("DNA", "LTR", "LINE", "SINE", "Rolling Circle", "Unclassified")

set_sf_factor <- function(dt) {
  sf_present <- intersect(SF_ORDER, unique(dt$superfam))
  extra       <- setdiff(unique(dt$superfam), SF_ORDER)
  dt[, superfam := factor(superfam, levels = c(sf_present, extra))]
}
phy_plt <- set_sf_factor(phy_plt)
rhi_plt <- set_sf_factor(rhi_plt)

# ---- Compute significance labels ----------------------------------------
sig_label <- function(p) {
  cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("***", "**", "*", "ns"), right = FALSE)
}

wilcox_summary <- function(dt) {
  dt[, .(
    n_with   = sum(has_6mA),
    n_without = sum(!has_6mA),
    med_with  = median(RNA_signal[has_6mA == TRUE],  na.rm = TRUE),
    med_without = median(RNA_signal[has_6mA == FALSE], na.rm = TRUE),
    p_val    = tryCatch(
      wilcox.test(RNA_signal[has_6mA == TRUE],
                  RNA_signal[has_6mA == FALSE],
                  exact = FALSE)$p.value,
      error = function(e) NA_real_)
  ), by = .(superfam, cond_label)]
}

phy_stat <- wilcox_summary(phy_plt)
rhi_stat <- wilcox_summary(rhi_plt)
phy_stat[, p_adj := p.adjust(p_val, "BH")]
rhi_stat[, p_adj := p.adjust(p_val, "BH")]
phy_stat[, sig   := sig_label(p_adj)]
rhi_stat[, sig   := sig_label(p_adj)]

# Save stats
fwrite(phy_stat, file.path(OUTDIR, "phyco_6mA_expression_stats.tsv"),  sep = "\t")
fwrite(rhi_stat, file.path(OUTDIR, "rhimi_6mA_expression_stats.tsv"),   sep = "\t")

cat("=== Phycomyces 6mA × RNA (sig results) ===\n")
print(phy_stat[p_adj < 0.05][order(superfam, cond_label)])
cat("=== Rhizopus 6mA × RNA (sig results) ===\n")
print(rhi_stat[p_adj < 0.05][order(superfam, cond_label)])

# ---- Plotting function --------------------------------------------------
make_violin <- function(dt, sig_dt, sp_label, y_cap = 10) {
  dt2 <- copy(dt)
  dt2[, RNA_capped := pmin(RNA_signal, y_cap)]
  dt2[, has_6mA_lbl := ifelse(has_6mA, "6mA+", "6mA\u2212")]
  dt2[, has_6mA_lbl := factor(has_6mA_lbl, levels = c("6mA\u2212","6mA+"))]

  # Y positions for significance labels (above each paired violin)
  top_y <- y_cap * 1.05

  ggplot(dt2, aes(x = has_6mA_lbl, y = RNA_capped, fill = has_6mA_lbl)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.75, linewidth = 0.3) +
    geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white",
                 color = "gray30", linewidth = 0.4) +
    geom_text(data = sig_dt[p_adj < 0.05],
              aes(x = 1.5, y = top_y, label = sig),
              inherit.aes = FALSE, size = 3.5, fontface = "bold") +
    facet_grid(cond_label ~ superfam, scales = "free_y") +
    scale_fill_manual(values = c("6mA\u2212" = "#AECDE8", "6mA+" = "#F4A261"),
                      name = NULL) +
    scale_y_continuous(name = "RNA signal (mean coverage)", limits = c(0, y_cap * 1.15)) +
    scale_x_discrete(name = NULL) +
    theme_bw(base_size = 10) +
    theme(
      strip.text.x    = element_text(size = 8, face = "bold"),
      strip.text.y    = element_text(size = 8),
      axis.text.x     = element_text(size = 8),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    ) +
    labs(title = paste0(sp_label, ": RNA expression in 6mA-marked vs unmarked TE copies"),
         subtitle = "Wilcoxon rank-sum test; BH-adjusted; boxes = median + IQR; y-axis capped")
}

# ---- Figure 1: Phycomyces violin (per superfamily) ----------------------
# Determine y cap from data
phy_cap <- quantile(phy_plt$RNA_signal, 0.99, na.rm = TRUE)
rhi_cap <- quantile(rhi_plt$RNA_signal, 0.99, na.rm = TRUE)

p_phy <- make_violin(phy_plt, phy_stat, "Phycomyces", y_cap = phy_cap)
p_rhi <- make_violin(rhi_plt, rhi_stat, "Rhizopus",   y_cap = rhi_cap)

phy_nsf <- uniqueN(phy_plt$superfam)
rhi_nsf <- uniqueN(rhi_plt$superfam)

CairoPDF(file.path(OUTDIR, "37b_phyco_6mA_expression.pdf"),
         width = max(8, phy_nsf * 2), height = 8)
print(p_phy); dev.off()
CairoPDF(file.path(OUTDIR, "37b_rhizopus_6mA_expression.pdf"),
         width = max(8, rhi_nsf * 2), height = 6)
print(p_rhi); dev.off()

# ---- Figure 2: Summary strip chart (fold change median) -----------------
# Median RNA in 6mA+ / median RNA in 6mA- per species × superfamily (dark condition only)
compare_conds <- c("Dark","MiceliumDark","Mycelium Dark")
phy_dark <- phy_stat[grepl("Dark|dark", cond_label, ignore.case = TRUE) & n_with >= MIN_WITH]
rhi_dark <- rhi_stat[grepl("Dark|dark", cond_label, ignore.case = TRUE) & n_with >= MIN_WITH]
phy_dark[, species := "Phycomyces"]
rhi_dark[, species := "Rhizopus"]
both_dark <- rbindlist(list(phy_dark, rhi_dark), fill = TRUE)
both_dark[, log2fc := log2((med_with + 0.001) / (med_without + 0.001))]

p_summary <- ggplot(both_dark[!is.na(log2fc)],
                    aes(x = superfam, y = log2fc, color = species, shape = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = c(Phycomyces = "#E15759", Rhizopus = "#59A14F")) +
  scale_shape_manual(values = c("***" = 17, "**" = 17, "*" = 16, "ns" = 1),
                     name = "Significance") +
  scale_y_continuous(name = expression(log[2](RNA["6mA+"]/RNA["6mA-"]))) +
  scale_x_discrete(name = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  labs(title = "6mA-marked TE copies: relative RNA expression vs unmarked (Dark condition)",
       subtitle = "Median RNA_signal ratio; filled = p.adj<0.05, open = ns")

CairoPDF(file.path(OUTDIR, "37b_summary_log2fc.pdf"), width = 9, height = 5)
print(p_summary); dev.off()
CairoPNG(file.path(OUTDIR, "37b_summary_log2fc.png"), width = 1200, height = 600, dpi = 130)
print(p_summary); dev.off()

# ---- Figure 3: Scatter – 6mA density vs RNA per TE copy -----------------
# Shows the continuous relationship between 6mA density (sites/kb) and RNA
# signal for individual TE copies, one panel per broad superfamily.
# Only one condition is shown (Dark / MiceliumDark) to keep the figure clean.

CAT_COLS <- c(
  "LTR"            = "#e76254",
  "LINE"           = "#1b7837",
  "DNA"            = "#4575b4",
  "SINE"           = "#f4a74d",
  "Rolling Circle" = "#9970ab",
  "Unclassified"   = "#888888"
)

make_scatter <- function(dt, sp_label, cond_sel, max_pts_per_sf = 4000) {
  dt2 <- dt[cond_label == cond_sel & RNA_signal >= 0 & !is.na(density_6mA)]
  if (nrow(dt2) == 0) {
    # fall back to first condition available
    cond_sel <- dt$cond_label[1]
    dt2 <- dt[cond_label == cond_sel & RNA_signal >= 0 & !is.na(density_6mA)]
  }
  dt2[, log2_rna   := log2(RNA_signal + 0.01)]
  dt2[, dens_log1p := log1p(density_6mA)]
  # keep only superfamilies with ≥ 20 copies and subsample heavy families
  dt2 <- dt2[superfam %in% dt2[, .N, by = superfam][N >= 20, superfam]]
  dt2 <- dt2[, .SD[sample(.N, min(.N, max_pts_per_sf))], by = superfam]

  n_sf <- uniqueN(dt2$superfam)
  cols <- CAT_COLS[levels(factor(dt2$superfam, levels = names(CAT_COLS)))]
  cols <- cols[!is.na(cols)]

  ggplot(dt2, aes(x = dens_log1p, y = log2_rna, colour = superfam)) +
    geom_point(alpha = 0.18, size = 0.7, show.legend = FALSE) +
    geom_smooth(method  = "loess", span = 0.8, se = TRUE,
                linewidth = 0.9, alpha = 0.25, show.legend = FALSE) +
    facet_wrap(~ superfam, scales = "free_y", ncol = min(n_sf, 3)) +
    scale_colour_manual(values = CAT_COLS) +
    scale_x_continuous(name = "6mA density  log1p (sites/kb)") +
    scale_y_continuous(name = expression(log[2](RNA + 0.01))) +
    theme_bw(base_size = 10) +
    theme(strip.text       = element_text(size = 9, face = "bold"),
          panel.grid.minor = element_blank()) +
    labs(title    = paste0(sp_label, ": RNA vs 6mA density per TE copy  (", cond_sel, ")"),
         subtitle = paste0("Each point = one TE locus; loess trend \u00b1 SE; ",
                           "max ", max_pts_per_sf, " loci/superfamily shown"))
}

# Detect correct condition label for dark mycelium
phy_cond_dark <- grep("Dark|dark", unique(phy$cond_label), value = TRUE)[1]
rhi_cond_dark <- grep("Dark|dark", unique(rhi$cond_label), value = TRUE)[1]

p_scat_phy <- make_scatter(phy_plt, "Phycomyces", phy_cond_dark)
p_scat_rhi <- make_scatter(rhi_plt, "Rhizopus",   rhi_cond_dark)

phy_nsf2 <- uniqueN(phy_plt[cond_label == phy_cond_dark & !is.na(density_6mA), superfam])
rhi_nsf2 <- uniqueN(rhi_plt[cond_label == rhi_cond_dark & !is.na(density_6mA), superfam])

CairoPDF(file.path(OUTDIR, "37b_phyco_scatter_6mA_vs_RNA.pdf"),
         width = max(6, ceiling(phy_nsf2 / 3) * 5), height = max(4, ceiling(phy_nsf2 / 3) * 4))
print(p_scat_phy); dev.off()

CairoPDF(file.path(OUTDIR, "37b_rhizopus_scatter_6mA_vs_RNA.pdf"),
         width = max(6, ceiling(rhi_nsf2 / 3) * 5), height = max(4, ceiling(rhi_nsf2 / 3) * 4))
print(p_scat_rhi); dev.off()

message("All figures saved to: ", OUTDIR)
