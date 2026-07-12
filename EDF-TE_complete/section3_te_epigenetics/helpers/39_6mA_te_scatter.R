#!/usr/bin/env Rscript
# 39: 6mA overview scatter + top-5 and top-30 family comparison
#
# (A) Scatter overview: all TE copies, x=log1p(6mA density), y=log2(RNA+0.01)
#     colored by broad superfamily – one for Phycomyces, one for Rhizopus
# (B) Combined boxplot: top-5 fine families with largest 6mA expression effect
#     Side-by-side panels for Phycomyces and Rhizopus
# (C) Extended top-30 version for supplementary

suppressPackageStartupMessages({
  .libPaths(c("~/R/libs", .libPaths()))
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(Cairo)
})

PHY_TSV <- "/legserv/Temp/Carlos/EDF-TE/results/TE_dynamics/20_chiapello_extension/data/te_phybl2_rna_6mA.tsv"
RHI_TSV <- "/legserv/Temp/Carlos/EDF-TE/results/TE_dynamics/37_6mA_expression/te_rhimi59_rna_6mA.tsv"
OUTDIR  <- "/legserv/Temp/Carlos/EDF-TE/results/TE_dynamics/39_6mA_scatter"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ── Fixed classify_sf (handles full paths AND bare fine names) ─────────────
classify_sf <- function(x) {
  ifelse(grepl("Gypsy|Copia|DIRS|BelPao|^LTR|LTR/", x, ignore.case = TRUE), "LTR",
  ifelse(grepl(paste0("^LINE|LINE/|^Retroposon",
                      "|^L1|^L2|^CR1|^RTE|^I-|^R1$|^R2",
                      "|^Rex|Jockey|^Tad1|^SLACS|^TRAS|^Hero"),
               x, ignore.case = TRUE), "LINE",
  ifelse(grepl("^SINE|SINE/|^Alu/|^Alu$|^B1$|^B2$|^B4$|^ID$|^MIR", x, ignore.case = TRUE), "SINE",
  ifelse(grepl(paste0("^DNA|DNA/|TcMar|CMC|hAT|Harbinger|Maverick|Kolobok|",
                      "MuLE|Transib|Polinton|Crypton|Zisupton|ISRm|Merlin|Academ"),
               x, ignore.case = TRUE), "DNA",
  ifelse(grepl("^RC|RC/|Helitron", x, ignore.case = TRUE), "Rolling Circle",
                                                            "Unclassified")))))
}

CAT_COLS <- c(LTR = "#e76254", LINE = "#1b7837", DNA = "#4575b4",
              SINE = "#f4a74d", "Rolling Circle" = "#9970ab", Unclassified = "#888888")
SF_ORDER <- names(CAT_COLS)

# ── Load & classify ────────────────────────────────────────────────────────
phy <- fread(PHY_TSV)
rhi <- fread(RHI_TSV)
phy[, superfam := classify_sf(superfam_fine)]
rhi[, superfam := classify_sf(superfam_fine)]
phy[, species  := "Phycomyces"]
rhi[, species  := "Rhizopus"]
phy[, cond_label := condition]
rhi[, cond_label := condition]

cat("Phycomyces: broad superfam counts\n"); print(phy[, .N, by=superfam][order(-N)])
cat("Rhizopus:  broad superfam counts\n"); print(rhi[, .N, by=superfam][order(-N)])

# ── Dark-condition subsets ─────────────────────────────────────────────────
phy_dark_lbl <- grep("Dark|dark", unique(phy$cond_label), value = TRUE)[1]
rhi_dark_lbl <- grep("Dark|dark", unique(rhi$cond_label), value = TRUE)[1]

phy_dark <- phy[cond_label == phy_dark_lbl & RNA_signal >= 0 & !is.na(density_6mA)]
rhi_dark <- rhi[cond_label == rhi_dark_lbl & RNA_signal >= 0 & !is.na(density_6mA)]

# ── (A) Overview scatter ───────────────────────────────────────────────────
# Compute shared axis limits across both species
compute_scatter_ranges <- function(dt_list, max_pts = 6000, seed = 42) {
  set.seed(seed)
  combined <- rbindlist(lapply(dt_list, function(dt) {
    dt[, .SD[sample(.N, min(.N, max_pts))], by = superfam][,
       .(log2_rna = log2(RNA_signal + 0.01), dens_log1p = log1p(density_6mA))]
  }))
  list(x = range(combined$dens_log1p, na.rm = TRUE),
       y = range(combined$log2_rna,   na.rm = TRUE))
}

make_overview_scatter <- function(dt, sp_label, cond_lbl, xlim, ylim,
                                  max_pts = 6000, seed = 42) {
  set.seed(seed)
  dt2 <- dt[, .SD[sample(.N, min(.N, max_pts))], by = superfam]
  dt2[, log2_rna   := log2(RNA_signal + 0.01)]
  dt2[, dens_log1p := log1p(density_6mA)]
  dt2[, superfam   := factor(superfam, levels = SF_ORDER)]

  ggplot(dt2, aes(x = dens_log1p, y = log2_rna, colour = superfam)) +
    geom_point(alpha = 0.22, size = 0.7) +
    scale_colour_manual(values = CAT_COLS, name = "Superfamily", drop = FALSE) +
    scale_x_continuous(name = "6mA density  log1p (sites/kb)", limits = xlim) +
    scale_y_continuous(name = expression(log[2](RNA + 0.01)), limits = ylim) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right", panel.grid.minor = element_blank()) +
    labs(title = paste0(sp_label, "  –  6mA density vs RNA expression per TE copy"),
         subtitle = paste0(cond_lbl, "  |  max ", max_pts, " copies/superfamily shown"))
}

ax <- compute_scatter_ranges(list(phy_dark, rhi_dark))
p_phy_ov <- make_overview_scatter(phy_dark, "Phycomyces", phy_dark_lbl, ax$x, ax$y)
p_rhi_ov <- make_overview_scatter(rhi_dark, "Rhizopus",   rhi_dark_lbl, ax$x, ax$y)

CairoPDF(file.path(OUTDIR, "39a_phyco_6mA_scatter.pdf"),     width = 8, height = 5)
print(p_phy_ov); dev.off()
CairoPDF(file.path(OUTDIR, "39a_rhizopus_6mA_scatter.pdf"),  width = 8, height = 5)
print(p_rhi_ov); dev.off()

combined_ov <- (p_phy_ov + p_rhi_ov) + plot_layout(guides = "collect") &
  theme(legend.position = "right")
CairoPDF(file.path(OUTDIR, "39a_6mA_scatter_combined.pdf"), width = 15, height = 5)
print(combined_ov); dev.off()
cat("(A) overview scatter done\n")

# ── (B) Top-N families by |log2FC|  ────────────────────────────────────────
# Use ALL conditions (one row per superfam_fine × cond), find max |log2FC|
MIN_TOTAL <- 10L   # minimum total copies per family (with + without 6mA)

compute_fc <- function(dt, sp_lbl) {
  dt2 <- dt[RNA_signal >= 0 & !is.na(density_6mA)]
  fc <- dt2[, .(
    n_with    = sum(has_6mA),
    n_without = sum(!has_6mA),
    med_with  = median(RNA_signal[has_6mA],  na.rm = TRUE),
    med_without = median(RNA_signal[!has_6mA], na.rm = TRUE)
  ), by = .(superfam_fine, cond_label)]
  fc[, n_total := n_with + n_without]
  fc <- fc[n_total >= MIN_TOTAL & n_with >= 1L & n_without >= 1L]
  fc[, log2fc := log2((med_with + 0.01) / (med_without + 0.01))]
  fc[, max_abs_fc := max(abs(log2fc), na.rm = TRUE), by = superfam_fine]
  fc[, species := sp_lbl]
  fc
}

phy_fc <- compute_fc(phy, "Phycomyces")
rhi_fc <- compute_fc(rhi, "Rhizopus")

top_n_families <- function(fc_dt, n) {
  ranked <- unique(fc_dt[, .(superfam_fine, max_abs_fc)])[order(-max_abs_fc)]
  ranked[1:min(n, .N), superfam_fine]
}

top5_phy  <- top_n_families(phy_fc, 5)
top5_rhi  <- top_n_families(rhi_fc, 5)
top30_phy <- top_n_families(phy_fc, 30)
top30_rhi <- top_n_families(rhi_fc, 30)

cat("Phycomyces top 5:\n"); print(top5_phy)
cat("Rhizopus top 5:\n");   print(top5_rhi)

# ── Boxplot function for top-N families ────────────────────────────────────
make_topn_boxplot <- function(dt, top_fams, sp_lbl, dark_lbl) {
  dt2 <- dt[superfam_fine %in% top_fams & cond_label == dark_lbl & RNA_signal >= 0]
  dt2[, superfam := classify_sf(superfam_fine)]
  dt2[, has_lbl  := ifelse(has_6mA, "6mA+", "6mA\u2212")]
  dt2[, has_lbl  := factor(has_lbl, levels = c("6mA\u2212", "6mA+"))]
  # order families by |FC|
  order_fams <- phy_fc[species == sp_lbl & cond_label == dark_lbl & superfam_fine %in% top_fams,
                       .(abs_fc = max(abs(log2fc))), by = superfam_fine][order(-abs_fc), superfam_fine]
  order_fams <- intersect(order_fams, unique(dt2$superfam_fine))
  if (length(order_fams) == 0) order_fams <- top_fams
  dt2[, superfam_fine := factor(superfam_fine, levels = order_fams)]
  cat_cols  <- setNames(CAT_COLS[classify_sf(levels(dt2$superfam_fine))],
                         levels(dt2$superfam_fine))
  ycap <- quantile(dt2$RNA_signal, 0.98, na.rm = TRUE)
  dt2[, rna_cap := pmin(RNA_signal, ycap)]

  ggplot(dt2, aes(x = superfam_fine, y = rna_cap, fill = superfam_fine)) +
    geom_boxplot(aes(linetype = has_lbl),
                 position = position_dodge(0.7), width = 0.55,
                 outlier.size = 0.3, outlier.alpha = 0.25, linewidth = 0.4) +
    scale_fill_manual(values = cat_cols, guide = "none") +
    scale_linetype_manual(values = c("6mA\u2212" = "dashed", "6mA+" = "solid"), name = NULL) +
    scale_x_discrete(name = NULL) +
    scale_y_continuous(name = "RNA signal (capped 98th pctile)") +
    coord_flip() +
    theme_bw(base_size = 10) +
    theme(legend.position = "top", panel.grid.major.y = element_blank()) +
    labs(title = sp_lbl,
         subtitle = paste0("Solid = 6mA+, dashed = 6mA\u2212  |  ", dark_lbl))
}

make_topn_boxplot_v2 <- function(dt, fc_dt, top_fams, sp_lbl, dark_lbl, ncols = 3) {
  # One facet per family; both 6mA+ and 6mA- as fill
  dt2 <- dt[superfam_fine %in% top_fams & cond_label == dark_lbl & RNA_signal >= 0]
  if (nrow(dt2) == 0) return(ggplot() + labs(title = paste(sp_lbl, "– no data")))
  dt2[, has_lbl := ifelse(has_6mA, "6mA+", "6mA\u2212")]
  dt2[, has_lbl := factor(has_lbl, levels = c("6mA\u2212", "6mA+"))]
  # order by FC
  ranked <- fc_dt[cond_label == dark_lbl & superfam_fine %in% top_fams,
                  .(abs_fc = max(abs(log2fc), na.rm = TRUE)), by = superfam_fine][order(-abs_fc)]
  fam_lvls <- intersect(ranked$superfam_fine, unique(dt2$superfam_fine))
  dt2[, superfam_fine := factor(superfam_fine, levels = fam_lvls)]
  ycap <- quantile(dt2$RNA_signal, 0.99, na.rm = TRUE)
  dt2[, rna_cap := pmin(RNA_signal, ycap)]

  ggplot(dt2, aes(x = has_lbl, y = rna_cap, fill = has_lbl)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2,
                 linewidth = 0.4, width = 0.55) +
    scale_fill_manual(values = c("6mA\u2212" = "#9ecae1", "6mA+" = "#fc8d59"), name = NULL) +
    facet_wrap(~ superfam_fine, scales = "free_y", ncol = ncols) +
    scale_y_continuous(name = "RNA signal") +
    scale_x_discrete(name = NULL) +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(size = 8, face = "bold"),
          legend.position = "bottom", panel.grid.minor = element_blank()) +
    labs(title = sp_lbl, subtitle = dark_lbl)
}

p5_phy <- make_topn_boxplot_v2(phy, phy_fc, top5_phy,  "Phycomyces", phy_dark_lbl, ncols = 5)
p5_rhi <- make_topn_boxplot_v2(rhi, rhi_fc, top5_rhi,  "Rhizopus",   rhi_dark_lbl, ncols = 5)

combined5 <- p5_phy / p5_rhi + plot_layout(heights = c(1, 1))
CairoPDF(file.path(OUTDIR, "39b_top5_family_boxplot.pdf"), width = 14, height = 8)
print(combined5); dev.off()
cat("(B) top-5 boxplot done\n")

# ── (C) Top-30 supplementary ───────────────────────────────────────────────
p30_phy <- make_topn_boxplot_v2(phy, phy_fc, top30_phy, "Phycomyces", phy_dark_lbl, ncols = 6)
p30_rhi <- make_topn_boxplot_v2(rhi, rhi_fc, top30_rhi, "Rhizopus",   rhi_dark_lbl, ncols = 6)

CairoPDF(file.path(OUTDIR, "39c_top30_phyco_supp.pdf"),
         width = 18, height = max(4, ceiling(length(top30_phy)/6) * 3))
print(p30_phy); dev.off()
CairoPDF(file.path(OUTDIR, "39c_top30_rhizopus_supp.pdf"),
         width = 18, height = max(4, ceiling(length(top30_rhi)/6) * 3))
print(p30_rhi); dev.off()
cat("(C) top-30 supplementary done\n")

# ── Save FC tables ─────────────────────────────────────────────────────────
fwrite(phy_fc[order(-max_abs_fc)], file.path(OUTDIR, "39_phyco_6mA_fc_by_family.tsv"), sep="\t")
fwrite(rhi_fc[order(-max_abs_fc)], file.path(OUTDIR, "39_rhizopus_6mA_fc_by_family.tsv"), sep="\t")

cat("\nDone. Outputs in:", OUTDIR, "\n")
