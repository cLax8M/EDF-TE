#!/usr/bin/env Rscript
# 45_mac_tss_expression.R
# Boxplot: TE expression stratified by presence/absence of 6mA MAC at the TSS
# Uses script 42 output (42_te_with_mac_at_tss.tsv) merged with RNA data
# Both species, all conditions, Wilcoxon test annotation

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
WDIR   <- "/legserv/Temp/Carlos/EDF-TE/results/TE_dynamics"
OUTDIR <- file.path(WDIR, "45_mac_tss_expression")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

MAC_TSV <- file.path(WDIR, "42_6mA_mac_at_tss/42_te_with_mac_at_tss.tsv")
PHY_RNA <- file.path(WDIR, "20_chiapello_extension/data/te_phybl2_rna_6mA.tsv")
RHI_RNA <- file.path(WDIR, "37_6mA_expression/te_rhimi59_rna_6mA.tsv")

# ── Classify superfamily ──────────────────────────────────────────────────────
classify_sf <- function(x) {
  ifelse(grepl("Gypsy|Copia|DIRS|BelPao|^LTR|LTR/", x, ignore.case = TRUE), "LTR",
  ifelse(grepl(paste0("^LINE|LINE/|^Retroposon|^L1|^L2|^CR1|^RTE|^I-|^R1$|^R2",
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

CAT_COLS <- c(LTR = "#e76254", LINE = "#1b7837", DNA = "#4575b4",
              SINE = "#f4a74d", "Rolling Circle" = "#9970ab",
              Unclassified = "#888888")

# ── Load MAC at TSS list ───────────────────────────────────────────────────────
cat("Loading MAC-at-TSS list (script 42 output)...\n")
mac <- fread(MAC_TSV)
# cols: species, chrom, start, end, name, strand, superfam
# Create position key
mac[, te_id_pos := paste(chrom, start, end, sep = "_")]
# Set of positions per species with MAC
phy_mac_pos <- mac[species == "Phycomyces", te_id_pos]
rhi_mac_pos <- mac[species == "Rhizopus",   te_id_pos]
cat(sprintf("  Phycomyces: %d TEs with MAC at TSS\n", length(phy_mac_pos)))
cat(sprintf("  Rhizopus:   %d TEs with MAC at TSS\n", length(rhi_mac_pos)))

# ── Helper to extract position key from te_id ──────────────────────────────────
# te_id: scaffold_N_start_end_extra  →  scaffold_N_start_end
parse_te_pos <- function(te_id) {
  sub("^(scaffold_[0-9]+_[0-9]+_[0-9]+)_.*$", "\\1", te_id)
}

# ── Load RNA data and annotate has_mac ─────────────────────────────────────────
cat("Loading RNA data...\n")
phy_rna <- fread(PHY_RNA)
phy_rna[, te_id_pos  := parse_te_pos(te_id)]
phy_rna[, has_mac    := te_id_pos %in% phy_mac_pos]
phy_rna[, species    := "Phycomyces blakesleeanus"]
phy_rna[, superfam_broad := classify_sf(superfam_fine)]
# Rename condition for display
phy_rna[, cond_label := sub("Micelium", "Mycelium ", condition)]

rhi_rna <- fread(RHI_RNA)
if (ncol(rhi_rna) >= 8) {
  setnames(rhi_rna, names(rhi_rna),
           c("te_id","condition","RNA_signal","n6mA","superfam_fine",
             "superfam_old","len","density_6mA","has_6mA")[seq_len(ncol(rhi_rna))])
}
rhi_rna[, te_id_pos  := parse_te_pos(te_id)]
rhi_rna[, has_mac    := te_id_pos %in% rhi_mac_pos]
rhi_rna[, species    := "Rhizopus microsporus"]
rhi_rna[, superfam_broad := classify_sf(superfam_fine)]
rhi_rna[, cond_label := condition]

# Keep useful cols
keep_cols <- c("te_id","te_id_pos","condition","cond_label","RNA_signal",
               "superfam_fine","superfam_broad","has_mac","species")
dat <- rbind(
  phy_rna[, ..keep_cols],
  rhi_rna[, intersect(names(rhi_rna), keep_cols), with = FALSE],
  fill = TRUE
)
dat[, log2_rna := log2(RNA_signal + 1)]
dat[, mac_label := fifelse(has_mac, "MAC at TSS", "No MAC")]
dat[, mac_label := factor(mac_label, levels = c("No MAC", "MAC at TSS"))]

cat(sprintf("  Total rows: %d\n", nrow(dat)))
cat("\nFraction with MAC at TSS per species:\n")
print(dat[, .(n_total = .N, n_mac = sum(has_mac),
              pct_mac = round(100*mean(has_mac),1)),
          by = .(species, condition)])

# ── Wilcoxon tests per species × condition ────────────────────────────────────
wtest <- dat[, {
  a <- log2_rna[has_mac  == TRUE]
  b <- log2_rna[has_mac  == FALSE]
  if (length(a) < 5 || length(b) < 5) {
    .(W_stat = NA_real_, p_val = NA_real_, n_mac = length(a), n_no_mac = length(b))
  } else {
    wt <- wilcox.test(a, b, exact = FALSE)
    .(W_stat = wt$statistic, p_val = wt$p.value,
      n_mac = length(a), n_no_mac = length(b))
  }
}, by = .(species, condition)]
wtest[, p_adj := p.adjust(p_val, method = "BH")]
wtest[, sig_label := fcase(
  p_adj < 0.001, "***",
  p_adj < 0.01,  "**",
  p_adj < 0.05,  "*",
  default = "ns"
)]
cat("\n=== Wilcoxon test results ===\n")
print(wtest)
fwrite(wtest, file.path(OUTDIR, "45_wilcoxon_results.tsv"), sep = "\t")

# ── Plot: overall (all conditions, both species) ───────────────────────────────
# Per-species summary for annotation
y_max_dat <- dat[, max(log2_rna, na.rm = TRUE), by = .(species, condition)]
wtest_annotate <- wtest[y_max_dat, on = .(species, condition)]
wtest_annotate[, y_annot := V1 * 1.05]
# Add cond_label for facet matching
wtest_annotate[, cond_label := fifelse(
  species == "Phycomyces blakesleeanus",
  sub("Micelium", "Mycelium ", condition),
  condition
)]

# n per group for labels
n_labels <- dat[, .(n = .N), by = .(species, cond_label, mac_label)]

# sample up to 300 points per group for jitter overlay
set.seed(42)
pts <- dat[, .SD[sample(.N, min(.N, 300))], by = .(species, cond_label, mac_label)]

p_main <- ggplot(dat, aes(x = mac_label, y = log2_rna, fill = mac_label)) +
  geom_violin(alpha = 0.55, linewidth = 0.3, trim = TRUE) +
  geom_point(data = pts, size = 0.15, alpha = 0.12, colour = "grey25",
             position = position_jitter(width = 0.18, height = 0, seed = 42)) +
  geom_boxplot(width = 0.25, outlier.shape = NA, linewidth = 0.4,
               fill = "white", alpha = 0.9) +
  geom_text(data = wtest_annotate,
            aes(x = 1.5, y = Inf, label = sig_label),
            inherit.aes = FALSE, size = 4, fontface = "bold", vjust = 1.5) +
  geom_text(data = n_labels,
            aes(x = mac_label, y = -Inf,
                label = paste0("n=", format(n, big.mark = ","))),
            inherit.aes = FALSE, size = 2.8, vjust = -0.4, colour = "grey30") +
  scale_fill_manual(values = c("No MAC" = "#4575b4", "MAC at TSS" = "#d73027"),
                    guide = "none") +
  facet_grid(species ~ cond_label, scales = "free_y") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = expression(log[2](RNA + 1)),
       title = "TE expression: MAC at TSS vs no MAC",
       subtitle = "Wilcoxon test (BH-adjusted): * p<0.05  ** p<0.01  *** p<0.001") +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(),
        strip.text.y = element_text(face = "italic"),
        strip.text.x = element_text(size = 9),
        plot.margin = margin(5, 5, 20, 5))

out_pdf <- file.path(OUTDIR, "45_mac_tss_expression.pdf")
out_png <- file.path(OUTDIR, "45_mac_tss_expression.png")
ggsave(out_pdf, p_main, width = 10, height = 7)
ggsave(out_png, p_main, width = 10, height = 7, dpi = 150)
cat(sprintf("Saved: %s\n", out_pdf))

# ── Bonus: by broad superfamily (Dark condition only) ─────────────────────────
dat_dark <- dat[condition %in% c("MiceliumDark", "Dark")]
# Remove Unclassified and small categories
dat_dark <- dat_dark[superfam_broad != "Unclassified"]

sf_order <- dat_dark[, .(n = .N), by = superfam_broad][order(-n), superfam_broad]
dat_dark[, sf_f := factor(superfam_broad, levels = sf_order)]

p_sf <- ggplot(dat_dark, aes(x = mac_label, y = log2_rna, fill = mac_label)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.3, alpha = 0.8) +
  scale_fill_manual(values = c("No MAC" = "#4575b4", "MAC at TSS" = "#d73027"),
                    guide = "none") +
  facet_grid(species ~ sf_f, scales = "free_y") +
  labs(x = NULL, y = expression(log[2](RNA + 1)),
       title = "TE expression by MAC at TSS status — by superfamily (Dark)") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
        panel.grid.major.x = element_blank(),
        strip.text.y = element_text(face = "italic"),
        strip.text.x = element_text(size = 8))

out_sf_pdf <- file.path(OUTDIR, "45_mac_tss_expression_by_superfam.pdf")
out_sf_png <- file.path(OUTDIR, "45_mac_tss_expression_by_superfam.png")
ggsave(out_sf_pdf, p_sf, width = 14, height = 6)
ggsave(out_sf_png, p_sf, width = 14, height = 6, dpi = 150)
cat(sprintf("Saved: %s\n", out_sf_pdf))

# ── Median RNA by mac group per species ───────────────────────────────────────
cat("\n=== Median log2(RNA+1) by MAC group ===\n")
print(dat[, .(med = median(log2_rna, na.rm=TRUE), n = .N),
          by = .(species, condition, mac_label)])
cat("Done.\n")
