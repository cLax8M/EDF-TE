#!/usr/bin/env Rscript
# 44b_mac_tss_families.R
# Repeat 44alt structure but split by MAC-at-TSS instead of 6mA,
# and select the 10 most frequent (most copies) families per species.
# Dark / MiceliumDark condition only.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

WDIR    <- "/legserv/Temp/Carlos/EDF-TE/results/TE_dynamics"
OUTDIR  <- file.path(WDIR, "44b_mac_tss_families")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

MAC_TSV <- file.path(WDIR, "42_6mA_mac_at_tss/42_te_with_mac_at_tss.tsv")
PHY_RNA <- file.path(WDIR, "20_chiapello_extension/data/te_phybl2_rna_6mA.tsv")
RHI_RNA <- file.path(WDIR, "37_6mA_expression/te_rhimi59_rna_6mA.tsv")

TOP_N_FAM <- 10
MIN_N     <- 5

# ── Helper: extract scaffold_N_start_end key from te_id ───────────────────────
parse_te_pos <- function(te_id) {
  sub("^(scaffold_[0-9]+_[0-9]+_[0-9]+)_.*$", "\\1", te_id)
}

# ── Load MAC at TSS ─────────────────────────────────────────────────────────
cat("Loading MAC-at-TSS data...\n")
mac <- fread(MAC_TSV)
mac[, te_id_pos := paste(chrom, start, end, sep = "_")]
phy_mac_pos <- mac[species == "Phycomyces", te_id_pos]
rhi_mac_pos <- mac[species == "Rhizopus",   te_id_pos]
cat(sprintf("  Phycomyces: %d TEs with MAC at TSS\n", length(phy_mac_pos)))
cat(sprintf("  Rhizopus:   %d TEs with MAC at TSS\n", length(rhi_mac_pos)))

# ── Load RNA data ────────────────────────────────────────────────────────────
cat("Loading RNA data...\n")
phy_rna <- fread(PHY_RNA)
phy_rna[, te_id_pos := parse_te_pos(te_id)]
phy_rna[, has_mac   := te_id_pos %in% phy_mac_pos]
phy_rna[, species   := "Phycomyces blakesleeanus"]

rhi_rna <- fread(RHI_RNA)
if (ncol(rhi_rna) >= 8) {
  setnames(rhi_rna, names(rhi_rna),
           c("te_id","condition","RNA_signal","n6mA","superfam_fine",
             "superfam_old","len","density_6mA","has_6mA")[seq_len(ncol(rhi_rna))])
}
rhi_rna[, te_id_pos := parse_te_pos(te_id)]
rhi_rna[, has_mac   := te_id_pos %in% rhi_mac_pos]
rhi_rna[, species   := "Rhizopus microsporus"]

# Keep Dark condition only
keep_cols <- c("te_id","te_id_pos","condition","RNA_signal",
               "superfam_fine","has_mac","species")
dat <- rbind(
  phy_rna[condition == "MiceliumDark",  ..keep_cols],
  rhi_rna[condition == "Dark",          ..keep_cols],
  fill = TRUE
)

# Remove structural/non-TE annotations and Unknown
dat <- dat[!grepl("tRNA|rRNA|snRNA|srpRNA|Simple_repeat|Low_complexity|Satellite|Unknown",
                  superfam_fine, ignore.case = TRUE)]

dat[, log2_rna := log2(RNA_signal + 1)]
dat[, mac_label := fifelse(has_mac, "MAC at TSS", "No MAC")]
dat[, mac_label := factor(mac_label, levels = c("No MAC", "MAC at TSS"))]

cat(sprintf("Total rows: %d\n", nrow(dat)))

# ── Select top 10 families by total copy count per species ────────────────────
all_counts <- dat[, .(n_copies = .N), by = .(species, superfam_fine)]
ok_all     <- dat[, .(n = .N), by = .(species, superfam_fine, mac_label)]
ok_wide    <- dcast(ok_all, species + superfam_fine ~ mac_label, value.var = "n", fill = 0L)
ok_wide    <- ok_wide[`No MAC` >= MIN_N & `MAC at TSS` >= MIN_N]
eligible   <- all_counts[ok_wide[, .(species, superfam_fine)],
                          on = .(species, superfam_fine), nomatch = NULL]
top_fams   <- eligible[order(species, -n_copies)][, head(.SD, TOP_N_FAM), by = species]
cat("\nTop", TOP_N_FAM, "families per species (after MIN_N filter, by copy count):\n")
print(top_fams)

plot_dat <- dat[top_fams[, .(species, superfam_fine)],
                on = .(species, superfam_fine), nomatch = NULL]
cat(sprintf("Families: %d total, %d rows\n", nrow(top_fams), nrow(plot_dat)))
plot_dat[, fam_f := factor(superfam_fine)]

# ── Wilcoxon per family × species ─────────────────────────────────────────────
wtest <- plot_dat[, {
  a <- log2_rna[mac_label == "MAC at TSS"]
  b <- log2_rna[mac_label == "No MAC"]
  if (length(a) >= MIN_N && length(b) >= MIN_N) {
    wt <- wilcox.test(a, b, exact = FALSE)
    .(p_val = wt$p.value, n_mac = length(a), n_no_mac = length(b),
      med_mac = median(a), med_no_mac = median(b))
  } else {
    .(p_val = NA_real_, n_mac = length(a), n_no_mac = length(b),
      med_mac = NA_real_, med_no_mac = NA_real_)
  }
}, by = .(species, superfam_fine)]
wtest[, p_adj := p.adjust(p_val, method = "BH"), by = species]
wtest[, sig_label := fcase(
  p_adj < 0.001, "***",
  p_adj < 0.01,  "**",
  p_adj < 0.05,  "*",
  default = "ns"
)]
cat("\n=== Wilcoxon results (MAC at TSS vs no MAC) ===\n")
print(wtest[order(species, p_adj)])
fwrite(wtest, file.path(OUTDIR, "44b_wilcoxon.tsv"), sep = "\t")

# n labels per group
n_labels <- plot_dat[, .(n = .N), by = .(species, superfam_fine, mac_label)]

# ── Colour palette ────────────────────────────────────────────────────────────
MAC_COLS <- c("No MAC" = "#4575b4", "MAC at TSS" = "#d73027")

# Map superfamily name -> TE type (order matters: specific first)
te_type_from_fam <- function(fam) {
  type <- rep("Unclassified", length(fam))
  type[grepl("Helitron",         fam, ignore.case=TRUE)] <- "Rolling Circle"
  type[grepl("Penelope|Pluto",   fam, ignore.case=TRUE)] <- "Penelope"
  type[grepl("\\bSINE\\b",       fam, ignore.case=TRUE)] <- "SINE"
  type[grepl("Gypsy|Copia|DIRS|Bel.Pao|Ngaro|ERV|BEL",
             fam, ignore.case=TRUE)] <- "LTR"
  type[grepl("LINE|\\bL1\\b|\\bL2\\b|RTE|I-Jockey|Jockey|TART|CR1|Ingi",
             fam, ignore.case=TRUE)] <- "LINE"
  type[grepl("CMC|TcMar|MULE|MuDR|EnSpm|PIF|Harbinger|Zisupton|Crypton|Maverick|Academ|Chapaev|Dada|Ginger|hAT|Kolobok|Mariner|Novosib|Polinton|Transib|ZaMe|Zator|ISRm",
             fam, ignore.case=TRUE)] <- "DNA"
  type
}

TYPE_COLS <- c(
  "DNA"            = "#4472C4",
  "LTR"            = "#E8736C",
  "LINE"           = "#3E9E8C",
  "SINE"           = "#F5C242",
  "Rolling Circle" = "#F07F00",
  "Penelope"       = "#9B3D8A",
  "Unclassified"   = "#AAAAAA"
)

# ── Build one panel per species ────────────────────────────────────────────────
make_panel <- function(sp_label) {
  d   <- plot_dat[species == sp_label]
  wt  <- wtest[species == sp_label]
  nl  <- n_labels[species == sp_label]
  if (nrow(d) == 0) return(NULL)

  # Order families by this species' copy count (descending)
  sp_levels <- top_fams[species == sp_label][order(-n_copies), superfam_fine]
  d[,  fam_f := factor(superfam_fine, levels = sp_levels)]
  nl[, fam_f := factor(superfam_fine, levels = sp_levels)]
  wt[, fam_f := factor(superfam_fine, levels = sp_levels)]

  # Sample points for jitter (max 400 per family × group)
  set.seed(42)
  pts <- d[, .SD[sample(.N, min(.N, 400))], by = .(superfam_fine, mac_label)]

  # X-axis label colours by TE type
  x_types <- te_type_from_fam(levels(d$fam_f))
  x_cols  <- TYPE_COLS[x_types]

  ggplot(d, aes(x = fam_f, y = log2_rna, fill = mac_label)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, linewidth = 0.35,
                 position = position_dodge(width = 0.75), width = 0.65) +
    geom_point(data = pts,
               aes(colour = mac_label, group = mac_label),
               alpha = 0.12, size = 0.35,
               position = position_jitterdodge(jitter.width = 0.2,
                                               dodge.width = 0.75),
               show.legend = FALSE) +
    geom_text(data = wt,
              aes(x = fam_f, y = Inf, label = sig_label),
              inherit.aes = FALSE, size = 3.5, fontface = "bold", vjust = 1.5) +
    geom_text(data = nl,
              aes(x = fam_f, y = -Inf,
                  label = paste0("n=", format(n, big.mark = ",")),
                  group = mac_label),
              position = position_dodge(width = 0.75),
              inherit.aes = FALSE, size = 2.3, vjust = -0.3, colour = "grey40") +
    scale_fill_manual(values = MAC_COLS, name = "MAC at TSS") +
    scale_colour_manual(values = MAC_COLS, guide = "none") +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = expression(log[2](RNA + 1)),
         title = sp_label,
         subtitle = sprintf("Top %d families by copy count | Dark condition | n >= %d per group",
                            TOP_N_FAM, MIN_N)) +
    theme_bw(base_size = 10) +
    theme(axis.text.x       = element_text(angle = 35, hjust = 1, size = 9,
                                              colour = x_cols, face = "bold"),
          panel.grid.major.x = element_blank(),
          plot.margin        = margin(5, 5, 18, 5),
          legend.position    = "right")
}

p_phy <- make_panel("Phycomyces blakesleeanus")
p_rhi <- make_panel("Rhizopus microsporus")

panels <- Filter(Negate(is.null), list(p_phy, p_rhi))
p_combined <- wrap_plots(panels, ncol = 1) +
  plot_annotation(
    title    = "TE expression: MAC at TSS vs no MAC - top 10 families per species",
    subtitle = "Dark/MiceliumDark condition | Wilcoxon BH-adjusted: * p<0.05  ** p<0.01  *** p<0.001"
  )

out_pdf <- file.path(OUTDIR, "44b_mac_tss_families.pdf")
out_png <- file.path(OUTDIR, "44b_mac_tss_families.png")
ggsave(out_pdf, p_combined, width = 12, height = 6)
ggsave(out_png, p_combined, width = 12, height = 6, dpi = 150)
cat(sprintf("Saved: %s\n", out_pdf))

# ── Summary TSV ───────────────────────────────────────────────────────────────
summ <- plot_dat[, .(
  n_total  = .N,
  n_mac    = sum(mac_label == "MAC at TSS"),
  n_no_mac = sum(mac_label == "No MAC"),
  med_mac    = median(log2_rna[mac_label == "MAC at TSS"],  na.rm = TRUE),
  med_no_mac = median(log2_rna[mac_label == "No MAC"],      na.rm = TRUE)
), by = .(species, superfam_fine)]
summ <- merge(summ, wtest[, .(species, superfam_fine, p_val, p_adj, sig_label)],
              by = c("species","superfam_fine"), all.x = TRUE)
fwrite(summ, file.path(OUTDIR, "44b_summary.tsv"), sep = "\t")
cat("Done.\n")
