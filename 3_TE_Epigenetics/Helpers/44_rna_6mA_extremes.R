#!/usr/bin/env Rscript
# 44_rna_6mA_extremes.R
# Boxplot: TE families with largest contrast between
#   Group A: high RNA expression (Q4) + 6mA present
#   Group B: low  RNA expression (Q1) + 6mA absent
# Both species, Dark condition, ≥10 copies per group, no tRNA families

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
WDIR   <- "/legserv/Temp/Carlos/EDF-TE/results/TE_dynamics"
OUTDIR <- file.path(WDIR, "44_rna_6mA_extremes")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

PHY_RNA <- file.path(WDIR, "20_chiapello_extension/data/te_phybl2_rna_6mA.tsv")
RHI_RNA <- file.path(WDIR, "37_6mA_expression/te_rhimi59_rna_6mA.tsv")

TOP_N_FAM <- 15
MIN_N     <- 10

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

# ── Load and prepare data per species ─────────────────────────────────────────
prep_species <- function(tsv_path, species_cond, species_label, n_col_fix = FALSE) {
  dt <- fread(tsv_path)
  if (n_col_fix && ncol(dt) >= 8) {
    setnames(dt, names(dt),
             c("te_id","condition","RNA_signal","n6mA","superfam_fine",
               "superfam_old","len","density_6mA","has_6mA")[seq_len(ncol(dt))])
  }
  # Standardise has_6mA to logical
  if (is.character(dt$has_6mA)) dt[, has_6mA := has_6mA == "TRUE"]
  dt <- dt[condition == species_cond]
  dt[, superfam_broad := classify_sf(superfam_fine)]
  # Exclude non-TE annotations
  dt <- dt[!grepl("tRNA|rRNA|snRNA|srpRNA|Simple_repeat|Low_complexity|Satellite",
                  superfam_fine, ignore.case = TRUE)]
  # RNA quartiles
  q25 <- quantile(dt$RNA_signal, 0.25, na.rm = TRUE)
  q75 <- quantile(dt$RNA_signal, 0.75, na.rm = TRUE)
  dt[, rna_group := fcase(
    RNA_signal <= q25, "Q1",
    RNA_signal >= q75, "Q4",
    default = NA_character_
  )]
  dt[, species := species_label]
  dt[]
}

cat("Loading Phycomyces (MiceliumDark)...\n")
phy <- prep_species(PHY_RNA, "MiceliumDark", "Phycomyces blakesleeanus")

cat("Loading Rhizopus (Dark)...\n")
rhi <- prep_species(RHI_RNA, "Dark", "Rhizopus microsporus", n_col_fix = TRUE)

dat <- rbind(phy, rhi, fill = TRUE)
cat(sprintf("Total rows: %d\n", nrow(dat)))

# ── Per-family contrast analysis ───────────────────────────────────────────────
# Group A: Q4 high RNA AND has_6mA = TRUE
# Group B: Q1 low  RNA AND has_6mA = FALSE
dat[, extremes_group := fcase(
  rna_group == "Q4" & has_6mA == TRUE,  "A: Q4 + 6mA",
  rna_group == "Q1" & has_6mA == FALSE, "B: Q1 - 6mA",
  default = NA_character_
)]

# Keep only extreme groups
dat_ex <- dat[!is.na(extremes_group)]

# Per-family per-species, compute n per group and median RNA
fam_stats <- dat_ex[, .(
  n       = .N,
  med_rna = median(log2(RNA_signal + 1), na.rm = TRUE)
), by = .(species, superfam_fine, superfam_broad, extremes_group)]

# Need both groups present and n >= MIN_N for each
fam_both <- fam_stats[n >= MIN_N, .(n_groups = .N), by = .(species, superfam_fine)]
fam_both <- fam_both[n_groups == 2, .(species, superfam_fine)]

fam_stats <- fam_stats[fam_both, on = .(species, superfam_fine)]

# Compute contrast: median(A) - median(B) per species × family
fam_wide <- dcast(fam_stats,
                  species + superfam_fine + superfam_broad ~ extremes_group,
                  value.var = "med_rna")
setnames(fam_wide, c("A: Q4 + 6mA", "B: Q1 - 6mA"), c("med_A", "med_B"))
fam_wide[, contrast := med_A - med_B]

# Top families per species by contrast
top_fams <- fam_wide[order(-contrast), head(.SD, TOP_N_FAM), by = species]
cat("\n=== Top families with highest A-B contrast ===\n")
print(top_fams[, .(species, superfam_fine, med_A, med_B, contrast)])

fwrite(top_fams, file.path(OUTDIR, "44_top_families_contrast.tsv"), sep = "\t")

# ── Plot ───────────────────────────────────────────────────────────────────────
# Get raw data for top families
top_fam_ids <- top_fams[, .(species, superfam_fine)]
plot_dat <- dat_ex[top_fam_ids, on = .(species, superfam_fine), nomatch = 0]
plot_dat <- plot_dat[!is.na(extremes_group)]
plot_dat[, log2_rna := log2(RNA_signal + 1)]

# For each species: order families by contrast
order_phy <- top_fams[species == "Phycomyces blakesleeanus"][order(-contrast), superfam_fine]
order_rhi <- top_fams[species == "Rhizopus microsporus"     ][order(-contrast), superfam_fine]

# Truncate long family names
plot_dat[, fam_label := sub("/rnd-.*", "", superfam_fine)]
plot_dat[, fam_label := sub("^DNA/", "", fam_label)]
plot_dat[, fam_label := sub("^LTR/", "", fam_label)]
plot_dat[, fam_label := sub("^RC/",  "", fam_label)]
plot_dat[, fam_label := sub("^LINE/","", fam_label)]

top_fams[, fam_label := sub("/rnd-.*", "", superfam_fine)]
top_fams[, fam_label := sub("^DNA/|^LTR/|^RC/|^LINE/", "", fam_label)]

GROUP_COLS2 <- c("A: Q4 + 6mA" = "#d73027", "B: Q1 - 6mA" = "#4575b4")

make_panel <- function(sp, fam_order) {
  d  <- plot_dat[species == sp]
  tf <- top_fams[species == sp]
  # Ordered factor for families
  ord  <- tf[match(fam_order, tf$superfam_fine), fam_label]
  d[, fam_f := factor(fam_label,
                      levels = rev(ord[!is.na(ord)]))]
  tf[, fam_f := factor(fam_label, levels = rev(ord[!is.na(ord)]))]

  broad_map <- unique(d[, .(fam_label, superfam_broad)])
  broad_map[, fam_f := factor(fam_label, levels = levels(d$fam_f))]

  # Colour points by broad category
  pt_col <- broad_map[order(fam_f), superfam_broad]
  names(pt_col) <- levels(d$fam_f)

  ggplot(d, aes(x = fam_f, y = log2_rna, fill = extremes_group)) +
    geom_jitter(aes(fill = extremes_group), alpha = 0.15, size = 0.35,
                position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
                show.legend = FALSE) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.4, alpha = 0.85,
                 position = position_dodge(0.75)) +
    coord_flip() +
    scale_fill_manual(values = GROUP_COLS2, name = NULL) +
    labs(x = NULL, y = expression(log[2](RNA + 1)),
         subtitle = sp) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          panel.grid.major.y = element_blank(),
          strip.text = element_text(size = 8))
}

p_phy <- make_panel("Phycomyces blakesleeanus", order_phy)
p_rhi <- make_panel("Rhizopus microsporus",     order_rhi)

p_combined <- (p_phy / p_rhi) +
  plot_annotation(
    title = "TE families: high expression + 6mA present vs low expression + 6mA absent",
    subtitle = sprintf("Dark condition | top %d families per species | n ≥ %d per group",
                       TOP_N_FAM, MIN_N),
    theme = theme(plot.title = element_text(size = 11, face = "bold"))
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

out_pdf <- file.path(OUTDIR, "44_rna_6mA_extremes.pdf")
out_png <- file.path(OUTDIR, "44_rna_6mA_extremes.png")
ggsave(out_pdf, p_combined, width = 10, height = 12)
ggsave(out_png, p_combined, width = 10, height = 12, dpi = 150)
cat(sprintf("Saved: %s\n", out_pdf))
cat("Done.\n")
