#!/usr/bin/env Rscript
# =============================================================================
# plot_circular_tree.R   — v4
#
# Changes vs v3:
#   - Tracks tighter (offset = 0.01, pwidth = 0.08)
#   - Each track uses a visually distinct colour palette
#   - Clade arc labels placed OUTSIDE all data tracks
#   - Tree topology rotated to standard fungal arrangement:
#       Holozoa | Microsporidia | early fungi | Chytrids+NU |
#       Glomeromycota | Mucoro | Mortierello | Zoopago | Kickx |
#       Entomoph | Basidiobolo | Basidiomycota | Ascomycota
#
# Usage:
#   Rscript plot_circular_tree.R [--tree HOclade]
# =============================================================================

.libPaths(c("~/R/library", .libPaths()))

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(viridisLite)
})

args         <- commandArgs(trailingOnly = TRUE)
use_HO_clade <- "--tree" %in% args &&
                args[which(args == "--tree") + 1] == "HOclade"

# ─── Paths ────────────────────────────────────────────────────────────────────
if (use_HO_clade) {
  TREE_FILE <- "/data/carlos/EDF-TE/SpeciesTree/08_guide_tree/guide_rooted_HOclade.treefile"
  OUT_PDF   <- "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures/circular_tree_HOclade_v4.pdf"
  cat("[MODE] HO-clade-rooted (test)\n")
} else {
  TREE_FILE <- "/data/carlos/EDF-TE/SpeciesTree/08_guide_tree/guide_rooted_Naegrub.treefile"
  OUT_PDF   <- "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures/circular_tree_v4.pdf"
  cat("[MODE] Naegrub-rooted (original)\n")
}
DATA_FILE <- "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figure_species_data.tsv"
dir.create(dirname(OUT_PDF), showWarnings = FALSE, recursive = TRUE)

# ─── Load ─────────────────────────────────────────────────────────────────────
tree <- read.tree(TREE_FILE)
dat  <- read.delim(DATA_FILE, stringsAsFactors = FALSE) %>%
  rename(label = tree_label) %>%
  filter(!is.na(label)) %>%
  select(label, everything())

cat("[INFO] Tree tips:", length(tree$tip.label),
    "  Matched:", sum(dat$label %in% tree$tip.label), "\n")

# ─── Clade colours ────────────────────────────────────────────────────────────
clade_cols <- c(
  Holozoa               = "#888888", Microsporidia    = "#d32f2f",
  Cryptomycota          = "#e64a19", Aphelida         = "#f57c00",
  Blastocladiomycota    = "#fbc02d", Chytridiomycota  = "#558b2f",
  Neocallimastigomycota = "#00897b", Olpidiomycota    = "#0277bd",
  Zoopagomycota         = "#4527a0", Kickxellomycotina= "#9c27b0",
  Entomophthoromycota   = "#c2185b", Basidiobolales   = "#6a1c1c",
  Mortierellomycota     = "#ef6c00", Mucoromycota     = "#f9a825",
  Glomeromycota         = "#2e7d32", Basidiomycota    = "#1565c0",
  Ascomycota            = "#6a1b9a", Unknown          = "#9e9e9e"
)

# ─── MRCA nodes (guide_rooted_Naegrub) ──────────────────────────────────────
clade_nodes_naegrub <- c(
  Ascomycota = 276, Basidiomycota = 288, Microsporidia = 249,
  Mucoromycota = 309, Mortierellomycota = 371, Chytridiomycota = 419,
  Blastocladiomycota = 455, Cryptomycota = 247, Zoopagomycota = 399,
  Kickxellomycotina = 407, Entomophthoromycota = 403, Glomeromycota = 303,
  Holozoa = 239, Aphelida = 246, Basidiobolales = 398,
  Olpidiomycota = 174, Neocallimastigomycota = 459
)

if (use_HO_clade) {
  prefix_map <- list(
    Ascomycota="^AS_", Basidiomycota="^BA_", Microsporidia="^MI_",
    Mucoromycota="^MU_", Mortierellomycota="^MO_", Chytridiomycota="^CH_",
    Blastocladiomycota="^BL_", Cryptomycota="^CR_", Zoopagomycota="^ZO_",
    Kickxellomycotina="^KX_", Entomophthoromycota="^EN_", Glomeromycota="^GL_",
    Holozoa="^HO_", Aphelida="^AP_", Basidiobolales="^BS_",
    Olpidiomycota="^OL_", Neocallimastigomycota="^NU_"
  )
  clade_nodes <- sapply(names(prefix_map), function(nm) {
    tips <- tree$tip.label[grepl(prefix_map[[nm]], tree$tip.label)]
    if (length(tips) > 1) as.integer(MRCA(tree, tips))
    else if (length(tips) == 1) which(tree$tip.label == tips)
    else NA_integer_
  })
} else {
  clade_nodes <- clade_nodes_naegrub
}

# ─── Helpers ──────────────────────────────────────────────────────────────────
mk_dot <- function(col, size_range = c(0.2, 6.5)) {
  dat %>%
    select(label, value = !!sym(col)) %>%
    filter(!is.na(value)) %>%
    mutate(dot_size = rescale(value, to = size_range))
}
mk_dot_log <- function(col, size_range = c(0.2, 6.5)) {
  dat %>%
    select(label, value = !!sym(col)) %>%
    filter(!is.na(value)) %>%
    mutate(dot_size = rescale(log10(pmax(value, 0.01)), to = size_range))
}

brks3 <- function(col) {
  x <- dat[[col]]
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(0, 0.5, 1))
  r <- range(x)
  if (r[1] == r[2]) return(c(r[1], r[1] + 1e-6, r[2] + 2e-6))
  c(r[1], (r[1] + r[2]) / 2, r[2])
}

brks3_log <- function(col) {
  x <- dat[[col]]
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) return(c(0.1, 1, 10))
  r <- range(x)
  if (r[1] == r[2]) return(c(r[1], r[1] * 1.0001, r[2] * 1.0002))
  c(r[1], 10 ^ ((log10(r[1]) + log10(r[2])) / 2), r[2])
}

FIRST_DOT_OFF <- 0.006   # keep the first track very close to the tips
DOT_OFF       <- 0.001   # minimal spacing between tracks
DOT_PW        <- 0.050   # track width (wider to fit larger dots)
BAR_OFF       <- 0.020   # extra separation from the %TE track to avoid overlap
cl      <- function(nm, ord)
  guide_legend(title.position = "top", title = nm, order = ord,
               title.theme = element_text(size = 4.8, face = "bold"),
               label.theme = element_text(size = 4.5),
               override.aes = list(size = c(2.0, 3.6, 5.0), alpha = 1))

# ─── Base tree ────────────────────────────────────────────────────────────────
p <- ggtree(tree, layout = "fan", open.angle = 35,
            linewidth = 0.20, color = "grey20") %<+% dat +
  geom_tippoint(aes(color = clade), size = 0.4, alpha = 0.85) +
  scale_color_manual(values = clade_cols, name = NULL,
    guide = guide_legend(ncol = 3, order = 1,
                         override.aes = list(size = 2.2), title = NULL)) +
  theme(legend.position = "bottom",
        legend.text   = element_text(size = 4.8),
        legend.key.size = unit(0.26, "cm"),
        plot.margin   = margin(15, 15, 15, 15))

# (no internal highlights — clade identity shown by tip points + outer arc ring)

# ─── Track 1: % multi-isoform  — lime green ───────────────────────────────────
d9 <- mk_dot("pct_multiisoform")
p <- p + new_scale_color() +
  geom_fruit(data = d9, geom = geom_point,
    mapping = aes(y = label, x = 1, color = value, size = I(dot_size)),
    offset = FIRST_DOT_OFF, pwidth = DOT_PW) +
  scale_color_gradient(low = "#f4ffe0", high = "#33691e",
    breaks = brks3("pct_multiisoform"),
    na.value = NA, guide = cl("% multi-\nisoform\ngenes", 2))

# ─── Track 2: % unknown genes  — steel blue ───────────────────────────────────
d1 <- mk_dot("pct_unknown")
p <- p + new_scale_color() +
  geom_fruit(data = d1, geom = geom_point,
    mapping = aes(y = label, x = 1, color = value, size = I(dot_size)),
    offset = DOT_OFF, pwidth = DOT_PW) +
  scale_color_gradient(low = "#dceefb", high = "#1a3a6b",
    breaks = brks3("pct_unknown"),
    na.value = "grey90", guide = cl("% unknown\ngenes", 3))

# ─── Track 3: % paralog genes  — deep crimson ────────────────────────────────
d2 <- mk_dot("pct_paralog_genes")
p <- p + new_scale_color() +
  geom_fruit(data = d2, geom = geom_point,
    mapping = aes(y = label, x = 1, color = value, size = I(dot_size)),
    offset = DOT_OFF, pwidth = DOT_PW) +
  scale_color_gradient(low = "#fde8e8", high = "#7f0000",
    breaks = brks3("pct_paralog_genes"),
    na.value = "grey90", guide = cl("% paralog\ngenes", 4))

# ─── Track 4: % OGs > fungal average  — forest green ─────────────────────────
d3 <- mk_dot("pct_OGs_above_fungal_avg")
p <- p + new_scale_color() +
  geom_fruit(data = d3, geom = geom_point,
    mapping = aes(y = label, x = 1, color = value, size = I(dot_size)),
    offset = DOT_OFF, pwidth = DOT_PW) +
  scale_color_gradient(low = "#e8f5e9", high = "#1b5e20",
    breaks = brks3("pct_OGs_above_fungal_avg"),
    na.value = "grey90", guide = cl("% OGs >\nfungal avg", 5))

# ─── Track 5: Genome size  — warm brown, log scale ───────────────────────────
d4 <- mk_dot_log("genome_size_mb", size_range = c(0.5, 3.5))
p <- p + new_scale_color() +
  geom_fruit(data = d4, geom = geom_point,
    mapping = aes(y = label, x = 1, color = value, size = I(dot_size)),
    offset = DOT_OFF, pwidth = DOT_PW) +
  scale_color_gradient(low = "#fde8d0", high = "#7b2400",
    breaks = brks3_log("genome_size_mb"),
    trans = "log10", na.value = "grey90",
    guide = cl("Genome\nsize (Mb)", 6))

# ─── Track 6: Gene density  — amber ──────────────────────────────────────────
d5 <- mk_dot_log("gene_density")
p <- p + new_scale_color() +
  geom_fruit(data = d5, geom = geom_point,
    mapping = aes(y = label, x = 1, color = value, size = I(dot_size)),
    offset = DOT_OFF, pwidth = DOT_PW) +
  scale_color_gradient(low = "#fff8e1", high = "#bf6000",
    breaks = brks3("gene_density"),
    na.value = "grey90", guide = cl("Gene density\n(genes/Mb)", 7))

# ─── Track 7: Avg introns / gene  — deep teal ────────────────────────────────
d6 <- mk_dot("Avg_Introns_All")
p <- p + new_scale_color() +
  geom_fruit(data = d6, geom = geom_point,
    mapping = aes(y = label, x = 1, color = value, size = I(dot_size)),
    offset = DOT_OFF, pwidth = DOT_PW) +
  scale_color_gradient(low = "#e0f4f4", high = "#004c4c",
    breaks = brks3("Avg_Introns_All"),
    na.value = "grey90", guide = cl("Avg introns\n/ gene", 8))

# ─── Track 8: % total TE  — violet/magenta ───────────────────────────────────
d7 <- mk_dot("pct_TE_total")
p <- p + new_scale_color() +
  geom_fruit(data = d7, geom = geom_point,
    mapping = aes(y = label, x = 1, color = value, size = I(dot_size)),
    offset = DOT_OFF, pwidth = DOT_PW) +
  scale_color_gradient(low = "#f8e8ff", high = "#5e0080",
    breaks = brks3("pct_TE_total"),
    na.value = "grey90", guide = cl("% total TE\n(genome)", 9))

# ─── Track 9: TE class proportion — normalised stacked bar ───────────────────
# Values normalised to 0-1 so bars stay within pwidth in circular layout
te_norm <- dat %>%
  mutate(total_cls = pct_ClassI + pct_ClassII) %>%
  filter(!is.na(total_cls), total_cls > 0) %>%
  mutate(`Class I`  = pct_ClassI  / total_cls,
         `Class II` = pct_ClassII / total_cls) %>%
  select(label, `Class I`, `Class II`) %>%
  pivot_longer(c("Class I","Class II"), names_to = "te_class", values_to = "frac")

p <- p + new_scale_fill() +
  geom_fruit(data = te_norm, geom = geom_bar,
    mapping     = aes(y = label, x = frac, fill = te_class),
    stat        = "identity", orientation = "y",
    offset      = BAR_OFF, pwidth = 0.040) +
  scale_fill_manual(
    values = c("Class I" = "#d97706", "Class II" = "#5b21b6"),
    name   = "TE class\n(prop. within TEs)",
    guide  = guide_legend(title.position = "top", order = 10,
               title.theme = element_text(size = 4.8, face = "bold"),
               label.theme = element_text(size = 4.5),
               override.aes = list(size = 3), nrow = 1))

# ─── Final theme ──────────────────────────────────────────────────────────────
p <- p + theme(
  legend.position   = "bottom",
  legend.box        = "horizontal",
  legend.text       = element_text(size = 4.5),
  legend.title      = element_text(size = 5.0, face = "bold"),
  legend.key.size   = unit(0.26, "cm"),
  legend.spacing.x  = unit(0.10, "cm"),
  plot.margin       = margin(20, 20, 20, 20),
  panel.border      = element_blank(),
  panel.background  = element_rect(fill = "transparent", colour = NA)
)

# ─── Save ─────────────────────────────────────────────────────────────────────
cat("[PLOT] Saving:", OUT_PDF, "\n")
ggsave(OUT_PDF, plot = p, width = 44, height = 46, units = "cm",
       device = "pdf", limitsize = FALSE)
cat("[OK] Done!\n")
