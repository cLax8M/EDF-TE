#!/usr/bin/env Rscript
# =============================================================================
# plot_rectangular_tree.R   — v1
#
# Rectangular phylogram (real branch lengths) with the same data tracks as
# the circular tree v4.  All species, proportional branches.
#
# Output: rectangular_tree_v1.pdf
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
})

# ─── Paths ────────────────────────────────────────────────────────────────────
TREE_FILE <- "/legserv/Temp/Carlos/EDF-TE/results/SpeciesTree/final_LG_G_rotated.nwk"
DATA_FILE <- "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figure_species_data.tsv"
OUT_PDF   <- "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures/rectangular_tree_v1.pdf"
dir.create(dirname(OUT_PDF), showWarnings = FALSE, recursive = TRUE)

# ─── Load ─────────────────────────────────────────────────────────────────────
options(expressions = 100000L)
tree <- read.tree(TREE_FILE)
dat  <- read.delim(DATA_FILE, stringsAsFactors = FALSE) %>%
  rename(label = tree_label) %>%
  mutate(label = sub("\\.aa$", "", label)) %>%
  filter(!is.na(label)) %>%
  select(label, everything())

# ─── Rotate MO below GL+MU ────────────────────────────────────────────────────
NE_SPECIES_circ <- c("CH_Anasp1","CH_Neocon1","CH_Neolan1","CH_Orpsp1",
                     "CH_Pecora1","CH_Caecom1","CH_Pirfi3")
get_clade_circ <- function(sp) {
  ifelse(sp %in% NE_SPECIES_circ, "NE", substr(sp, 1L, 2L))
}
get_clade_rep <- function(clade) tree$tip.label[get_clade_circ(tree$tip.label)==clade][1]
mo_rep <- get_clade_rep("MO"); gl_rep <- get_clade_rep("GL"); mu_rep <- get_clade_rep("MU")
if (!any(is.na(c(mo_rep, gl_rep, mu_rep)))) {
  tree <- ape::rotate(tree, getMRCA(tree, c(mo_rep, gl_rep, mu_rep)))
  cat("[INFO] Rotated MO below GL+MU\n")
}

# ─── Force CR together ────────────────────────────────────────────────────────
if ("CR_Rozal1" %in% tree$tip.label && "CR_Parsa1" %in% tree$tip.label) {
  tree <- ape::drop.tip(tree, "CR_Rozal1")
  where_parsa <- which(tree$tip.label == "CR_Parsa1")
  tree <- phytools::bind.tip(tree, tip.label = "CR_Rozal1",
                             where = where_parsa, position = 0)
  cat("[INFO] Forced CR_Rozal1 sister to CR_Parsa1 (display only)\n")
}

cat("[INFO] Tree tips:", length(tree$tip.label),
    "  Matched:", sum(dat$label %in% tree$tip.label), "\n")

# ─── Clade colours (canonical palette) ────────────────────────────────────────
clade_cols <- c(
  Ascomycota            = "#e8887a", Basidiomycota        = "#4472c4",
  Mucoromycota          = "#e07840", Mortierellomycota    = "#f4a8b0",
  Glomeromycota         = "#8b5e48", Kickxellomycotina    = "#f0a820",
  Zoopagomycota         = "#c8dc70", Entomophthoromycota  = "#e01010",
  Basidiobolales        = "#d0b4e0", Olpidiomycota        = "#3668b0",
  Chytridiomycota       = "#e0d040", Nucleariida          = "#1a2850",
  Blastocladiomycota    = "#7030a0", Microsporidia        = "#3aab98",
  Cryptomycota          = "#6aafe6", Aphelida             = "#5c9e3a",
  Holozoa               = "#b4d0e8", Unknown              = "#9e9e9e"
)

# ─── Helpers ──────────────────────────────────────────────────────────────────
mk_dot <- function(col, size_range = c(0.15, 4.0)) {
  dat %>%
    select(label, value = !!sym(col)) %>%
    filter(!is.na(value)) %>%
    mutate(dot_size = rescale(value, to = size_range))
}
mk_dot_log <- function(col, size_range = c(0.15, 4.0)) {
  dat %>%
    select(label, value = !!sym(col)) %>%
    filter(!is.na(value)) %>%
    mutate(dot_size = rescale(log10(pmax(value, 0.01)), to = size_range))
}
brks3 <- function(col) {
  x <- dat[[col]]; x <- x[is.finite(x)]
  if (length(x)==0) return(c(0,0.5,1))
  r <- range(x)
  if (r[1]==r[2]) return(c(r[1],r[1]+1e-6,r[2]+2e-6))
  c(r[1],(r[1]+r[2])/2,r[2])
}
brks3_log <- function(col) {
  x <- dat[[col]]; x <- x[is.finite(x) & x>0]
  if (length(x)==0) return(c(0.1,1,10))
  r <- range(x)
  if (r[1]==r[2]) return(c(r[1],r[1]*1.0001,r[2]*1.0002))
  c(r[1], 10^((log10(r[1])+log10(r[2]))/2), r[2])
}

FIRST_DOT_OFF <- 0.15
DOT_OFF       <- 0.02
DOT_PW        <- 0.06
BAR_OFF       <- 0.06

cl <- function(nm, ord)
  guide_legend(title.position="top", title=nm, order=ord,
               title.theme=element_text(size=4.8, face="bold"),
               label.theme=element_text(size=4.5),
               override.aes=list(size=c(2.0,3.6,5.0), alpha=1))

# ─── Base tree (rectangular, real branch lengths) ─────────────────────────────
p <- ggtree(tree, layout = "rectangular", ladderize = FALSE,
            linewidth = 0.18, color = "grey20") %<+% dat +
  geom_tippoint(aes(color = clade), size = 0.7, alpha = 0.85) +
  scale_color_manual(values = clade_cols, name = NULL,
    guide = guide_legend(ncol = 2, order = 1,
                         override.aes = list(size = 2.5), title = NULL)) +
  theme_tree2() +
  theme(legend.position  = "right",
        legend.text      = element_text(size = 5),
        legend.key.size  = unit(0.28, "cm"),
        axis.text.x      = element_text(size = 5),
        plot.margin      = margin(10, 5, 10, 10))

# ─── Track 1: % multi-isoform ─────────────────────────────────────────────────
d9 <- mk_dot("pct_multiisoform")
p <- p + new_scale_color() +
  geom_fruit(data=d9, geom=geom_point,
    mapping=aes(y=label, x=1, color=value, size=I(dot_size)),
    offset=FIRST_DOT_OFF, pwidth=DOT_PW) +
  scale_color_gradient(low="#f4ffe0", high="#33691e",
    breaks=brks3("pct_multiisoform"), na.value=NA,
    guide=cl("% multi-\nisoform\ngenes",2))

# ─── Track 2: % unknown genes ─────────────────────────────────────────────────
d1 <- mk_dot("pct_unknown")
p <- p + new_scale_color() +
  geom_fruit(data=d1, geom=geom_point,
    mapping=aes(y=label, x=1, color=value, size=I(dot_size)),
    offset=DOT_OFF, pwidth=DOT_PW) +
  scale_color_gradient(low="#dceefb", high="#1a3a6b",
    breaks=brks3("pct_unknown"), na.value="grey90",
    guide=cl("% unknown\ngenes",3))

# ─── Track 3: % paralog genes ─────────────────────────────────────────────────
d2 <- mk_dot("pct_paralog_genes")
p <- p + new_scale_color() +
  geom_fruit(data=d2, geom=geom_point,
    mapping=aes(y=label, x=1, color=value, size=I(dot_size)),
    offset=DOT_OFF, pwidth=DOT_PW) +
  scale_color_gradient(low="#fde8e8", high="#7f0000",
    breaks=brks3("pct_paralog_genes"), na.value="grey90",
    guide=cl("% paralog\ngenes",4))

# ─── Track 4: % OGs > fungal average ─────────────────────────────────────────
d3 <- mk_dot("pct_OGs_above_fungal_avg")
p <- p + new_scale_color() +
  geom_fruit(data=d3, geom=geom_point,
    mapping=aes(y=label, x=1, color=value, size=I(dot_size)),
    offset=DOT_OFF, pwidth=DOT_PW) +
  scale_color_gradient(low="#e8f5e9", high="#1b5e20",
    breaks=brks3("pct_OGs_above_fungal_avg"), na.value="grey90",
    guide=cl("% OGs >\nfungal avg",5))

# ─── Track 5: Genome size (log) ───────────────────────────────────────────────
d4 <- mk_dot_log("genome_size_mb", size_range=c(0.3,3.0))
p <- p + new_scale_color() +
  geom_fruit(data=d4, geom=geom_point,
    mapping=aes(y=label, x=1, color=value, size=I(dot_size)),
    offset=DOT_OFF, pwidth=DOT_PW) +
  scale_color_gradient(low="#fde8d0", high="#7b2400",
    breaks=brks3_log("genome_size_mb"), trans="log10", na.value="grey90",
    guide=cl("Genome\nsize (Mb)",6))

# ─── Track 6: Gene density ────────────────────────────────────────────────────
d5 <- mk_dot_log("gene_density")
p <- p + new_scale_color() +
  geom_fruit(data=d5, geom=geom_point,
    mapping=aes(y=label, x=1, color=value, size=I(dot_size)),
    offset=DOT_OFF, pwidth=DOT_PW) +
  scale_color_gradient(low="#fff8e1", high="#bf6000",
    breaks=brks3("gene_density"), na.value="grey90",
    guide=cl("Gene density\n(genes/Mb)",7))

# ─── Track 7: Avg introns / gene ──────────────────────────────────────────────
d6 <- mk_dot("Avg_Introns_All")
p <- p + new_scale_color() +
  geom_fruit(data=d6, geom=geom_point,
    mapping=aes(y=label, x=1, color=value, size=I(dot_size)),
    offset=DOT_OFF, pwidth=DOT_PW) +
  scale_color_gradient(low="#e0f4f4", high="#004c4c",
    breaks=brks3("Avg_Introns_All"), na.value="grey90",
    guide=cl("Avg introns\n/ gene",8))

# ─── Track 8: % total TE ──────────────────────────────────────────────────────
d7 <- mk_dot("pct_TE_total")
p <- p + new_scale_color() +
  geom_fruit(data=d7, geom=geom_point,
    mapping=aes(y=label, x=1, color=value, size=I(dot_size)),
    offset=DOT_OFF, pwidth=DOT_PW) +
  scale_color_gradient(low="#f8e8ff", high="#5e0080",
    breaks=brks3("pct_TE_total"), na.value="grey90",
    guide=cl("% total TE\n(genome)",9))

# ─── Track 9: TE class proportion (stacked bar) ───────────────────────────────
te_norm <- dat %>%
  mutate(total_cls = pct_ClassI + pct_ClassII) %>%
  filter(!is.na(total_cls), total_cls > 0) %>%
  mutate(`Class I`  = pct_ClassI  / total_cls,
         `Class II` = pct_ClassII / total_cls) %>%
  select(label, `Class I`, `Class II`) %>%
  pivot_longer(c("Class I","Class II"), names_to="te_class", values_to="frac")

p <- p + new_scale_fill() +
  geom_fruit(data=te_norm, geom=geom_bar,
    mapping=aes(y=label, x=frac, fill=te_class),
    stat="identity", orientation="y",
    offset=BAR_OFF, pwidth=0.045) +
  scale_fill_manual(
    values=c("Class I"="#d97706","Class II"="#5b21b6"),
    name="TE class\n(prop. within TEs)",
    guide=guide_legend(title.position="top", order=10,
             title.theme=element_text(size=4.8, face="bold"),
             label.theme=element_text(size=4.5),
             override.aes=list(size=3), nrow=1))

# ─── Save ─────────────────────────────────────────────────────────────────────
cat("[PLOT] Saving:", OUT_PDF, "\n")
ggsave(OUT_PDF, plot=p, width=60, height=120, units="cm",
       device="pdf", limitsize=FALSE)
cat("[OK] Done!\n")
