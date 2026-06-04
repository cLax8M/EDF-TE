#!/usr/bin/env Rscript
# fig_B_clade_boxplot.R  v2 — Proportion of OGs with co-option per clade

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(scales)
})

OG_FILE <- "/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures/species_og_cooption_stats.tsv"
OUT_DIR <- "/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures"

CLADE_COLORS <- c(
  AS="#e8887a", BA="#4472c4", MU="#e07840", MO="#f4a8b0", GL="#8b5e48",
  KX="#f0a820", ZO="#c8dc70", EN="#e01010", BS="#d0b4e0", OL="#3668b0",
  CH="#e0d040", NE="#9060b0", BL="#7030a0", MI="#3aab98", CR="#a8c8e8",
  AP="#5c9e3a", NU="#1a2850", HO="#b4d0e8"
)
CLADE_NAMES <- c(
  AS="Ascomycota", BA="Basidiomycota", MU="Mucoromycotina",
  MO="Mortierellomycotina", GL="Glomeromycotina", KX="Kickxellomycotina",
  ZO="Zoopagomycotina", EN="Entomophthoromycotina", BS="Basidiobolomycota",
  OL="Olpidiomycota", CH="Chytridiomycetes", NE="Neocallimastigomycetes",
  BL="Blastocladiomycota", MI="Microsporidia", CR="Cryptomycota",
  AP="Aphelidiomycota", NU="Nucleariida", HO="Holozoa+Amoebozoa"
)
NE_SPECIES <- c("CH_Anasp1","CH_Neocon1","CH_Neolan1","CH_Orpsp1",
                "CH_Pecora1","CH_Caecom1","CH_Pirfi3")
# Phylogenetic order from the rotated tree (HO = outgroup, AS/BA at the crown end)
CLADE_ORDER <- c("HO","NU","AP","CR","MI","BL","OL","CH","NE",
                 "BS","EN","ZO","KX","GL","MO","MU","BA","AS")

og <- read.table(OG_FILE, header=TRUE, sep="\t", stringsAsFactors=FALSE)
og$clade <- ifelse(og$species %in% NE_SPECIES, "NE", substr(og$species, 1, 2))
og <- og[og$clade %in% names(CLADE_COLORS), ]
og$clade_label <- CLADE_NAMES[og$clade]

# Require ≥3 species per clade
clade_n <- table(og$clade)
og <- og[og$clade %in% names(clade_n[clade_n >= 3]), ]

# Order by phylogenetic position (tree order: outgroup → crown; Asco/Basidio at right)
ordered_codes  <- CLADE_ORDER[CLADE_ORDER %in% unique(og$clade)]
ordered_labels <- CLADE_NAMES[ordered_codes]
og$clade_label <- factor(og$clade_label, levels=ordered_labels)

# Colour vector aligned to factor levels
label2code <- setNames(names(CLADE_NAMES), CLADE_NAMES)
fill_vals <- CLADE_COLORS[label2code[levels(og$clade_label)]]
fill_vals[is.na(fill_vals)] <- "grey60"
names(fill_vals) <- levels(og$clade_label)

p <- ggplot(og, aes(x=clade_label, y=pct_og_coopted, fill=clade_label)) +
  geom_boxplot(outlier.shape=NA, alpha=0.72, width=0.55, color="grey30", linewidth=0.5) +
  geom_jitter(aes(color=clade_label), width=0.16, size=1.8, alpha=0.75) +
  scale_fill_manual(values=fill_vals,  guide="none") +
  scale_color_manual(values=fill_vals, guide="none") +
  coord_cartesian(ylim=c(0, 5)) +
  labs(x=NULL,
       y="% orthogroups with \u22651 co-opted gene",
       title="Proportion of gene families with TE co-option per clade",
       subtitle="One dot = one species. Clades with \u22653 species. Ordered by phylogenetic position (outgroup \u2192 crown).") +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10, face="bold",
                                 color=fill_vals[levels(og$clade_label)]),
        panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold"),
        plot.margin=margin(10, 10, 10, 60))

ggsave(file.path(OUT_DIR, "figB_clade_og_proportion.pdf"), p, width=12, height=7)
ggsave(file.path(OUT_DIR, "figB_clade_og_proportion.png"), p, width=12, height=7, dpi=150)
cat("Saved figB_clade_og_proportion.pdf/.png\n")
