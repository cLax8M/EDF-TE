#!/usr/bin/env Rscript
# fig_A_tree_cooption_v3.R
# Co-option tree with SUBTREE-EXCLUSIVE OG counts at internal nodes.
#
# At each internal node N, the annotation shows the number of OG families
# (orthogroups with co-opted TE domains) that are present in at least one
# descendant clade of N and ABSENT from every clade outside N's subtree.
# This answers: "how many co-option innovations are exclusive to this lineage
# vs. shared with distant outgroups?"
#
# This differs from Dollo parsimony (v2) which shows gains at the EXACT
# ancestral node (+1 at Dikarya = 1 OG found in BOTH AS AND BA).
# v3 shows cumulative lineage-restricted innovations (227 at Dikarya =
# 227 OG families in AS and/or BA but in no other clade).
#
# Tips: observed OG family count (black, above tip label).
# Internal nodes: subtree-exclusive count (blue bubble + label).
# Root node: suppressed (trivially = all OGs).

suppressPackageStartupMessages({
  library(ape); library(ggtree); library(ggplot2)
  library(dplyr); library(phytools); library(scales)
})

TREE_FILE  <- "/legserv/Temp/Carlos/EDF-TE/results/SpeciesTree/final_LG_G_rotated.nwk"
OG_FILE    <- "/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures/species_og_cooption_stats.tsv"
STEP9_FILE <- "/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures/step9_mixed_fate_orthogroups.tsv"
OUT_DIR    <- "/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures"

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
get_clade <- function(sp) ifelse(sp %in% NE_SPECIES, "NE", substr(sp, 1, 2))

# ─── Load & prepare ──────────────────────────────────────────────────────────
tree <- read.newick(TREE_FILE)
og   <- read.table(OG_FILE, header=TRUE, sep="\t", stringsAsFactors=FALSE)
og$clade <- get_clade(og$species)
og <- og[og$clade %in% names(CLADE_COLORS), ]

clade_sum <- og %>%
  group_by(clade) %>%
  summarise(n_sp=n(),
            median_pct = median(pct_og_coopted, na.rm=TRUE),
            median_n   = median(n_og_with_coopt, na.rm=TRUE),
            .groups="drop")

# ─── Collapse tree ────────────────────────────────────────────────────────────
clade_reps <- og %>%
  group_by(clade) %>%
  arrange(abs(pct_og_coopted - median(pct_og_coopted, na.rm=TRUE))) %>%
  slice(1) %>% ungroup()

tr_c <- drop.tip(tree, setdiff(tree$tip.label, clade_reps$species))
tr_c$tip.label <- clade_reps$clade[match(tr_c$tip.label, clade_reps$species)]
cat("Collapsed tree:", length(tr_c$tip.label), "tips\n")

if (all(c("MO","GL","MU") %in% tr_c$tip.label)) {
  tr_c <- ape::rotate(tr_c, getMRCA(tr_c, c("MO","GL","MU")))
}

n_tip <- length(tr_c$tip.label); n_int <- tr_c$Nnode; n_all <- n_tip + n_int
root_nd <- n_tip + 1L

# ─── Step9: per-clade OG presence ────────────────────────────────────────────
step9 <- read.table(STEP9_FILE, header=TRUE, sep="\t", stringsAsFactors=FALSE)
step9 <- step9[step9$n_coopted_genes_in_species_og > 0, ]
step9$clade <- get_clade(step9$species)
step9 <- step9[step9$clade %in% tr_c$tip.label, ]

og_clades <- step9 %>%
  group_by(orthogroup) %>%
  summarise(clades = list(unique(clade)), .groups="drop")

cat("[INFO]", nrow(og_clades), "OGs with co-option across", length(tr_c$tip.label), "collapsed clades\n")

# ─── Precompute descendant tip sets ──────────────────────────────────────────
get_desc_tips <- function(nd) {
  desc <- integer(0)
  visit <- function(n2) {
    ch <- tr_c$edge[tr_c$edge[,1L]==n2, 2L]
    if (length(ch)==0L) { desc <<- c(desc, n2); return() }
    for (ci in ch) visit(ci)
  }
  visit(nd); desc
}
desc_tips_list <- lapply(1L:n_all, get_desc_tips)

# ─── Subtree-exclusive counts ─────────────────────────────────────────────────
# For each OG, precompute the tip indices that have it
cat("[INFO] Computing subtree-exclusive OG counts...\n")
og_tip_idx <- lapply(seq_len(nrow(og_clades)), function(i) {
  idx <- match(og_clades$clades[[i]], tr_c$tip.label)
  as.integer(idx[!is.na(idx)])
})

# subtree_excl[nd] = number of OGs where ALL clades are within desc_tips_list[[nd]]
subtree_excl <- vapply(1L:n_all, function(nd) {
  dt <- desc_tips_list[[nd]]
  sum(vapply(og_tip_idx, function(idx) length(idx)>0 && all(idx %in% dt), logical(1)))
}, integer(1))

# Suppress root (trivially = all OGs, not informative)
subtree_excl[root_nd] <- NA_integer_

cat("[INFO] Subtree-exclusive counts computed.\n")
cat("[INFO] Dikarya (AS+BA):", subtree_excl[getMRCA(tr_c, c("AS","BA"))], "\n")

# ─── Tip observed counts ──────────────────────────────────────────────────────
n_og_tip_distinct <- vapply(tr_c$tip.label, function(cl) {
  length(unique(step9$orthogroup[step9$clade == cl]))
}, integer(1))

is_tip  <- c(rep(TRUE, n_tip), rep(FALSE, n_int))
cl_node <- c(tr_c$tip.label, rep(NA_character_, n_int))

# ─── Bubble size ──────────────────────────────────────────────────────────────
# Tips: based on observed count; internal: based on subtree_excl
internal_excl <- subtree_excl[(n_tip+1):n_all]
internal_excl[is.na(internal_excl)] <- 0L
node_number <- c(n_og_tip_distinct, internal_excl)
bub_size <- rescale(sqrt(pmax(node_number, 1)), to=c(3, 20))

# ─── Bubble fill ──────────────────────────────────────────────────────────────
bub_fill <- ifelse(is_tip, "#cccccc", "#5b9bd5")
# Tips with ANY exclusive OGs get a slightly warmer tint
tip_has_excl <- n_og_tip_distinct > 0
bub_fill[1:n_tip][tip_has_excl] <- "#a8c0d8"

# ─── Labels ───────────────────────────────────────────────────────────────────
# Internal nodes: exclusive count (blue); suppress root
lbl_excl <- c(
  rep("", n_tip),
  ifelse(is.na(subtree_excl[(n_tip+1):n_all]), "",
    format(subtree_excl[(n_tip+1):n_all], big.mark=","))
)
# Tips: observed count (black)
lbl_obs <- c(format(n_og_tip_distinct, big.mark=","), rep("", n_int))

node_df <- data.frame(
  node       = seq_len(n_all),
  is_tip     = is_tip,
  fill_color = bub_fill,
  bub_size   = bub_size,
  lbl_excl   = lbl_excl,
  lbl_obs    = lbl_obs,
  lbl_hjust  = ifelse(is_tip, -0.35, -0.15),
  clade_code = cl_node,
  clade_name = CLADE_NAMES[cl_node],
  stringsAsFactors = FALSE
)

# ─── Plot ─────────────────────────────────────────────────────────────────────
p <- ggtree(tr_c, layout="rectangular", branch.length="none", ladderize=FALSE,
            linewidth=1.1, color="grey20") %<+% node_df +
  # Bubbles
  geom_point2(aes(size=bub_size, fill=fill_color),
              shape=21, color="grey50", stroke=0.5, alpha=0.60) +
  scale_fill_identity() + scale_size_identity() +
  # Internal node exclusive labels (blue)
  geom_text2(aes(label=lbl_excl, hjust=lbl_hjust),
             vjust=-0.55, size=2.5, color="#1a5fa8", fontface="bold") +
  # Tip observed labels (dark)
  geom_text2(aes(label=lbl_obs, hjust=lbl_hjust),
             vjust=-0.55, size=2.5, color="grey15", fontface="bold") +
  # Clade name labels
  geom_tiplab(aes(label=clade_name, color=clade_code),
              offset=0.42, size=3.8, fontface="bold") +
  scale_color_manual(values=CLADE_COLORS, guide="none") +
  xlim(NA, max(node.depth(tr_c)) + 3.2) +
  theme_tree2() +
  theme(legend.position="none",
        plot.title    = element_text(size=11, face="bold"),
        plot.subtitle = element_text(size=7.5, color="grey45"),
        axis.text.x   = element_blank(), axis.ticks=element_blank(),
        plot.margin   = margin(10, 50, 10, 10)) +
  labs(title="Co-option of TE domains: lineage-restricted innovations",
       subtitle=paste0(
         "Tips: observed co-opted OG families (black).\n",
         "Internal nodes (blue): OG families exclusive to that subtree, ",
         "i.e., present in ≥1 descendant clade and absent from all outgroups.\n",
         "Root suppressed (trivially equals total OGs). Bubble size ∝ exclusive count."))

ggsave(file.path(OUT_DIR, "figA_tree_cooption_v3.pdf"), p, width=12, height=10)
ggsave(file.path(OUT_DIR, "figA_tree_cooption_v3.png"), p, width=12, height=10, dpi=150)
cat("Saved figA_tree_cooption_v3.pdf/.png\n")
warnings()
