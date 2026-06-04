#!/usr/bin/env Rscript
# fig_A_tree_cooption_v2.1.R  — OG co-option: Dollo parsimony gains + losses shown
# v2.1: FUNGI-ONLY — Holozoa+Amoebozoa (HO) excluded from tree and Dollo computation.
# HO is dropped after tree collapse so inode losses inflated by HO-private OGs are removed.
# Gains (blue, +N) AND losses (red, -M) shown at ALL nodes (internal + tips).
# Dollo parsimony: single gain at MRCA; each loss assigned to the root of the
# highest subtree where the OG is completely absent (minimises loss events).

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
cat("Clades:\n"); print(clade_sum %>% arrange(desc(median_pct)))

# ─── Collapse tree ────────────────────────────────────────────────────────────
clade_reps <- og %>%
  group_by(clade) %>%
  arrange(abs(pct_og_coopted - median(pct_og_coopted, na.rm=TRUE))) %>%
  slice(1) %>% ungroup()

tr_c <- drop.tip(tree, setdiff(tree$tip.label, clade_reps$species))
tr_c$tip.label <- clade_reps$clade[match(tr_c$tip.label, clade_reps$species)]
# v2.1: drop HO (Holozoa+Amoebozoa outgroup) — fungi-only analysis
if ("HO" %in% tr_c$tip.label) tr_c <- drop.tip(tr_c, "HO")
cat("Collapsed tree (v2.1, no HO):", length(tr_c$tip.label), "tips\n")
# NWK is pre-rotated; ladderize=FALSE in ggtree preserves the order directly.
# Rotate MO to appear below GL and MU (flip children of their MRCA node).
if (all(c("MO","GL","MU") %in% tr_c$tip.label)) {
  tr_c <- ape::rotate(tr_c, getMRCA(tr_c, c("MO","GL","MU")))
}
# ─── fastAnc on pct_og_coopted (for internal node annotation only) ────────────
tip_vec <- setNames(clade_sum$median_pct[match(tr_c$tip.label, clade_sum$clade)], tr_c$tip.label)
tip_vec[is.na(tip_vec)] <- 0
anc_vals <- fastAnc(tr_c, tip_vec, vars=FALSE, CI=FALSE)

n_tip <- length(tr_c$tip.label); n_int <- tr_c$Nnode; n_all <- n_tip + n_int

# ─── Dollo parsimony: OG gains AND losses per node ───────────────────────────
cat("[INFO] Computing Dollo-parsimony OG gains+losses on collapsed tree...\n")

step9 <- read.table(STEP9_FILE, header=TRUE, sep="\t", stringsAsFactors=FALSE)
step9 <- step9[step9$n_coopted_genes_in_species_og > 0, ]
step9$clade <- get_clade(step9$species)
step9 <- step9[step9$clade %in% tr_c$tip.label, ]

og_clades <- step9 %>%
  group_by(orthogroup) %>%
  summarise(clades = list(unique(clade)), .groups="drop")

cat("[INFO]", nrow(og_clades), "OGs with co-option across", length(tr_c$tip.label), "collapsed clades\n")

gain_counts <- integer(n_all)
loss_counts <- integer(n_all)

# Helper: get all descendant tips of a node
get_desc_tips <- function(nd) {
  desc <- integer(0)
  visit <- function(n2) {
    ch <- tr_c$edge[tr_c$edge[,1L]==n2, 2L]
    if (length(ch)==0L) { desc <<- c(desc, n2); return() }
    for (ci in ch) visit(ci)
  }
  visit(nd); desc
}

# Pre-compute descendant tip sets for all nodes
desc_tips_list <- vector("list", n_all)
for (nd in 1L:n_all) desc_tips_list[[nd]] <- get_desc_tips(nd)

# Pre-compute immediate children for each node
node_children <- vector("list", n_all)
for (nd in 1L:n_all) node_children[[nd]] <- tr_c$edge[tr_c$edge[,1L]==nd, 2L]

# Dollo parsimony loss assignment:
# Each loss is placed at the root of the highest subtree fully lacking the OG.
# This minimises loss events (vs counting each absent tip separately).
find_dollo_losses <- function(present_tips, mrca_nd) {
  result <- integer(0)
  check <- function(nd) {
    desc <- desc_tips_list[[nd]]
    if (length(intersect(desc, present_tips)) == 0L) {
      result <<- c(result, nd); return()   # whole subtree absent -> 1 loss here
    }
    for (ch in node_children[[nd]]) check(ch)   # some present -> go deeper
  }
  for (ch in node_children[[mrca_nd]]) check(ch)
  result
}

find_mrca_set <- function(clade_set) {
  tip_idx <- match(clade_set, tr_c$tip.label)
  tip_idx <- tip_idx[!is.na(tip_idx)]
  if (length(tip_idx) == 0) return(NA_integer_)
  if (length(tip_idx) == 1) return(tip_idx[1])
  getMRCA(tr_c, tip_idx)
}

for (i in seq_len(nrow(og_clades))) {
  og_tips <- og_clades$clades[[i]]  # tip codes that HAVE this OG
  tip_idx <- match(og_tips, tr_c$tip.label)
  tip_idx <- tip_idx[!is.na(tip_idx)]
  if (length(tip_idx) == 0) next

  mrca_nd <- if (length(tip_idx)==1) tip_idx[1] else getMRCA(tr_c, tip_idx)
  if (is.na(mrca_nd) || mrca_nd < 1 || mrca_nd > n_all) next

  # Gain at MRCA
  gain_counts[mrca_nd] <- gain_counts[mrca_nd] + 1L

  # Losses: Dollo parsimony - assign each loss to the root of the highest absent subtree
  for (ln in find_dollo_losses(tip_idx, mrca_nd))
    loss_counts[ln] <- loss_counts[ln] + 1L
}

cat("[INFO] Total gains across all nodes:", sum(gain_counts), "\n")
cat("[INFO] Total losses across all nodes:", sum(loss_counts), "\n")

# ─── Node annotation table ───────────────────────────────────────────────────
n_og_tip_distinct <- vapply(tr_c$tip.label, function(cl) {
  length(unique(step9$orthogroup[step9$clade == cl]))
}, integer(1))
is_tip  <- c(rep(TRUE, n_tip), rep(FALSE, n_int))
cl_node <- c(tr_c$tip.label, rep(NA_character_, n_int))

# ── Bubble size: proportional to gains (internal) or observed (tips) ──────────
node_number <- c(n_og_tip_distinct, gain_counts[(n_tip+1):n_all])
bub_size <- rescale(sqrt(pmax(node_number, 1)), to=c(3, 20))

# ── Bubble fill: blue = gains, red-pink = losses only, grey = neither ──────────
bub_fill <- ifelse(is_tip,
  ifelse(loss_counts > 0, "#e8887a", "#cccccc"),
  ifelse(gain_counts > 0, "#5b9bd5",
    ifelse(loss_counts > 0, "#e8887a", "#cccccc")))

# ── Labels ─────────────────────────────────────────────────────────────────────
# Internal nodes: "+N" in blue
lbl_gain <- c(
  rep("", n_tip),
  ifelse(gain_counts[(n_tip+1):n_all] > 0,
         paste0("+", format(gain_counts[(n_tip+1):n_all], big.mark=",")),
         "")
)
# Tip observed count (black) and losses (red) at ALL nodes
lbl_obs      <- c(format(n_og_tip_distinct, big.mark=","), rep("", n_int))
loss_fmt     <- ifelse(loss_counts > 0,
                       paste0("\u2212", format(loss_counts, big.mark=",")),
                       "")
# Split: tip losses use below-tip position; internal node losses to the left of node
lbl_loss_tip <- ifelse(is_tip,  loss_fmt, "")
lbl_loss_int <- ifelse(!is_tip, loss_fmt, "")

node_df <- data.frame(
  node         = seq_len(n_all),
  is_tip       = is_tip,
  fill_color   = bub_fill,
  bub_size     = bub_size,
  lbl_gain     = lbl_gain,     # blue "+N" at internal nodes
  lbl_obs      = lbl_obs,      # black observed count at tips
  lbl_loss_tip = lbl_loss_tip, # red "-M" below tip nodes
  lbl_loss_int = lbl_loss_int, # red "-M" left of internal nodes
  lbl_hjust    = ifelse(is_tip, -0.35, -0.15),
  clade_code   = cl_node,
  clade_name   = CLADE_NAMES[cl_node],
  stringsAsFactors = FALSE
)

# ─── Plot ─────────────────────────────────────────────────────────────────────
p <- ggtree(tr_c, layout="rectangular", branch.length="none", ladderize=FALSE,
            linewidth=1.1, color="grey20") %<+% node_df +
  # Bubbles
  geom_point2(aes(size=bub_size, fill=fill_color),
              shape=21, color="grey50", stroke=0.5, alpha=0.60) +
  scale_fill_identity() + scale_size_identity() +
  # Internal node gain labels (blue, above center)
  geom_text2(aes(label=lbl_gain, hjust=lbl_hjust),
             vjust=-0.55, size=2.5, color="#1a5fa8", fontface="bold") +
  # Tip observed labels (dark, above center) — tips: observed count only, no gain/loss
  geom_text2(aes(label=lbl_obs, hjust=lbl_hjust),
             vjust=-0.55, size=2.5, color="grey15", fontface="bold") +
  # Internal node loss labels (red, to the LEFT of branching point = on parent branch)
  geom_text2(aes(label=lbl_loss_int),
             hjust=1.15, vjust=0.5, size=2.5, color="#c0392b", fontface="bold") +
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
  labs(title="Co-option of TE domains across the fungal tree of life",
       subtitle=paste0(
         "Fungi-only (v2.1: HO outgroup excluded). ",
         "Tips: observed OG families (black, above) / Dollo losses at highest absent clade (red, below).\n",
         "Internal nodes: +N gains (blue, above) / \u2212M Dollo losses at ancestral lineages (red, below)."))

ggsave(file.path(OUT_DIR, "figA_tree_cooption_v2.1.pdf"), p, width=12, height=10)
ggsave(file.path(OUT_DIR, "figA_tree_cooption_v2.1.png"), p, width=12, height=10, dpi=150)
cat("Saved figA_tree_cooption_v2.1.pdf/.png\n")
