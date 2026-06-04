#!/usr/bin/env Rscript
# fig_A_tree_cooption.R  v4 — OG co-option: Dollo parsimony gains at ancestors

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
cat("Collapsed tree:", length(tr_c$tip.label), "tips\n")
# NWK is pre-rotated; ladderize=FALSE in ggtree preserves the order directly.# Rotate MO to appear below GL and MU (flip children of their MRCA node).
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

  # Losses: any descendant tips of MRCA that DON'T have this OG
  mrca_desc <- desc_tips_list[[mrca_nd]]
  lost_tips <- setdiff(mrca_desc, tip_idx)
  for (lt in lost_tips) loss_counts[lt] <- loss_counts[lt] + 1L
}

cat("[INFO] Total gains across all nodes:", sum(gain_counts), "\n")
cat("[INFO] Total losses across all nodes:", sum(loss_counts), "\n")

# ─── Node annotation table ───────────────────────────────────────────────────
n_og_tip_distinct <- vapply(tr_c$tip.label, function(cl) {
  length(unique(step9$orthogroup[step9$clade == cl]))
}, integer(1))
is_tip  <- c(rep(TRUE, n_tip), rep(FALSE, n_int))
cl_node <- c(tr_c$tip.label, rep(NA_character_, n_int))

# Number shown: tips = n_og_distinct; internal = gain count
node_number <- c(n_og_tip_distinct, gain_counts[(n_tip+1):n_all])
loss_number <- loss_counts  # losses at every node (mostly at tips)

# Bubble size proportional to the displayed number
bub_size <- rescale(sqrt(pmax(node_number, 1)), to=c(3, 20))

# Labels: tips show OG count; internal nodes show "+N gains / -M losses"
lbl_tip      <- as.character(n_og_tip_distinct)
lbl_internal <- ifelse(
  loss_counts[(n_tip+1):n_all] > 0,
  paste0("+", gain_counts[(n_tip+1):n_all], " / -", loss_counts[(n_tip+1):n_all]),
  paste0("+", gain_counts[(n_tip+1):n_all])
)
lbl <- c(lbl_tip, lbl_internal)
lbl[!is_tip & gain_counts[(n_tip+1):n_all] == 0 & loss_counts[(n_tip+1):n_all] == 0] <- ""

node_df <- data.frame(
  node=seq_len(n_all), is_tip=is_tip,
  fill_color = "#cccccc",   # all same grey, semi-transparent
  bub_size=bub_size, node_label=lbl,
  lbl_hjust=ifelse(is_tip, -0.45, -0.15),
  clade_code=cl_node, clade_name=CLADE_NAMES[cl_node],
  stringsAsFactors=FALSE
)

# ─── Plot ─────────────────────────────────────────────────────────────────────
p <- ggtree(tr_c, layout="rectangular", branch.length="none", ladderize=FALSE,
            linewidth=1.1, color="grey20") %<+% node_df +
  geom_point2(aes(size=bub_size, fill=fill_color),
              shape=21, color="grey50", stroke=0.5, alpha=0.55) +
  scale_fill_identity() + scale_size_identity() +
  geom_text2(aes(label=node_label, hjust=lbl_hjust), vjust=0.5,
             size=2.7, color="grey10", fontface="bold") +
  geom_tiplab(aes(label=clade_name, color=clade_code),
              offset=0.38, size=3.8, fontface="bold") +
  scale_color_manual(values=CLADE_COLORS, guide="none") +
  xlim(NA, max(node.depth(tr_c)) + 3.2) +
  theme_tree2() +
  theme(legend.position="none",
        plot.title=element_text(size=11, face="bold"),
        plot.subtitle=element_text(size=7.5, color="grey45"),
        axis.text.x=element_blank(), axis.ticks=element_blank(),
        plot.margin=margin(10, 50, 10, 10)) +
  labs(title="Co-option of TE domains across the eukaryote tree of life",
       subtitle=paste0(
         "Tips: distinct OG families with co-option (circle size \u221d count).\n",
         "Internal nodes: +N gains / -M losses (Dollo parsimony)."))

ggsave(file.path(OUT_DIR, "figA_tree_cooption.pdf"), p, width=12, height=10)
ggsave(file.path(OUT_DIR, "figA_tree_cooption.png"), p, width=12, height=10, dpi=150)
cat("Saved figA_tree_cooption.pdf/.png\n")
