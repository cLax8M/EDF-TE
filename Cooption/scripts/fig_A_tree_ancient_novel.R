#!/usr/bin/env Rscript
# fig_A_tree_ancient_novel.R  — Ancient OG retention and Fungal novelty
# Visual style mirrors figA_tree_cooption_v2 (cladogram, tip bubbles, compact
# node labels).
#
# TIPS  : bubble sized by total observed co-opted OGs (= ancient + novel),
#         coloured by clade; ABOVE = total OG count (black, like v2);
#         BELOW = ancestral HO OGs still retained in this clade (red).
#
# INTERNAL NODES: +N (blue, above)  = novel OG families first gained on this
#                                      stem (Dollo MRCA exactly at this node)
#                 -M (red,  left)   = ancestral HO OGs lost on this stem
#                                      (= ancient_abs[parent] − ancient_abs[node])

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

# ── Load data ─────────────────────────────────────────────────────────────────
tree  <- read.newick(TREE_FILE)
og    <- read.table(OG_FILE,    header=TRUE, sep="\t", stringsAsFactors=FALSE)
step9 <- read.table(STEP9_FILE, header=TRUE, sep="\t", stringsAsFactors=FALSE)

og$clade    <- get_clade(og$species)
step9$clade <- get_clade(step9$species)

# Keep only species with data in co-option
og <- og[og$clade %in% names(CLADE_COLORS), ]

# ── HO ancestral OG set ───────────────────────────────────────────────────────
ho_ogs <- unique(step9$orthogroup[
  grepl("^HO_", step9$species) & step9$n_coopted_genes_in_species_og > 0
])
cat("[INFO] HO ancestral OGs:", length(ho_ogs), "\n")

# ── Per-clade OG sets (including HO for the tree) ────────────────────────────
step9_active <- step9[step9$n_coopted_genes_in_species_og > 0, ]

clade_ogs <- step9_active %>%
  group_by(clade) %>%
  summarise(ogs = list(unique(orthogroup)), .groups="drop")
clade_ogs_map <- setNames(clade_ogs$ogs, clade_ogs$clade)

# ── Collapse tree to clade representatives ───────────────────────────────────
clade_reps <- og %>%
  group_by(clade) %>%
  arrange(abs(pct_og_coopted - median(pct_og_coopted, na.rm=TRUE))) %>%
  slice(1) %>% ungroup()

tr_c <- drop.tip(tree, setdiff(tree$tip.label, clade_reps$species))
tr_c$tip.label <- clade_reps$clade[match(tr_c$tip.label, clade_reps$species)]

# Rotate for aesthetics
if (all(c("MO","GL","MU") %in% tr_c$tip.label))
  tr_c <- ape::rotate(tr_c, getMRCA(tr_c, c("MO","GL","MU")))

cat("Collapsed tree tips:", paste(tr_c$tip.label, collapse=", "), "\n")

n_tip <- length(tr_c$tip.label)
n_all <- n_tip + tr_c$Nnode

# ── Helper: all descendant nodes (internal + tip) of a node ──────────────────
get_all_descendants <- function(nd) {
  # Returns all nodes in the subtree rooted at nd (including nd itself)
  desc <- c(nd)
  visit <- function(n2) {
    ch <- tr_c$edge[tr_c$edge[,1L] == n2, 2L]
    if (length(ch) == 0L) return()
    desc <<- c(desc, ch)
    for (ci in ch) visit(ci)
  }
  visit(nd); desc
}

# Helper: descendant TIPS only
get_desc_tips <- function(nd) {
  desc <- integer(0)
  visit <- function(n2) {
    ch <- tr_c$edge[tr_c$edge[,1L] == n2, 2L]
    if (length(ch) == 0L) { desc <<- c(desc, n2); return() }
    for (ci in ch) visit(ci)
  }
  visit(nd); desc
}

# Tip index by label
tip_idx <- setNames(seq_len(n_tip), tr_c$tip.label)

# ── Tip total OG counts (= ancient + novel = total observed per clade) ────────
n_og_tip <- vapply(tr_c$tip.label, function(cl)
  length(unique(step9_active$orthogroup[step9_active$clade == cl])),
  integer(1))

# ── ANCIENT counts: HO OGs present in ≥1 descendant tip (DECREASES root→tips) ─
# ancient_abs[X] = |HO OGs ∩ union(OG sets of all descendant tips of X)|
# Pre-compute descendant tip lists (reused below)
desc_tips_list <- lapply(seq_len(n_all), get_desc_tips)

ancient_abs <- integer(n_all)
for (nd in seq_len(n_all)) {
  tip_labs       <- tr_c$tip.label[desc_tips_list[[nd]]]
  tip_labs_fungi <- tip_labs[tip_labs != "HO"]   # ancient = OGs retained in FUNGI, not HO itself
  all_ogs  <- unique(unlist(lapply(tip_labs_fungi, function(cl) clade_ogs_map[[cl]])))
  ancient_abs[nd] <- sum(all_ogs %in% ho_ogs, na.rm=TRUE)
}

# ── ANCIENT LOSSES per stem branch ────────────────────────────────────────────
# ancient_loss[X] = ancient_abs[parent(X)] − ancient_abs[X]  (≥ 0)
parent_of <- function(nd) {
  idx <- which(tr_c$edge[,2L] == nd)
  if (length(idx)==0L) NA_integer_ else tr_c$edge[idx, 1L]
}
ancient_loss <- integer(n_all)
for (nd in seq_len(n_all)) {
  pa <- parent_of(nd)
  if (!is.na(pa)) ancient_loss[nd] <- max(0L, ancient_abs[pa] - ancient_abs[nd])
}

cat("[INFO] Root ancient OGs:", ancient_abs[n_tip+1L], "\n")
cat("[INFO] Tip ancient (retained):", paste(tr_c$tip.label, ancient_abs[seq_len(n_tip)], sep="=", collapse=", "), "\n")

# ── NOVEL OG FAMILIES: MRCA-exact count per stem branch ───────────────────────
# For each non-HO OG, find its MRCA in the collapsed fungal tree.
# novel_mrca[X] = count of non-HO OGs whose MRCA is exactly node X.
# This mirrors v2's +gains logic (single Dollo gain at MRCA).

step9_fungi <- step9_active[!grepl("^HO_", step9_active$species), ]
tip_idx_map <- setNames(seq_len(n_tip), tr_c$tip.label)

og2clades <- step9_fungi %>%
  filter(!orthogroup %in% ho_ogs) %>%
  group_by(orthogroup) %>%
  summarise(clades = list(intersect(unique(clade), tr_c$tip.label)), .groups="drop") %>%
  filter(lengths(clades) > 0)

cat("[INFO] Non-HO fungal OGs:", nrow(og2clades), "\n")

novel_mrca <- integer(n_all)
for (i in seq_len(nrow(og2clades))) {
  ti <- tip_idx_map[og2clades$clades[[i]]]
  ti <- ti[!is.na(ti)]
  if (length(ti) == 0L) next
  M <- if (length(ti) == 1L) ti else getMRCA(tr_c, ti)
  if (is.null(M) || is.na(M)) next
  novel_mrca[M] <- novel_mrca[M] + 1L
}

cat("[INFO] Internal novel MRCA counts:", sum(novel_mrca[(n_tip+1L):n_all]), "\n")
cat("[INFO] Tip novel MRCA counts:     ", sum(novel_mrca[seq_len(n_tip)]), "\n")

# ── Build node annotation data frame ──────────────────────────────────────────
is_tip_v     <- seq_len(n_all) <= n_tip
n_int        <- tr_c$Nnode
clade_code_v <- c(tr_c$tip.label, rep(NA_character_, n_int))

# Bubble: sized by total OGs at tips (same scale as figA_v2), tiny dot at internals
bub_tip_total <- c(n_og_tip, rep(0L, n_int))
bub_size <- ifelse(is_tip_v,
  rescale(sqrt(pmax(bub_tip_total, 1)), to=c(3, 20)),
  1.5)

# Bubble fill: clade colour at tips (like v2), light grey at internal nodes
bub_fill <- ifelse(is_tip_v,
  CLADE_COLORS[clade_code_v],
  "#dddddd")
bub_fill[is.na(bub_fill)] <- "#cccccc"

# Label: total OGs at tips only (black, above bubble)
lbl_total <- c(format(n_og_tip, big.mark=","), rep("", n_int))

# Label: ancestral HO OGs retained in descendants — at internal nodes only (red)
# ancient_abs[node] = union of HO OGs across all tip descendants → decreases root→tips
lbl_int_anc <- ifelse(
  !is_tip_v & ancient_abs > 0,
  as.character(ancient_abs),
  "")

# Label: novel OG families first co-opted at exactly this node — internal nodes (blue, "+N")
lbl_int_novel <- c(
  rep("", n_tip),
  ifelse(novel_mrca[(n_tip+1L):n_all] > 0,
         paste0("+", novel_mrca[(n_tip+1L):n_all]), "")
)

# Cumulative novel OGs from root to each internal node (green, increasing root→tips)
# BFS from root, summing novel_mrca along the path
cumul_novel <- novel_mrca
root_node   <- n_tip + 1L
cumul_novel[root_node] <- novel_mrca[root_node]
bfs_queue <- root_node
while (length(bfs_queue) > 0) {
  curr <- bfs_queue[1]; bfs_queue <- bfs_queue[-1]
  children <- tr_c$edge[tr_c$edge[, 1] == curr, 2]
  for (ch in children) {
    cumul_novel[ch] <- cumul_novel[curr] + novel_mrca[ch]
    if (ch > n_tip) bfs_queue <- c(bfs_queue, ch)
  }
}
lbl_int_cumul <- ifelse(!is_tip_v & cumul_novel > 0, as.character(cumul_novel), "")

# Separate size scales per metric so ancient OG bubbles are large and prominent,
# independent of novel/cumulative counts
max_anc <- max(c(ancient_abs[!is_tip_v], 1L))
max_nov <- max(c(novel_mrca,             1L))   # includes tips (non-zero tips ok)
max_cum <- max(c(cumul_novel[!is_tip_v], 1L))

INT_MAX_ANC <- 14    # mm — ancient retained OGs (most prominent)
INT_MAX_NOV <- 10    # mm — novel OGs first gained at this node
INT_MAX_CUM <-  8    # mm — cumulative novel (root→node)
INT_Y       <- 0.82  # data-unit offset between bubble rows

bub_int_anc <- ifelse(!is_tip_v & ancient_abs > 0,
                      pmax(sqrt(ancient_abs / max_anc) * INT_MAX_ANC, 0.5), 0)
bub_int_nov <- ifelse(!is_tip_v & novel_mrca > 0,
                      pmax(sqrt(novel_mrca  / max_nov) * INT_MAX_NOV, 0.5), 0)
bub_int_cum <- ifelse(!is_tip_v & cumul_novel > 0,
                      pmax(sqrt(cumul_novel / max_cum) * INT_MAX_CUM, 0.5), 0)

lbl_hjust <- ifelse(is_tip_v, -0.35, -0.15)

node_df <- data.frame(
  node           = seq_len(n_all),
  is_tip         = is_tip_v,
  clade_code     = clade_code_v,
  clade_name     = CLADE_NAMES[clade_code_v],
  fill_color     = bub_fill,
  bub_size       = bub_size,
  bub_int_anc    = bub_int_anc,
  bub_int_nov    = bub_int_nov,
  bub_int_cum    = bub_int_cum,
  lbl_total      = lbl_total,
  lbl_int_anc    = lbl_int_anc,
  lbl_int_novel  = lbl_int_novel,
  lbl_int_cumul  = lbl_int_cumul,
  lbl_hjust      = lbl_hjust,
  stringsAsFactors = FALSE
)

# Suppress ancient_abs at the new root (HO+fungi ancestor): the ancestral-OG
# tracking is only meaningful one node *after* the HO split (i.e., within the
# fungal subtree). The root value equals the fungi-MRCA value anyway.
root_idx <- which(node_df$node == (n_tip + 1L))
node_df$bub_int_anc[root_idx] <- 0
node_df$lbl_int_anc[root_idx] <- ""

# ── Plot ───────────────────────────────────────────────────────────────────────
p <- ggtree(tr_c, layout="rectangular", branch.length="none", ladderize=FALSE,
            linewidth=1.1, color="grey20") %<+% node_df +
  # Tip bubbles (large, coloured by clade)
  geom_point2(aes(size=ifelse(is_tip, bub_size, 0), fill=fill_color),
              shape=21, color="grey50", stroke=0.5, alpha=0.60) +
  scale_fill_identity() + scale_size_identity() +
  # Total distinct OGs at tips (black, above bubble)
  geom_text2(aes(label=lbl_total, hjust=lbl_hjust),
             vjust=-0.55, size=2.5, color="grey15", fontface="bold") +
  # Internal — blue bubble: novel OGs first gained here (above, +INT_Y)
  geom_point2(aes(size=bub_int_nov), fill="#1a5fa8",
              shape=21, color="white", stroke=0.3, alpha=0.80,
              position=position_nudge(y=INT_Y)) +
  geom_text2(aes(label=lbl_int_novel), hjust=0.5,
             vjust=-1.1, size=1.9, color="#1a5fa8", fontface="bold",
             position=position_nudge(y=INT_Y)) +
  # Internal — red bubble: ancestral HO OGs retained (centre)
  geom_point2(aes(size=bub_int_anc), fill="#c0392b",
              shape=21, color="white", stroke=0.3, alpha=0.80) +
  geom_text2(aes(label=lbl_int_anc), hjust=0.5,
             vjust=-1.1, size=1.9, color="#c0392b", fontface="bold") +
  # Internal — green bubble: cumulative novel OGs root→node (below, -INT_Y)
  geom_point2(aes(size=bub_int_cum), fill="#27ae60",
              shape=21, color="white", stroke=0.3, alpha=0.80,
              position=position_nudge(y=-INT_Y)) +
  geom_text2(aes(label=lbl_int_cumul), hjust=0.5,
             vjust=-1.1, size=1.9, color="#27ae60", fontface="bold",
             position=position_nudge(y=-INT_Y)) +
  # Clade name labels — coloured by clade
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
  labs(title="Ancient co-opted OG retention and Fungal novelty across the tree of life",
       subtitle=paste0(
         "Tips: total distinct co-opted OG families observed in that clade (black number, bubble size \u221d total).\n",
         "Internal nodes \u2014 red: ancestral HO OGs retained (decreasing) / ",
         "green: cumulative novel OGs root\u2192node (increasing) / ",
         "blue: novel OGs first gained at this exact node (+N). Bubble size \u221d sqrt(count)."))

# ── Save ──────────────────────────────────────────────────────────────────────
for (ext in c("pdf", "png")) {
  out <- file.path(OUT_DIR, paste0("figA_ancient_novel.", ext))
  ggsave(out, p, width=12, height=10, dpi=150)
  cat("[OK] Saved", out, "\n")
}

# TSV summary
node_df$node_label <- ifelse(node_df$is_tip,
  tr_c$tip.label[node_df$node], paste0("inode_", node_df$node))
node_df$ancient_abs  <- ancient_abs
node_df$ancient_loss <- ancient_loss
node_df$novel_mrca   <- novel_mrca
write.table(
  node_df[, c("node_label","is_tip","ancient_abs","ancient_loss","novel_mrca")],
  file.path(OUT_DIR, "figA_ancient_novel_counts.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE)
cat("[OK] Saved figA_ancient_novel_counts.tsv\n")
