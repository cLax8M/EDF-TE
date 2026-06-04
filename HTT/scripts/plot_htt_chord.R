#!/usr/bin/env Rscript
# =============================================================================
# plot_htt_chord.R  — v2
#
# HTT chord diagram (Circos layout):
#   - Nearly full-circle fan tree (open.angle = 2) as a thin OUTER RING
#   - 223 individual candidate arcs drawn (one per species pair)
#   - Within-fungi arcs = dark red; HO-involving (cross-kingdom) = blue
# =============================================================================

.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tibble)
})

TREE_NWK <- "/legserv/Temp/Carlos/EDF-TE/results/SpeciesTree/final_LG_G_rotated.nwk"
HTT_TSV  <- "/legserv/Temp/Carlos/EDF-TE/results/HTT/06_htt_candidates/htt_candidates/htt_by_clade_pair.tsv"
OUT_PDF  <- "/legserv/Temp/Carlos/EDF-TE/results/HTT/figures/htt_chord_v2.pdf"
dir.create(dirname(OUT_PDF), showWarnings = FALSE, recursive = TRUE)

# ── Read ALL pairwise transfers ───────────────────────────────────────────────
htt_raw <- read.table(HTT_TSV, header = TRUE, sep = "\t")
# Normalise pair order (alphabetical), merge symmetric duplicates
htt_all <- htt_raw %>%
  mutate(
    cladeA = pmin(qclade, sclade),
    cladeB = pmax(qclade, sclade)
  ) %>%
  group_by(cladeA, cladeB) %>%
  summarise(n = sum(htt_count), .groups = "drop") %>%
  arrange(desc(n))
cat("[INFO] Clade pairs:", nrow(htt_all), " | Total candidates:", sum(htt_all$n), "\n")

# ── Clade colors ──────────────────────────────────────────────────────────────
CLADE_COL <- c(
  HO = "#7B241C", CH = "#E74C3C", BL = "#E67E22", AS = "#1E8449",
  BA = "#2471A3", BS = "#1A5276", MO = "#7D3C98", MU = "#BB8FCE",
  GL = "#D4AC0D", KX = "#F39C12", EN = "#D35400", ZO = "#85C1E9",
  NU = "#1C2833", AP = "#27AE60",
  MI = "#5DADE2", CR = "#76D7C4", OL = "#F7DC6F"  # Microsporidia, Rozellomycota, Olpidium
)
CLADE_NAMES <- c(
  HO = "Holozoa",           CH = "Chytridiomycota",    BL = "Blastocladiomycota",
  BA = "Basidiobolales",    MO = "Mortierellomycota",   MU = "Mucoromycotina",
  BS = "Basidiomycota",     AS = "Ascomycota",
  GL = "Glomeromycotina",   KX = "Kickxellomycotina",   EN = "Entomophthoromycota",
  ZO = "Zoopagomycota",     NU = "Nuclearida",          AP = "Aphelida",
  MI = "Microsporidia",     CR = "Rozellomycota",       OL = "Olpidiomycota"
)
# NOTE: LABEL_CLADES will be set to all tree clades after clade_sec is built

# ── Load tree + same rotations as main script ─────────────────────────────────
tree <- read.tree(TREE_NWK)

NE_SPECIES <- c("CH_Anasp1","CH_Neocon1","CH_Neolan1","CH_Orpsp1",
                "CH_Pecora1","CH_Caecom1","CH_Pirfi3")
get_pfx <- function(sp) {
  if (sp %in% NE_SPECIES) return("CH")
  p <- substr(sp, 1, 2)
  if (p == "NE") return("CH")
  p
}
get_rep <- function(p) { hits <- tree$tip.label[sapply(tree$tip.label, get_pfx) == p]; hits[1] }

mo <- get_rep("MO"); gl <- get_rep("GL"); mu <- get_rep("MU")
if (!anyNA(c(mo, gl, mu)))
  tree <- ape::rotate(tree, getMRCA(tree, c(mo, gl, mu)))
if (all(c("CR_Rozal1","CR_Parsa1") %in% tree$tip.label)) {
  idx <- which(tree$tip.label == "CR_Rozal1")
  bl  <- tree$edge.length[tree$edge[,2] == idx]
  tree <- ape::drop.tip(tree, "CR_Rozal1")
  tree <- phytools::bind.tip(tree, "CR_Rozal1",
    where = which(tree$tip.label == "CR_Parsa1"), position = 0, edge.length = bl)
}
n_tips <- length(tree$tip.label)
cat("[INFO] Tips:", n_tips, "\n")

# ── Fan layout → extract angles and depths for every node ─────────────────────
p_fan <- ggtree(tree, layout = "fan", open.angle = 2,
                ladderize = FALSE)
nd    <- p_fan$data           # x = depth (node rank), y = tipindex, angle = degrees

# Strict linear y→angle mapping from actual tip positions
tips_nd <- nd %>% filter(isTip)
y_lo <- min(tips_nd$y); y_hi <- max(tips_nd$y)
a_lo <- min(tips_nd$angle); a_hi <- max(tips_nd$angle)
y2deg <- function(y) a_lo + (y - y_lo) / (y_hi - y_lo) * (a_hi - a_lo)
y2rad <- function(y) y2deg(y) * pi / 180

# Node depth and y-value look-up vectors (indexed by ape node ID)
n_total <- n_tips + tree$Nnode
nd_by_node <- nd[match(seq_len(n_total), nd$node), ]  # row i = ape node i
nd_x <- nd_by_node$x                                  # depth
nd_y <- nd_by_node$y                                  # tipindex / angle proxy

# ── Rescale depths: outer tree with real branch lengths (non-ultrametric) ──────────
# Use max-tip-distance per node (bottom-up traversal) so ALL tips land at
# the same inner radius regardless of tip-to-root distance variation.
OUTER_FRAC <- 0.12
nd_x_max   <- max(nd_x, na.rm = TRUE)   # ggtree coord scale (largest root-to-node)

# max_tip_dist[i] = max cumulative branch length from node i to any descendant leaf
edges        <- tree$edge          # columns: [parent, child] — needed here and below
max_tip_dist <- rep(0.0, n_total)
el <- tree$edge.length
for (k in nrow(edges):1) {
  pa <- edges[k, 1]; ch <- edges[k, 2]
  cand <- max_tip_dist[ch] + el[k]
  if (cand > max_tip_dist[pa]) max_tip_dist[pa] <- cand
}
max_fwd <- max_tip_dist[n_tips + 1L]    # = max root-to-tip distance
# tips (max_tip_dist=0) → nd_x_max (inner circle); root → nd_x_max*(1+OUTER_FRAC)
nd_x_s <- nd_x_max * (1 + OUTER_FRAC * max_tip_dist / max_fwd)

# ── Compute tip y-range for every node (bottom-up on ape edge matrix) ─────────
# ape edges are in preorder → reversing gives children before parents
# (edges already defined above)
tip_min <- rep(Inf,  n_total)
tip_max <- rep(-Inf, n_total)
for (i in seq_len(n_tips)) { tip_min[i] <- nd_y[i]; tip_max[i] <- nd_y[i] }
for (i in nrow(edges):1) {
  p_nd <- edges[i,1]; c_nd <- edges[i,2]
  if (tip_min[c_nd] < tip_min[p_nd]) tip_min[p_nd] <- tip_min[c_nd]
  if (tip_max[c_nd] > tip_max[p_nd]) tip_max[p_nd] <- tip_max[c_nd]
}

# ── Build tree branches in Cartesian space ────────────────────────────────────
# Each edge: one radial segment at the child's angle (using scaled depth)
rad_list <- lapply(seq_len(nrow(edges)), function(i) {
  pn <- edges[i,1]; cn <- edges[i,2]
  a  <- y2rad(nd_y[cn])
  tibble(x1 = nd_x_s[pn]*cos(a), y1 = nd_x_s[pn]*sin(a),
         x2 = nd_x_s[cn]*cos(a), y2 = nd_x_s[cn]*sin(a))
})
radial_df <- bind_rows(rad_list)

# Each internal node: shoulder arc at scaled depth
int_nodes <- unique(edges[,1])
arc_list <- lapply(int_nodes, function(ni) {
  r  <- nd_x_s[ni]
  a1 <- y2rad(tip_min[ni]);  a2 <- y2rad(tip_max[ni])
  if (is.infinite(a1) || is.infinite(a2) || a1 == a2) return(NULL)
  n_pts <- max(3L, as.integer(abs(a2 - a1) * 60 / pi) + 2L)
  angs  <- seq(a1, a2, length.out = n_pts)
  tibble(xc = r * cos(angs), yc = r * sin(angs), node_id = ni)
})
arc_tree_df <- bind_rows(arc_list)
cat("[INFO] Radial segs:", nrow(radial_df),
    "| Shoulder arcs:", n_distinct(arc_tree_df$node_id), "\n")

# ── Tip → clade assignment ────────────────────────────────────────────────────
tips_nd <- tips_nd %>%
  mutate(
    pfx = sapply(label, function(lb) {
      p <- get_pfx(lb)
      if (p %in% names(CLADE_COL)) p else NA_character_
    }),
    xc  = nd_x_max * cos(angle * pi/180),   # all tips at same inner radius
    yc  = nd_x_max * sin(angle * pi/180)
  )

# ── Geometry parameters ───────────────────────────────────────────────────────
max_d   <- nd_x_max              # nominal radius (= tip radius = inner circle edge)
R1      <- max_d * 0.99          # inner sector ring edge (at tips)
R2      <- max_d * (1 + OUTER_FRAC + 0.02)  # outer sector ring edge (just past root)
R_chord <- max_d * 0.97          # chord endpoints (just inside tips)
R_lbl   <- max_d * (1 + OUTER_FRAC + 0.16)  # clade labels outside the tree
ctrl_r  <- max_d * 0.10          # Bezier control point (near center)
n_max   <- max(htt_all$n)     # for normalising arc width/alpha

# ── Scale bar (branch-length scale, substitutions per site) ──────────────────
bar_per_sub <- nd_x_max * OUTER_FRAC / max_fwd     # plot units per 1 sub/site
SCALE_SUB   <- max(0.5, round(0.08 * nd_x_max / bar_per_sub / 0.5) * 0.5)
bar_len     <- SCALE_SUB * bar_per_sub
sb_y        <- -(nd_x_max * (1 + OUTER_FRAC) + nd_x_max * 0.25)
sb_seg_df   <- tibble(x = -bar_len/2, xend = bar_len/2, y = sb_y, yend = sb_y)
sb_tick_df  <- bind_rows(
  tibble(x = -bar_len/2, xend = -bar_len/2,
         y = sb_y - nd_x_max*0.015, yend = sb_y + nd_x_max*0.015),
  tibble(x =  bar_len/2, xend =  bar_len/2,
         y = sb_y - nd_x_max*0.015, yend = sb_y + nd_x_max*0.015))
sb_lbl_df   <- tibble(x = 0, y = sb_y - nd_x_max*0.08,
                      label = sprintf("%.4g subs/site", SCALE_SUB))

# ── Clade sector positions ────────────────────────────────────────────────────
clade_sec <- tips_nd %>%
  filter(!is.na(pfx)) %>%
  group_by(pfx) %>%
  summarise(
    a_min  = min(angle), a_max = max(angle),
    a_mid  = mean(angle), n_tips = n(),
    .groups = "drop"
  )

# ── Clade ring polygons ───────────────────────────────────────────────────────
GAP <- 0.3   # half-degree gap between adjacent sectors
ring_df <- bind_rows(lapply(seq_len(nrow(clade_sec)), function(i) {
  r  <- clade_sec[i,]
  a1 <- (r$a_min + GAP) * pi/180
  a2 <- (r$a_max - GAP) * pi/180
  if (a1 >= a2) return(NULL)
  n  <- max(4L, as.integer(abs(a2 - a1) * 60 / pi) + 2L)
  ag <- seq(a1, a2, length.out = n)
  tibble(
    xc  = c(R1 * cos(ag), rev(R2 * cos(ag)), R1 * cos(a1)),
    yc  = c(R1 * sin(ag), rev(R2 * sin(ag)), R1 * sin(a1)),
    pfx = r$pfx
  )
}))



# ── Clade labels (ALL clades in tree ring) ───────────────────────────────────
# Label every clade that has a colored sector (including Microsporidia etc.)
LABEL_CLADES <- unique(clade_sec$pfx)
lbl_df <- clade_sec %>%
  filter(pfx %in% LABEL_CLADES) %>%
  mutate(
    # Small clades pre-staggered radially to reduce ggrepel displacement
    # (smaller stagger now: R_lbl is already well outside the outer tree ring)
    r_mult = case_when(n_tips <= 3  ~ 1.12,
                       n_tips <= 7  ~ 1.06,
                       TRUE         ~ 1.00),
    # Place label at R_lbl in the radial direction of the sector midpoint
    lx    = (R_lbl * r_mult) * cos(a_mid * pi/180),
    ly    = (R_lbl * r_mult) * sin(a_mid * pi/180),
    # Anchor point: outer ring midpoint (for ggrepel segment)
    ax    = R2    * cos(a_mid * pi/180),
    ay    = R2    * sin(a_mid * pi/180),
    txt   = ifelse(pfx %in% names(CLADE_NAMES), CLADE_NAMES[pfx], pfx),
    rot   = ifelse(a_mid > 90 & a_mid < 270, a_mid - 180, a_mid),
    hjust = ifelse(a_mid > 90 & a_mid < 270, 1, 0)
  )

# ── HTT individual chord Bezier arcs (one per candidate pair) ────────────────
CAND_TSV  <- sub("htt_by_clade_pair", "htt_candidates", HTT_TSV)
htt_cand  <- read.table(CAND_TSV, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)

angmid_deg <- function(a1, a2) {
  d <- ((a2 - a1 + 180) %% 360) - 180
  (a1 + d/2) %% 360
}

# Genome → tip angle lookup
tip_angle <- setNames(tips_nd$angle, tips_nd$label)
get_ang <- function(genome, clade) {
  if (genome %in% names(tip_angle)) return(tip_angle[[genome]])
  sec <- clade_sec[clade_sec$pfx == clade, ]
  if (nrow(sec) > 0) return(sec$a_mid) else return(NA_real_)
}

chord_raw <- lapply(seq_len(nrow(htt_cand)), function(i) {
  qg <- htt_cand$qgenome[i]; sg <- htt_cand$sgenome[i]
  qc <- htt_cand$qclade[i];  sc <- htt_cand$sclade[i]
  a1d <- get_ang(qg, qc); a2d <- get_ang(sg, sc)
  if (is.na(a1d) || is.na(a2d)) { cat("[WARN] genome not found:", qg, sg, "\n"); return(NULL) }
  amd <- angmid_deg(a1d, a2d)
  a1  <- a1d * pi/180; a2 <- a2d * pi/180; am <- amd * pi/180
  ex1 <- R_chord * cos(a1); ey1 <- R_chord * sin(a1)
  ex2 <- R_chord * cos(a2); ey2 <- R_chord * sin(a2)
  cx  <- ctrl_r  * cos(am); cy  <- ctrl_r  * sin(am)
  t   <- seq(0, 1, length.out = 80)
  col <- if ("HO" %in% c(qc, sc)) "#3B8FD4" else "#555555"
  tibble(
    xc        = (1-t)^2*ex1 + 2*t*(1-t)*cx + t^2*ex2,
    yc        = (1-t)^2*ey1 + 2*t*(1-t)*cy + t^2*ey2,
    pair      = htt_cand$pair_id[i],
    arc_color = col
  )
})
chord_df <- bind_rows(chord_raw)
chord_r  <- chord_df %>% filter(arc_color == "#555555")
chord_b  <- chord_df %>% filter(arc_color == "#3B8FD4")
cat("[HTT] Chord arcs drawn:", n_distinct(chord_df$pair), "\n")

# ── Build the plot (pure Cartesian — no coord_polar) ──────────────────────────
gg <- ggplot() +
  # ── HTT arcs drawn FIRST (bottom layer — masked by opaque sector ring)
  geom_path(data = chord_r,
            aes(x = xc, y = yc, group = pair),
            color = "#555555", linewidth = 0.16, alpha = 0.28, lineend = "round") +
  geom_path(data = chord_b,
            aes(x = xc, y = yc, group = pair),
            color = "#3B8FD4", linewidth = 0.55, alpha = 0.75, lineend = "round") +
  # ── Tree branches
  geom_path(data = arc_tree_df,
            aes(x = xc, y = yc, group = node_id),
            color = "grey20", linewidth = 0.28, lineend = "round") +
  geom_segment(data = radial_df,
               aes(x = x1, y = y1, xend = x2, yend = y2),
               color = "grey20", linewidth = 0.28, lineend = "round") +
  scale_color_manual(values = CLADE_COL, guide = "none") +
  # ── Tip points on top
  geom_point(data = tips_nd %>% filter(!is.na(pfx)),
             aes(x = xc, y = yc, color = pfx), size = 0.5) +
  # ── Clade labels: ggrepel to avoid overlap
  ggrepel::geom_text_repel(
            data = lbl_df,
            aes(x = lx, y = ly, label = txt, angle = rot),
            size = 2.5, fontface = "bold", color = "grey10",
            force = 1.5, force_pull = 1.2, max.overlaps = Inf,
            box.padding = 0.20, point.padding = 0.1,
            segment.color = "grey55", segment.size = 0.25, segment.linetype = 1,
            min.segment.length = 0.05, direction = "both", seed = 42) +
  # ── Scale bar
  geom_segment(data = sb_seg_df, aes(x=x, y=y, xend=xend, yend=yend),
               inherit.aes=FALSE, color="grey20", linewidth=0.5, lineend="butt") +
  geom_segment(data = sb_tick_df, aes(x=x, y=y, xend=xend, yend=yend),
               inherit.aes=FALSE, color="grey20", linewidth=0.5, lineend="butt") +
  geom_text(data = sb_lbl_df, aes(x=x, y=y, label=label),
            inherit.aes=FALSE, size=2.0, color="grey20", hjust=0.5, vjust=1) +
  coord_fixed(clip = "off") +
  expand_limits(x = c(-max_d * 2.40, max_d * 2.40),
                y = c(-max_d * 2.40, max_d * 2.40)) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5, unit = "mm"))

cat("[PLOT] Saving:", OUT_PDF, "\n")
ggsave(OUT_PDF, gg, width = 27, height = 27, units = "cm",
       device = "pdf", limitsize = FALSE)
cat("[OK] Done!\n")
