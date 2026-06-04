#!/usr/bin/env Rscript
# =============================================================================
# plot_circular_tree_htt.R
#
# Extends circular_tree_v4.5 with HTT chord overlay:
#   - Bezier arcs between clade midpoints, width ∝ log(n HTT candidates)
#   - Color scheme: warm (#b85c00) semi-transparent chords
#   - Only pairs with n >= 8 HTT candidates shown
#
# Output:
#   /legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figures/
#     circular_tree_v4.5_htt.pdf
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
  library(purrr)
  library(tibble)
})

# ─── Paths ────────────────────────────────────────────────────────────────────
TREE_FILE <- "/legserv/Temp/Carlos/EDF-TE/results/SpeciesTree/final_LG_G_rotated.nwk"
DATA_FILE <- "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/figure_species_data.tsv"
OUT_PDF   <- "/legserv/Temp/Carlos/EDF-TE/results/HTT/figures/circular_tree_v4.5_htt.pdf"
dir.create(dirname(OUT_PDF), showWarnings = FALSE, recursive = TRUE)

# ─── Load ─────────────────────────────────────────────────────────────────────
options(expressions = 100000L)
tree <- read.tree(TREE_FILE)
dat  <- read.delim(DATA_FILE, stringsAsFactors = FALSE) %>%
  rename(label = tree_label) %>%
  mutate(label = sub("\\.aa$", "", label)) %>%
  filter(!is.na(label)) %>%
  select(label, everything())

# ─── Same rotations as v4.5 ───────────────────────────────────────────────────
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
if ("CR_Rozal1" %in% tree$tip.label && "CR_Parsa1" %in% tree$tip.label) {
  rozal_idx <- which(tree$tip.label == "CR_Rozal1")
  rozal_bl  <- tree$edge.length[tree$edge[, 2] == rozal_idx]
  tree <- ape::drop.tip(tree, "CR_Rozal1")
  where_parsa <- which(tree$tip.label == "CR_Parsa1")
  tree <- phytools::bind.tip(tree, tip.label = "CR_Rozal1",
                             where = where_parsa, position = 0,
                             edge.length = rozal_bl)
  cat("[INFO] Forced CR_Rozal1 sister to CR_Parsa1\n")
}
cat("[INFO] Tree tips:", length(tree$tip.label), "\n")

# ─── Clade colours ────────────────────────────────────────────────────────────
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

# (Data-track helpers removed — clean layout with arcs only)

# ─── Base tree ────────────────────────────────────────────────────────────────
p <- ggtree(tree, layout = "fan", open.angle = 35, ladderize = FALSE,
            linewidth = 0.25, color = "grey25") %<+% dat +
  geom_tippoint(aes(color = clade), size = 1.5, alpha = 0.90) +
  scale_color_manual(values = clade_cols, name = NULL,
    guide = guide_legend(ncol = 3, order = 1,
                         override.aes = list(size = 3.0), title = NULL)) +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 5.5),
        legend.key.size  = unit(0.30, "cm"),
        plot.margin      = margin(25, 25, 25, 25),
        panel.border     = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA))

# ─── ggtree fan coordinate system ────────────────────────────────────────────
# In ggtree layout="fan", p$data stores PRE-TRANSFORMATION values:
#   p$data$x = depth (branch length from root) — polar radius at render time
#   p$data$y = tip index (1 to ntips)          — mapped to angle at render time
#   p$data$angle = visual angle in degrees (for label rotation)
# ggtree applies the fan transformation (x*cos(angle), x*sin(angle)) to ALL
# layers including external geom_* at render time via its coord class.
# Therefore arc coordinates MUST use (depth, tipindex) space, just like p$data.
# arc at outer ring = (arc_r, y_idx),  arc dipping toward center = (ctrl_r, y_mid)

all_tips   <- p$data %>% filter(isTip)
angle_min  <- min(all_tips$angle)   # degrees, first tip in fan
angle_max  <- max(all_tips$angle)   # degrees, last tip in fan
ntips_all  <- nrow(all_tips)
angle_to_yidx <- function(ang_deg) {
  1 + (ang_deg - angle_min) / (angle_max - angle_min) * (ntips_all - 1)
}
arc_r  <- max(all_tips$x) * 1.08   # arc endpoints just OUTSIDE outermost tree tips
ctrl_r <- arc_r * 0.12             # control point near centre → classic deep chord arcs

# ─── Extract HTT clade tip coordinates ───────────────────────────────────────
tip_coords <- p$data %>%
  filter(isTip) %>%
  mutate(pfx = case_when(
    grepl("^HO_", label) ~ "HO",
    grepl("^AP_", label) ~ "AP",
    grepl("^NU_", label) ~ "NU",
    grepl("^CH_", label) ~ "CH",   # includes NE_ Neocallimastigomycota
    grepl("^BL_", label) ~ "BL",
    grepl("^ZO_", label) ~ "ZO",
    grepl("^EN_", label) ~ "EN",
    grepl("^BA_", label) ~ "BA",   # Basidiobolales
    grepl("^KX_", label) ~ "KX",
    grepl("^GL_", label) ~ "GL",
    grepl("^MO_", label) ~ "MO",
    grepl("^MU_", label) ~ "MU",
    grepl("^BS_", label) ~ "BS",   # Basidiomycota
    grepl("^AS_", label) ~ "AS",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(pfx))

clade_ang <- tip_coords %>%
  group_by(pfx) %>%
  summarize(mean_angle_deg = mean(angle), .groups = "drop") %>%
  mutate(
    ang_rad = mean_angle_deg * pi / 180,
    y_idx   = angle_to_yidx(mean_angle_deg)
  )

cat("[HTT] Clades found:", paste(sort(clade_ang$pfx), collapse = ", "), "\n")
cat("[HTT] arc_r:", round(arc_r, 3), "| ctrl_r:", round(ctrl_r, 3), "\n")

# ─── Clade label annotations (placed just outside the arc ring) ───────────────
# Label only the clades involved in HTT transfers (so the arcs are interpretable)
htt_clades <- c("HO", "CH", "BL", "BA", "MO", "MU", "BS", "AS")
clade_label_names <- c(
  HO = "Holozoa", CH = "Chytridiomycota", BL = "Blastoclad.",
  BA = "Basidiobolales", MO = "Mortierellomycota", MU = "Mucoromycota",
  BS = "Basidiomycota",  AS = "Ascomycota"
)
clade_labels_df <- clade_ang %>%
  filter(pfx %in% htt_clades) %>%
  mutate(
    label_text = clade_label_names[pfx],
    # place text at 1.22 x arc_r (just beyond arc ring)
    lx = arc_r * 1.22,
    ly = y_idx,
    # rotate to follow the fan: angle + 90 for readability toward outside
    rot_angle = mean_angle_deg
  )

# ─── Bezier arc builder in ggtree (depth, tipindex) coordinate space ───────────
# Endpoints: (arc_r, y_idx1) and (arc_r, y_idx2)  — just outside tree tips.
# Control point: (ctrl_r, y_mid)  — near centre → classic deep chord arc.
# ggtree renders each (depth, tipindex) point as:
#   screen = (depth · cos(θ(tipindex)), depth · sin(θ(tipindex)))
make_arc <- function(pfx1, pfx2, n_htt, in_romeijn, n_pts = 150) {
  y1 <- clade_ang$y_idx[clade_ang$pfx == pfx1]
  y2 <- clade_ang$y_idx[clade_ang$pfx == pfx2]
  if (length(y1) == 0 || length(y2) == 0) return(NULL)

  y_mid <- (y1 + y2) / 2   # direct midpoint (no wrap-around: all pairs within 180°)

  t <- seq(0, 1, length.out = n_pts)
  tibble(
    x         = (1-t)^2 * arc_r + 2*t*(1-t) * ctrl_r + t^2 * arc_r,
    y         = (1-t)^2 * y1    + 2*t*(1-t) * y_mid   + t^2 * y2,
    pair      = paste0(pfx1, "-", pfx2),
    n         = n_htt,
    arc_color = ifelse(in_romeijn, "#9b2000", "#1a3a8a"),
    lwd       = 0.5 + 2.5 * log(n_htt) / log(33)
  )
}

# HTT pairs (n >= 8 candidates, from our analysis)
# concordant_romeijn: TRUE if both lineages covered by Romeijn et al. 2026
#   (Mucoromycotina + Chytrids + Blastocladiomycota + Basidio + Asco)
# FALSE for CH-HO (Holozoa is non-fungal, not in Romeijn)
htt_pairs <- tribble(
  ~c1,  ~c2,   ~n,  ~in_romeijn,
  "MO", "MU",  33,  TRUE,
  "CH", "MU",  19,  TRUE,
  "BS", "MU",  19,  TRUE,
  "BL", "MU",  15,  TRUE,
  "CH", "HO",  14,  FALSE,   # HO = Holozoa (non-fungal, not in Romeijn)
  "AS", "MU",  14,  TRUE,
  "BS", "MO",  12,  TRUE,
  "AS", "CH",  10,  TRUE,
  "AS", "MO",   8,  TRUE,
  "BA", "MO",   8,  TRUE
)

arc_df <- bind_rows(
  lapply(seq_len(nrow(htt_pairs)), function(i)
    make_arc(htt_pairs$c1[i], htt_pairs$c2[i], htt_pairs$n[i], htt_pairs$in_romeijn[i]))
)

cat("[HTT] Arc rows:", nrow(arc_df), "for", n_distinct(arc_df$pair), "pairs\n")

cat("[HTT] Arc widths (n=8/19/33):",
    round(0.5 + 2.5*log(c(8,19,33))/log(33), 2), "\n")

# ─── Add clade labels just outside the arc ring ────────────────────────────
p <- p +
  geom_text(
    data        = clade_labels_df,
    mapping     = aes(x = lx, y = ly, label = label_text, angle = rot_angle),
    size        = 2.8,
    hjust       = 0.5,
    color       = "grey20",
    fontface    = "bold",
    inherit.aes = FALSE
  )

# ─── Draw HTT chord arcs: concordant (red) first, novel (blue) on top ────────
arc_concordant <- arc_df %>% filter(arc_color == "#9b2000")
arc_novel      <- arc_df %>% filter(arc_color == "#1a3a8a")
p <- p +
  geom_path(
    data        = arc_concordant,
    mapping     = aes(x = x, y = y, group = pair, linewidth = I(lwd)),
    color       = "#9b2000",
    alpha       = 0.65,
    lineend     = "round",
    inherit.aes = FALSE
  ) +
  geom_path(
    data        = arc_novel,
    mapping     = aes(x = x, y = y, group = pair, linewidth = I(lwd)),
    color       = "#1a3a8a",
    alpha       = 0.80,
    lineend     = "round",
    inherit.aes = FALSE
  ) +
  # Legend text near centre of tree (ctrl_r zone)
  annotate("text",
           x = ctrl_r * 0.6,
           y = mean(range(clade_ang$y_idx)),
           label = paste0(
             "HTT arcs (n >= 8 candidates)\n",
             "Red = concordant with Romeijn et al. 2026\n",
             "Blue = novel (CH-HO, cross-kingdom)\n",
             "Width proportional to log(n)\n",
             "MO-MU:33  CH/BS-MU:19  BL-MU:15\n",
             "CH-HO:14  AS-MU:14  BS-MO:12\n",
             "AS-CH:10  AS-MO:8  BA-MO:8"
           ),
           size = 2.5, hjust = 0.5, vjust = 0.5, lineheight = 1.4,
           color = "grey30")

# ─── Save ─────────────────────────────────────────────────────────────────────
cat("[PLOT] Saving:", OUT_PDF, "\n")
ggsave(OUT_PDF, plot = p, width = 28, height = 30, units = "cm",
       device = "pdf", limitsize = FALSE)
cat("[OK] Done!\n")
