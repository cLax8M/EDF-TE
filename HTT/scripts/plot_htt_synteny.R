#!/usr/bin/env Rscript
# HTT LINE locus synteny snapshot
# HO_Monbrevi (NW_001865040.1) vs MO_MoelGBA34 (scaffold_76)
# ZYMTRI_2023_FAMILY_141, 574 bp, 80.37% identity

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

OUT <- "/legserv/Temp/Carlos/EDF-TE/results/HTT/figures/htt_line_synteny.pdf"

WIN_LEN <- 20575L  # bp, same for both windows

# ---- TE colour palette -------------------------------------------------------
te_colours <- c(
  "Simple_repeat"   = "#CCCCCC",
  "Low_complexity"  = "#AAAAAA",
  "LINE_focal"      = "#E84B37",   # ZYMTRI_2023_FAMILY_141 (the HTT element)
  "LINE_other"      = "#F4A460",   # other LINE families
  "Unknown"         = "#98C379",
  "LTR"             = "#C678DD",
  "DNA"             = "#56B6C2"
)

# helper: map raw type string → colour key
te_col_key <- function(type, focal) {
  ifelse(focal, "LINE_focal",
  ifelse(type %in% c("Simple_repeat"),         "Simple_repeat",
  ifelse(type %in% c("Low_complexity"),        "Low_complexity",
  ifelse(grepl("^LINE", type),                 "LINE_other",
  ifelse(grepl("^Unknown", type),              "Unknown",
  ifelse(grepl("^LTR",  type),                 "LTR",
  ifelse(grepl("^DNA",  type),                 "DNA",
                                               "Unknown")))))))
}

# ---- HO_Monbrevi TEs ---------------------------------------------------------
# Window: NW_001865040.1:1290160-1310734 (coords relative to window start)
ho_tes <- tibble::tribble(
  ~start, ~end,  ~type,            ~name,            ~focal,
  560,    594,   "Simple_repeat",  "(CTCACA)N",       FALSE,
  666,    763,   "Simple_repeat",  "(TC)N",           FALSE,
  764,    783,   "Simple_repeat",  "(CTCTTT)N",       FALSE,
  846,    864,   "Simple_repeat",  "(GCT)N",          FALSE,
 3165,   3207,   "Simple_repeat",  "(GTGC)N",         FALSE,
 7225,   7243,   "Simple_repeat",  "(G)N",            FALSE,
10001,  10574,   "LINE",           "ZYMTRI_141",      TRUE,
11888,  11906,   "Simple_repeat",  "(CT)N",           FALSE,
12172,  12188,   "Simple_repeat",  "(A)N",            FALSE,
14639,  15972,   "LINE",           "ZYMTRI_532",      FALSE,
15758,  15829,   "Unknown",        "ZYMTRI_407",      FALSE,
16753,  16800,   "Simple_repeat",  "(AC)N",           FALSE,
17807,  17846,   "Simple_repeat",  "(TCT)N",          FALSE,
18803,  18824,   "Simple_repeat",  "(GCG)N",          FALSE
) %>% mutate(col_key = te_col_key(type, focal))

# ---- MO_MoelGBA34 TEs --------------------------------------------------------
# Window: scaffold_76:25087-45661 (coords relative to window start)
mo_tes <- tibble::tribble(
  ~start, ~end,   ~type,              ~name,              ~focal,
     148,    197,  "Low_complexity",  "GA-RICH",           FALSE,
    1362,   1391,  "Simple_repeat",   "(TTGTTGT)N",        FALSE,
    1730,   1812,  "Simple_repeat",   "(CT)N",             FALSE,
    2961,   2982,  "Simple_repeat",   "(AGAG)N",           FALSE,
    4306,   4364,  "Simple_repeat",   "(TCT)N",            FALSE,
    4757,   4800,  "Simple_repeat",   "(ACA)N",            FALSE,
    4907,   4940,  "Simple_repeat",   "(GGT)N",            FALSE,
    5002,   5170,  "Unknown",         "ACOL_JOCKEY_ELE1",  FALSE,
    5171,   5266,  "LINE/I-Jockey",   "AG-JOCK-1",         FALSE,
    5267,   5450,  "LINE/R1",         "RTAG4",             FALSE,
    6843,   6862,  "Simple_repeat",   "(AGC)N",            FALSE,
    7989,   8072,  "LINE/R1",         "RTAG4",             FALSE,
    8091,   8134,  "Simple_repeat",   "(CAA)N",            FALSE,
    8231,   8276,  "LINE/R1",         "RTAG4",             FALSE,
    9586,   9690,  "Low_complexity",  "A-RICH",            FALSE,
   10001,  10574,  "LINE",            "ZYMTRI_141",        TRUE,
   10958,  11038,  "Unknown",         "ZYMTRI_1492",       FALSE,
   12007,  12040,  "Simple_repeat",   "(TGTCGGT)N",        FALSE,
   12041,  12056,  "Simple_repeat",   "(GTT)N",            FALSE,
   12098,  12114,  "Simple_repeat",   "(A)N",              FALSE,
   13658,  13715,  "Simple_repeat",   "(CAACAG)N",         FALSE,
   14116,  14204,  "LTR/ERVK",        "ERV2-2B_BT-I",      FALSE,
   16130,  16271,  "DNA/hAT-Charlie", "HAT-N4_ALY",        FALSE,
   16448,  16493,  "Simple_repeat",   "(GTG)N",            FALSE,
   18039,  18071,  "Low_complexity",  "A-RICH",            FALSE,
   18454,  18524,  "Simple_repeat",   "(TGGAGG)N",         FALSE,
   18525,  18617,  "LINE/I-Jockey",   "JOCKEY-6_DK",       FALSE,
   18799,  18922,  "DNA/Maverick",    "POLINTON-2A_NV",    FALSE,
   19377,  19423,  "Simple_repeat",   "(TGT)N",            FALSE,
   20164,  20220,  "Low_complexity",  "A-RICH",            FALSE,
   20310,  20340,  "Low_complexity",  "A-RICH",            FALSE
) %>% mutate(col_key = te_col_key(type, focal))

# ---- nucmer alignment --------------------------------------------------------
aln_q_start  <- 9991L; aln_q_end  <- 10746L
aln_s_start  <- 9991L; aln_s_end  <- 10746L
aln_pct_id   <- 80.37

# ---- plot helpers ------------------------------------------------------------
TRACK_Y      <- 0.5   # centre of genomic track
TRACK_H      <- 0.35  # half-height of track box
TE_H         <- 0.28  # half-height of TE boxes (taller than track)

build_track_plot <- function(tes, win_len, species_label, chrom_label,
                              flip = FALSE, show_x = TRUE) {
  # backbone
  p <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = win_len, ymin = TRACK_Y - TRACK_H,
                  ymax = TRACK_Y + TRACK_H),
              fill = "#E8E8E8", colour = NA)

  # TE boxes
  p <- p + geom_rect(
    data = tes,
    aes(xmin = start, xmax = end,
        ymin = TRACK_Y - TE_H, ymax = TRACK_Y + TE_H,
        fill = col_key),
    colour = NA
  ) +
  # focal outline
  geom_rect(
    data = filter(tes, focal),
    aes(xmin = start - 50, xmax = end + 50,
        ymin = TRACK_Y - TE_H - 0.04, ymax = TRACK_Y + TE_H + 0.04),
    fill = NA, colour = "#E84B37", linewidth = 0.6
  ) +
  scale_fill_manual(
    values = te_colours,
    labels = c(
      "LINE_focal"     = "LINE (ZYMTRI_141, HTT)",
      "LINE_other"     = "LINE (other)",
      "Simple_repeat"  = "Simple repeat",
      "Low_complexity" = "Low complexity",
      "Unknown"        = "Unknown",
      "LTR"            = "LTR retrotransposon",
      "DNA"            = "DNA transposon"
    ),
    name = "TE class"
  ) +
  coord_cartesian(xlim = c(-3000, win_len + 200), ylim = c(0, 1), clip = "off") +
  theme_void() +
  theme(
    legend.position  = "none",
    plot.margin      = margin(0, 0, 0, 0)
  )

  # strand arrows showing direction
  focal_row <- filter(tes, focal)
  mid <- (focal_row$start + focal_row$end) / 2
  p <- p + annotate("segment",
    x = mid - 180, xend = mid + 180,
    y = TRACK_Y, yend = TRACK_Y,
    colour = "white", linewidth = 0.5,
    arrow = arrow(length = unit(0.06, "inches"), type = "open")
  )

  # species + chrom label
  p <- p + annotate("text",
    x = -400, y = TRACK_Y,
    label = paste0(species_label, "\n", chrom_label),
    hjust = 1, vjust = 0.5, size = 2.5, fontface = "bold", colour = "#333333"
  )

  # x-axis (kb)
  if (show_x) {
    kb_breaks <- seq(0, win_len, by = 5000)
    kb_labels <- paste0(kb_breaks / 1000, " kb")
    p <- p + annotate("segment",
      x = kb_breaks, xend = kb_breaks,
      y = TRACK_Y - TRACK_H - 0.02, yend = TRACK_Y - TRACK_H - 0.06,
      colour = "#555555", linewidth = 0.25
    ) + annotate("text",
      x = kb_breaks, y = TRACK_Y - TRACK_H - 0.12,
      label = kb_labels, size = 1.8, colour = "#555555"
    ) + annotate("segment",
      x = 0, xend = win_len,
      y = TRACK_Y - TRACK_H - 0.02, yend = TRACK_Y - TRACK_H - 0.02,
      colour = "#555555", linewidth = 0.25
    )
  }

  p
}

# ---- ribbon panel ------------------------------------------------------------
build_ribbon_plot <- function(q_s, q_e, s_s, s_e, pct_id, win_len) {
  # The ribbon x coords match the window coordinates of each track
  xs <- c(q_s, q_e, s_e, s_s, q_s)
  ys <- c(1,   1,   0,   0,   1)   # top = HO track bottom, bottom = MO track top

  lbl_x <- (q_s + q_e) / 2
  lbl_y <- 0.5

  ggplot() +
    geom_polygon(aes(x = xs, y = ys),
                 fill = "#E84B3750", colour = "#E84B3790", linewidth = 0.3) +
    annotate("text",
      x = lbl_x, y = lbl_y,
      label = paste0(round(pct_id, 1), "% id"),
      size = 2.6, colour = "#C0392B", fontface = "italic"
    ) +
coord_cartesian(xlim = c(-3000, win_len + 200), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

# ---- legend panel ------------------------------------------------------------
build_legend <- function() {
  # Two-row legend; names(te_colours): Simple_repeat, Low_complexity, LINE_focal, LINE_other, Unknown, LTR, DNA
  leg_df <- data.frame(
    col_key = names(te_colours),
    label   = c(
      "Simple repeat",          "Low complexity",
      "LINE (ZYMTRI_141, HTT)", "LINE (other families)",
      "Unknown / unclassified", "LTR retrotransposon",
      "DNA transposon"
    ),
    x = c(0.0, 3.0, 6.0,  0.0, 3.0, 6.0, 9.0),
    y = c(0.75, 0.75, 0.75, 0.28, 0.28, 0.28, 0.28)
  )

  ggplot(leg_df) +
    geom_rect(aes(xmin = x, xmax = x + 0.28,
                  ymin = y - 0.18, ymax = y + 0.18, fill = col_key),
              colour = "#888888", linewidth = 0.2) +
    geom_text(aes(x = x + 0.40, y = y, label = label),
              hjust = 0, size = 2.3) +
    scale_fill_manual(values = te_colours) +
    coord_cartesian(xlim = c(-0.2, 12.5), ylim = c(0, 1)) +
    theme_void() +
    theme(legend.position = "none")
}

# ---- assemble ----------------------------------------------------------------
p_ho  <- build_track_plot(ho_tes, WIN_LEN,
                          "Monosiga brevicollis\n(Holozoa)",
                          "NW_001865040.1",
                          show_x = FALSE)

p_rib <- build_ribbon_plot(aln_q_start, aln_q_end,
                           aln_s_start, aln_s_end,
                           aln_pct_id, WIN_LEN)

p_mo  <- build_track_plot(mo_tes, WIN_LEN,
                          "Mortierella elongata GBA34\n(Mortierellomycotina)",
                          "scaffold_76",
                          show_x = TRUE)

p_leg <- build_legend()

# Combine with patchwork
layout <- "
AA
BB
CC
DD
"
fig <- p_ho / p_rib / p_mo / p_leg +
  plot_layout(heights = c(1.2, 0.6, 1.2, 0.8)) +
  plot_annotation(
    title = "HTT LINE locus: Monosiga brevicollis ↔ Mortierella elongata GBA34",
    subtitle = "ZYMTRI_2023_FAMILY_141 | 574 bp | 80.4% nucleotide identity | no flanking synteny",
    theme = theme(
      plot.title    = element_text(size = 9, face = "bold"),
      plot.subtitle = element_text(size = 7.5, colour = "#555555"),
      plot.margin   = margin(6, 12, 4, 70)
    )
  )

ggsave(OUT, fig, width = 10, height = 5.5, device = "pdf")
cat("[HTT-synteny] Figure saved to:", OUT, "\n")
