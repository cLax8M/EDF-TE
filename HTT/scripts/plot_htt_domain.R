#!/usr/bin/env Rscript
# HTT LINE domain architecture diagram
# ZYMTRI_2023_FAMILY_141 — 574 bp consensus (LINE element)
# Found in Monosiga brevicollis and Mortierella elongata GBA34

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

OUT <- "/legserv/Temp/Carlos/EDF-TE/results/HTT/figures/htt_line_domain.pdf"

# ---- Background on the element -----------------------------------------------
# ZYMTRI_2023_FAMILY_141 consensus = 574 bp
# Annotated as class: LINE (non-LTR retrotransposon) by Earl Grey / RepeatMasker
# Monosiga brevicollis: NW_001865040.1:1300161-1300734, kimura80 = 0.247 (strand -)
# M. elongata GBA34:   scaffold_76:35088-35661,         kimura80 = 0.175 (strand -)
#
# Domain analysis:
#   hmmscan (E <= 10) against:
#     - Kapitonov et al. LINE-specific HMM database  → NO HITS
#     - REXdb plant v4.0 + metazoa v3.1              → NO HITS
#   TEsorter classification (filtered_tes.tsv):   LINE (no specific subfamily)
#   RepeatMasker classification:                  LINE (no specific subfamily)
#
# Two reading frames have 0 premature stop codons over 574 bp (191 aa):
#   plus-strand frame 3: SNVHGHLLDASRVEGLNVAH... (no stop)
#   minus-strand frame 1: SATMPQEVLEVTEKFMRD... (no stop)  ← same orientation as annotation
# Conclusion: highly diverged LINE; protein sequence too diverged for REXdb/Kapitonov assignment.
#             Classification as LINE based on nucleotide similarity to library.
#             Element is transcribed from minus strand (annotation); intact putative ORF.

CON_LEN <- 574L   # consensus bp

# ---- Domains -----------------------------------------------------------------
# LINEs (non-LTR retrotransposons) canonical ORF2 structure:
#   5' UTR | [ORF1: gag-like] | [ORF2: APE - RT - RNaseH] | 3' UTR (poly-A tail)
# This 574 bp fragment covers only the RT-containing region (partial ORF2).
# No in-frame stop codons in two frames → intact translatable ORF.

domains <- tibble::tribble(
  ~start, ~end, ~label,             ~fill,     ~yc, ~h,
  1,      574,  "ORF (no stops)",   "#F0C08A",  0.5, 0.22,  # ORF bar
  1,      574,  "RT-related\ndomain\n(putative)",  "#E84B37",  0.5, 0.15   # domain
)

# ---- Per-species instances ---------------------------------------------------
instances <- tibble::tribble(
  ~species,                     ~clade,              ~kimura, ~start, ~end, ~y,
  "Monosiga brevicollis",       "Holozoa",            0.247,     1,    574, 1.8,
  "Mortierella elongata GBA34", "Mortierellomycotina",0.175,     1,    574, 1.2
)

# ---- Colour scale for kimura divergence --------------------------------------
kimura_pal <- function(k) {
  # low divergence = darker red, high = lighter orange
  r <- floor(210 + (1 - k/0.5) * 45)
  grDevices::rgb(
    pmin(255, pmax(0, 228 - k * 200)),
    pmin(255, pmax(0, 75  - k * 50)),
    55,
    maxColorValue = 255
  )
}

# ---- Element overview panel --------------------------------------------------
p_overview <- ggplot() +
  # grey backbone (full consensus)
  geom_rect(aes(xmin = 0, xmax = CON_LEN,
                ymin = 0.42, ymax = 0.58),
            fill = "#DDDDDD", colour = NA) +
  # ORF backing
  geom_rect(aes(xmin = 1, xmax = CON_LEN,
                ymin = 0.40, ymax = 0.60),
            fill = "#F4C97A", colour = "#C8993A", linewidth = 0.4) +
  # Domain bar (labeled as LINE class; no specific subfamily identified)
  geom_rect(aes(xmin = 1, xmax = CON_LEN,
                ymin = 0.44, ymax = 0.56),
            fill = "#E84B37", colour = "#C0392B", linewidth = 0.3) +
  # Domain label
  annotate("text", x = CON_LEN / 2, y = 0.5,
           label = "LINE (no specific subfamily by protein hmmscan)", colour = "white",
           size = 2.6, fontface = "bold") +
  # Element name label
  annotate("text", x = CON_LEN / 2, y = 0.74,
           label = "ZYMTRI_2023_FAMILY_141  |  LINE (non-LTR retrotransposon, unassigned lineage)",
           size = 3.0, fontface = "bold", colour = "#333333") +
  annotate("text", x = CON_LEN / 2, y = 0.25,
           label = "574 bp consensus | intact ORF (0 stops in minus-strand frame 1, 191 aa)\nhmmscan vs REXdb+Kapitonov LINE HMMs: no significant hit (highly diverged)",
           size = 2.3, colour = "#666666", lineheight = 0.9) +
  # ORF label
  annotate("text", x = -20, y = 0.5,
           label = "ORF", size = 2.6, hjust = 1, colour = "#5C3D0A") +
  # bp scale
  annotate("segment", x = 0, xend = CON_LEN, y = 0.12, yend = 0.12,
           colour = "#888888", linewidth = 0.3) +
  annotate("segment",
           x = c(0, 100, 200, 300, 400, 500, 574),
           xend = c(0, 100, 200, 300, 400, 500, 574),
           y = 0.12, yend = 0.08, colour = "#888888", linewidth = 0.25) +
  annotate("text",
           x = c(0, 100, 200, 300, 400, 500, 574),
           y = 0.04,
           label = c("0", "100", "200", "300", "400", "500", "574 bp"),
           size = 2.2, colour = "#666666") +
  coord_cartesian(xlim = c(-60, CON_LEN + 20), ylim = c(0, 1)) +
  theme_void() +
  theme(plot.margin = margin(4, 4, 0, 4))

# ---- Per-species panel -------------------------------------------------------
p_instances <- ggplot() +
  # HO instance
  geom_rect(data = filter(instances, clade == "Holozoa"),
    aes(xmin = start, xmax = end, ymin = y - 0.18, ymax = y + 0.18),
    fill = kimura_pal(0.247), colour = "#C0392B", linewidth = 0.4) +
  annotate("text", x = -20, y = 1.8,
    label = paste0("Monosiga brevicollis (m = 0.247)"),
    hjust = 1, size = 2.5, fontface = "italic", colour = "#333333") +
  annotate("text", x = CON_LEN + 10, y = 1.8,
    label = "Holozoa", hjust = 0, size = 2.3, colour = "#555555") +
  # MO instance
  geom_rect(data = filter(instances, clade == "Mortierellomycotina"),
    aes(xmin = start, xmax = end, ymin = y - 0.18, ymax = y + 0.18),
    fill = kimura_pal(0.175), colour = "#C0392B", linewidth = 0.4) +
  annotate("text", x = -20, y = 1.2,
    label = paste0("Mortierella elongata GBA34 (m = 0.175)"),
    hjust = 1, size = 2.5, fontface = "italic", colour = "#333333") +
  annotate("text", x = CON_LEN + 10, y = 1.2,
    label = "Mortierellomycotina", hjust = 0, size = 2.3, colour = "#555555") +
  # Kimura annotation
  annotate("segment", x = CON_LEN / 2, xend = CON_LEN / 2,
           y = 1.0, yend = 2.0, colour = "#AAAAAA", linewidth = 0.2,
           linetype = "dashed") +
  annotate("text", x = CON_LEN / 2, y = 0.75,
    label = "Sequence divergence from consensus (Kimura-80 distance)",
    size = 2.3, colour = "#888888") +
  coord_cartesian(xlim = c(-280, CON_LEN + 130), ylim = c(0.5, 2.3)) +
  theme_void() +
  theme(plot.margin = margin(2, 4, 4, 4))

# ---- HTT annotation ----------------------------------------------------------
p_htt <- ggplot() +
  geom_rect(aes(xmin = 0.05, xmax = 0.95, ymin = 0.2, ymax = 0.8),
            fill = "#EBF5FB", colour = "#2980B9", linewidth = 0.5) +
  annotate("text", x = 0.5, y = 0.65,
    label = "HTT event: Holozoa  ↔  Mortierellomycotina",
    size = 3, fontface = "bold", colour = "#1A5276") +
  annotate("text", x = 0.5, y = 0.38,
    label = paste0(
      "Nucleotide identity: 80.4%  |  ",
      "No flanking synteny (HTT score = 1.0)  |  ",
      "Single validated event (step 09)"
    ),
    size = 2.5, colour = "#2471A3") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_void() +
  theme(plot.margin = margin(2, 4, 2, 4))

# ---- Assemble ----------------------------------------------------------------
fig <- p_overview / p_instances / p_htt +
  plot_layout(heights = c(2.5, 1.8, 0.8)) +
  plot_annotation(
    title = "LINE TE architecture — ZYMTRI_2023_FAMILY_141",
    theme = theme(
      plot.title  = element_text(size = 10, face = "bold"),
      plot.margin = margin(6, 6, 4, 6)
    )
  )

ggsave(OUT, fig, width = 8, height = 5, device = "pdf")
cat("[HTT-domain] Figure saved to:", OUT, "\n")
