#!/usr/bin/env Rscript
# HTT LINE locus – synteny snapshot (v2, gene-arrow style)
# HO_Monbrevi (NW_001865040.1:1290160-1310734)  vs
# MO_MoelGBA34 (scaffold_76:25087-45661)
# TE: ZYMTRI_2023_FAMILY_141 (LINE), 574 bp, 80.37 % identity
# Domain: RepeatMasker/TEsorter → LINE (no more specific assignment possible)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(grid)
})

OUT <- "/legserv/Temp/Carlos/EDF-TE/results/HTT/figures/htt_line_synteny_v2.pdf"

# Species for each MO track (window start coords)
# MO_MoelGBA34  : scaffold_76:25087       (win start 25087, LINE rel 10001-10574)
# MO_MocaN2610  : scaffold_92:14128       (win start 14128, LINE rel 10001-10574)
# MO_MozyN2592  : scaffold_191:20855      (win start 20855, LINE rel 10001-10574)

WIN_LEN  <- 20575L          # bp (same for both windows)
GENE_Y   <- 0.50            # y-centre of gene track
GENE_H   <- 0.16            # half-height of gene body
TIP_H    <- 0.25            # half-height at arrow tip
TIP_W    <- 500             # bp occupied by the arrowhead
TE_Y     <- 0.50
TE_H     <- 0.22            # taller than gene body

# ── colour palette ──────────────────────────────────────────────────────────
COLS <- list(
  gene_def   = "#BDBDBD",   # default gene (gray)
  gene_host  = "#7EC8C8",   # gene hosting the LINE (teal)
  line_focal = "#E84B37",   # ZYMTRI_141 LINE
  line_other = "#F5A470",   # other LINEs
  ltr        = "#C678DD",
  dna        = "#56B6C2",
  simple     = "#DDDDDD",
  low_cplx   = "#C0C0C0",
  unknown    = "#A8D5A2",
  ribbon     = "#E84B3740"  # semi-transparent ribbon fill
)

# ── identity → color mapping ─────────────────────────────────────────────
# pct: 0-100; returns hex color (no alpha). Monochromatic: 50%=pale blue → 100%=dark blue
ident_color <- function(pct) {
  ramp <- colorRampPalette(c("#C8E6FA", "#014F86"))(100)
  idx  <- max(1L, min(100L, as.integer((pct - 50L) * 2L)))
  ramp[idx]
}

# ── helper: build chevron polygon coords for one gene ──────────────────────
gene_polygon <- function(start, end, strand, y = GENE_Y, h = GENE_H, tip = TIP_H, tw = TIP_W, id = "g") {
  span <- end - start
  if (span <= 0) return(NULL)
  tw2 <- min(tw, span * 0.4)  # never let tip > 40% of gene
  if (strand == "+") {
    body <- end - tw2
    xs <- c(start, body,  body,         end,   body,         body,  start)
    ys <- c(y - h, y - h, y - tip, y,   y + tip, y + h, y + h)
  } else {
    body <- start + tw2
    xs <- c(end,   body,  body,         start, body,         body,  end)
    ys <- c(y + h, y + h, y + tip, y,   y - tip, y - h, y - h)
  }
  data.frame(x = xs, y = ys, id = id)
}

# ── helper: build a simple TE rectangle ─────────────────────────────────────
te_rect <- function(start, end, col_key, label = NA) {
  data.frame(start = start, end = end, col_key = col_key, label = label)
}

# ── HO_Monbrevi genes (relative to window start = 1290160) ──────────────────
#  (from NCBI GFF3; MONBRDRAFT_17631 contains the LINE)
#  MONBRDRAFT_17631 exon structure (window coords, minus strand):
#    Exon1 (5'): 11599-11708 (UTR)  |  Intron: 11029-11598  |  Exon2 (3'): 9633-11028
#    CDS within Exon2: 9991-10857  <- LINE (10001-10574) is WITHIN the CDS!
ho_genes <- list(
  list(start = 1119,  end = 4547,  strand = "+", name = "MONBRDRAFT_31239", host = FALSE,
       exons = list(c(1120,1281),c(1355,1391),c(1482,1680),c(1783,1900),c(1978,2608),
                    c(2738,2987),c(3077,3161),c(3363,3584),c(3735,3856),c(3995,4110),
                    c(4204,4364),c(4471,4548))),
  list(start = 4868,  end = 5924,  strand = "-", name = "MONBRDRAFT_14285", host = FALSE,
       exons = list(c(4869,4973),c(5090,5158),c(5240,5329),c(5435,5605),c(5728,5925))),
  list(start = 6040,  end = 6735,  strand = "-", name = "MONBRDRAFT_5820",  host = FALSE,
       exons = list(c(6041,6130),c(6293,6736))),
  list(start = 7405,  end = 8438,  strand = "-", name = "MONBRDRAFT_5821",  host = FALSE,
       exons = list(c(7406,8003),c(8132,8439))),
  list(start = 9632,  end = 11707, strand = "-", name = "MONBRDRAFT_17631", host = TRUE,
       exons = list(c(9633, 11028), c(11599, 11707)),
       utrs  = list(c(9633, 9990), c(10858, 11028), c(11599, 11707)),
       cds   = list(c(9991, 10857))),                   # CDS: contains the LINE!
  list(start = 12344, end = 14253, strand = "-", name = "MONBRDRAFT_22954", host = FALSE,
       exons = list(c(12345,13877),c(13920,14038),c(14146,14254))),
  list(start = 14545, end = 16078, strand = "-", name = "TUA",              host = FALSE,
       exons = list(c(14546,16079))),
  list(start = 16978, end = WIN_LEN, strand = "+", name = "MONBRDRAFT_35916", host = FALSE,
       exons = list(c(16979,17113),c(17628,17726),c(17966,18352),c(18758,18841),
                    c(19068,19439),c(19712,19815),c(19893,20253),c(20363,20575)))
)

# ── HO TEs (relative coords) ────────────────────────────────────────────────
ho_tes <- data.frame(
  start   = c(560,  666,  764,  846,  3165, 7225, 10001, 11888,12172,14639,15758,16753,17807,18803),
  end     = c(594,  763,  783,  864,  3207, 7243, 10574, 11906,12188,15972,15829,16800,17846,18824),
  col_key = c(rep("simple", 6), "line_focal", "simple","simple","line_other","unknown","simple","simple","simple"),
  label   = c(rep(NA, 6), "ZYMTRI_141", rep(NA, 7))
)

# ── MO_MoelGBA34 genes (relative to window start = 25087) ───────────────────
#  BGZ88_011182 contains the LINE
#  BGZ88_011182 exon structure (4 exons, all CDS, minus strand):
#    Exon1: 9991-11056 | Intron: 11057-11173 | Exon2: 11174-11253 | Int: 11254-11412
#    Exon3: 11413-11455 | Intron: 11456-12160 | Exon4: 12161-12173
#    LINE (10001-10574) is WITHIN Exon1/CDS!
mo_genes <- list(
  list(start = 2223,  end = 4047,  strand = "+", name = "BGZ88_011180",    host = FALSE,
       exons = list(c(2224,2300),c(2436,2713),c(3143,3198),c(3626,4048))),
  list(start = 4721,  end = 9536,  strand = "+", name = "BGZ88_011181",    host = FALSE,
       exons = list(c(4722,5648),c(5782,5803),c(5906,9537))),
  list(start = 9990,  end = 12173, strand = "-", name = "BGZ88_011182",    host = TRUE,
       exons = list(c(9991,11056), c(11174,11253), c(11413,11455), c(12161,12173)),
       utrs  = list(),
       cds   = list(c(9991,11056), c(11174,11253), c(11413,11455), c(12161,12173))),
  list(start = 12864, end = 15654, strand = "+", name = "BGZ88_011183",    host = FALSE,
       exons = list(c(12865,12998),c(13106,13140),c(13246,13351),c(13457,13523),
                    c(13613,13957),c(14056,14229),c(14339,14642),c(14741,15219),
                    c(15356,15442),c(15554,15655)))
)

# ── MO TEs ──────────────────────────────────────────────────────────────────
mo_tes <- data.frame(
  start   = c(148,  1362, 1730, 2961, 4306, 4757, 4907, 5002, 5171, 5267, 6843, 7989, 8231, 9586, 10001,10958,12007,13658,14116,16130,18454,18525,18799,19377,20164,20310),
  end     = c(197,  1391, 1812, 2982, 4364, 4800, 4940, 5170, 5266, 5450, 6862, 8072, 8276, 9690, 10574,11038,12040,13715,14204,16271,18524,18617,18922,19423,20220,20340),
  col_key = c("low_cplx","simple","simple","simple","simple","simple","simple","unknown","line_other","line_other","simple","line_other","line_other","low_cplx","line_focal","unknown","simple","simple","ltr","dna","simple","line_other","dna","simple","low_cplx","low_cplx"),
  label   = c(rep(NA,14),"ZYMTRI_141",rep(NA,11))
)

# ── MO_MocaN2610 genes (relative to window start = 14128 on scaffold_92) ────
#  BGZ47_011432 (- strand) is the host gene; LINE (10001-10574) within its CDS
moca_genes <- list(
  list(start = 2232,  end = 4016,  strand = "+", name = "BGZ47_011430",    host = FALSE,
       exons = list(c(2232,2308),c(2440,2717),c(3118,3173),c(3594,4016))),
  list(start = 4808,  end = 9546,  strand = "+", name = "BGZ47_011431",    host = FALSE,
       exons = list(c(4808,5575),c(5684,5714),c(5819,9546))),
  list(start = 9990,  end = 11509, strand = "-", name = "BGZ47_011432",    host = TRUE,
       exons = list(c(9990,11055),c(11182,11261),c(11409,11444),c(11504,11509)),
       utrs  = list(),
       cds   = list(c(9990,11055),c(11182,11261),c(11409,11444),c(11504,11509))),
  list(start = 14604, end = 16919, strand = "+", name = "BGZ47_011433",    host = FALSE,
       exons = list(c(14604,14607),c(14665,14814),c(15101,15266),c(15375,15409),
                    c(15522,15627),c(15729,15795),c(15880,16227),c(16324,16497),
                    c(16616,16919)))
)

moca_tes <- data.frame(
  start   = c(5027, 5132, 7834, 9127, 10001, 12545, 15926),
  end     = c(5131, 5212, 8037, 9228, 10574, 12583, 16057),
  col_key = c("line_other","line_other","line_other","line_other","line_focal","dna","dna"),
  label   = c(rep(NA,4),"ZYMTRI_141",rep(NA,2))
)

# ── MO_MozyN2592 genes (relative to window start = 20855 on scaffold_191) ───
#  EC991_003374 (- strand) is the host gene; LINE (10001-10574) within its CDS
mozy_genes <- list(
  list(start = 2522,  end = 4326,  strand = "+", name = "EC991_003372",    host = FALSE,
       exons = list(c(2522,2631),c(2723,2994),c(3436,3485),c(3901,4326))),
  list(start = 4758,  end = 9503,  strand = "+", name = "EC991_003373",    host = FALSE,
       exons = list(c(4758,4763),c(5244,5759),c(5955,7910),c(8367,9503))),
  list(start = 9990,  end = 11651, strand = "-", name = "EC991_003374",    host = TRUE,
       exons = list(c(9990,11055),c(11193,11272),c(11415,11457),c(11617,11651)),
       utrs  = list(),
       cds   = list(c(9990,11055),c(11193,11272),c(11415,11457),c(11617,11651))),
  list(start = 12756, end = 15279, strand = "+", name = "EC991_003375",    host = FALSE,
       exons = list(c(12756,12889),c(12996,13030),c(13132,13237),c(13487,13553),
                    c(13659,14006),c(14101,14274),c(14386,14686),c(14795,15279)))
)

mozy_tes <- data.frame(
  start   = c(5145, 10001, 13673, 13691, 18733),
  end     = c(5263, 10574, 13820, 13778, 18867),
  col_key = c("dna","line_focal","dna","line_other","ltr"),
  label   = c(NA,"ZYMTRI_141",NA,NA,NA)
)

# ── nucmer alignment coords ──────────────────────────────────────────────────
ALN_Q_S <- 9991L;  ALN_Q_E <- 10746L
ALN_S_S <- 9991L;  ALN_S_E <- 10746L

# Identity of Monosiga LINE vs each MO species (nt %)
ALN_PCTS <- c(
  "L. elongata GBA34"  = 80.4,
  "L. camargensis N2610" = 80.6,
  "L. zychae N2592"    = 79.6
)
# Host-gene protein aa identity: Monosiga vs Linnemannia (mean over 3 species)
HO_MO_HOST_PCTID <- 56.2

# ── build track: gene arrows (exon/intron) + TE boxes ───────────────────────
build_track <- function(genes_list, tes_df, species_label, chrom_label,
                        xlim_l = -3200, xlim_r = WIN_LEN + 300,
                        show_xaxis = FALSE, compact = FALSE) {

  # filter TEs: keep only full elements (remove simple repeats / low complexity)
  tes_df <- tes_df[!tes_df$col_key %in% c("simple", "low_cplx", "unknown"), ]

  # Baseline (chromosome line)
  p <- ggplot() +
    annotate("segment", x = 0, xend = WIN_LEN, y = GENE_Y, yend = GENE_Y,
             colour = "#AAAAAA", linewidth = 0.25)

  for (i in seq_along(genes_list)) {
    g <- genes_list[[i]]
    gs <- max(0L, g$start); ge <- min(WIN_LEN, g$end)
    if (gs >= ge) next

    gcol  <- if (g$host) COLS$gene_host else COLS$gene_def
    gbord <- if (g$host) "#4DAAAA"      else "#999999"
    bbcol <- if (g$host) "#4DAAAA"      else "#BBBBBB"   # intron backbone

    if (is.null(g$exons) || length(g$exons) == 0) {
      # fallback: single chevron arrow
      poly <- gene_polygon(gs, ge, g$strand, id = paste0("g", i))
      p <- p + geom_polygon(data = poly, aes(x = x, y = y, group = id),
                            fill = gcol, colour = gbord, linewidth = 0.22)
    } else {
      exons_s <- do.call(rbind, lapply(g$exons, function(e) data.frame(s = e[1], e = e[2])))
      exons_s <- exons_s[order(exons_s$s), ]
      n_ex    <- nrow(exons_s)
      # terminal exon (gets arrowhead): last on + strand, first on - strand
      term_j  <- if (g$strand == "+") n_ex else 1L

      # intron backbone
      p <- p + annotate("segment", x = gs, xend = ge, y = GENE_Y, yend = GENE_Y,
                        colour = bbcol, linewidth = 0.5)

      for (j in seq_len(n_ex)) {
        exs <- max(0L, as.integer(exons_s$s[j]))
        exe <- min(WIN_LEN, as.integer(exons_s$e[j]))
        if (exs >= exe) next

        if (j == term_j) {
          # terminal exon: chevron arrow
          poly <- gene_polygon(exs, exe, g$strand, id = paste0("g", i, "_e", j))
          p <- p + geom_polygon(data = poly, aes(x = x, y = y, group = id),
                                fill = gcol, colour = gbord, linewidth = 0.22)
        } else {
          # non-terminal exon: plain rectangle
          p <- p + annotate("rect", xmin = exs, xmax = exe,
                            ymin = GENE_Y - GENE_H, ymax = GENE_Y + GENE_H,
                            fill = gcol, color = gbord, linewidth = 0.22)
        }
      }

      # host genes: overlay CDS (taller, solid teal) and UTR (shorter, light teal)
      if (g$host) {
        if (!is.null(g$utrs) && length(g$utrs) > 0) {
          for (u in g$utrs) {
            us <- max(0L, u[1]); ue <- min(WIN_LEN, u[2])
            if (us >= ue) next
            p <- p + annotate("rect", xmin = us, xmax = ue,
                              ymin = GENE_Y - GENE_H, ymax = GENE_Y + GENE_H,
                              fill = "#B5E5E5", color = "#4DAAAA", linewidth = 0.22)
          }
        }
        if (!is.null(g$cds) && length(g$cds) > 0) {
          for (cd in g$cds) {
            cs <- max(0L, cd[1]); ce <- min(WIN_LEN, cd[2])
            if (cs >= ce) next
            p <- p + annotate("rect", xmin = cs, xmax = ce,
                              ymin = GENE_Y - TIP_H * 0.82, ymax = GENE_Y + TIP_H * 0.82,
                              fill = COLS$gene_host, color = "#4DAAAA", linewidth = 0.22)
          }
        }
      }
    }
  }

  # Gene labels (non-host: angled above) — suppressed in compact mode
  if (!compact) {
    lbl_df <- do.call(rbind, lapply(genes_list, function(g) {
      if (is.na(g$name) || g$host) return(NULL)
      mid <- (max(0L, g$start) + min(WIN_LEN, g$end)) / 2
      data.frame(x = mid, y = GENE_Y + TIP_H + 0.08, label = g$name)
    }))
    if (!is.null(lbl_df) && nrow(lbl_df) > 0)
      p <- p + geom_text(data = lbl_df, aes(x = x, y = y, label = label),
          size = 1.9, fontface = "italic", colour = "#555555",
          angle = 35, hjust = 0, vjust = 0)
  }

  # Host gene labels (teal, below)
  host_lbl <- do.call(rbind, lapply(genes_list, function(g) {
    if (!g$host) return(NULL)
    mid <- (max(0L, g$start) + min(WIN_LEN, g$end)) / 2
    data.frame(x = mid, y = GENE_Y - TIP_H - 0.09, label = g$name)
  }))
  if (!is.null(host_lbl) && nrow(host_lbl) > 0)
    p <- p + geom_text(data = host_lbl, aes(x = x, y = y, label = label),
        size = 2.2, fontface = "italic", colour = "#2E8B8B", vjust = 1)

  # TE boxes
  te_cols <- c(
    line_focal = COLS$line_focal, line_other = COLS$line_other,
    ltr = COLS$ltr, dna = COLS$dna
  )
  p <- p + geom_rect(data = tes_df,
    aes(xmin = start, xmax = end, ymin = TE_Y - TE_H, ymax = TE_Y + TE_H, fill = col_key),
    colour = NA) +
  geom_text(data = tes_df[!is.na(tes_df$label), ],
    aes(x = (start + end)/2, y = TE_Y + TE_H + 0.07, label = label),
    size = 2.3, fontface = "bold.italic", colour = COLS$line_focal, vjust = 0) +
  scale_fill_manual(values = te_cols, na.value = COLS$gene_def) +
  coord_cartesian(xlim = c(xlim_l, xlim_r), ylim = c(0, 1), clip = "off") +
  theme_void() + theme(legend.position = "none", plot.margin = margin(2,2,2,2))

  # Species label on left
  lbl_x <- xlim_l + 80
  p <- p +
    annotate("text", x = lbl_x, y = GENE_Y + 0.06, label = species_label,
             hjust = 0, vjust = 0.5, size = 2.7, fontface = "bold.italic", colour = "#1A1A1A") +
    annotate("text", x = lbl_x, y = GENE_Y - 0.14, label = chrom_label,
             hjust = 0, vjust = 0.5, size = 2.0, colour = "#666666")

  if (show_xaxis) {
    kb_breaks <- seq(0, WIN_LEN, by = 5000)
    p <- p +
      annotate("segment", x = kb_breaks, xend = kb_breaks,
               y = 0.07, yend = 0.02, colour = "#888888", linewidth = 0.25) +
      annotate("text", x = kb_breaks, y = -0.04,
               label = paste0(kb_breaks/1000, " kb"), size = 1.8, colour = "#777777") +
      annotate("segment", x = 0, xend = WIN_LEN,
               y = 0.07, yend = 0.07, colour = "#888888", linewidth = 0.2)
  }

  p
}

# ── ribbon panel (only TE alignment) ────────────────────────────────────────
build_ribbon <- function(xlim_l = -3200, xlim_r = WIN_LEN + 300) {
  # Ribbon polygon connecting the two tracks (top = 1, bottom = 0)
  xs <- c(ALN_Q_S, ALN_Q_E, ALN_S_E, ALN_S_S, ALN_Q_S)
  ys <- c(1,       1,       0,       0,       1)

  mid_x <- (ALN_Q_S + ALN_Q_E) / 2

  # Build multi-identity label from ALN_PCTS
  ident_labels <- paste(
    sprintf("%s: %.1f%%", names(ALN_PCTS), ALN_PCTS),
    collapse = "    "
  )

  line_col   <- ident_color(HO_MO_HOST_PCTID)   # protein aa identity, same scale as MO-MO ribbons
  fill_col   <- paste0(line_col, "40")
  border_col <- paste0(line_col, "90")

  ggplot() +
    geom_polygon(aes(x = xs, y = ys),
                 fill = fill_col, colour = border_col, linewidth = 0.3) +
    # TE identity labels for all 3 MO species
    annotate("text", x = mid_x, y = 0.55,
             label = ident_labels,
             size = 2.3, colour = "#C0392B", fontface = "italic", hjust = 0.5) +
    # "No flanking synteny" label in grey on the far left
    annotate("text", x = xlim_l + 300, y = 0.5,
             label = "no flanking\nsynteny",
             size = 2.0, colour = "#AAAAAA", hjust = 0, lineheight = 0.9) +
    coord_cartesian(xlim = c(xlim_l, xlim_r), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 2, 0, 2))
}

# ── ribbon panel (gene-to-gene synteny between two species) ─────────────────
# pct_ids: numeric vector of protein % aa identity, one per gene pair (same order as genes)
build_ribbon_synteny <- function(from_genes, to_genes, pct_ids,
                                  xlim_l = -3200, xlim_r = WIN_LEN + 300) {
  p <- ggplot()

  for (i in seq_along(from_genes)) {
    fg <- from_genes[[i]]
    tg <- to_genes[[i]]
    fs <- max(0L, fg$start); fe <- min(WIN_LEN, fg$end)
    ts <- max(0L, tg$start); te <- min(WIN_LEN, tg$end)
    if (fs >= fe || ts >= te) next
    base_col <- ident_color(pct_ids[i])
    fll <- paste0(base_col, "50")
    brd <- paste0(base_col, "A0")
    poly <- data.frame(x = c(fs, fe, te, ts), y = c(1, 1, 0, 0))
    p <- p + geom_polygon(data = poly, aes(x = x, y = y),
                          fill = fll, colour = brd, linewidth = 0.22)
    # % identity label in centre of each trapezoid
    mx <- (fs + fe + ts + te) / 4
    p <- p + annotate("text", x = mx, y = 0.5,
                      label = sprintf("%.0f%%", pct_ids[i]),
                      size = 1.7, colour = "white", fontface = "bold")
  }

  # LINE ribbon on top (always at 10001-10574 in both species)
  host_i   <- which(sapply(from_genes, function(g) isTRUE(g$host)))
  line_pct <- if (length(host_i) > 0) pct_ids[host_i[1]] else mean(pct_ids)
  lc       <- ident_color(line_pct)
  line_poly <- data.frame(x = c(10001L, 10574L, 10574L, 10001L), y = c(1, 1, 0, 0))
  p <- p +
    geom_polygon(data = line_poly, aes(x = x, y = y),
                 fill = paste0(lc, "60"), colour = paste0(lc, "CC"), linewidth = 0.3) +
    coord_cartesian(xlim = c(xlim_l, xlim_r), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 2, 0, 2))
  p
}

# ── legend panel ────────────────────────────────────────────────────────────
build_legend <- function() {
  # Build a small identity gradient swatch (5 steps from 60% to 95%)
  ident_steps <- seq(60, 95, length.out = 5)
  ident_df <- data.frame(
    xi  = seq(6.05, 7.65, length.out = 5),
    col = sapply(ident_steps, ident_color)
  )
  items <- data.frame(
    col_key = c("gene_def", "gene_utr", "gene_host", "line_focal", "line_other", "ltr", "dna"),
    label   = c("Gene (exons shown)",
                "UTR exon",
                "CDS exon / host gene (contains LINE)",
                "LINE (transferred, ZYMTRI_141)",
                "LINE (other)", "LTR retrotransposon", "DNA transposon"),
    x = c(0, 1.7, 3.4, 6.0, 8.0, 9.8, 11.5),
    y = rep(0.55, 7)
  )
  col_vals <- c(
    gene_def   = COLS$gene_def,   gene_utr = "#B5E5E5",
    gene_host  = COLS$gene_host,
    line_focal = COLS$line_focal, line_other = COLS$line_other,
    ltr = COLS$ltr, dna = COLS$dna
  )
  ggplot(items) +
    geom_rect(aes(xmin = x, xmax = x + 0.25, ymin = 0.35, ymax = 0.75, fill = col_key),
              colour = "#888888", linewidth = 0.18) +
    geom_text(aes(x = x + 0.37, y = 0.55, label = label), hjust = 0, size = 2.2) +
    # identity gradient swatch (60%–95%)
    geom_rect(data = ident_df,
              aes(xmin = xi, xmax = xi + 0.32, ymin = 0.0, ymax = 0.30, fill = I(col)),
              colour = NA) +
    annotate("text", x = 5.9,  y = 0.15, label = "55%", hjust = 1,   size = 1.8, colour = "#555555") +
    annotate("text", x = 8.05, y = 0.15, label = "95%", hjust = 0,   size = 1.8, colour = "#555555") +
    annotate("text", x = 6.97, y = 0.15, label = "ribbon = % aa identity", hjust = 0.5, size = 1.8, colour = "#555555") +
    scale_fill_manual(values = col_vals) +
    coord_cartesian(xlim = c(-0.1, 13.5), ylim = c(0, 1)) +
    theme_void() + theme(legend.position = "none")
}

# ── scale bar panel ──────────────────────────────────────────────────────────
build_scalebar <- function(xlim_l = -3200, xlim_r = WIN_LEN + 300) {
  # 10 kb scale bar placed near right end
  sb_end   <- WIN_LEN - 500
  sb_start <- sb_end - 10000
  ggplot() +
    annotate("segment", x = sb_start, xend = sb_end, y = 0.5, yend = 0.5,
             colour = "#222222", linewidth = 0.8) +
    annotate("segment", x = c(sb_start, sb_end), xend = c(sb_start, sb_end),
             y = 0.35, yend = 0.65, colour = "#222222", linewidth = 0.6) +
    annotate("text", x = (sb_start + sb_end)/2, y = 0.2,
             label = "10 kb", size = 2.6, colour = "#222222") +
    coord_cartesian(xlim = c(xlim_l, xlim_r), ylim = c(0, 1)) +
    theme_void() + theme(plot.margin = margin(0, 2, 0, 2))
}

# ── assemble ─────────────────────────────────────────────────────────────────
XL <- -4800L
XR <- WIN_LEN + 400L

p_ho  <- build_track(ho_genes, ho_tes,
           "Monosiga brevicollis",
           "NW_001865040.1:1,290,160\u20131,310,734  (+)",
           xlim_l = XL, xlim_r = XR, show_xaxis = FALSE)

p_rib <- build_ribbon(xlim_l = XL, xlim_r = XR)

p_mo  <- build_track(mo_genes, mo_tes,
           "Linnemannia elongata GBA34",
           "scaffold_76:25,087\u201345,661  (+)",
           xlim_l = XL, xlim_r = XR, show_xaxis = FALSE)

# Protein aa identity (pairwise local BLOSUM62 alignment) per gene pair
# Order: gene1(+), gene2(+), host-gene(-), gene4(+)
pct_mo_moca   <- c(90.0, 68.8, 88.8, 86.2)   # L. elongata GBA34 vs L. camargensis N2610
pct_moca_mozy <- c(72.5, 61.2, 86.2, 75.0)   # L. camargensis N2610 vs L. zychae N2592

p_rib2 <- build_ribbon_synteny(mo_genes, moca_genes, pct_mo_moca,   xlim_l = XL, xlim_r = XR)

p_moca <- build_track(moca_genes, moca_tes,
           "Linnemannia camargensis N2610",
           "scaffold_92:14,128–34,703  (+)",
           xlim_l = XL, xlim_r = XR, show_xaxis = FALSE, compact = FALSE)

p_rib3 <- build_ribbon_synteny(moca_genes, mozy_genes, pct_moca_mozy, xlim_l = XL, xlim_r = XR)

p_mozy <- build_track(mozy_genes, mozy_tes,
           "Linnemannia zychae N2592",
           "scaffold_191:20,855\u201341,430  (+)",
           xlim_l = XL, xlim_r = XR, show_xaxis = FALSE, compact = FALSE)

p_sb  <- build_scalebar(xlim_l = XL, xlim_r = XR)
p_leg <- build_legend()

fig <- (p_ho / p_rib / p_mo / p_rib2 / p_moca / p_rib3 / p_mozy / p_sb / p_leg) +
  plot_layout(heights = c(2.2, 0.8, 2.0, 0.5, 2.0, 0.5, 2.0, 0.5, 0.7)) +
  plot_annotation(
    title    = "HTT event: LINE (ZYMTRI_2023_FAMILY_141) shared between Holozoa and Mortierellomycotina",
    subtitle = paste0("574 bp | LINE family in CDS of host genes | TEsorter/RepeatMasker class: LINE | no flanking synteny (HTT score = 1.0)\n",
                      "Three Linnemannia species share the LINE at the same syntenic locus (~80% identity to Monosiga brevicollis)"),
    theme = theme(
      plot.title    = element_text(size = 9, face = "bold",   margin = margin(b = 2)),
      plot.subtitle = element_text(size = 7.5, colour = "#555555"),
      plot.margin   = margin(6, 6, 4, 8)
    )
  )

ggsave(OUT, fig, width = 11, height = 11.0, device = "pdf")
cat("[HTT-synteny-v2] Saved:", OUT, "\n")
