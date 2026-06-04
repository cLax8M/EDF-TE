#!/usr/bin/env Rscript
# fig_H_plot_phylo.R
# Visualize co-opted gene phylogenetic trees (one per architecture)
# Colored by clade; handles multiple architectures automatically.

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(cowplot)
})

OUT_DIR <- "/legserv/Temp/Carlos/EDF-TE/results/Co-option/figures"

CLADE_LABELS <- c(
  AP="Aphelida", AS="Ascomycota", BA="Basidiomycota",
  BL="Blastocladiomycota", BS="Basidiobolomycota", CH="Chytridiomycota",
  CR="Cryptomycota", EN="Entomophthoromycota", GL="Glomeromycota",
  HO="Holozoa", KX="Kickxellomycotina", MI="Microsporidia",
  MO="Mortierellomycota", MU="Mucoromycotina", NU="Nuclearida",
  OL="Olpidiomycota", ZO="Zoopagomycota"
)
CLADE_COLORS <- c(
  HO="#e07b54", MI="#b0b0b0", NU="#d4b8e0", AP="#5ba3d9",
  GL="#2e7d32", MO="#81c784", MU="#c8e6c9", ZO="#ff8a65",
  EN="#ffb74d", BS="#f06292", KX="#ba68c8", BL="#64b5f6",
  CH="#4dd0e1", CR="#fff176", OL="#a5d6a7", AS="#ff7043",
  BA="#8d6e63"
)

# Build named color vector by clade label (for ggplot scale_color_manual)
named_colors <- setNames(
  CLADE_COLORS[names(CLADE_LABELS)],
  CLADE_LABELS[names(CLADE_LABELS)]
)

# ── Domain architecture schemas ────────────────────────────────────────────────
# Each entry: list of rows, each row = list(label, domains)
# Domain: name, x0, x1, fill ("te"=TE-derived, "host"=host, "bg"=backbone)
# x scale: 0–100 (relative positions, schematic)
COL_TE   <- "#c0392b"    # deep red   – TE-derived domain
COL_HOST <- "#2980b9"    # steel blue – host (non-TE) domain
COL_BG   <- "#d9d9d9"    # light grey – gene backbone bar

DOMAIN_SCHEMAS <- list(
  Chromo_Helitron = list(
    title = "Chromo + Helitron architecture",
    subtitle = "Chromo domain co-opted from Helitron transposons",
    rows = list(
      list(label = "Helitron TE",
           domains = list(
             list(nm="", x0=0, x1=100, fill=COL_BG),        # backbone
             list(nm="HEL1",  x0=2,  x1=30, fill=COL_HOST),
             list(nm="HEL2",  x0=33, x1=60, fill=COL_HOST),
             list(nm="HA2",   x0=63, x1=78, fill=COL_HOST),
             list(nm="Chromo",x0=82, x1=98, fill=COL_TE)
           )),
      list(label = "Co-opted gene",
           domains = list(
             list(nm="", x0=0, x1=100, fill=COL_BG),        # backbone
             list(nm="SNF2_N", x0=2,  x1=38, fill=COL_HOST),
             list(nm="DEAD",   x0=41, x1=68, fill=COL_HOST),
             list(nm="Chromo", x0=82, x1=98, fill=COL_TE)
           ))
    )
  ),
  DUF4217_LINE = list(
    title = "DUF4217 + LINE architecture",
    subtitle = "DUF4217 domain co-opted from LINE retrotransposons",
    rows = list(
      list(label = "LINE TE",
           domains = list(
             list(nm="", x0=0, x1=100, fill=COL_BG),
             list(nm="ORF1",   x0=2,  x1=28, fill=COL_HOST),
             list(nm="RVT_1",  x0=31, x1=60, fill=COL_HOST),
             list(nm="RT_RNaseH", x0=63, x1=80, fill=COL_HOST),
             list(nm="DUF4217",x0=83, x1=98, fill=COL_TE)
           )),
      list(label = "Co-opted gene",
           domains = list(
             list(nm="", x0=0, x1=100, fill=COL_BG),
             list(nm="host",   x0=2,  x1=30, fill=COL_HOST),
             list(nm="DUF4217",x0=33, x1=60, fill=COL_TE),
             list(nm="host",   x0=63, x1=98, fill=COL_HOST)
           ))
    )
  ),
  APE_LINE = list(
    title = "AP endonuclease + LINE architecture",
    subtitle = "AP endonuclease (EN) domain co-opted from LINE retrotransposons",
    rows = list(
      list(label = "LINE TE",
           domains = list(
             list(nm="", x0=0, x1=100, fill=COL_BG),
             list(nm="AP_endo", x0=2,  x1=25, fill=COL_TE),
             list(nm="RVT_1",   x0=28, x1=60, fill=COL_HOST),
             list(nm="RT_RNaseH",x0=63,x1=80, fill=COL_HOST)
           )),
      list(label = "Co-opted gene",
           domains = list(
             list(nm="", x0=0, x1=100, fill=COL_BG),
             list(nm="host",   x0=2,  x1=30, fill=COL_HOST),
             list(nm="AP_endo",x0=33, x1=65, fill=COL_TE),
             list(nm="host",   x0=68, x1=98, fill=COL_HOST)
           ))
    )
  ),
  ExoDNase_LINE = list(
    title = "Exodeoxyribonuclease III + LINE architecture",
    subtitle = "ExoDNase III (xth/Endo IV) domain co-opted from LINE retrotransposons",
    rows = list(
      list(label = "LINE TE",
           domains = list(
             list(nm="", x0=0, x1=100, fill=COL_BG),
             list(nm="ExoDNase",x0=2,  x1=30, fill=COL_TE),
             list(nm="RVT_1",  x0=33, x1=63, fill=COL_HOST),
             list(nm="RNaseH",  x0=66, x1=82, fill=COL_HOST)
           )),
      list(label = "Co-opted gene",
           domains = list(
             list(nm="", x0=0, x1=100, fill=COL_BG),
             list(nm="host",    x0=2,  x1=28, fill=COL_HOST),
             list(nm="ExoDNase",x0=31, x1=65, fill=COL_TE),
             list(nm="host",    x0=68, x1=98, fill=COL_HOST)
           ))
    )
  ),
  LRR_LINE = list(
    title = "Leucine-rich repeat + LINE architecture",
    subtitle = "LRR domain gene fused with LINE-derived domain",
    rows = list(
      list(label = "LINE TE",
           domains = list(
             list(nm="", x0=0, x1=100, fill=COL_BG),
             list(nm="ORF1",  x0=2,  x1=25, fill=COL_HOST),
             list(nm="RVT_1", x0=28, x1=58, fill=COL_HOST),
             list(nm="LINE_dom",x0=62,x1=80, fill=COL_TE)
           )),
      list(label = "Co-opted gene",
           domains = list(
             list(nm="", x0=0, x1=100, fill=COL_BG),
             list(nm="LRR_1", x0=2,  x1=20, fill=COL_HOST),
             list(nm="LRR_2", x0=22, x1=40, fill=COL_HOST),
             list(nm="LRR_3", x0=42, x1=60, fill=COL_HOST),
             list(nm="LINE_dom",x0=63,x1=80, fill=COL_TE)
           ))
    )
  ),
  RhoGAP_LINE = list(
    title = "RhoGAP + LINE architecture",
    subtitle = "RhoGAP domain co-opted from LINE retrotransposons",
    rows = list(
      list(label = "LINE TE",
           domains = list(
             list(nm="", x0=0, x1=100, fill=COL_BG),
             list(nm="GAP_dom",x0=2,  x1=30, fill=COL_TE),
             list(nm="RVT_1",  x0=33, x1=63, fill=COL_HOST),
             list(nm="RNaseH",  x0=66, x1=82, fill=COL_HOST)
           )),
      list(label = "Co-opted gene",
           domains = list(
             list(nm="", x0=0, x1=100, fill=COL_BG),
             list(nm="host",   x0=2,  x1=28, fill=COL_HOST),
             list(nm="RhoGAP", x0=31, x1=68, fill=COL_TE),
             list(nm="host",   x0=71, x1=98, fill=COL_HOST)
           ))
    )
  ),
  DUF4208_Helitron = list(
    title = "DUF4208 + Helitron/Gypsy architecture",
    subtitle = "DUF4208 co-opted in Helitron-Gypsy chimera; Chromo domain from Chromovirus-type Gypsy",
    rows = list(
      list(label = "Helitron+Gypsy TE",
           domains = list(
             list(nm="",        x0=0,  x1=100, fill=COL_BG),
             list(nm="SNF2_N",  x0=2,  x1=26,  fill=COL_HOST),
             list(nm="Helicase",x0=29, x1=56,  fill=COL_HOST),
             list(nm="Chromo",  x0=59, x1=76,  fill=COL_TE),
             list(nm="DUF4208", x0=80, x1=98,  fill=COL_TE)
           )),
      list(label = "Co-opted gene",
           domains = list(
             list(nm="",        x0=0,  x1=100, fill=COL_BG),
             list(nm="SNF2_N",  x0=2,  x1=26,  fill=COL_HOST),
             list(nm="Helicase",x0=29, x1=56,  fill=COL_HOST),
             list(nm="Chromo",  x0=59, x1=76,  fill=COL_TE),
             list(nm="DUF4208", x0=80, x1=98,  fill=COL_TE)
           ))
    )
  )
)

make_domain_panel <- function(arch_label) {
  schema <- DOMAIN_SCHEMAS[[arch_label]]
  if (is.null(schema)) return(NULL)

  n_rows <- length(schema$rows)
  bar_h  <- 0.35   # height of each domain bar
  row_y  <- seq(n_rows, 1) * 1.0  # y positions (top row = highest)

  rects <- data.frame()
  texts <- data.frame()
  row_labels <- data.frame()

  for (i in seq_along(schema$rows)) {
    row  <- schema$rows[[i]]
    y_c  <- row_y[i]
    y0   <- y_c - bar_h / 2
    y1   <- y_c + bar_h / 2
    row_labels <- rbind(row_labels, data.frame(x=-2, y=y_c, label=row$label,
                                                stringsAsFactors=FALSE))
    for (d in row$domains) {
      is_bg <- d$nm == ""
      rects <- rbind(rects, data.frame(
        xmin=d$x0, xmax=d$x1,
        ymin=if(is_bg) y_c-0.08 else y0,
        ymax=if(is_bg) y_c+0.08 else y1,
        fill=d$fill,
        stringsAsFactors=FALSE
      ))
      if (d$nm != "") {
        texts <- rbind(texts, data.frame(
          x    = (d$x0 + d$x1) / 2,
          y    = y_c,
          label= d$nm,
          color= ifelse(d$fill == COL_BG, "black", "white"),
          stringsAsFactors=FALSE
        ))
      }
    }
  }

  p <- ggplot() +
    geom_rect(data=rects,
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=fill),
              color="white", linewidth=0.4) +
    scale_fill_identity() +
    geom_text(data=texts,
              aes(x=x, y=y, label=label, color=color),
              size=3.0, fontface="bold") +
    scale_color_identity() +
    geom_text(data=row_labels,
              aes(x=x, y=y, label=label),
              hjust=1, size=3.2, fontface="italic", color="grey20") +
    # Legend patches
    annotate("rect", xmin=2, xmax=10, ymin=-0.3, ymax=0.1,
             fill=COL_TE,   color="white", linewidth=0.3) +
    annotate("text", x=11, y=-0.1, label="TE-derived domain",
             hjust=0, size=2.8, color=COL_TE) +
    annotate("rect", xmin=40, xmax=48, ymin=-0.3, ymax=0.1,
             fill=COL_HOST, color="white", linewidth=0.3) +
    annotate("text", x=49, y=-0.1, label="Host domain",
             hjust=0, size=2.8, color=COL_HOST) +
    labs(title=schema$title, subtitle=schema$subtitle) +
    theme_void() +
    theme(
      plot.title    = element_text(size=9,  face="bold",   hjust=0),
      plot.subtitle = element_text(size=7.5, color="grey40", hjust=0),
      plot.margin   = margin(6, 6, 10, 80)
    ) +
    coord_cartesian(xlim=c(-30, 105), ylim=c(-0.6, n_rows + 0.6), clip="off")

  p
}

plot_arch_tree <- function(arch_label) {
  # IQ-TREE writes PREFIX.treefile
  tree_file <- file.path(OUT_DIR, paste0("figH_", arch_label, "_seqs_iqtree.treefile"))
  cand_file <- file.path(OUT_DIR, paste0("figH_", arch_label, "_candidates.tsv"))

  if (!file.exists(tree_file)) {
    cat("[SKIP] No tree file for", arch_label, ":", tree_file, "\n")
    return(invisible(NULL))
  }
  if (!file.exists(cand_file)) {
    cat("[SKIP] No candidates file for", arch_label, "\n")
    return(invisible(NULL))
  }

  tree  <- read.tree(tree_file)
  cands <- read.table(cand_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

  # Tip labels from python: species|clade|gene_id
  tip_meta <- data.frame(
    label   = tree$tip.label,
    species = str_extract(tree$tip.label, "^[^|]+"),
    clade   = str_extract(tree$tip.label, "(?<=\\|)[A-Z]{2}(?=\\|)"),
    gene_id = str_extract(tree$tip.label, "[^|]+$"),
    stringsAsFactors = FALSE
  )
  tip_meta$clade_label <- CLADE_LABELS[tip_meta$clade]
  tip_meta$clade_label[is.na(tip_meta$clade_label)] <- tip_meta$clade[is.na(tip_meta$clade_label)]

  # Short label for tip: genus first letter + species suffix
  tip_meta$short_sp <- sub("^([A-Z])[a-z]+_", "\\1.", tip_meta$species)

  n_tips <- length(tree$tip.label)
  show_labels <- n_tips <= 200

  p <- ggtree(tree, layout="fan", open.angle=15, size=0.25, color="grey40") %<+% tip_meta +
    geom_tippoint(aes(color=clade_label), size=ifelse(n_tips > 100, 1.5, 2.5), alpha=0.85)

  if (show_labels) {
    p <- p + geom_tiplab2(
      aes(label=paste0(short_sp, " | ", gene_id), color=clade_label),
      size=2.0, offset=0.005, fontface="italic"
    )
  }

  p <- p +
    scale_color_manual(values=named_colors, na.value="grey60", name="Clade") +
    geom_nodelab(
      aes(label=ifelse(!is.na(as.numeric(label)) & as.numeric(label) >= 70, label, "")),
      size=2.0, hjust=-0.1, color="grey40"
    ) +
    ggtitle(gsub("_", " + ", arch_label)) +
    theme_tree2() +
    theme(
      legend.position = "right",
      plot.title = element_text(size=12, face="bold")
    )

  # Fan layout: expand radially so tip labels don't clip
  if (show_labels) {
    p <- p + hexpand(0.4, direction = 1)
  }

  # Square canvas for fan layout
  fig_sq <- max(12, min(n_tips * 0.2, 60))
  fig_h  <- fig_sq
  fig_w  <- fig_sq

  pdf_out <- file.path(OUT_DIR, paste0("figH_", arch_label, "_tree.pdf"))
  png_out <- file.path(OUT_DIR, paste0("figH_", arch_label, "_tree.png"))
  ggsave(pdf_out, p, width=fig_w, height=fig_h, limitsize=FALSE)
  ggsave(png_out, p, width=fig_w, height=fig_h, dpi=150, limitsize=FALSE)
  cat("[OK]", arch_label, "->", pdf_out, "\n")

  # Domain architecture panel (separate PDF)
  dom_panel <- make_domain_panel(arch_label)
  if (!is.null(dom_panel)) {
    dom_pdf <- file.path(OUT_DIR, paste0("figH_", arch_label, "_domains.pdf"))
    dom_png <- file.path(OUT_DIR, paste0("figH_", arch_label, "_domains.png"))
    ggsave(dom_pdf, dom_panel, width=6, height=2.5)
    ggsave(dom_png, dom_panel, width=6, height=2.5, dpi=200)
    cat("[OK] Domain panel ->", dom_pdf, "\n")
  }
}

# Architectures to SKIP (known false-positive enrichment)
SKIP_ARCHS <- c("Helicase_Helitron")  # Helitron Rep subdomains ≠ co-opted host helicases

# Find all architectures with IQ-TREE output
tree_files <- list.files(OUT_DIR, pattern="figH_.*_seqs_iqtree\\.treefile$", full.names=FALSE)
  tree_files <- tree_files[!grepl("_seqs_seqs_", tree_files)]  # exclude old double-seqs files
arch_labels <- gsub("figH_(.*)_seqs_iqtree\\.treefile", "\\1", tree_files)
arch_labels <- arch_labels[!arch_labels %in% SKIP_ARCHS]

if (length(arch_labels) == 0) {
  cat("[ERROR] No IQ-TREE .treefile files found in", OUT_DIR, "\n")
  cat("  Expected pattern: figH_ARCH_seqs_aligned_iqtree.treefile\n")
  quit(status=1)
}

cat("[INFO] Plotting", length(arch_labels), "architectures:", paste(arch_labels, collapse=", "), "\n")
for (al in arch_labels) plot_arch_tree(al)
cat("[DONE]\n")

