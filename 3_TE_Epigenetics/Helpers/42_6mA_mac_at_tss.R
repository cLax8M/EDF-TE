#!/usr/bin/env Rscript
# 42_6mA_mac_at_tss.R
# Count TE copies that have a 6mA island (MAC) at their TSS region
#   TSS window : [TSS - 50 bp, TSS + 150 bp]  (strand-aware)
#   Min MAC len: 20 bp
# Species: Phycomyces blaslakeleeanus (Phybl2)
#          Rhizopus microsporus (Rhimi59)
# MAC union across all available conditions per species.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
  library(Cairo)
  library(patchwork)
})

setDTthreads(4)

# ── parameters ──────────────────────────────────────────────────────────────
OUTDIR         <- "42_6mA_mac_at_tss"
MIN_ISLAND_LEN <- 20L    # minimum MAC length (bp)
TSS_UP         <- 50L    # upstream of TSS
TSS_DOWN       <- 150L   # downstream of TSS

ISLAND_DIR <- paste0("/legserv/Temp/Carlos/EDF-TE/data/",
                     "Genome_Assemblies_Conditions/Epigenetics/6mA")

PHY_TE_PATH  <- paste0("/legserv/Temp/Carlos/EDF-TE/results/TE_dynamics/",
                        "20_chiapello_extension/data/TE_Phybl2.bed")
RHI_TE_PATH  <- paste0("/legserv/Temp/Carlos/EDF-TE/data/",
                        "TEs Mucoromycotina/MU_Rhimi59.filteredRepeats.bed")

CAT_COLS <- c(LTR="#e76254", LINE="#1b7837", DNA="#4575b4",
              SINE="#f4a74d", "Rolling Circle"="#9970ab", Unclassified="#888888")

# ── classify_sf ──────────────────────────────────────────────────────────────
classify_sf <- function(x) {
  ifelse(grepl("Gypsy|Copia|DIRS|BelPao|^LTR|LTR/", x, ignore.case=TRUE), "LTR",
  ifelse(grepl(paste0("^LINE|LINE/|^Retroposon|^L1|^L2|^CR1|^RTE|^I-|^R1$|^R2",
                      "|^Rex|Jockey|^Tad1|^SLACS|^TRAS|^Hero"),
               x, ignore.case=TRUE), "LINE",
  ifelse(grepl("^SINE|SINE/|^Alu/|^Alu$|^B1$|^B2$|^B4$|^ID$|^MIR",
               x, ignore.case=TRUE), "SINE",
  ifelse(grepl(paste0("^DNA|DNA/|TcMar|CMC|hAT|Harbinger|Maverick|Kolobok|",
                      "MuLE|Transib|Polinton|Crypton|Zisupton|ISRm|Merlin|Academ"),
               x, ignore.case=TRUE), "DNA",
  ifelse(grepl("^RC|RC/|Helitron", x, ignore.case=TRUE), "Rolling Circle",
                                                          "Unclassified")))))
}

# ── helper: load 6mA island files and filter by length ──────────────────────
load_islands <- function(paths, min_len=MIN_ISLAND_LEN) {
  lst <- lapply(paths, fread, sep="\t", header=TRUE)
  dt  <- rbindlist(lst, fill=TRUE)
  dt  <- dt[length >= min_len, .(chrom, start, end, length, num_mods, id)]
  # Deduplicate identical regions (same MAC called in multiple conditions)
  unique(dt, by=c("chrom","start","end"))
}

# ── core analysis per species ─────────────────────────────────────────────────
analyze_species <- function(label, te_path, island_paths) {
  message("\n=== ", label, " ===")

  # Load TEs
  te <- fread(te_path,
              col.names=c("chrom","start","end","name","score","strand"),
              header=FALSE)
  te[, te_idx    := .I]
  te[, superfam  := classify_sf(name)]
  message("  TE copies: ", nrow(te))

  # Load 6mA islands (union across all conditions)
  islands <- load_islands(island_paths)
  message("  MAC files: ", paste(basename(island_paths), collapse=", "))
  message("  MACs (len>=", MIN_ISLAND_LEN, "): ", nrow(islands))

  # ── Strand-aware TSS window ─────────────────────────────────────────────
  # + strand: [start - TSS_UP,  start + TSS_DOWN]
  # - strand: [end   - TSS_DOWN, end  + TSS_UP  ]
  # (BED is 0-based half-open; keep in BED coords, then convert to 1-based for GRanges)
  tss_start <- ifelse(te$strand == "+",
                      pmax(0L, te$start - TSS_UP),
                      pmax(0L, te$end   - TSS_DOWN))
  tss_end   <- ifelse(te$strand == "+",
                      te$start + TSS_DOWN,
                      te$end   + TSS_UP)

  # Build GRanges (BED 0-based → GRanges 1-based: start+1)
  gr_tss <- GRanges(
    seqnames = te$chrom,
    ranges   = IRanges(start=tss_start + 1L, end=tss_end),
    te_idx   = te$te_idx
  )
  gr_mac <- GRanges(
    seqnames = islands$chrom,
    ranges   = IRanges(start=islands$start + 1L, end=islands$end)
  )

  # Overlap: which TSS windows contain a MAC?
  hits <- findOverlaps(gr_tss, gr_mac)
  te_with_mac <- unique(queryHits(hits))
  te[, has_mac_at_tss := te_idx %in% te_with_mac]

  message("  TE copies with MAC at TSS: ", sum(te$has_mac_at_tss),
          " / ", nrow(te),
          sprintf("  (%.1f%%)", 100*mean(te$has_mac_at_tss)))

  te[, species := label]
  te
}

# ── Run analysis ─────────────────────────────────────────────────────────────
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

phy_island_files <- sort(list.files(ISLAND_DIR,
                                    pattern="Phybl2.*islands.*\\.tab$",
                                    full.names=TRUE))
rhi_island_files <- sort(list.files(ISLAND_DIR,
                                    pattern="Rhimi59.*islands.*\\.tab$",
                                    full.names=TRUE))

phy <- analyze_species("Phycomyces", PHY_TE_PATH,  phy_island_files)
rhi <- analyze_species("Rhizopus",   RHI_TE_PATH,  rhi_island_files)

all_te <- rbind(phy, rhi)

# ── Summary tables ───────────────────────────────────────────────────────────
# Per-species overall
overall <- all_te[, .(
  n_total    = .N,
  n_with_mac = sum(has_mac_at_tss),
  pct        = round(100 * mean(has_mac_at_tss), 1)
), by=species]
message("\n--- Overall ---")
print(overall)

# Per species × superfamily
by_sf <- all_te[, .(
  n_total    = .N,
  n_with_mac = sum(has_mac_at_tss),
  pct        = round(100 * mean(has_mac_at_tss), 1)
), by=.(species, superfam)][order(species, -n_with_mac)]
message("\n--- By superfamily ---")
print(by_sf)

# Save tables
fwrite(overall, file.path(OUTDIR, "42_overall.tsv"), sep="\t")
fwrite(by_sf,   file.path(OUTDIR, "42_by_superfam.tsv"), sep="\t")
fwrite(all_te[has_mac_at_tss == TRUE,
              .(species, chrom, start, end, name, strand, superfam)],
       file.path(OUTDIR, "42_te_with_mac_at_tss.tsv"), sep="\t")

# ── Plot ─────────────────────────────────────────────────────────────────────
plot_dt <- copy(by_sf)[superfam != "Unclassified"]
plot_dt[, superfam := factor(superfam,
                              levels=c("LTR","LINE","DNA","SINE","Rolling Circle"))]
plot_dt[, label := sprintf("%d/%d", n_with_mac, n_total)]
plot_dt[, species := factor(species, levels=c("Phycomyces","Rhizopus"))]

p_bar <- ggplot(plot_dt, aes(x=superfam, y=pct, fill=superfam, label=label)) +
  geom_col(width=0.7, colour="white", linewidth=0.3) +
  geom_text(vjust=-0.25, size=2.8, colour="grey30") +
  scale_fill_manual(values=CAT_COLS, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.15))) +
  facet_wrap(~species, ncol=2) +
  labs(
    x    = NULL,
    y    = "% TE copies with MAC at TSS",
    title= sprintf("6mA islands (MAC) within -%d / +%d bp of TE TSS  (min MAC %d bp)",
                   TSS_UP, TSS_DOWN, MIN_ISLAND_LEN)
  ) +
  theme_bw(base_size=11) +
  theme(
    axis.text.x       = element_text(angle=30, hjust=1),
    strip.background  = element_rect(fill="grey92"),
    panel.grid.major.x= element_blank()
  )

CairoPDF(file.path(OUTDIR, "42_mac_at_tss.pdf"), width=9, height=5)
print(p_bar)
dev.off()
message("\nDone. Results in: ", OUTDIR)
