#!/usr/bin/env Rscript
# Test: do MU_ species with DNMT1 + DIM2 have higher genomic TE content?

library(ggplot2)

# Paths
te_file   <- "/legserv/Temp/Carlos/EDF-TE/results/Section2_TElandscape/D01_te_genomics_merged_provisional.tsv"
dnmt_file <- "/legserv/Temp/Carlos/EDF-TE/results/TE_dynamics/dnmt_machinery/data/dnmt_matrix_final.tsv"
out_dir   <- "/legserv/Temp/Carlos/EDF-TE/results/TE_dynamics/dnmt_machinery/"

# Read data
te   <- read.table(te_file,   header=TRUE, sep="\t", stringsAsFactors=FALSE)
dnmt <- read.table(dnmt_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Filter to MU_ clade
te_mu   <- te[grepl("^MU_", te$species), ]
dnmt_mu <- dnmt[grepl("^MU_", dnmt$species), ]

cat("MU species in TE file:  ", nrow(te_mu), "\n")
cat("MU species in DNMT file:", nrow(dnmt_mu), "\n")

# Merge
df <- merge(te_mu[, c("species","pct_TE_total","genome_size_TE")],
            dnmt_mu[, c("species","DNMT1_n","DIM2_n")],
            by="species", all=FALSE)

cat("MU species in both files:", nrow(df), "\n\n")

# Classify: has_both = DNMT1_n >= 1 AND DIM2_n >= 1
df$group <- ifelse(df$DNMT1_n >= 1 & df$DIM2_n >= 1, "DNMT1+DIM2", "Other MU")
df$group <- factor(df$group, levels=c("Other MU","DNMT1+DIM2"))

# Summary statistics
cat("=== Group sizes ===\n")
print(table(df$group))

cat("\n=== Median % TE total by group ===\n")
print(tapply(df$pct_TE_total, df$group, median))

cat("\n=== Mean % TE total by group ===\n")
print(tapply(df$pct_TE_total, df$group, mean))

cat("\n=== Range by group ===\n")
for (g in levels(df$group)) {
  cat(g, ": min=", round(min(df[df$group==g,"pct_TE_total"]),2),
         " max=", round(max(df[df$group==g,"pct_TE_total"]),2),"\n")
}

cat("\n=== DNMT1+DIM2 species ===\n")
d <- df[df$group=="DNMT1+DIM2", c("species","pct_TE_total","DNMT1_n","DIM2_n")]
d <- d[order(-d$pct_TE_total), ]
print(d)

# Wilcoxon rank-sum test (one-sided: DNMT1+DIM2 > Other)
wt <- wilcox.test(pct_TE_total ~ group, data=df,
                  alternative="less",   # "less" because "Other MU" is level 1
                  exact=FALSE)
cat("\n=== Wilcoxon rank-sum test (DNMT1+DIM2 > Other MU) ===\n")
cat("W =", wt$statistic, "\n")
cat("p-value =", wt$p.value, "\n")

# Two-sided for completeness
wt2 <- wilcox.test(pct_TE_total ~ group, data=df, exact=FALSE)
cat("two-sided p-value =", wt2$p.value, "\n")

# Effect size: rank-biserial correlation
n1 <- sum(df$group=="Other MU")
n2 <- sum(df$group=="DNMT1+DIM2")
r  <- 1 - (2*wt$statistic)/(n1*n2)
cat("Rank-biserial r =", round(r,3), "(effect size)\n")

# Plot
p <- ggplot(df, aes(x=group, y=pct_TE_total, fill=group)) +
  geom_boxplot(alpha=0.6, outlier.shape=NA) +
  geom_jitter(aes(color=group), width=0.15, size=2.5, alpha=0.85) +
  scale_fill_manual(values=c("Other MU"="#7fc8e8","DNMT1+DIM2"="#e87f7f")) +
  scale_color_manual(values=c("Other MU"="#3a8fb5","DNMT1+DIM2"="#c0392b")) +
  labs(title="Genomic TE content in Mucoromycotina",
       subtitle=paste0("Wilcoxon p = ", formatC(wt$p.value, format="e", digits=2),
                       " (one-sided),  n(Other MU)=", n1, ", n(DNMT1+DIM2)=", n2),
       x=NULL, y="% Genome covered by TEs",
       fill=NULL, color=NULL) +
  theme_classic(base_size=13) +
  theme(legend.position="none")

# Label DNMT1+DIM2 points
p <- p + ggrepel::geom_text_repel(
  data=df[df$group=="DNMT1+DIM2",],
  aes(label=species), size=3, color="#8B1A1A",
  box.padding=0.4, max.overlaps=20
)

out_png <- file.path(out_dir, "dnmt1dim2_vs_te_content_MU.png")
out_pdf <- file.path(out_dir, "dnmt1dim2_vs_te_content_MU.pdf")
ggsave(out_png, p, width=6, height=6, dpi=150)
ggsave(out_pdf, p, width=6, height=6)
cat("\nPlot saved to:", out_png, "\n")
cat("Plot saved to:", out_pdf, "\n")
