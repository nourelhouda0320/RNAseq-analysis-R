############################################################
# RNA-seq Analysis in R
# Mouse Mammary Gland RNA-seq Dataset
#
# Main packages:
# - edgeR
# - limma
# - Glimma
# - gplots
# - RColorBrewer
# - org.Mm.eg.db
############################################################


############################
# 1. Load packages
############################

library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(org.Mm.eg.db)


############################
# 2. Create output folders
############################

if (!dir.exists("results")) {
  dir.create("results")
}

if (!dir.exists("results/figures")) {
  dir.create("results/figures", recursive = TRUE)
}


############################
# 3. Import sample metadata
############################

sampleinfo <- read.delim(
  "data/SampleInfo_Corrected.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Convert categorical variables to factors
sampleinfo$CellType <- factor(sampleinfo$CellType)
sampleinfo$Status <- factor(sampleinfo$Status)


############################
# 4. Import RNA-seq count data
############################

counts <- read.delim(
  "data/GSE60450_Lactation-GenewiseCounts.txt",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Keep gene identifiers as row names
rownames(counts) <- counts$EntrezGeneID

# Remove gene annotation columns
counts <- counts[, -(1:2)]

# Standardize sample names
colnames(counts) <- substr(colnames(counts), 1, 7)


############################
# 5. Check sample names
############################

cat("\nSample names in count matrix:\n")
print(colnames(counts))

cat("\nSample names in metadata:\n")
print(sampleinfo$SampleName)

if (!all(colnames(counts) == sampleinfo$SampleName)) {
  stop(
    "ERROR: Sample names in the count matrix do not match the sample metadata."
  )
}

cat("\nSample names successfully matched.\n")


############################
# 6. Calculate CPM
############################

cpm_values <- cpm(counts)


############################
# 7. Expression filtering
############################

# Keep genes with CPM > 0.5 in at least 2 samples
keep <- rowSums(cpm_values > 0.5) >= 2

filtered_counts <- counts[keep, ]

cat("\nNumber of genes before filtering:", nrow(counts), "\n")
cat("Number of genes after filtering:", nrow(filtered_counts), "\n")

# Save filtered count matrix
write.csv(
  filtered_counts,
  "results/filtered_counts.csv"
)


############################
# 8. CPM versus raw counts
############################

pdf(
  "results/figures/CPM_vs_raw_counts.pdf",
  width = 8,
  height = 6
)

plot(
  counts[, 1],
  cpm_values[, 1],
  log = "xy",
  pch = 16,
  cex = 0.4,
  xlab = "Raw counts",
  ylab = "CPM",
  main = "CPM versus Raw Counts"
)

dev.off()


############################
# 9. Create DGEList object
############################

dge <- DGEList(
  counts = filtered_counts,
  samples = sampleinfo
)


############################
# 10. Library size visualization
############################

pdf(
  "results/figures/library_size.pdf",
  width = 8,
  height = 6
)

barplot(
  dge$samples$lib.size / 1e6,
  names.arg = sampleinfo$SampleName,
  las = 2,
  ylab = "Library size (millions)",
  main = "Library Sizes"
)

dev.off()


############################
# 11. Log2-CPM calculation
############################

logCPM <- cpm(
  dge,
  log = TRUE,
  prior.count = 2
)

# Save logCPM matrix
write.csv(
  logCPM,
  "results/logCPM_matrix.csv"
)


############################
# 12. Log2-CPM distribution
############################

pdf(
  "results/figures/logCPM_boxplot.pdf",
  width = 10,
  height = 6
)

boxplot(
  logCPM,
  las = 2,
  main = "Log2-CPM Distribution",
  ylab = "Log2-CPM"
)

dev.off()


############################
# 13. MDS analysis
############################

# Colors based on cell type
celltype_colors <- as.numeric(sampleinfo$CellType)

pdf(
  "results/figures/MDS_CellType.pdf",
  width = 8,
  height = 6
)

plotMDS(
  dge,
  col = celltype_colors,
  main = "MDS Plot - Cell Type"
)

legend(
  "topright",
  legend = levels(sampleinfo$CellType),
  col = seq_along(levels(sampleinfo$CellType)),
  pch = 16,
  title = "Cell Type"
)

dev.off()


############################
# 14. MDS analysis by status
############################

status_colors <- as.numeric(sampleinfo$Status)

pdf(
  "results/figures/MDS_Status.pdf",
  width = 8,
  height = 6
)

plotMDS(
  dge,
  col = status_colors,
  main = "MDS Plot - Biological Status"
)

legend(
  "topright",
  legend = levels(sampleinfo$Status),
  col = seq_along(levels(sampleinfo$Status)),
  pch = 16,
  title = "Status"
)

dev.off()


############################
# 15. Identify highly variable genes
############################

gene_variance <- apply(
  logCPM,
  1,
  var
)

# Select the 500 most variable genes
top_n <- min(500, length(gene_variance))

top_variable_genes <- names(
  sort(
    gene_variance,
    decreasing = TRUE
  )[1:top_n]
)

high_var_matrix <- logCPM[
  top_variable_genes,
  ]


############################
# 16. Save highly variable genes
############################

high_var_results <- data.frame(
  EntrezGeneID = top_variable_genes,
  Variance = gene_variance[top_variable_genes]
)

write.csv(
  high_var_results,
  "results/highly_variable_genes.csv",
  row.names = FALSE
)


############################
# 17. Heatmap
############################

pdf(
  "results/figures/High_var_genes_heatmap.pdf",
  width = 10,
  height = 12
)

heatmap.2(
  high_var_matrix,
  scale = "row",
  trace = "none",
  dendrogram = "both",
  key = TRUE,
  density.info = "none",
  col = colorRampPalette(
    rev(brewer.pal(9, "RdBu"))
  )(100),
  margins = c(8, 8),
  main = "Top Highly Variable Genes"
)

dev.off()


############################
# 18. Save heatmap as PNG
############################

png(
  "results/figures/High_var_genes_heatmap.png",
  width = 1200,
  height = 1400,
  res = 150
)

heatmap.2(
  high_var_matrix,
  scale = "row",
  trace = "none",
  dendrogram = "both",
  key = TRUE,
  density.info = "none",
  col = colorRampPalette(
    rev(brewer.pal(9, "RdBu"))
  )(100),
  margins = c(8, 8),
  main = "Top Highly Variable Genes"
)

dev.off()


############################
# 19. Analysis summary
############################

cat("\n")
cat("============================================\n")
cat("RNA-seq analysis completed successfully.\n")
cat("============================================\n")
cat("Genes before filtering :", nrow(counts), "\n")
cat("Genes after filtering  :", nrow(filtered_counts), "\n")
cat("Highly variable genes :", length(top_variable_genes), "\n")
cat("\nResults saved in: results/\n")
cat("Figures saved in: results/figures/\n")
cat("============================================\n")
