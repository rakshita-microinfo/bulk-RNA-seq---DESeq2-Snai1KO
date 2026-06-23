#------------------------------------------------------------------------------------------------------------
# 1. Script Header
# Differential Expression Analysis of Snai1 Knockout Breast Cancer Cells
# Dataset: E-MTAB-5244 (EBI Expression Atlas)
# Tool: DESeq2
# Author: Rakshita Govind
# Date: June 2026
# Description: Bulk RNA-seq analysis comparing Snai1 KO vs wild type in HS578 triple-negative cancer cells
#------------------------------------------------------------------------------------------------------------
# 2. Load Libraries
#-------------------------------------------------------------------------------------------------------------
# Note: Run this once in your R console to install packages before running:
# install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "ggplot2", "pheatmap", "EnhancedVolcano", "org.Hs.eg.db", "AnnotationDbi"))
  
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(AnnotationDbi)
#-----------------------------------------------------------------------------
# 3. LOAD DATA
#-----------------------------------------------------------------------------
# Set your working directory to the folder containing your data files
# setwd("path/to/your/folder")  # adjust this to your local path

# Load the count matrix
# Row names = Ensembl gene IDs, Columns = sample ERR accession numbers
counts <- read.delim("counts.tsv",
                     header    = TRUE,
                     sep       = "\t",
                     row.names = 1)

# Remove Gene.Name column (second identifier column, not needed)
counts <- counts[, -1]

# Load the experiment design file
metadata <- read.delim("experiment-design.tsv",
                       header    = TRUE,
                       sep       = "\t",
                       row.names = 1)

# Extract only the condition column we need
metadata_clean <- data.frame(
  condition = metadata$Factor.Value.genotype.,
  row.names = rownames(metadata)
)

# Verify sample names match between counts and metadata
stopifnot(all(colnames(counts) == rownames(metadata_clean)))

#--------------------------------------------------------------
# 4. Exploratory Data Analysis
#----------------------------------------------------------------

# Dimensions of count matrix
cat("Number of genes:", nrow(counts), "\n")
cat("Number of samples:", ncol(counts), "\n")
# Sample names
cat("Sample names:\n")
print(colnames(counts))

# Condition groups
cat("Sample conditions:\n")
print(metadata_clean)

# Check no negative values exist in count matrix
cat("Minimum count value (should be 0):", min(counts), "\n")
# -----------------------------------------------------------------------------
# 5. DESEQ2 ANALYSIS
# -----------------------------------------------------------------------------

# Build DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata_clean,
  design    = ~ condition
)

# Set wild type as reference level
dds$condition <- relevel(dds$condition, ref = "wild type genotype")

# Filter low count genes (removes noise)
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

# Run DESeq2
dds <- DESeq(dds)

# -----------------------------------------------------------------------------
# 6. EXTRACT RESULTS
# -----------------------------------------------------------------------------

# Extract results: Snai1 KO vs wild type
res <- results(dds,
               contrast = c("condition", "Snai1 knockout", "wild type genotype"),
               alpha    = 0.05)

# Sort by adjusted p-value (most significant first)
res <- res[order(res$padj), ]

# Summary of results
summary(res)

# Filter significant DEGs (padj < 0.05, fold change > 2x)
sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

cat("Total significant DEGs:", nrow(sig), "\n")
cat("Upregulated in Snai1 KO:", sum(sig$log2FoldChange > 0), "\n")
cat("Downregulated in Snai1 KO:", sum(sig$log2FoldChange < 0), "\n")

# Convert results to data frame and add gene symbols
res_df <- as.data.frame(res)

res_df$symbol <- mapIds(org.Hs.eg.db,
                        keys      = rownames(res_df),
                        column    = "SYMBOL",
                        keytype   = "ENSEMBL",
                        multiVals = "first")

# Save full results to CSV
write.csv(res_df,
          file      = "DESeq2_Snai1KO_results.csv",
          row.names = TRUE)

cat("Results saved to DESeq2_Snai1KO_results.csv\n")

# -----------------------------------------------------------------------------
# 7. VISUALIZATIONS
# -----------------------------------------------------------------------------

# Variance stabilizing transformation (required for PCA and heatmap)
vsd <- vst(dds, blind = FALSE)

# --- Plot 1: PCA ---
pca_plot <- plotPCA(vsd, intgroup = "condition") +
  theme_minimal() +
  ggtitle("PCA: Snai1 KO vs Wild Type — HS578T Breast Cancer Cells") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = name), vjust = -1, size = 3)

print(pca_plot)
ggsave("PCA_Snai1KO.png", plot = pca_plot, width = 6, height = 5, dpi = 300)
cat("PCA plot saved\n")

# --- Plot 2: Volcano Plot ---
volcano_plot <- EnhancedVolcano(as.data.frame(res),
                                lab            = res_df$symbol,
                                x              = "log2FoldChange",
                                y              = "padj",
                                title          = "Snai1 KO vs Wild Type",
                                subtitle       = "Triple-negative breast cancer cells — HS578T",
                                pCutoff        = 0.05,
                                FCcutoff       = 1,
                                pointSize      = 2,
                                labSize        = 3,
                                col            = c("grey70", "grey70", "grey70", "red3"),
                                legendPosition = "bottom")

print(volcano_plot)
ggsave("Volcano_Snai1KO.png", plot = volcano_plot, width = 8, height = 7, dpi = 300)
cat("Volcano plot saved\n")

# --- Plot 3: Heatmap ---
# Prepare data for heatmap
top30         <- rownames(res)[1:30]
mat           <- assay(vsd)[top30, ]
mat           <- mat - rowMeans(mat)
rownames(mat) <- res_df[top30, "symbol"]
anno          <- as.data.frame(colData(vsd)[, "condition", drop = FALSE])

# Draw and save heatmap
pheatmap(mat,
         annotation_col = anno,
         show_rownames  = TRUE,
         fontsize_row   = 8,
         color          = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main           = "Top 30 DEGs — Snai1 KO vs Wild Type")

# Export: use Export > Save as Image in RStudio Plots panel
# Save as: Heatmap_Snai1KO.png, width 800, height 900

# -----------------------------------------------------------------------------
# SESSION INFO
# -----------------------------------------------------------------------------
# Always include this at the end of any bioinformatics script
# It records exactly which package versions were used
# Essential for reproducibility

sessionInfo()
