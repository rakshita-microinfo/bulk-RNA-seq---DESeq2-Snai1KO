# Bulk RNA-seq Differential Expression Analysis — Snai1 Knockout in Breast Cancer

This project investigates the transcriptional role of **Snai1** in triple-negative breast cancer cells. Using publicly available RNA-seq data, I performed differential expression analysis comparing **Snai1 knockout** vs **wild type** HS578T cells to identify genes regulated by Snai1 and their potential role in cancer invasion and epithelial-mesenchymal transition (EMT).

## Dataset

The data used in this project was retrieved from the 
**EBI Expression Atlas** (https://www.ebi.ac.uk/gxa/home).
- **Accession:** E-MTAB-5244
- **Organism:** *Homo sapiens*
- **Cell line:** HS578T (triple-negative breast cancer)
- **Samples:** 6 total — 3 wild type controls, 3 Snai1 knockout

## Tools and Methods

**Language:** R (version 4.5.0)
**Packages:**
- `DESeq2` — differential expression analysis
- `ggplot2` — data visualization
- `pheatmap` — heatmap generation
- `EnhancedVolcano` — volcano plot generation
- `org.Hs.eg.db` — Ensembl ID to gene symbol annotation
- `AnnotationDbi` — annotation database interface

**Analysis Workflow:**
1. **Data Loading** — Loaded raw count matrix (genes × samples) and 
   experimental design file (sample metadata) from EBI Expression Atlas
2. **Data Cleaning** — Removed secondary gene identifier column from 
   count matrix and extracted relevant condition column from metadata
3. **Quality Verification** — Confirmed sample names matched between 
   count matrix and metadata using `stopifnot()`
4. **Exploratory Data Analysis** — Checked matrix dimensions, sample 
   conditions, and absence of negative values
5. **DESeq2 Analysis** — Built DESeq2 object, set wild type as reference, 
   filtered low count genes (minimum 10 reads), and ran differential 
   expression analysis
6. **Results Extraction** — Filtered significant DEGs (padj < 0.05, 
   |log2FoldChange| > 1), mapped Ensembl IDs to gene symbols, 
   and exported results to CSV
7. **Visualization** — Generated PCA plot, volcano plot, and heatmap 
   of top 30 DEGs
8. **Biological Interpretation** — Evaluated top DEGs in the context 
   of Snai1-driven EMT and cancer invasion

## Key Results

In total, **1,323 genes** were significantly differentially expressed upon Snai1 
knockout, of which **645 were upregulated** and **678 were downregulated**.

PCA showed that PC1 explained **86% of total variance**, indicating that the KO vs 
wild type difference is the dominant source of variation in the dataset. PC2 explained 
only 6%, confirming minimal variation between replicates and high reproducibility.

Top upregulated genes included **THY1** and **CXCL12**, both implicated in EMT and 
cancer cell signalling. Top downregulated genes included **PDPN** and **PRKCQ**, 
associated with cancer invasion and cell motility. The results suggest that Snai1 
regulates a complex transcriptional network in HS578T cells beyond simple EMT 
activation, with both pro- and anti-invasive gene programs responding to its loss

## Output Files

| `Snai1KO_DESeq2_analysis.R` | Full analysis script |
| `DESeq2_Snai1KO_results.csv` | Complete DESeq2 results table |
| `PCA_Snai1KO.png` | PCA plot — sample quality control |
| `Volcano_Snai1KO.png` | Volcano plot — DEG overview |
| `Heatmap_Snai1KO.png` | Heatmap — top 30 DEGs |


## Author

**Rakshita Govind**
Computational Biologist | Environmental Microbiomics | Bioinformatics
[LinkedIn](https://linkedin.com/in/rakshita-govind) | [GitHub](https://github.com/rakshita-microinfo)

