# PCA of RNA-seq data
R workflow for performing Principal Component Analysis (PCA) on CPM (Counts Per Million) expression data from RNA-seq experiments.
This script visualizes sample clustering across multiple conditions and timepoints, helping assess sample quality, batch effects, and biological grouping before downstream analysis.

# Overview
This R Markdown script takes a pre-computed CPM matrix, applies log2 transformation and scaling, runs PCA, and produces a publication-ready PCA plot saved as a high-resolution TIFF file.

Requirements \
Input File \
CPM Matrix \
A tab-delimited file containing CPM-normalized gene expression values with:

Rows — Gene IDs (non-zero genes only) \
Columns — Samples \
First 5 annotation columns are removed before analysis

# Pipeline Steps
CPM Matrix (non-zero genes) > Remove annotation columns (columns 1–5) > Log2 transformation: log2(CPM + 1) > Scaling across samples: scale() > PCA using prcomp() on transposed matrix > Extract PC scores + % variance explained > Plot PC1 vs PC2 with condition labels (ggplot2) > Export as high-resolution TIFF (1200 dpi)

# Notes
Only non-zero genes are included in the CPM matrix to reduce noise and improve PCA resolution. \
The +1 pseudocount before log transformation is critical to avoid undefined values for zero-count genes. \
Scaling with scale() before PCA ensures no single highly-expressed gene dominates the principal components. \
This PCA is intended as a quality control and exploratory step — run it before DESeq2 differential expression analysis to confirm expected sample grouping. \
If samples within the same condition do not cluster together, this may indicate batch effects or labeling issues.

# R code to perform PCA

```{r}
#Load necessary packages
library(readxl)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(dplyr)
```

```{r}
PCA_B157 <- read.table(
  "path to input.txt",
  header = TRUE,
  row.names = "GeneId",
  sep = "\t",
  check.names = FALSE
)

# Remove annotation columns (adjust if needed)
PCA_B157 <- PCA_B157[, -c(1:5)]
```

```{r}
PCA_sample <- data.frame(
  condition = c("M16", "M16", "M16", "M24", "M24", "M24", "MY", "MY", "MY",
                "W16", "W16", "W16", "W24", "W24", "W24", "WY", "WY", "WY"),
  row.names = colnames(PCA_B157)
)

PCA_cpm_log <- log2(PCA_B157 + 1)  # Log-transform to handle continuous values
PCA_cpm_scaled <- scale(PCA_cpm_log)  # Scale across samples
```

```{r}
pca_resultCPM <- prcomp(t(PCA_cpm_scaled), center = TRUE, scale. = TRUE)  # Transpose for PCA

#Extract PCA scores and calculate explained variance
pca_scores <- as.data.frame(pca_resultCPM$x)
percentVar_CPM <- round(100 * (pca_resultCPM$sdev^2 / sum(pca_resultCPM$sdev^2)), 4)

#Add condition labels to PCA data frame
pca_scores$condition <- PCA_sample$condition
```

```{r}
custom_colors <- c("M16" = "#C2A68C",
                   "M24" = "#FDB7EA",
                   "MY" = "#97A87A",
                   "W16" = "#FFC107",
                   "W24" = "#5A9CB5",
                   "WY" = "#D25353")

plot_PCA_CPM <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 10) +
  labs(
    title = "PCA of CPM",
    x = paste0("PC1 (", percentVar_CPM[1], "%)"),
    y = paste0("PC2 (", percentVar_CPM[2], "%)")
  ) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 24, face = "bold")
  )

print(plot_PCA_CPM)
```

```{r}
ggsave(filename = "PCA_CPM.tiff", path = "path to output folder", plot = plot_PCA_CPM, width = 12, height = 10, units = "in", dpi = 1200, device = "tiff", compression = "lzw")
```

