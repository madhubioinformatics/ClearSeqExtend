# ClearSeqExtend

ClearSeqExtend is an R package that integrates SeuratExtend for enhanced single-cell RNA sequencing (scRNA-seq) analysis, providing ambient RNA correction, PCA, UMAP, clustering, and RNA velocity analysis.

1. Introduction

Welcome to the ClearSeqExtend tutorial! This package integrates ambient RNA correction, clustering, RNA velocity, and more for enhanced scRNA-seq analysis using Seurat and SeuratExtend. In this tutorial, we will walk you through the installation, usage, and practical examples to get started with ClearSeqExtend.

2. Installation

To install the ClearSeqExtend package, follow these steps:

```r

# Install the devtools package if you don't have it
install.packages("devtools")

# Install ClearSeqExtend from GitHub
devtools::install_github("madhubioinformatics/ClearSeqExtend")

3. Loading the Package

After installation, load the ClearSeqExtend package along with other dependencies:

```r

# Load the required libraries
library(Seurat)
library(SeuratExtend)
library(ClearSeqExtend)

4. Example Dataset

In this example, we will use the 10X Genomics PBMC dataset. You can download the dataset using the following command:

```r

# Download the dataset (if not already downloaded)
download.file("https://cf.10xgenomics.com/samples/cell-exp/6.1.2/10k_PBMC_3p_nextgem_Chromium_X/10k_PBMC_3p_nextgem_Chromium_X_raw_feature_bc_matrix.h5", 
              destfile = "10k_PBMC_3p_nextgem_Chromium_X_raw_feature_bc_matrix.h5")

# Load the dataset into Seurat
pbmc_data <- Read10X_h5("10k_PBMC_3p_nextgem_Chromium_X_raw_feature_bc_matrix.h5")

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = pbmc_data, project = "PBMC")

5. Ambient RNA Correction and SeuratExtend Integration

Next, we will perform ambient RNA correction using ClearSeqExtend and integrate SeuratExtend functionalities for normalization and dimensionality reduction.

```r

# Perform ambient RNA correction and integrate SeuratExtend
seurat_obj <- clearseq_with_extend(seurat_obj, threshold = 10)

# Check the contents of the Seurat object
seurat_obj
This function performs ambient RNA correction, normalization, PCA, and UMAP.

6. Clustering and Visualization

Once the ambient RNA has been corrected and the data normalized, you can perform clustering and visualize the UMAP plot.

```r

# Find neighbors and clusters
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Visualize the clusters using UMAP
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

7. RNA Velocity (Optional)

If you have a .loom file for RNA velocity analysis, you can integrate RNA velocity using ClearSeqExtend.

```r

# RNA velocity using loom file (replace with your loom path)
loom_path <- "path/to/your/loom_file.loom"
velocity_obj <- run_velocity(seurat_obj, loom_path = loom_path)

8. Export Results

You can save the processed Seurat object or export the data in various formats such as .loom or .h5ad.

```r

# Save Seurat object
saveRDS(seurat_obj, file = "processed_pbmc_seurat.rds")

# Export to loom file
ExportToLoom(seurat_obj, filename = "pbmc_seurat.loom")

9. Conclusion

Congratulations! You have successfully run the ClearSeqExtend pipeline for ambient RNA correction, clustering, and RNA velocity analysis. This tutorial is a starting point for more advanced analyses.



