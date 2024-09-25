ClearSeqExtend: An Enhanced Toolkit for Ambient RNA Correction and scRNA-seq Analysis

ClearSeqExtend is an R package designed for robust ambient RNA correction in single-cell RNA-seq (scRNA-seq) data, integrating advanced features from SeuratExtend for comprehensive analysis. This package allows users to correct for ambient RNA contamination, perform clustering, and generate high-quality visualizations using SeuratExtend functionalities.

Installation
To install the ClearSeqExtend package from GitHub, follow these steps:

r
Copy code
# Install devtools if not already installed
install.packages("devtools")

# Install ClearSeqExtend from GitHub
devtools::install_github("madhubioinformatics/ClearSeqExtend")

if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}
remotes::install_github("huayc09/SeuratExtend")
Loading the Package
After installation, load ClearSeqExtend, Seurat, and SeuratExtend packages:

r
Copy code
library(Seurat)
library(SeuratExtend)
library(ClearSeqExtend)
Usage Example
Below is a step-by-step guide on how to use ClearSeqExtend for ambient RNA correction and advanced scRNA-seq analysis.

1. Load the Dataset and Create a Seurat Object

r

# Load your dataset (e.g., PBMC dataset from 10X Genomics)
pbmc_data <- Read10X_h5("path_to_data/10k_PBMC_3p_nextgem_Chromium_X_raw_feature_bc_matrix.h5")

# Create a Seurat object from the data
seurat_obj <- CreateSeuratObject(counts = pbmc_data, project = "PBMC")
2. Perform Ambient RNA Correction

Correct for ambient RNA contamination using the clearseq_with_extend function:

r

# Perform ambient RNA correction with a UMI threshold of 10
seurat_obj <- clearseq_with_extend(seurat_obj, threshold = 10)
3. Normalize the Data and Identify Variable Features

r

# Normalize the data using SeuratExtend's advanced normalization
seurat_obj <- enhanced_normalization(seurat_obj)

# Identify highly variable features in the dataset
seurat_obj <- find_variable_features_extend(seurat_obj)
4. Perform PCA and UMAP for Dimensionality Reduction

r

# Run PCA for dimensionality reduction using SeuratExtend
seurat_obj <- run_pca_extend(seurat_obj)

# Run UMAP for further dimensionality reduction
seurat_obj <- run_umap_extend(seurat_obj, dims = 1:10)
5. Cluster the Cells and Visualize with UMAP

r

# Perform clustering on the Seurat object
seurat_obj <- cluster_cells_extend(seurat_obj, resolution = 0.5)

# Generate UMAP plot to visualize clusters
umap_plot <- umap_plot_extend(seurat_obj)

# Display the UMAP plot
print(umap_plot)
