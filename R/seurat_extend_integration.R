# R/seurat_extend_integration.R

#' Enhanced Normalization Using SeuratExtend
#'
#' Applies enhanced normalization on the Seurat object using SeuratExtend.
#' @param seurat_obj A Seurat object containing scRNA-seq data.
#' @return A normalized Seurat object.
#' @export
enhanced_normalization <- function(seurat_obj) {
  # Use SeuratExtend's normalization method (e.g., replacing NormalizeData with an advanced method)
  seurat_obj <- SeuratExtend::EnhancedNormalizeData(seurat_obj)
  
  return(seurat_obj)
}

#' Find Variable Features Using SeuratExtend
#'
#' Identifies highly variable features in the Seurat object using SeuratExtend.
#' @param seurat_obj A Seurat object containing scRNA-seq data.
#' @return A Seurat object with variable features identified.
#' @export
find_variable_features_extend <- function(seurat_obj) {
  # Use SeuratExtend's method for finding variable features
  seurat_obj <- SeuratExtend::EnhancedFindVariableFeatures(seurat_obj)
  
  return(seurat_obj)
}


#' Run PCA Using SeuratExtend
#'
#' Performs dimensionality reduction via PCA using SeuratExtend's advanced PCA method.
#' @param seurat_obj A Seurat object containing scRNA-seq data.
#' @param features A vector of features to use for PCA.
#' @return A Seurat object with PCA applied.
#' @export
run_pca_extend <- function(seurat_obj, features = NULL) {
  if (is.null(features)) {
    features <- VariableFeatures(seurat_obj)
  }
  
  # Perform PCA using SeuratExtend's advanced PCA method
  seurat_obj <- SeuratExtend::RunExtendedPCA(seurat_obj, features = features)
  
  return(seurat_obj)
}

#' Run UMAP Using SeuratExtend
#'
#' Performs UMAP dimensionality reduction using SeuratExtend's advanced UMAP method.
#' @param seurat_obj A Seurat object containing scRNA-seq data.
#' @param dims Number of dimensions to use for UMAP.
#' @return A Seurat object with UMAP applied.
#' @export
run_umap_extend <- function(seurat_obj, dims = 1:10) {
  # Perform UMAP using SeuratExtend's advanced UMAP method
  seurat_obj <- SeuratExtend::RunExtendedUMAP(seurat_obj, dims = dims)
  
  return(seurat_obj)
}

#' Perform Clustering Using SeuratExtend
#'
#' Applies clustering on the Seurat object using SeuratExtend's advanced clustering methods.
#' @param seurat_obj A Seurat object containing scRNA-seq data.
#' @param resolution Resolution parameter for clustering (default: 0.5).
#' @return A Seurat object with clusters identified.
#' @export
cluster_cells_extend <- function(seurat_obj, resolution = 0.5) {
  # Use SeuratExtend's advanced clustering method
  seurat_obj <- SeuratExtend::EnhancedFindClusters(seurat_obj, resolution = resolution)
  
  return(seurat_obj)
}

#' Generate UMAP Plot Using SeuratExtend
#'
#' Generates a UMAP plot using SeuratExtend's advanced visualization method.
#' @param seurat_obj A Seurat object containing scRNA-seq data.
#' @param group.by Group cells by specified metadata (default: "seurat_clusters").
#' @return A UMAP plot.
#' @export
umap_plot_extend <- function(seurat_obj, group.by = "seurat_clusters") {
  # Generate a UMAP plot using SeuratExtend
  plot <- SeuratExtend::EnhancedDimPlot(seurat_obj, reduction = "umap", group.by = group.by, label = TRUE)
  
  return(plot)
}

