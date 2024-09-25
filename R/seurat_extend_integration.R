#' Enhanced Normalization Using SeuratExtend
#'
#' Applies enhanced normalization on the Seurat object using SeuratExtend.
#' @param seurat_obj A Seurat object containing scRNA-seq data.
#' @return A normalized Seurat object.
#' @export
enhanced_normalization <- function(seurat_obj) {
  # Apply SeuratExtend's normalization
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
  seurat_obj <- SeuratExtend::EnhancedFindVariableFeatures(seurat_obj)
  return(seurat_obj)
}
#' Run PCA Using SeuratExtend
#'
#' Performs PCA using SeuratExtend.
#' @param seurat_obj A Seurat object.
#' @return A Seurat object with PCA results.
#' @export
run_pca_extend <- function(seurat_obj) {
  seurat_obj <- SeuratExtend::RunExtendedPCA(seurat_obj)
  return(seurat_obj)
}
#' Run UMAP Using SeuratExtend
#'
#' Runs UMAP on the Seurat object using SeuratExtend's advanced UMAP method.
#' @param seurat_obj A Seurat object.
#' @return A Seurat object with UMAP results.
#' @export
run_umap_extend <- function(seurat_obj) {
  seurat_obj <- SeuratExtend::RunExtendedUMAP(seurat_obj)
  return(seurat_obj)
}
#' Cluster Cells Using SeuratExtend
#'
#' Performs clustering on the Seurat object using SeuratExtend's method.
#' @param seurat_obj A Seurat object.
#' @param resolution Resolution parameter for clustering.
#' @return A Seurat object with clusters identified.
#' @export
cluster_cells_extend <- function(seurat_obj, resolution = 0.5) {
  seurat_obj <- SeuratExtend::EnhancedFindClusters(seurat_obj, resolution = resolution)
  return(seurat_obj)
}
#' Generate UMAP Plot Using SeuratExtend
#'
#' Generates a UMAP plot for visualization using SeuratExtend.
#' @param seurat_obj A Seurat object.
#' @return A UMAP plot.
#' @export
umap_plot_extend <- function(seurat_obj) {
  plot <- SeuratExtend::EnhancedDimPlot(seurat_obj, reduction = "umap")
  return(plot)
}
