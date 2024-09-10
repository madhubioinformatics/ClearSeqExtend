# R/seurat_extend_integration.R

#' SeuratExtend Integration: Normalization and Dimensionality Reduction
#'
#' This function performs normalization, PCA, and UMAP using SeuratExtend methods.
#'
#' @param seurat_obj A Seurat object containing scRNA-seq data.
#' @return A Seurat object with enhanced analysis using SeuratExtend.
#' @export
integrate_seurat_extend <- function(seurat_obj) {
    # Normalize data using SeuratExtend
    seurat_obj <- SeuratExtend::EnhancedNormalizeData(seurat_obj)
    
    # Identify variable features using SeuratExtend
    seurat_obj <- SeuratExtend::EnhancedFindVariableFeatures(seurat_obj)
    
    # Perform PCA and UMAP using SeuratExtend methods
    seurat_obj <- SeuratExtend::RunExtendedPCA(seurat_obj, features = VariableFeatures(seurat_obj))
    seurat_obj <- SeuratExtend::RunExtendedUMAP(seurat_obj, dims = 1:10)
    
    return(seurat_obj)
}
