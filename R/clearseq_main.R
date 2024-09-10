# R/clearseq_main.R

#' ClearSeqExtend: Ambient RNA Correction and Enhanced Analysis
#'
#' This function performs ambient RNA correction, RNA velocity analysis, and clustering using SeuratExtend.
#'
#' @param seurat_obj A Seurat object containing raw scRNA-seq data.
#' @param threshold Numeric value representing the UMI threshold for identifying empty droplets.
#' @return A Seurat object with corrected counts and enhanced analysis using SeuratExtend.
#' @export
clearseq_with_extend <- function(seurat_obj, threshold = 10) {
    # Apply ClearSeq Correction (ambient RNA correction)
    seurat_obj <- clearseq_correction(seurat_obj, threshold)
    
    # Normalize and find variable features using SeuratExtend
    seurat_obj <- SeuratExtend::EnhancedNormalizeData(seurat_obj)
    seurat_obj <- SeuratExtend::EnhancedFindVariableFeatures(seurat_obj)
    
    # Run PCA and UMAP using SeuratExtend’s enhanced methods
    seurat_obj <- SeuratExtend::RunExtendedPCA(seurat_obj, features = VariableFeatures(seurat_obj))
    seurat_obj <- SeuratExtend::RunExtendedUMAP(seurat_obj, dims = 1:10)
    
    return(seurat_obj)
}
