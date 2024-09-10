# R/utils.R

#' Ambient RNA Correction with Bayesian Inference
#'
#' This function applies Bayesian inference to correct for ambient RNA contamination.
#'
#' @param data_matrix A matrix of raw RNA counts.
#' @param prior A numeric value representing the prior probability of contamination.
#' @return A corrected matrix with ambient RNA removed.
#' @export
ambient_correction_bayesian <- function(data_matrix, prior = 0.5) {
    ambient_profile <- colMeans(data_matrix) * prior
    corrected_matrix <- sweep(data_matrix, 2, ambient_profile, FUN = "-")
    corrected_matrix[corrected_matrix < 0] <- 0
    return(corrected_matrix)
}

#' ClearSeq Ambient RNA Correction
#'
#' Identifies empty droplets and removes ambient RNA contamination.
#'
#' @param seurat_obj A Seurat object.
#' @param threshold UMI threshold to identify empty droplets.
#' @return A Seurat object with ambient RNA-corrected counts.
#' @export
clearseq_correction <- function(seurat_obj, threshold = 10) {
    total_umi <- colSums(GetAssayData(seurat_obj, slot = "counts"))
    empty_droplets <- total_umi < threshold
    ambient_profile <- rowMeans(GetAssayData(seurat_obj, slot = "counts")[, empty_droplets])
    corrected_counts <- sweep(GetAssayData(seurat_obj, slot = "counts"), 1, ambient_profile, FUN = "-")
    corrected_counts[corrected_counts < 0] <- 0  # No negative values allowed
    
    # Add corrected counts to a new assay in the Seurat object
    corrected_assay <- CreateAssayObject(counts = corrected_counts)
    seurat_obj[["corrected"]] <- corrected_assay
    
    return(seurat_obj)
}
