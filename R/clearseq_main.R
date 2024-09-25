#' Perform Ambient RNA Correction and SeuratExtend Integration
#'
#' This function applies ambient RNA correction to a Seurat object and integrates SeuratExtend functionalities.
#' @param seurat_obj A Seurat object containing the scRNA-seq data.
#' @param threshold UMI count threshold for identifying ambient RNA contamination.
#' @return A Seurat object with ambient RNA corrected and SeuratExtend methods applied.
#' @export
clearseq_with_extend <- function(seurat_obj, threshold = 10) {
  # Ensure the counts matrix is in a sparse matrix format
  if (!inherits(GetAssayData(seurat_obj, slot = "counts"), "dgCMatrix")) {
    stop("The counts data must be a sparse matrix.")
  }
  
  # Apply ambient RNA correction
  total_umi <- Matrix::colSums(GetAssayData(seurat_obj, slot = "counts"))
  empty_droplets <- total_umi < threshold
  ambient_profile <- Matrix::rowMeans(GetAssayData(seurat_obj, slot = "counts")[, empty_droplets])
  corrected_counts <- sweep(GetAssayData(seurat_obj, slot = "counts"), 1, ambient_profile, FUN = "-")
  corrected_counts[corrected_counts < 0] <- 0
  
  # Add corrected counts as a new assay in Seurat object
  corrected_assay <- CreateAssayObject(counts = corrected_counts)
  seurat_obj[["corrected"]] <- corrected_assay
  
  return(seurat_obj)
}
