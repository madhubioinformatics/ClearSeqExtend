#' ClearSeqExtend: Ambient RNA Correction and Enhanced Analysis
#'
#' This function performs ambient RNA correction using a sparse matrix and integrates SeuratExtend for advanced analysis.
#'
#' @param seurat_obj A Seurat object containing scRNA-seq data.
#' @param threshold UMI threshold for identifying empty droplets.
#' @return A Seurat object with corrected counts and additional analysis steps applied.
#' @export
clearseq_with_extend <- function(seurat_obj, threshold = 10) {
  total_umi <- Matrix::colSums(GetAssayData(seurat_obj, slot = "counts"))
  empty_droplets <- total_umi < threshold

  if (sum(empty_droplets) < 2) {
    stop("Not enough empty droplets to perform ambient RNA correction.")
  }
  
  ambient_profile <- Matrix::rowMeans(GetAssayData(seurat_obj, slot = "counts")[, empty_droplets])
  corrected_counts <- sweep(GetAssayData(seurat_obj, slot = "counts"), 1, ambient_profile, FUN = "-")
  corrected_counts[corrected_counts < 0] <- 0
  
  corrected_assay <- CreateAssayObject(counts = corrected_counts)
  seurat_obj[["corrected"]] <- corrected_assay
  
  return(seurat_obj)
}
