# R/utils.R

#' Bayesian Ambient RNA Correction
#'
#' This function applies a Bayesian correction for ambient RNA contamination using sparse matrix operations.
#' @param data_matrix A sparse matrix of raw RNA counts.
#' @param prior A numeric value representing the prior probability of contamination.
#' @return A corrected sparse matrix with ambient RNA removed.
#' @export
ambient_correction_bayesian <- function(data_matrix, prior = 0.5) {
  # Estimate ambient RNA profile
  ambient_profile <- Matrix::colMeans(data_matrix) * prior
  
  # Correct the matrix by subtracting the ambient profile
  corrected_matrix <- sweep(data_matrix, 2, ambient_profile, FUN = "-")
  
  # Ensure no negative values
  corrected_matrix[corrected_matrix < 0] <- 0
  
  return(corrected_matrix)
}


#' Filter Cells Based on UMI Counts
#'
#' Filters out cells with total UMI counts below a given threshold.
#' @param seurat_obj A Seurat object.
#' @param threshold Numeric value for minimum UMI count per cell.
#' @return A filtered Seurat object with cells that meet the threshold.
#' @export
filter_cells_by_umi <- function(seurat_obj, threshold = 10) {
  # Filter cells based on UMI count threshold
  seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[colSums(GetAssayData(seurat_obj, slot = "counts")) >= threshold])
  
  return(seurat_obj)
}


#' Estimate Ambient RNA Profile
#'
#' This function calculates the ambient RNA profile based on cells with UMI counts below a threshold.
#' @param data_matrix A sparse matrix of raw RNA counts.
#' @param low_umi_cells A logical vector indicating which cells have low UMI counts.
#' @return A numeric vector representing the ambient RNA profile.
#' @export
estimate_ambient_profile <- function(data_matrix, low_umi_cells) {
  # Estimate the ambient RNA profile using low UMI cells
  ambient_profile <- Matrix::rowMeans(data_matrix[, low_umi_cells])
  
  return(ambient_profile)
}


#' Correct RNA Counts Based on Ambient RNA Profile
#'
#' Corrects RNA counts by subtracting the estimated ambient RNA profile from each cell.
#' @param data_matrix A sparse matrix of raw RNA counts.
#' @param ambient_profile A numeric vector representing the ambient RNA profile.
#' @return A sparse matrix with corrected RNA counts.
#' @export
correct_counts <- function(data_matrix, ambient_profile) {
  # Correct each cell's RNA counts by subtracting the ambient RNA profile
  corrected_counts <- sweep(data_matrix, 1, ambient_profile, FUN = "-")
  
  # Ensure no negative values
  corrected_counts[corrected_counts < 0] <- 0
  
  return(corrected_counts)
}

#' Ensure Sparse Matrix Format
#'
#' Converts a dense matrix to a sparse matrix, if needed.
#' @param matrix A matrix object.
#' @return A sparse matrix.
#' @export
ensure_sparse_matrix <- function(matrix) {
  if (!inherits(matrix, "dgCMatrix")) {
    matrix <- as(matrix, "dgCMatrix")
  }
  return(matrix)
}

