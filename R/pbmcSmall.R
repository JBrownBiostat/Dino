#' Subset of 500 peripheral blood mononuclear cells (PBMCs) from a healthy donor
#'
#' This dataset derives from the "3k PBMCs from a Healthy Donor" public dataset
#' from 10X Genomics.
#'
#' @docType data
#'
#' @usage data(pbmcSmall)
#'
#' @format An object of class \code{"dgCMatrix"}.
#'
#' @keywords datasets
#'
#' @source \href{https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k}{3k PBMCs from a Healthy Donor}
#'
#' @examples
#' data(pbmcSmall)
#' str(pbmcSmall)
#'
#' @importClassesFrom Matrix dgCMatrix
"pbmcSmall"
