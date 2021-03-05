#' Normalize scRNAseq data
#'
#' @description \code{Dino} removes cell-to-cell variation in observed
#' counts due to the effects of sequencing depth from single-cell mRNA
#' sequencing experiments. \code{Dino} was particularly designed with UMI
#' based protocols in mind, but is applicable to non-UMI based chemistries
#' in the library preparation stage of sequencing.
#'
#' @usage \code{Dino(counts, nCores = 1, prec = 3, seed = 1, minNZ = 10,
#'   nSubGene = 1e4, nSubCell = 1e4, depth = NULL, slope = NULL,
#'   minSlope = 1/2, maxSlope = 2, clusterSlope = T,
#'   returnMeta = F,
#'   emPar = list(maxIter = 100, tol = 0.1, conPar = 15, maxK = 100), ...)}
#'
#' @param counts A numeric matrix object of expression counts - usually in
#'   dgCMatrix format for memory efficiency. Column names denote cells
#'   (samples or droplets) and row names denote genes.
#' @param nCores A non-negative integer scalar denoting the number of cores
#'   which should be used. Setting nCores to 0 uses all cores as determined by
#'   running \code{parallel::detectCores()}
#' @param prec A positive integer denoting the number of decimals to which to
#'   round depth (if estimated internally via \code{depth = NULL}) and
#'   normalized counts for computational efficiency.
#' @param seed A numeric value denoting the random seed for replicability of
#'   results. \code{Dino} automatically records the current random seed prior
#'   to resetting to \emph{seed}, and resets to the recorded seed prior to
#'   exiting. \emph{Note}, results may not be reproducible if more than 1 core
#'   is used when running \code{Dino}.
#' @param minNZ A positive integer denoting the minimum number of non-zero
#'   counts for a gene to be normalized by the Dino algorithm. It is
#'   recommended to pre-filter the \emph{counts} matrix such that all genes
#'   meet this threshold. Otherwise, genes with fewer than \emph{minNZ}
#'   non-zeros will be scaled by depth for normalization.
#' @param nSubGene A positive integer denoting the number of genes to subset
#'   for calculation of \emph{slope}.
#' @param nSubCell A positive integer denoting the number of samples to subset
#'   for calculation of \emph{slope} and the EM algorithm.
#' @param depth A numeric vector of length equal to the columns of counts.
#'   \emph{depth} denotes a median-centered, log-scale measure of cell-wise
#'   sequencing depth. \code{Dino} defaults to defining depth as the
#'   (within-cell) sum of counts across genes, followed by a log and
#'   median-centering transformation.
#' @param slope A numeric scalar denoting the count-depth relationship on
#'   the log-log scale. Typical values are close to 1 (implying a unit increase
#'   in depth corresponds to a unit increase in expected counts on the log-log
#'   scale), but may be higher, particularly in the case of non-UMI
#'   protocols. \code{Dino} defaults to estimating \emph{slope} internally.
#' @param minSlope A numeric scalar denoting the minimum slope. Fitted slopes
#'   below this value will return a warning and be set to 1
#' @param maxSlope A numeric scalar denoting the maximum slope. Fitted slopes
#'   above this value will return a warning and be set to 1
#' @param clusterSlope A logical indicating whether cells should be
#'   pre-clustered prior to calculation of slope. Under the default where
#'   cells are pre-clustered, cluster is used as a factor in the regression.
#' @param returnMeta A logical indicating whether metadata (sequencing depth
#'   and slope) should be returned.
#' @param emPar A list of parameters to send to the EM algorithm.
#'   \emph{maxIter} denotes the maximum number of model updates. \emph{tol}
#'   denotes the cutoff threshold for reductions in the log likelihood
#'   function. \emph{conPar} denotes the concentration parameter for the
#'   resampling. \emph{conPar = 1} implies full resampling from the fitted
#'   distribution. As \emph{conPar} increases, the normalized expression
#'   converges to the scale-factor normalized values. \emph{maxK} denotes the
#'   maximum number of mixture components in the mixture model.
#' @param ... Additional parameters to pass to \code{Scran::quickCluster}.
#'
#' @return \code{Dino} by default returns a matrix of normalized expression
#'   with identical dimensions as \emph{counts}. If \emph{returnMeta=T}, then
#'   \code{Dino} returns a list of normalized expression, sequencing depth, and
#'   slope.
#'
#' @export
#'
#' @importFrom Matrix Matrix
#' @importFrom Matrix colSums
#' @importFrom Matrix rowSums
#' @importFrom Matrix t
Dino <- function(counts, nCores = 1, prec = 3, seed = 1, minNZ = 10,
                 nSubGene = 1e4, nSubCell = 1e4, depth = NULL, slope = NULL,
                 minSlope = 1/2, maxSlope = 2, clusterSlope = T,
                 returnMeta = F,
                 emPar = list(maxIter = 100, tol = 0.1, conPar = 15, maxK = 100),
                 ...) {
  ## Perform argument checks
  checkOut <- check_DinoIn(counts, nCores, prec, seed, minNZ,
                           nSubGene, nSubCell, depth, slope, returnMeta)
  counts <- t(checkOut$counts)
  nCores <- checkOut$nCores
  prec <- checkOut$prec
  seed <- checkOut$seed
  set.seed(seed)
  oldSeed <- checkOut$oldSeed
  minNZ <- checkOut$minNZ
  nSubGene <- checkOut$nSubGene
  nSubCell <- checkOut$nSubCell
  depth <- checkOut$depth
  slope <- checkOut$slope
  returnMeta <- checkOut$returnMeta
  rm(checkOut)


  ## Set depth
  if(is.null(depth)) {
    message("Computing sequencing depth")
    depth <- log(rowSums(counts))
    depth <- depth - median(depth)
    depth <- round(depth, prec)
  }
  depthRep <- calcDepthRep(depth)


  ## Calculate slope
  if(is.null(slope)) {
    message("Calculating regression slope")
    slope <- calcSlope(counts, depth, nSubGene, nSubCell, clusterSlope,
                       nCores, ...)
    if(slope < minSlope) {
      warning(paste0("Fitted slope (", round(slope, 3), ") below minSlope (",
                     round(minSlope, 3), ") -- setting slope = 1"))
      slope = 1
    }
    if(slope > maxSlope) {
      warning(paste0("Fitted slope (", round(slope, 3), ") above maxSlope (",
                     round(maxSlope, 3), ") -- setting slope = 1"))
      slope = 1
    }
  }


  ## Resample counts
  if(nrow(counts) > nSubCell) {
    subInd <- sample.int(nrow(counts), nSubCell)
    depthRep <- calcDepthRep(depth[subInd])
  } else {
    subInd <- 1:nrow(counts)
  }
  dinoGenes <- which(colSums(counts[subInd, ] > 0) >= minNZ)
  libGenes <- which(!(seq_len(ncol(counts[subInd, ])) %in% dinoGenes))
  if(length(libGenes) > 0){
    warning(paste0("Some genes have expression non-zero expression below ",
                   "'minNZ' when subsampled\nand will be normalized via
                   scale-factor"))
    libCounts <- counts[, libGenes, drop = FALSE]
    countList <- splitGenes(counts[, dinoGenes, drop = FALSE])
  } else {
    countList <- splitGenes(counts[, dinoGenes, drop = FALSE])
  }
  prll <- setPar(nCores, countList)

  message("Normalizing counts by resampling")
  if(length(libGenes) > 0) {
    normCounts.lib <- libCounts / exp(depth)
    colnames(normCounts.lib) <- colnames(counts)[libGenes]
    normCounts.Dino <- parResampCounts(countList, depth, depthRep, slope,
                                       prec, subInd, prll, emPar)
    colnames(normCounts.Dino) <- colnames(counts)[dinoGenes]
    normCounts <- cbind(normCounts.Dino, normCounts.lib)
    normCounts <- normCounts[, colnames(counts)]
  } else {
    normCounts <- parResampCounts(countList, depth, depthRep, slope,
                                  prec, subInd, prll, emPar)
  }

  colnames(normCounts) <- colnames(counts)
  rownames(normCounts) <- rownames(counts)
  normCounts <- t(normCounts)

  message("Done")
  set.seed(oldSeed)
  if(returnMeta == F) {
    return(normCounts)
  } else {
    return(list(
      normMat = normCounts,
      depth = depth,
      slope = slope
    ))
  }
}


#' Create Seurat object from Dino normalized data
#'
#' @description \code{SeuratFromDino} is a wrapper simplifying the export of
#'   \code{Dino} normalized counts to a \emph{Seurat} object for secondary
#'   analysis.
#'
#' @usage \code{SeuratFromDino(counts, doNorm = T, ...)}
#'
#' @param counts A numeric matrix of count data, either raw (eg. UMIs) or
#'   normalized expression.
#' @param doNorm A logical indicating whether to normalize the input
#'   \emph{counts} data before exporting results to a \emph{Seurat} object.
#'   By default, it is assumed that the contents of \emph{counts} raw
#'   expression which should be normalized.
#' @param doLog A logical indicating whether normalized counts should be log
#'   transformed with a psuedocount of 1 prior to export.
#' @param ... Additional parameters to pass to \code{Dino}
#'
#' @export
#'
#' @importFrom Seurat CreateAssayObject
#' @importFrom Seurat CreateSeuratObject
SeuratFromDino <- function(counts, doNorm = T, doLog = T, ...) {
  if(doNorm) {
    DinoOut <- Dino(counts, ...)
    if(is.list(DinoOut)) {
      counts <- DinoOut$normMat
    } else {
      counts <- DinoOut
    }
  }
  if(doLog) {
    normAssay <- CreateAssayObject(data = log1p(counts))
  } else {
    normAssay <- CreateAssayObject(data = counts)
  }
  return(CreateSeuratObject(normAssay))
}
