#' Check Dino input
#'
#' @description \code{check_DinoIn} is a helper function to verify the
#'   correct format of the arguments of \code{Dino}.
#'
#' @usage \code{counts, nCores, prec, seed, minNZ, nSubGene, nSubCell, depth,
#'   slope, returnMeta)}
#'
#' @importFrom Matrix Matrix
check_DinoIn <- function(counts, nCores, prec, seed, minNZ,
                         nSubGene, nSubCell, depth, slope, returnMeta) {
  ## Checks
  counts <- checkCounts(counts)
  minNZ <- checkMinNZ(minNZ, counts)
  checkRet <- checkPrecDepth(prec, depth, counts)
  prec <- checkRet$prec
  depth <- checkRet$depth
  nSubGene <- checkNSubGene(nSubGene)
  nSubCell <- checkNSubCell(nSubCell)
  nCores <- checkNCores(nCores)
  seedRet <- checkSeed(seed)
  seed <- seedRet$seed
  oldSeed <- seedRet$oldSeed
  sope <- checkSlope(slope)
  returnMeta <- checkRetM(returnMeta)
  return(list(counts = counts, nCores = nCores, prec = prec,
              seed = seed, oldSeed = oldSeed, minNZ = minNZ,
              nSubGene = nSubGene, nSubCell = nSubCell, depth = depth,
              slope = slope, returnMeta = returnMeta))
}


#' Check \emph{counts} input
#'
#' @description \code{checkCounts} is a helper function to verify the
#'   correct format of the \emph{counts} argument
#'
#' @usage \code{checkCounts(counts)}
#'
#' @importFrom Matrix Matrix
#' @importFrom Matrix colSums
checkCounts <- function(counts) {
  if(is.null(dim(counts))) {
    stop("'counts' not a matrix: dim(counts) returns null")
  }
  if(class(counts) != "dgCMatrix") {
    message("converting 'counts' to dgCMatrix")
    counts <- Matrix(counts, sparse = TRUE)
  }
  if(any(round(counts) != counts)) {
    message("rounding 'counts'")
    counts <- Matrix(round(counts), sparse = TRUE)
  }
  if(any(counts < 0)) {
    stop("'counts' contains negative values")
  }
  if(any(colSums(counts) == 0)) {
    stop("'counts' contains samples with all-zero genes")
  }
  return(counts)
}


#' Check \emph{minNZ} input
#'
#' @description \code{checkMinNZ} is a helper function to verify the
#'   correct format of the \emph{minNZ} argument
#'
#' @usage \code{checkMinNZ(minNZ, counts)}
#'
#' @importFrom Matrix Matrix
#' @importFrom Matrix rowSums
checkMinNZ <- function(minNZ, counts) {
  if(!is.numeric(minNZ)) {
    stop("'minNZ' is not a number")
  }
  if(!is.null(dim(minNZ))) {
    stop("'minNZ' has dimension greater than 1")
  }
  if(length(minNZ) > 1) {
    stop("'minNZ' doesn't have length 1")
  }
  if(minNZ != round(minNZ)) {
    message("rounding 'minNZ'")
    minNZ <- round(minNZ)
  }
  if(minNZ < 2) {
    stop("'minNZ' less than 2, choose a larger value")
  }
  if(minNZ < 10) {
    warning("'minNZ' less than 10, this is not recommended")
  }
  if(!any(rowSums(counts > 0) >= minNZ)) {
    stop("all genes have fewer than 'minNZ' non-zero values")
  }
  return(minNZ)
}


#' Check \emph{prec} and \emph{depth} inputs
#'
#' @description \code{checkPrecDepth} is a helper function to verify the
#'   correct format of the \emph{prec} and \emph{depth} arguments
#'
#' @usage \code{checkPrecDepth(depth, counts)}
#'
#' @importFrom Matrix Matrix
checkPrecDepth <- function(prec, depth, counts) {
  if(!is.numeric(prec)) {
    stop("'prec' is not a number")
  }
  if(!is.null(dim(prec))) {
    stop("'prec' has dimension greater than 1")
  }
  if(length(prec) > 1) {
    stop("'prec' doesn't have length 1")
  }
  if(prec != round(prec)) {
    message("rounding 'prec'")
    prec <- round(prec)
  }

  if(!is.null(depth))  {
    if(!is.numeric(depth)) {
      stop("'depth' is not of type numeric")
    }
    if(!is.null(dim(depth))) {
      stop("'depth' has dimension greater than 1")
    }
    if(length(depth) != ncol(counts)) {
      stop("'depth' has length not equal to the number of samples in 'counts'")
    }
  }

  return(list(
    prec = prec,
    depth = depth
  ))
}


#' Check \emph{nSubGene} inputs
#'
#' @description \code{checkNSubGene} is a helper function to verify the
#'   correct format of the \emph{nSubGene} argument.
#'
#' @usage \code{checkNSubGene(nSubGene)}
checkNSubGene <- function(nSubGene) {
  if(!is.numeric(nSubGene)) {
    stop("'nSubGene' is not of type numeric")
  }
  if(!is.null(dim(nSubGene))) {
    stop("'nSubGene' has dimension greater than 1")
  }
  if(nSubGene < 100) {
    warning("'nSubGene' less than 100, increase threshold before re-running")
  }
  if(nSubGene != round(nSubGene)) {
    message("rounding 'nSubGene'")
    nSubGene <- round(nSubGene)
  }
  return(nSubGene)
}


#' Check \emph{nSubCell} inputs
#'
#' @description \code{checkNSubCell} is a helper function to verify the
#'   correct format of the \emph{nSubCell} argument.
#'
#' @usage \code{checkNSubCell(nSubCell)}
checkNSubCell <- function(nSubCell) {
  if(!is.numeric(nSubCell)) {
    stop("'nSubCell' is not of type numeric")
  }
  if(!is.null(dim(nSubCell))) {
    stop("'nSubCell' has dimension greater than 1")
  }
  if(nSubCell < 100) {
    warning("'nSubCell' less than 100, increase threshold before re-running")
  }
  if(nSubCell != round(nSubCell)) {
    message("rounding 'nSubCell'")
    nSubCell <- round(nSubCell)
  }
  return(nSubCell)
}


#' Check \emph{nCores} input
#'
#' @description \code{checkNCores} is a helper function to verify the
#'   correct format of the \emph{nCores} argument
#'
#' @usage \code{checkNCores(nCores)}
#'
#' @importFrom parallel detectCores
checkNCores <- function(nCores) {
  if(!is.numeric(nCores)) {
    stop("'nCores' is not a number")
  }
  if(!is.null(dim(nCores))) {
    stop("'nCores' has dimension greater than 1")
  }
  if(length(nCores) > 1) {
    stop("'nCores' doesn't have length 1")
  }
  if(nCores != round(nCores)) {
    message("rounding 'nCores'")
    nCores <- round(nCores)
  }
  if(nCores < 0) {
    stop("'nCores' is less than 0")
  }
  if(nCores == 0) {
    message("using all available cores")
    nCores <- detectCores()
  }
  if(nCores > detectCores()) {
    warning("'nCores' is greater than the number of available cores, reducing
            to ", detectCores())
    nCores <- detectCores()
  }
  return(nCores)
}


#' Check \emph{seed} inputs
#'
#' @description \code{checkSeed} is a helper function to verify the
#'   correct format of the \emph{seed} argument and adjust the working seed
#'   for replicability if \emph{seed} is not \code{NULL}.
#'
#' @usage \code{checkSeed(seed)}
checkSeed <- function(seed) {
  if(is.null(seed)) {
    return(NULL)
  }
  if(!is.numeric(seed)) {
    stop("'seed' is not of type numeric")
  }
  if(!is.null(dim(seed))) {
    stop("'seed' has dimension greater than 1")
  }
  if(seed != round(seed)) {
    stop("'seed' is not a valid integer")
  }
  oldSeed <- .Random.seed
  return(list(oldSeed = oldSeed,
              seed = seed))
}


#' Check \emph{slope} inputs
#'
#' @description \code{checkSlope} is a helper function to verify the
#'   correct format of the \emph{slope} argument.
#'
#' @usage \code{checkSlope(slope)}
checkSlope <- function(slope) {
  if(is.null(slope)) {
    return(NULL)
  }
  if(!is.numeric(slope)) {
    stop("'slope' is not of type numeric")
  }
  if(!is.null(dim(slope))) {
    stop("'slope' has dimension greater than 1")
  }
  if(slope < 0.25 | slope > 4) {
    warning("'slope' is outside of typical bounds for sequencing experiments")
  }
  return(slope)
}


#' Check \emph{returnMeta} inputs
#'
#' @description \code{returnMeta} is a helper function to verify the
#'   correct format of the \emph{returnMeta} argument.
#'
#' @usage \code{checkRetM(returnMeta)}
checkRetM <- function(returnMeta) {
  if(!is.logical(returnMeta)) {
    stop("'returnMeta' is not of type logical")
  }
  if(length(returnMeta) > 1) {
    stop("'returnMeta' has length greater than 1")
  }
  return(returnMeta)
}
