#' Calculate replications of depth
#'
#' @description \code{calcDepthRep} is a helper function to calculate the
#'   replications in the depth vector. Procedure is close to that of
#'   \code{table}, but maintains the same rounding scheme as other functions.
#'
#' @usage \code{calcDepthRep(depth)}
calcDepthRep <- function(depth) {
  depth <- sort(depth)
  dupDepth <- duplicated(depth)
  freq <- c(which(!dupDepth), length(depth) + 1)
  return(data.frame(
    depth = depth[!dupDepth],
    rep = freq[-1] - freq[-length(freq)]
  ))
}


#' Setup parallel environment
#'
#' @description \code{setPar} is a helper function to initiate the parallel
#'   environment
#'
#' @usage \code{setPar(nCores, taskList)}
#'
#' @importFrom BiocParallel SnowParam
#' @importFrom BiocParallel MulticoreParam
#' @importFrom BiocParallel register
setPar <- function(nCores, taskList = NULL) {
  if(is.null(taskList)) {
    nTask = 0L
  } else {
    nTask <- ceiling(sqrt(length(taskList) / nCores)) * nCores
  }
  if (.Platform$OS.type == "windows") {
    prll <- SnowParam(workers = nCores, tasks = nTask, progressbar = T)
    register(BPPARAM = prll, default = TRUE)
  } else {
    prll <- MulticoreParam(workers = nCores, tasks = nTask, progressbar = T)
    register(BPPARAM = prll, default = TRUE)
  }
  return(prll)
}


#' Split sparse matrix into list of columns
#'
#' @description \code{splitGenes} is a helper function to split a sparse count
#'   matrix into a list of gene sets for parallel computation.
#'
#' @usage \code(splitGenes(counts, nCol = 40))
#'
#' @importFrom Matrix Matrix
splitGenes <- function(counts, nCol = 40) {
  if(nCol >= ncol(counts)) {
    return(list(counts))
  }
  m <- floor(ncol(counts) / nCol)
  ret <- lapply(seq_len(m), function(i) {
    counts[, (nCol * i - (nCol - 1)):(nCol * i), drop = FALSE]
  })
  if(m * nCol < ncol(counts)) {
    ret[[m + 1]] <- counts[, (nCol * m + 1):ncol(counts), drop = FALSE]
  }
  return(ret)
}
