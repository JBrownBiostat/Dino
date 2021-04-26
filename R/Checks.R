# check_DinoIn is a helper function to verify the correct format of the
# arguments of Dino
#
#' @importFrom Matrix Matrix
check_DinoIn <- function(counts, nCores, prec, minNZ,
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
    sope <- checkSlope(slope)
    returnMeta <- checkRetM(returnMeta)
    return(list(counts = counts, nCores = nCores, prec = prec,
                minNZ = minNZ, nSubGene = nSubGene, nSubCell = nSubCell,
                depth = depth, slope = slope, returnMeta = returnMeta))
}


# checkCounts is a helper function to verify the correct format of the counts
# argument
#
#' @importFrom Matrix Matrix
#' @importFrom Matrix colSums
#' @importFrom methods is
checkCounts <- function(counts) {
    if(is.null(dim(counts))) {
        stop("'counts' not a matrix; dim returns null")
    }
    if(!is(counts, "dgCMatrix")) {
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


# checkMinNZ is a helper function to verify the correct format of the minNZ
# argument
#
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


# checkPrecDepth is a helper function to verify the correct format of the prec
# and depth arguments
#
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
            stop("'depth' has length not equal to the number of samples")
        }
    }

    return(list(
        prec = prec,
        depth = depth
    ))
}


# checkNSubGene is a helper function to verify the nSubGene argument
checkNSubGene <- function(nSubGene) {
    if(!is.numeric(nSubGene)) {
        stop("'nSubGene' is not of type numeric")
    }
    if(!is.null(dim(nSubGene))) {
        stop("'nSubGene' has dimension greater than 1")
    }
    if(nSubGene < 100) {
        warning("'nSubGene' less than 100, raise threshold before re-running")
    }
    if(nSubGene != round(nSubGene)) {
        message("rounding 'nSubGene'")
        nSubGene <- round(nSubGene)
    }
    return(nSubGene)
}


# checkNSubCell is a helper function to verify the correct format of the
# nSubCell argument
checkNSubCell <- function(nSubCell) {
    if(!is.numeric(nSubCell)) {
        stop("'nSubCell' is not of type numeric")
    }
    if(!is.null(dim(nSubCell))) {
        stop("'nSubCell' has dimension greater than 1")
    }
    if(nSubCell < 100) {
        warning("'nSubCell' less than 100, raise threshold before re-running")
    }
    if(nSubCell != round(nSubCell)) {
        message("rounding 'nSubCell'")
        nSubCell <- round(nSubCell)
    }
    return(nSubCell)
}


# checkNCores is a helper function to verify the correct format of the nCores
# argument
#
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
        warning(cat("'nCores' is greater than the number of available ",
                    "cores, reducing to ", detectCores()))
        nCores <- detectCores()
    }
    return(nCores)
}


# checkSlope is a helper function to verify the correct format of the slope
# argument
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
        warning("'slope' is outside of usual bounds for sequencing experiments")
    }
    return(slope)
}


# returnMeta is a helper function to verify the correct format of the
# returnMeta argument
checkRetM <- function(returnMeta) {
    if(!is.logical(returnMeta)) {
        stop("'returnMeta' is not of type logical")
    }
    if(length(returnMeta) > 1) {
        stop("'returnMeta' has length greater than 1")
    }
    return(returnMeta)
}
