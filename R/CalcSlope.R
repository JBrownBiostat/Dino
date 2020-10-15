#' Calculate expression-depth slope
#'
#' @description \code{calcSlope} is a helper function to calculate the
#'   expression-depth slope
#'
#' @usage \code{calcSlope(counts, depth, nSubCell, nSubGene,
#'   clusterSlope, nCores, ...)}
#'
#' @importFrom Matrix Matrix
#' @importFrom Matrix colSums
#' @importFrom scran quickCluster
#' @importFrom BiocSingular IrlbaParam
calcSlope <- function(counts, depth, nSubGene, nSubCell, clusterSlope,
                      nCores, ...) {
  if(nrow(counts) > nSubCell) {
    subInd <- sample.int(nrow(counts), nSubCell)
    subCounts <- counts[subInd, ]
    subDepth <- depth[subInd]
  } else {
    subCounts <- counts
    subDepth <- depth
  }
  prll <- setPar(nCores)
  if(clusterSlope) {
    clustVec <- quickCluster(x = t(subCounts), BSPARAM = IrlbaParam(),
                             BPPARAM = prll, ...)
  } else {
    clustVec <- rep(1, nSubCell)
  }

  subGenes <- colSums(subCounts > 0) / nrow(subCounts) >= 1e-2
  geneMu <- log(apply(subCounts, 2, function(x) {
    exp(mean(log(x + 1))) - 1
  }))
  muDens <- density(geneMu)
  muFunc <- approxfun(muDens$x, muDens$y)
  geneSamp <- sample(which(subGenes),
                     min(sum(subGenes), nSubGene),
                     prob = 1 / muFunc(geneMu[subGenes]))
  datList <- splitGenes(subCounts[, geneSamp], nCol = 40)
  clustTab <- table(clustVec)
  prll <- setPar(nCores, datList)
  slopeVec <- parPoisSlope(datList, subDepth, clustVec, clustTab,
                           prll, pctNZ = 0.01)

  return(mean(slopeVec))
}


#' Calculate regression slopes
#'
#' @description \code{poisSlope} is a helper function to compute poisson
#'   regression slopes
#'
#' @usage \code(poisSlope(y, depth, clustVec, clustTab, pctNZ = 0.01))
poisSlope <- function(y, depth, clustVec, clustTab, pctNZ = 0.01) {
  clustNZ <- table(clustVec[y > 0])
  regInd <- clustVec %in% names(clustNZ)[clustNZ > 1]
  y <- y[regInd]; depth <- depth[regInd]; clustVec <- clustVec[regInd]
  clustVec <- factor(clustVec, levels = unique(clustVec))
  if(length(levels(clustVec)) >= 2) {
    return(glm(y ~ depth + clustVec, family = "poisson")$coefficients[2])
  } else {
    return(glm(y ~ depth, family = "poisson")$coefficients[2])
  }
}


#' Calculate parallel regression slopes
#'
#' @description \code{parPoisSlope} is a helper function to compute poisson
#'   regression slopes
#'
#' @usage \code(parPoisSlope(datList, depth, clustVec, clustTab,
#'   prll, pctNZ = 0.01))
#'
#' @importFrom Matrix Matrix
parPoisSlope <- function(datList, depth, clustVec, clustTab, prll,
                         pctNZ = 0.01) {
  slopeVec <- bplapply(datList,
                       function(subDat, depth, clustVec, clustTab, pctNZ) {
                         subSlopes <- rep(1, ncol(subDat))
                         for(i in 1:ncol(subDat)) {
                           subSlopes[i] <- poisSlope(subDat[, i], depth,
                                                     clustVec, clustTab,
                                                     pctNZ)
                         }
                         return(subSlopes)
                       }, depth = depth, clustVec = clustVec,
                       clustTab = clustTab, pctNZ = pctNZ, BPPARAM = prll)
  do.call(c, slopeVec)
}
