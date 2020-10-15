#' Resample counts
#'
#' @description \code{resampCounts} is a helper function to normalize counts
#'   by resampling from an estimated expression distribution conditional on
#'   the observed counts.
#'
#' @usage \code{resampCounts(countMat, depth, slope, prec, subInd, emPar)}
#'
#' @importFrom Matrix Matrix
resampCounts <- function(countMat, depth, depthRep, slope,
                         prec, subInd, emPar) {
  subNorm <- apply(countMat,
                   2,
                   function(y, depth, depthRep, slope, prec, subInd) {
    # Calculate eCDF
    outDev <- calcDev(y[subInd], depth[subInd], depthRep, slope)
    cdfList <- list(
      eCDF_fun = approxfun(
        x = outDev$quant,
        y = outDev$pct,
        yleft = 0, yright = 1
      ),
      eQuant_fun = approxfun(
        x = outDev$pct,
        y = outDev$quant,
        yleft = -Inf, yright = max(outDev$quant)
      )
    )

    # Run EM
    minD <- min(depth)
    if(minD < 0) {
      fitSlope <- min(slope, (exp(minD - 3) - 1) / (exp(minD) - 1))
    } else {
      fitSlope <- max(slope, (exp(minD - 3) - 1) / (exp(minD) - 1))
    }
    parList <- initEM.func(
      y = y[subInd],
      slope = fitSlope,
      subCDF = cdfList,
      depth = depth[subInd],
      stepInit = 5
    )

    nUpdate <- emPar$maxIter
    deltaStop <- emPar$tol
    iter <- 0

    while(iter < nUpdate) {
      iter <- iter + 1

      parList$d <- -parList$gTilde$gTilde +
        parList$S %*% parList$g$g
      parList$alpha <- -1

      gamInd <- parList$gamList$gamVec + parList$alpha * parList$d[1:parList$J]
      if(!all(gamInd >= -75)) {
        if(parList$J > 2) {
          parList <- parListUpdate.func(parList)
          parList$lLik <- lLik.func(parList)
          parList$S <- parList$S * 0
        } else {break}
      }

      betaInd <- parList$betaVec +
                      parList$alpha * parList$d[(parList$J + 1):(2 * parList$J)]
      if(!all(betaInd > -75)) {
        if(parList$J > 2) {
          parList <- parListUpdate.func(parList, doBeta = T)
          parList$lLik <- lLik.func(parList)
          parList$S <- parList$S * 0
        } else {break}
      }
      if(any(gamInd > 75) | any(betaInd > 75) | max(gamInd) - min(gamInd) > log(2 * length(subInd))) {
        parList$S <- parList$S * 0
      }

      parList <- calcStep.func(parList)

      dTheta <- parList$parUpdate$parVec -
        c(parList$gamList$gamVec, parList$betaVec)
      dg <- parList$gUpdate$g - parList$g$g
      parList$dgTilde <- gTildeUpdate.func(parList)
      dgTilde <- parList$dgTilde$gTilde - parList$gTilde$gTilde

      dThetaStar <- -dgTilde + parList$S %*% dg
      parList$S <- parList$S +
        c(1 + (t(dg) %*% dThetaStar) / (t(dg) %*% dTheta)) *
        (dTheta %*% t(dTheta)) / c(t(dg) %*% dTheta) -
        (dThetaStar %*% t(dTheta) + t(dThetaStar %*% t(dTheta))) /
        c(t(dg) %*% dTheta)

      parList$pList <- parList$pListUpdate
      parList$gamList <- parList$parUpdate$gamList
      parList$betaVec <- parList$parUpdate$betaVec
      parList$g <- parList$gUpdate
      parList$gTilde <- parList$dgTilde
      if(parList$lLikUpdate - parList$lLik <= deltaStop) {
        break
      }
      parList$lLik <- parList$lLikUpdate
    }

    # Update model for all points
    parList$y <- y
    parList$depth <- exp(depth) * fitSlope - fitSlope + 1
    parList$pList <- lpMat.func(parList)

    tauPar <- parList$gamList$gamVec - max(parList$gamList$gamVec)
    tauPar <- exp(tauPar)
    tauPar <- tauPar / sum(tauPar)
    tauPar <- pmax(tauPar, .Machine$double.xmin * 1e1)
    tauPar <- tauPar / sum(tauPar)
    lamPar <- pmax(exp(parList$betaVec), .Machine$double.xmin * 1e1)
    parList$pList$spMat <- parList$pList$spMat * exp(parList$pList$lpMax)

    # Resample counts
    distVec <- apply(parList$pList$spMat[-parList$J, , drop = F],
                     2,
                     function(pVec) {
      which.min(runif(1) > c(cumsum(pVec), 1))
    })

    bwVec <- sapply(1:25, function(k) {
      gamSamp <- lamPar[sample.int(parList$J,
                                   parList$J,
                                   replace = T,
                                   prob = tauPar)]
      tryCatch({
        suppressWarnings(bw.ucv(gamSamp))
      }, error = function(e) {
        tryCatch({
          bw.nrd(gamSamp)
        }, error = function(e) {
          1 / sqrt(parList$J)
        })
      }
      )
    })
    bw <- mean(bwVec)

    shapeVec <- lamPar / bw
    concVec <- emPar$conPar

    normVec <- rgamma(n = length(parList$y),
                      shape = parList$y * concVec + shapeVec[distVec],
                      scale = 1 / (parList$depth * concVec + 1 / bw))
    normVec <- round(normVec, prec)

    return(normVec)
  }, depth = depth, depthRep = depthRep, slope = slope,
  prec = prec, subInd = subInd)
  subNorm <- Matrix(subNorm, sparse = T)
  return(subNorm)
}


#' Parallelize count resampling
#'
#' @description \code{parResampCounts} is a helper function to parallelize the
#'   resampling algorithm.
#'
#' @usage \code{parResampCounts(countList, depth, depthRep,
#'   slope, prec, subInd, prll, emPar)}
#'
#' @importFrom BiocParallel bplapply
#' @importFrom Matrix Matrix
parResampCounts <- function(countList, depth, depthRep, slope,
                            prec, subInd, prll, emPar) {
  normList <- bplapply(countList,
                       function(subCounts, depth, depthRep,
                                slope, prec, subInd) {
                         resampCounts(subCounts, depth, depthRep,
                                      slope, prec, subInd, emPar)
                       }, depth = depth, depthRep = depthRep, slope = slope,
                       prec = prec, subInd = subInd, BPPARAM = prll)
  while(length(normList) > 2) {
    for(i in 1:floor(length(normList) / 2)) {
      if(length(normList) >= 2 * i) {
        normList[[i]] <- cbind(normList[[2 * i - 1]], normList[[2 * i]])
      } else {
        normList[[i]] <- normList[[2 * i - 1]]
      }
    }
    if(length(normList) > 2 * i) {
      i <- i + 1
      normList[[i]] <- normList[[2 * i - 1]]
    }
    while(length(normList) > i) {
      normList[[i + 1]] <- NULL
    }
  }
  if(length(normList) == 2) {
    normList <- cbind(normList[[1]], normList[[2]])
  }
  return(normList)
}


#' Calculate EM log-likelihood
#'
#' @description \code{lLik.func} is a helper function to compute the EM
#'   log-likelihood up to an addative constant
#'
#' @usage \code(lLik.func(parList))
#'
#' @importFrom matrixStats rowMaxs
lLik.func <- function(parList) {
  M <- outer(parList$y, parList$betaVec)
  M <- M - outer(parList$depth, exp(parList$betaVec))
  M <- M + outer(rep(1, length(parList$depth)), parList$gamList$gamVec)
  S <- rowMaxs(M)
  return(sum(log(rowSums(exp(M - S))) + S) -
           length(parList$y) * log(parList$gamList$gamDot))
}


#' Calculate Full EM log-likelihood
#'
#' @description \code{lLikFull.func} is a helper function to compute the EM
#'   log-likelihood up to an addative constant
#'
#' @usage \code(lLikFull.func(parList))
#'
#' @importFrom matrixStats rowMaxs
lLikFull.func <- function(parList) {
  M <- outer(parList$y, parList$betaVec)
  M <- M - outer(parList$depth, exp(parList$betaVec))
  M <- M + outer(rep(1, length(parList$depth)), parList$gamList$gamVec)
  S <- rowMaxs(M)
  return(sum(log(rowSums(exp(M - S))) + S) -
           length(parList$y) * log(parList$gamList$gamDot) +
           sum(parList$y * log(parList$depth) - lfactorial(parList$y)))
}


#' Calculate EM log-likelihood at update pars
#'
#' @description \code{lLikUpdate.func} is a helper function to compute the EM
#'   log-likelihood up to an addative constant at update pars
#'
#' @usage \code(lLikUpdate.func(parList))
#'
#' @importFrom matrixStats rowMaxs
lLikUpdate.func <- function(parList) {
  M <- outer(parList$y, parList$parUpdate$betaVec)
  M <- M - outer(parList$depth, exp(parList$parUpdate$betaVec))
  M <- M + outer(rep(1, length(parList$depth)), parList$parUpdate$gamList$gamVec)
  S <- rowMaxs(M)
  return(sum(log(rowSums(exp(M - S))) + S) -
           length(parList$y) * log(parList$parUpdate$gamList$gamDot))
}


#' Calculate EM membership log-probabilities
#'
#' @description \code{lpMat.func} is a helper function to compute the EM
#'   log membership probabilities
#'
#' @usage \code(lpMat.func(parList))
#'
#' @importFrom matrixStats rowMaxs
#' @importFrom matrixStats colMaxs
lpMat.func <- function(parList) {
  lpMat <- outer(parList$betaVec, parList$y)
  lpMat <- lpMat - outer(exp(parList$betaVec), parList$depth)
  lpMat <- lpMat + parList$gamList$gamVec
  lpMat <- lpMat - outer(rep(1, parList$J), colMaxs(lpMat))
  lpMat <- lpMat - outer(rep(1, parList$J), log(colSums(exp(lpMat))))
  lpMax <- rowMaxs(lpMat)
  return(list(
    spMat = exp(lpMat - lpMax),
    lpMax = lpMax
  ))
}


#' Calculate EM membership log-probabilities at update pars
#'
#' @description \code{lpMatUpdate.func} is a helper function to compute the EM
#'   log membership probabilities at update pars
#'
#' @usage \code(lpMatUpdate.func(parList))
#'
#' @importFrom matrixStats rowMaxs
#' @importFrom matrixStats colMaxs
lpMatUpdate.func <- function(parList) {
  lpMat <- outer(parList$parUpdate$betaVec, parList$y)
  lpMat <- lpMat - outer(exp(parList$parUpdate$betaVec), parList$depth)
  lpMat <- lpMat + parList$parUpdate$gamList$gamVec
  lpMat <- lpMat - outer(rep(1, parList$J), colMaxs(lpMat))
  lpMat <- lpMat - outer(rep(1, parList$J), log(colSums(exp(lpMat))))
  lpMax <- rowMaxs(lpMat)
  return(list(
    spMat = exp(lpMat - lpMax),
    lpMax = lpMax
  ))
}


#' Calculate EM gamma updates
#'
#' @description \code{gamUpdate.func} is a helper function to compute the EM
#'   update to the gamma parameter
#'
#' @usage \code(gamUpdate.func(parList))
gamUpdate.func <- function(parList) {
  gamUpdate <- log(rowSums(parList$pList$spMat)) +
    parList$pList$lpMax
  return(gamUpdate - gamUpdate[parList$zInd])
}


#' Calculate EM gamma updates at update pars
#'
#' @description \code{gamUpdate2.func} is a helper function to compute the EM
#'   update to the gamma parameter at update pars
#'
#' @usage \code(gamUpdate2.func(parList))
gamUpdate2.func <- function(parList) {
  gamUpdate <- log(rowSums(parList$pListUpdate$spMat)) +
    parList$pListUpdate$lpMax
  return(gamUpdate - gamUpdate[parList$zInd])
}


#' Calculate EM beta updates
#'
#' @description \code{betaUpdate.func} is a helper function to compute the EM
#'   update to the beta parameter
#'
#' @usage \code(betaUpdate.func(parList))
betaUpdate.func <- function(parList) {
  betaUpdate <- log(pmax(parList$pList$spMat %*% parList$y, sqrt(.Machine$double.xmin))) -
    log(pmax(parList$pList$spMat %*% parList$depth, sqrt(.Machine$double.xmin)))
  return(c(betaUpdate))
}


#' Calculate EM beta updates at update pars
#'
#' @description \code{betaUpdate2.func} is a helper function to compute the EM
#'   update to the beta parameter at update pars
#'
#' @usage \code(betaUpdate2.func(parList))
betaUpdate2.func <- function(parList) {
  betaUpdate <- log(pmax(parList$pListUpdate$spMat %*% parList$y, sqrt(.Machine$double.xmin))) -
    log(pmax(parList$pListUpdate$spMat %*% parList$depth, sqrt(.Machine$double.xmin)))
  return(c(betaUpdate))
}


#' Calculate EM parameter updates
#'
#' @description \code{gTilde.func} is a helper function to compute the EM
#'   parameter update deltas
#'
#' @usage \code(gTilde.func(parList))
gTilde.func <- function(parList) {
  gamUpdate <- gamUpdate.func(parList)
  betaUpdate <- betaUpdate.func(parList)
  return(list(
    gamTilde = gamUpdate - parList$gamList$gamVec,
    betaTilde = betaUpdate - parList$betaVec,
    gTilde = c(gamUpdate - parList$gamList$gamVec,
               betaUpdate - parList$betaVec)
  ))
}


#' Calculate EM parameter updates at update pars
#'
#' @description \code{gTildeUpdate.func} is a helper function to compute the EM
#'   parameter update deltas at update pars
#'
#' @usage \code(gTildeUpdate.func(parList))
gTildeUpdate.func <- function(parList) {
  gamUpdate <- gamUpdate2.func(parList)
  betaUpdate <- betaUpdate2.func(parList)
  return(list(
    gamTilde = gamUpdate - parList$parUpdate$gamList$gamVec,
    betaTilde = betaUpdate - parList$parUpdate$betaVec,
    gTilde = c(gamUpdate - parList$parUpdate$gamList$gamVec,
               betaUpdate - parList$parUpdate$betaVec)
  ))
}


#' Calculate likelihood gradient in gamma
#'
#' @description \code{gGamma.func} is a helper function to compute the
#'   likelihood gradient in gamma
#'
#' @usage \code(gGamma.func(parList))
gGamma.func <- function(parList) {
  M <- rowSums(parList$pList$spMat) * exp(parList$pList$lpMax)
  M <- M - sum(M) * exp(parList$gamList$gamVec - log(parList$gamList$gamDot))
  M[parList$zInd] <- 0
  return(M)
}


#' Calculate likelihood gradient in gamma at update pars
#'
#' @description \code{gGamma2.func} is a helper function to compute the
#'   likelihood gradient in gamma at update pars
#'
#' @usage \code(gGamma2.func(parList))
gGamma2.func <- function(parList) {
  M <- rowSums(parList$pListUpdate$spMat) * exp(parList$pListUpdate$lpMax)
  M <- M - sum(M) * exp(parList$parUpdate$gamList$gamVec -
                          log(parList$parUpdate$gamList$gamDot))
  M[parList$zInd] <- 0
  return(M)
}


#' Calculate likelihood gradient in beta
#'
#' @description \code{gBeta.func} is a helper function to compute the
#'   likelihood gradient in beta
#'
#' @usage \code(gBeta.func(parList))
gBeta.func <- function(parList) {
  M <- parList$pList$spMat %*% parList$y -
    rowSums(parList$pList$spMat *
              outer(exp(parList$betaVec), parList$depth))
  return(c(exp(parList$pList$lpMax) * M))
}


#' Calculate likelihood gradient in beta at update pars
#'
#' @description \code{gBeta2.func} is a helper function to compute the
#'   likelihood gradient in beta at update pars
#'
#' @usage \code(gBeta2.func(parList))
gBeta2.func <- function(parList) {
  M <- parList$pListUpdate$spMat %*% parList$y -
    rowSums(parList$pListUpdate$spMat *
              outer(exp(parList$parUpdate$betaVec), parList$depth))
  return(c(exp(parList$pListUpdate$lpMax) * M))
}


#' Calculate EM parameter gradients
#'
#' @description \code{parGrad.func} is a helper function to compute the EM
#'   parameter gradients
#'
#' @usage \code(parGrad.func(parList))
parGrad.func <- function(parList) {
  gGam <- gGamma.func(parList)
  gBeta <- gBeta.func(parList)
  return(list(
    gGam = gGam,
    gBeta = gBeta,
    g = c(gGam, gBeta)
  ))
}


#' Calculate EM parameter gradients at update pars
#'
#' @description \code{parGradUpdate.func} is a helper function to compute the EM
#'   parameter gradients at update pars
#'
#' @usage \code(parGradUpdate.func(parList))
parGradUpdate.func <- function(parList) {
  gGam <- gGamma2.func(parList)
  gBeta <- gBeta2.func(parList)
  return(list(
    gGam = gGam,
    gBeta = gBeta,
    g = c(gGam, gBeta)
  ))
}


#' Calculate initial EM step length
#'
#' @description \code{alpha.func} is a helper function to compute the
#'   initial EM step length
#'
#' @usage \code{alpha.func(parList)}
alpha.func <- function(parList) {
  return(c(sign(t(parList$d) %*% parList$g$g)))
}


#' Calculate EM parameter updates
#'
#' @description \code{parUpdate.func} is a helper function to compute the EM
#'   parameter updates
#'
#' @usage \code(parUpdate.func(parList))
parUpdate.func <- function(parList) {
  parVec <- c(parList$gamList$gamVec, parList$betaVec) +
    parList$alpha * parList$d
  return(list(
    gamList = list(
      gamVec = parVec[1:parList$J],
      gamDot = sum(exp(parVec[1:parList$J]))
    ),
    betaVec = parVec[(parList$J + 1):(2 * parList$J)],
    parVec = parVec
  ))
}


#' Initialize EM steps
#'
#' @description \code{initEM.func} is a helper function to compute the EM
#'   parameter updates
#'
#' @usage \code(initEM.func(y ,slope, subCDF, depth, stepInit))
initEM.func <- function(y, slope, subCDF, depth, stepInit) {
  minD <- min(depth)
  if(minD < 0) {
    slope <- min(slope, (exp(minD - 3) - 1) / (exp(minD) - 1))
  } else {
    slope <- max(slope, (exp(minD - 3) - 1) / (exp(minD) - 1))
  }
  parList <- list(y = y,
                  depth = exp(depth) * slope - slope + 1)
  parList$J <- min(
    200,
    max(
      2,
      ceiling(sqrt(sum(parList$y > 0)))
    )
  )
  parList$gamList <- list(
    gamVec = rep(0, parList$J),
    gamDot = sum(exp(rep(0, parList$J)))
  )
  parList$zInd <- which.max(parList$gamList$gamVec)
  lBound <- max(1e-4, subCDF$eCDF_fun(-100))
  parList$betaVec <- subCDF$eQuant_fun(
    seq(lBound,
        1 - 1e-4,
        length = parList$J)
  )

  for(i in 1:stepInit) {
    parList$pList <- lpMat.func(parList)
    parList$gTilde <- gTilde.func(parList)
    parList$gamList$gamVec <- parList$gamList$gamVec + parList$gTilde$gamTilde
    parList$zInd <- which.max(parList$gamList$gamVec)
    parList$gamList$gamDot <- sum(exp(parList$gamList$gamVec))
    parList$betaVec <- parList$betaVec + parList$gTilde$betaTilde
  }

  parList$zInd <- which.max(parList$gamList$gamVec)

  parList$pList <- lpMat.func(parList)
  parList$g <- parGrad.func(parList)
  parList$gTilde <- gTilde.func(parList)
  parList$S <- matrix(0, 2 * parList$J, 2 * parList$J)

  parList$lLik <- lLik.func(parList)

  return(parList)
}


#' Update parList by removing parameter
#'
#' @description \code{parListUpdate.func} is a helper function to update
#'   parList when a parameter is to be removed
#'
#' @usage \code(parListUpdate.func(parList, doBeta = F))
parListUpdate.func <- function(parList, doBeta = F) {
  if(doBeta) {
    remInd <- which.min(parList$betaVec)
  } else {
    remInd <- which.min(parList$gamList$gamVec)
  }
  parList$J <- parList$J - 1
  parList$gamList$gamVec <- parList$gamList$gamVec[-remInd]
  parList$zInd <- which.max(parList$gamList$gamVec)
  parList$gamList$gamVec <- parList$gamList$gamVec -
    parList$gamList$gamVec[parList$zInd]
  parList$gamList$gamDot <- sum(exp(parList$gamList$gamVec))
  parList$betaVec <- parList$betaVec[-remInd]
  parList$pList <- lpMat.func(parList)
  parList$g <- parGrad.func(parList)
  parList$gTilde <- gTilde.func(parList)
  parList$S <- matrix(0, 2 * parList$J, 2 * parList$J)
  return(parList)
}

#' Calculate EM step direction/length
#'
#' @description \code{calcStep.func} is a helper function to calculate the
#'   step direction/length for the modified EM iteration
#'
#' @usage \code(calcStep.func(parList))
calcStep.func <- function(parList) {
  c1 <- 1e-3; c2 <- 0.8

  # step direction and initial step length
  parList$d <- -parList$gTilde$gTilde +
    parList$S %*% parList$g$g
  parList$alpha <- -1
  if(any(abs(parList$d) > 250)) {
    parList$alpha <- -250 / max(abs(parList$d))
  }

  # check initial Wolfe conditions
  parList$parUpdate <- parUpdate.func(parList)
  parList$lLikUpdate <- lLikUpdate.func(parList)
  w1 <- parList$lLikUpdate >= parList$lLik +
    c1 * parList$alpha * t(parList$d) %*% parList$g$g
  parList$pListUpdate <- lpMatUpdate.func(parList)
  parList$gUpdate <- parGradUpdate.func(parList)
  w2 <- abs(t(parList$d) %*% parList$gUpdate$g) <=
    c2 * abs(t(parList$d) %*% parList$g$g)

  if(w1 & w2) {
    # Return after conditions met
    return(parList)
  } else if(t(parList$g$g) %*% parList$S %*% parList$g$g > 0) {
    # Reset S matrix (indefinite) and take standard EM step
    parList$S <- parList$S * 0
    parList$D <- -parList$gTilde$gTilde
    parList$parUpdate <- parUpdate.func(parList)
    parList$lLikUpdate <- lLikUpdate.func(parList)
    parList$pListUpdate <- lpMatUpdate.func(parList)
    parList$gUpdate <- parGradUpdate.func(parList)
    parList$alpha <- -1
  } else {
    # Line search for step length
    alphaVec <- c(0, parList$alpha)
    alphaMin <- -100
    iter <- 1
    while(TRUE) {
      if(!w1 | (parList$lLikUpdate <= parList$lLik & iter > 1)) {
        parList$alpha <- zoom.func(parList, alphaVec[iter:(iter + 1)])
        break
      } else if (w2) {
        parList$alpha <- alphaVec[(iter + 1)]
        break
      } else if (t(parList$d) %*% parList$gUpdate$g >= 0) {
        parList$alpha <- zoom.func(parList, alphaVec[(iter + 1):iter])
        break
      } else {
        if(parList$alpha == alphaMin) {break}
        iter <- iter + 1
        alphaVec[iter + 1] <- max(alphaVec[iter] * 2, alphaMin)
        parList$alpha <- alphaVec[iter + 1]
        parList$parUpdate <- parUpdate.func(parList)
        parList$lLikUpdate <- lLikUpdate.func(parList)
        w1 <- parList$lLikUpdate >= parList$lLik +
          c1 * parList$alpha * t(parList$d) %*% parList$g$g
        parList$pListUpdate <- lpMatUpdate.func(parList)
        parList$gUpdate <- parGradUpdate.func(parList)
        w2 <- abs(t(parList$d) %*% parList$gUpdate$g) <=
          c2 * abs(t(parList$d) %*% parList$g$g)
      }
    }
  }

  if(is.null(parList$alpha)) {
    # Reset S matrix (indefinite) and take standard EM step
    parList$S <- parList$S * 0
    parList$D <- -parList$gTilde$gTilde
    parList$parUpdate <- parUpdate.func(parList)
    parList$lLikUpdate <- lLikUpdate.func(parList)
    parList$pListUpdate <- lpMatUpdate.func(parList)
    parList$gUpdate <- parGradUpdate.func(parList)
    parList$alpha <- -1
  }

  parList$parUpdate <- parUpdate.func(parList)
  parList$lLikUpdate <- lLikUpdate.func(parList)
  parList$pListUpdate <- lpMatUpdate.func(parList)
  parList$gUpdate <- parGradUpdate.func(parList)

  return(parList)
}


#' Find acceptable step length from bounds
#'
#' @description \code{zoom.func} is a helper function to calculate the
#'   acceptable step length from bounds
#'
#' @usage \code(zoom.func(parList, alphaBounds))
zoom.func <- function(parList, alphaBounds) {
  c1 <- 1e-3; c2 <- 0.8

  if(alphaBounds[1] != 0) {
    parList$alpha <- alphaBounds[1]
    parList$parUpdate <- parUpdate.func(parList)
    lLikLow <- lLikUpdate.func(parList)
  } else {
    lLikLow <- parList$lLik
  }
  iter <- 0
  while(T) {
    iter <- iter + 1
    if(iter == 10) {
      return(NULL)
    }
    parList$alpha <- mean(alphaBounds)
    parList$parUpdate <- parUpdate.func(parList)
    parList$lLikUpdate <- lLikUpdate.func(parList)
    w1 <- parList$lLikUpdate >= parList$lLik +
      c1 * parList$alpha * t(parList$d) %*% parList$g$g
    parList$pListUpdate <- lpMatUpdate.func(parList)
    parList$gUpdate <- parGradUpdate.func(parList)
    w2 <- abs(t(parList$d) %*% parList$gUpdate$g) <=
      c2 * abs(t(parList$d) %*% parList$g$g)
    if(!w1 | parList$lLikUpdate <= lLikLow) {
      alphaBounds[2] <- parList$alpha
    } else {
      if(w2) {
        return(parList$alpha)
      } else if(t(parList$d) %*% parList$gUpdate$g * diff(alphaBounds) <= 0) {
        alphaBounds[2] <- alphaBounds[1]
      }
      alphaBounds[1] <- parList$alpha
      lLikLow <- parList$lLikUpdate
    }
  }
}
