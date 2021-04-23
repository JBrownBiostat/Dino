# resampCounts is a helper function to normalize counts by resampling from an
# estimated expression distribution conditional on the observed counts
#
#' @importFrom Matrix Matrix
#' @importFrom stats approxfun
#' @importFrom stats bw.ucv
#' @importFrom stats bw.nrd
#' @importFrom stats ppois
#' @importFrom stats qgamma
#' @importFrom stats pgamma
#' @importFrom stats runif
#' @importFrom stats rgamma
resampCounts <- function(countMat, depth, depthRep, slope,
                         prec, subInd, doRQS, emPar) {
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
    parList <- initEM_func(
        y = y[subInd],
        slope = fitSlope,
        subCDF = cdfList,
        depth = depth[subInd],
        stepInit = 5,
        emPar
    )

    nUpdate <- emPar$maxIter
    deltaStop <- emPar$tol
    iter <- 0

    while(iter < nUpdate) {
        iter <- iter + 1

        parList$d <- -parList$gTilde$gTilde +
            parList$S %*% parList$g$g
        parList$alpha <- -1

        gamInd <- parList$gamList$gamVec +
            parList$alpha * parList$d[seq_len(parList$J)]
        if(!all(gamInd >= -75)) {
            if(parList$J > 2) {
                parList <- parListUpdate_func(parList)
                parList$lLik <- lLik_func(parList)
                parList$S <- parList$S * 0
            } else {break}
        }

        betaInd <- parList$betaVec +
            parList$alpha * parList$d[(parList$J + 1):(2 * parList$J)]
        if(!all(betaInd > -75)) {
            if(parList$J > 2) {
                parList <- parListUpdate_func(parList, doBeta = TRUE)
                parList$lLik <- lLik_func(parList)
                parList$S <- parList$S * 0
            } else {break}
        }
        if(any(gamInd > 75) | any(betaInd > 75) |
           max(gamInd) - min(gamInd) > log(2 * length(subInd))) {
            parList$S <- parList$S * 0
        }

        parList <- calcStep_func(parList)

        dTheta <- parList$parUpdate$parVec -
            c(parList$gamList$gamVec, parList$betaVec)
        dg <- parList$gUpdate$g - parList$g$g
        parList$dgTilde <- gTildeUpdate_func(parList)
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

    if(doRQS) {
        lamPar <- pmax(exp(parList$betaVec), .Machine$double.xmin * 1e1)
        bw <- tryCatch({
            suppressWarnings(bw.ucv(lamPar))
        }, error = function(e) {
            tryCatch({
                bw.nrd(lamPar)
            }, error = function(e) {
                1
            })
        })
        phi <- min(bw, 1)
        tau <- exp(parList$gamList$gamVec - max(parList$gamList$gamVec))
        tau <- tau / sum(tau)

        # Update model for all points
        parList$y <- y
        parList$depth <- exp(depth) * fitSlope - fitSlope + 1
        lowP <- rep(0, length(y))
        highP <- rep(0, length(y))
        nzInd <- y > 0
        for(i in seq_len(parList$J)) {
            lowP[nzInd] <- lowP[nzInd] + tau[i] *
                ppois(y[nzInd] - 1, lamPar[i] * parList$depth[nzInd])
            highP <- highP + tau[i] * ppois(y, lamPar[i] * parList$depth)
        }

        # Estimate mixture CDF
        lBound <- max(10 ^ (-prec),
                      qgamma(min(lowP), min(lamPar) / phi, scale = phi))
        hBound <- max(qgamma(max(highP), max(lamPar) / phi, scale = phi),
                      lBound + 2e-1)
        qSeq <- seq(lBound,
                    hBound,
                    min(1e-1))
        pSeq <- rep(0, length(qSeq))
        for(i in seq_len(parList$J)) {
            pSeq <- pSeq + tau[i] * pgamma(qSeq, lamPar[i] / phi, scale = phi)
        }
        cdfFun <- approxfun(pSeq, qSeq, yleft = 0, yright = max(qSeq))
        normVec <- cdfFun(runif(length(y), lowP, highP))
        normVec <- round(normVec, prec)

        return(normVec)
    } else {
        lamPar <- pmax(exp(parList$betaVec), .Machine$double.xmin * 1e1)
        bw <- tryCatch({
            suppressWarnings(bw.ucv(lamPar))
        }, error = function(e) {
            tryCatch({
                bw.nrd(lamPar)
            }, error = function(e) {
                1
            })
        })
        phi <- min(bw, 1)

        # Update model for all points
        parList$y <- y
        parList$depth <- exp(depth) * fitSlope - fitSlope + 1
        parList$pList <- lpMat_func(parList)
        parList$pList$spMat <- parList$pList$spMat * exp(parList$pList$lpMax)
        lamPar <- pmax(exp(parList$betaVec), .Machine$double.xmin * 1e1)

        # Resample counts
        distVec <- apply(parList$pList$spMat[-parList$J, , drop = FALSE], 2,
                         function(pVec) {
                             which.min(runif(1) > c(cumsum(pVec), 1))
                         })

        shapeVec <- lamPar / phi
        concVec <- emPar$conPar

        normVec <- rgamma(n = length(parList$y),
                          shape = parList$y * concVec + shapeVec[distVec],
                          scale = 1 / (parList$depth * concVec + 1 / phi))
        normVec <- round(normVec, prec)

        return(normVec)
        }
    }, depth = depth, depthRep = depthRep, slope = slope,
    prec = prec, subInd = subInd)
    subNorm <- Matrix(subNorm, sparse = TRUE)
    return(subNorm)
}


# parResampCounts is a helper function to parallelize the  resampling algorithm
#
#' @importFrom BiocParallel bplapply
#' @importFrom Matrix Matrix
parResampCounts <- function(countList, depth, depthRep, slope,
                            prec, subInd, prll, doRQS, emPar) {
    normList <- bplapply(countList,
                         function(subCounts, depth, depthRep,
                                  slope, prec, subInd) {
                             resampCounts(subCounts, depth, depthRep,
                                          slope, prec, subInd, doRQS, emPar)
                             }, depth = depth, depthRep = depthRep,
                         slope = slope, prec = prec, subInd = subInd,
                         BPPARAM = prll)
    while(length(normList) > 2) {
        for(i in seq_len(floor(length(normList) / 2))) {
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


# lLik_func is a helper function to compute the EM log-likelihood up to an
# additive constant
#
#' @importFrom matrixStats rowMaxs
lLik_func <- function(parList) {
    M <- outer(parList$y, parList$betaVec)
    M <- M - outer(parList$depth, exp(parList$betaVec))
    M <- M + outer(rep(1, length(parList$depth)), parList$gamList$gamVec)
    S <- rowMaxs(M)
    return(sum(log(rowSums(exp(M - S))) + S) -
               length(parList$y) * log(parList$gamList$gamDot))
}


# lLikFull_func is a helper function to compute the EM log-likelihood up to an
# additive constant
#
#' @importFrom matrixStats rowMaxs
lLikFull_func <- function(parList) {
    M <- outer(parList$y, parList$betaVec)
    M <- M - outer(parList$depth, exp(parList$betaVec))
    M <- M + outer(rep(1, length(parList$depth)), parList$gamList$gamVec)
    S <- rowMaxs(M)
    return(sum(log(rowSums(exp(M - S))) + S) -
               length(parList$y) * log(parList$gamList$gamDot) +
               sum(parList$y * log(parList$depth) - lfactorial(parList$y)))
}


# lLikUpdate_func is a helper function to compute the EM log-likelihood up to
# an addative constant at update pars
#
#' @importFrom matrixStats rowMaxs
lLikUpdate_func <- function(parList) {
    M <- outer(parList$y, parList$parUpdate$betaVec)
    M <- M - outer(parList$depth, exp(parList$parUpdate$betaVec))
    M <- M + outer(rep(1, length(parList$depth)),
                   parList$parUpdate$gamList$gamVec)
    S <- rowMaxs(M)
    return(sum(log(rowSums(exp(M - S))) + S) -
               length(parList$y) * log(parList$parUpdate$gamList$gamDot))
}


# lpMat_func is a helper function to compute the EM log membership probabilities
#
#' @importFrom matrixStats rowMaxs
#' @importFrom matrixStats colMaxs
lpMat_func <- function(parList) {
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


# lpMatUpdate_func is a helper function to compute the EM log membership
# probabilities at update pars
#
#' @importFrom matrixStats rowMaxs
#' @importFrom matrixStats colMaxs
lpMatUpdate_func <- function(parList) {
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


# gamUpdate_func is a helper function to compute the EM update to the gamma
# parameter
gamUpdate_func <- function(parList) {
    gamUpdate <- log(rowSums(parList$pList$spMat)) +
        parList$pList$lpMax
    return(gamUpdate - gamUpdate[parList$zInd])
}


# gamUpdate2_func is a helper function to compute the EM update to the gamma
# parameter at update pars
gamUpdate2_func <- function(parList) {
    gamUpdate <- log(rowSums(parList$pListUpdate$spMat)) +
        parList$pListUpdate$lpMax
    return(gamUpdate - gamUpdate[parList$zInd])
}


# betaUpdate_func is a helper function to compute the EM update to the beta
# parameter
betaUpdate_func <- function(parList) {
    betaUpdate <- log(pmax(parList$pList$spMat %*% parList$y,
                           sqrt(.Machine$double.xmin))) -
        log(pmax(parList$pList$spMat %*% parList$depth,
                 sqrt(.Machine$double.xmin)))
    return(c(betaUpdate))
}


# betaUpdate2_func is a helper function to compute the EM update to the beta
# parameter at update pars
betaUpdate2_func <- function(parList) {
    betaUpdate <- log(pmax(parList$pListUpdate$spMat %*% parList$y,
                           sqrt(.Machine$double.xmin))) -
        log(pmax(parList$pListUpdate$spMat %*% parList$depth,
                 sqrt(.Machine$double.xmin)))
    return(c(betaUpdate))
}


# gTilde_func is a helper function to compute the EM parameter update deltas
gTilde_func <- function(parList) {
    gamUpdate <- gamUpdate_func(parList)
    betaUpdate <- betaUpdate_func(parList)
    return(list(
        gamTilde = gamUpdate - parList$gamList$gamVec,
        betaTilde = betaUpdate - parList$betaVec,
        gTilde = c(gamUpdate - parList$gamList$gamVec,
                   betaUpdate - parList$betaVec)
    ))
}


# gTildeUpdate_func is a helper function to compute the EM parameter update
# deltas at update pars
gTildeUpdate_func <- function(parList) {
    gamUpdate <- gamUpdate2_func(parList)
    betaUpdate <- betaUpdate2_func(parList)
    return(list(
        gamTilde = gamUpdate - parList$parUpdate$gamList$gamVec,
        betaTilde = betaUpdate - parList$parUpdate$betaVec,
        gTilde = c(gamUpdate - parList$parUpdate$gamList$gamVec,
                   betaUpdate - parList$parUpdate$betaVec)
    ))
}


# gGamma_func is a helper function to compute the likelihood gradient in gamma
gGamma_func <- function(parList) {
    M <- rowSums(parList$pList$spMat) * exp(parList$pList$lpMax)
    M <- M - sum(M) * exp(parList$gamList$gamVec - log(parList$gamList$gamDot))
    M[parList$zInd] <- 0
    return(M)
}


# gGamma2_func is a helper function to compute the likelihood gradient in
# gamma at update pars
gGamma2_func <- function(parList) {
    M <- rowSums(parList$pListUpdate$spMat) * exp(parList$pListUpdate$lpMax)
    M <- M - sum(M) * exp(parList$parUpdate$gamList$gamVec -
                              log(parList$parUpdate$gamList$gamDot))
    M[parList$zInd] <- 0
    return(M)
}


# gBeta_func is a helper function to compute the likelihood gradient in beta
gBeta_func <- function(parList) {
    M <- parList$pList$spMat %*% parList$y -
        rowSums(parList$pList$spMat *
                    outer(exp(parList$betaVec), parList$depth))
    return(c(exp(parList$pList$lpMax) * M))
}


# gBeta2_func is a helper function to compute the likelihood gradient in beta
# at update pars
gBeta2_func <- function(parList) {
    M <- parList$pListUpdate$spMat %*% parList$y -
        rowSums(parList$pListUpdate$spMat *
                    outer(exp(parList$parUpdate$betaVec), parList$depth))
    return(c(exp(parList$pListUpdate$lpMax) * M))
}


# parGrad_func is a helper function to compute the EM parameter gradients
parGrad_func <- function(parList) {
    gGam <- gGamma_func(parList)
    gBeta <- gBeta_func(parList)
    return(list(
        gGam = gGam,
        gBeta = gBeta,
        g = c(gGam, gBeta)
    ))
}


# parGradUpdate_func is a helper function to compute the EM parameter
# gradients at update pars
parGradUpdate_func <- function(parList) {
    gGam <- gGamma2_func(parList)
    gBeta <- gBeta2_func(parList)
    return(list(
        gGam = gGam,
        gBeta = gBeta,
        g = c(gGam, gBeta)
    ))
}


# alpha_func is a helper function to compute the initial EM step length
alpha_func <- function(parList) {
    return(c(sign(t(parList$d) %*% parList$g$g)))
}


# parUpdate_func is a helper function to compute the EM parameter updates
parUpdate_func <- function(parList) {
    parVec <- c(parList$gamList$gamVec, parList$betaVec) +
        parList$alpha * parList$d
    return(list(
        gamList = list(
            gamVec = parVec[seq_len(parList$J)],
            gamDot = sum(exp(parVec[seq_len(parList$J)]))
            ),
        betaVec = parVec[(parList$J + 1):(2 * parList$J)],
        parVec = parVec
    ))
}


# initEM_func is a helper function to compute the EM parameter updates
initEM_func <- function(y, slope, subCDF, depth, stepInit, emPar) {
    minD <- min(depth)
    if(minD < 0) {
        slope <- min(slope, (exp(minD - 3) - 1) / (exp(minD) - 1))
    } else {
        slope <- max(slope, (exp(minD - 3) - 1) / (exp(minD) - 1))
    }
    parList <- list(y = y,
                    depth = exp(depth) * slope - slope + 1)
    parList$J <- min(
        emPar$maxK,
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

    for(i in seq_len(stepInit)) {
        parList$pList <- lpMat_func(parList)
        parList$gTilde <- gTilde_func(parList)
        parList$gamList$gamVec <- parList$gamList$gamVec +
            parList$gTilde$gamTilde
        parList$zInd <- which.max(parList$gamList$gamVec)
        parList$gamList$gamDot <- sum(exp(parList$gamList$gamVec))
        parList$betaVec <- parList$betaVec + parList$gTilde$betaTilde
    }

    parList$zInd <- which.max(parList$gamList$gamVec)

    parList$pList <- lpMat_func(parList)
    parList$g <- parGrad_func(parList)
    parList$gTilde <- gTilde_func(parList)
    parList$S <- matrix(0, 2 * parList$J, 2 * parList$J)

    parList$lLik <- lLik_func(parList)

    return(parList)
}


# parListUpdate_func is a helper function to update parList when a parameter
# is to be removed
parListUpdate_func <- function(parList, doBeta = FALSE) {
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
    parList$pList <- lpMat_func(parList)
    parList$g <- parGrad_func(parList)
    parList$gTilde <- gTilde_func(parList)
    parList$S <- matrix(0, 2 * parList$J, 2 * parList$J)
    return(parList)
}

# calcStep_func is a helper function to calculate the step direction/length
# for the modified EM iteration
calcStep_func <- function(parList) {
    c1 <- 1e-3; c2 <- 0.8

    # step direction and initial step length
    parList$d <- -parList$gTilde$gTilde +
        parList$S %*% parList$g$g
    parList$alpha <- -1
    if(any(abs(parList$d) > 250)) {
        parList$alpha <- -250 / max(abs(parList$d))
    }

    # check initial Wolfe conditions
    parList$parUpdate <- parUpdate_func(parList)
    parList$lLikUpdate <- lLikUpdate_func(parList)
    w1 <- parList$lLikUpdate >= parList$lLik +
        c1 * parList$alpha * t(parList$d) %*% parList$g$g
    parList$pListUpdate <- lpMatUpdate_func(parList)
    parList$gUpdate <- parGradUpdate_func(parList)
    w2 <- abs(t(parList$d) %*% parList$gUpdate$g) <=
        c2 * abs(t(parList$d) %*% parList$g$g)

    if(w1 & w2) {
        # Return after conditions met
        return(parList)
    } else if(t(parList$g$g) %*% parList$S %*% parList$g$g > 0) {
        # Reset S matrix (indefinite) and take standard EM step
        parList$S <- parList$S * 0
        parList$D <- -parList$gTilde$gTilde
        parList$parUpdate <- parUpdate_func(parList)
        parList$lLikUpdate <- lLikUpdate_func(parList)
        parList$pListUpdate <- lpMatUpdate_func(parList)
        parList$gUpdate <- parGradUpdate_func(parList)
        parList$alpha <- -1
    } else {
        # Line search for step length
        alphaVec <- c(0, parList$alpha)
        alphaMin <- -100
        iter <- 1
        while(TRUE) {
            if(!w1 | (parList$lLikUpdate <= parList$lLik & iter > 1)) {
                parList$alpha <- zoom_func(parList, alphaVec[iter:(iter + 1)])
                break
            } else if (w2) {
                parList$alpha <- alphaVec[(iter + 1)]
                break
            } else if (t(parList$d) %*% parList$gUpdate$g >= 0) {
                parList$alpha <- zoom_func(parList, alphaVec[(iter + 1):iter])
                break
            } else {
                if(parList$alpha == alphaMin) {break}
                iter <- iter + 1
                alphaVec[iter + 1] <- max(alphaVec[iter] * 2, alphaMin)
                parList$alpha <- alphaVec[iter + 1]
                parList$parUpdate <- parUpdate_func(parList)
                parList$lLikUpdate <- lLikUpdate_func(parList)
                w1 <- parList$lLikUpdate >= parList$lLik +
                    c1 * parList$alpha * t(parList$d) %*% parList$g$g
                parList$pListUpdate <- lpMatUpdate_func(parList)
                parList$gUpdate <- parGradUpdate_func(parList)
                w2 <- abs(t(parList$d) %*% parList$gUpdate$g) <=
                    c2 * abs(t(parList$d) %*% parList$g$g)
            }
        }
    }

    if(is.null(parList$alpha)) {
        # Reset S matrix (indefinite) and take standard EM step
        parList$S <- parList$S * 0
        parList$D <- -parList$gTilde$gTilde
        parList$parUpdate <- parUpdate_func(parList)
        parList$lLikUpdate <- lLikUpdate_func(parList)
        parList$pListUpdate <- lpMatUpdate_func(parList)
        parList$gUpdate <- parGradUpdate_func(parList)
        parList$alpha <- -1
    }

    parList$parUpdate <- parUpdate_func(parList)
    parList$lLikUpdate <- lLikUpdate_func(parList)
    parList$pListUpdate <- lpMatUpdate_func(parList)
    parList$gUpdate <- parGradUpdate_func(parList)

    return(parList)
}


# zoom_func is a helper function to calculate the acceptable step length
# from bounds
zoom_func <- function(parList, alphaBounds) {
    c1 <- 1e-3; c2 <- 0.8

    if(alphaBounds[1] != 0) {
        parList$alpha <- alphaBounds[1]
        parList$parUpdate <- parUpdate_func(parList)
        lLikLow <- lLikUpdate_func(parList)
    } else {
        lLikLow <- parList$lLik
    }
    iter <- 0
    while(TRUE) {
        iter <- iter + 1
        if(iter == 10) {
            return(NULL)
        }
        parList$alpha <- mean(alphaBounds)
        parList$parUpdate <- parUpdate_func(parList)
        parList$lLikUpdate <- lLikUpdate_func(parList)
        w1 <- parList$lLikUpdate >= parList$lLik +
            c1 * parList$alpha * t(parList$d) %*% parList$g$g
        parList$pListUpdate <- lpMatUpdate_func(parList)
        parList$gUpdate <- parGradUpdate_func(parList)
        w2 <- abs(t(parList$d) %*% parList$gUpdate$g) <=
            c2 * abs(t(parList$d) %*% parList$g$g)
        if(!w1 | parList$lLikUpdate <= lLikLow) {
            alphaBounds[2] <- parList$alpha
        } else {
            if(w2) {
                return(parList$alpha)
            } else if(t(parList$d) %*% parList$gUpdate$g * diff(alphaBounds) <=
                      0) {
                alphaBounds[2] <- alphaBounds[1]
            }
            alphaBounds[1] <- parList$alpha
            lLikLow <- parList$lLikUpdate
        }
    }
}

