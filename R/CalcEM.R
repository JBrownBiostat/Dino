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
    subNorm <- apply(countMat, 2,
                function(y, depth, depthRep, slope, prec, subInd) {
        # Initialize EM parameters
        parList <- emParInit_func(y, depth, depthRep, slope, subInd, emPar)
        fitSlope <- parList$fitSlope
        parList <- parList$parList

        # Fit EM parameters
        parList <- fitEM_func(parList, emPar, subInd)

        # Resample normalized expression values
        if(doRQS) {
            return(rqsSample_func(parList, y, depth, fitSlope, prec))
        } else {
            return(posteriorSample_func(
                parList, y, depth, fitSlope, prec, emPar
            ))
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
                    function(subCounts, depth, depthRep, slope, prec, subInd) {
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
    M <- outer(em_y(parList), em_betaVec(parList))
    M <- M - outer(em_depth(parList), exp(em_betaVec(parList)))
    M <- M +
        outer(rep(1, length(em_depth(parList))), em_gamList(parList)$gamVec)
    S <- rowMaxs(M)
    return(sum(log(rowSums(exp(M - S))) + S) -
        length(em_y(parList)) * log(em_gamList(parList)$gamDot))
}


# lLikFull_func is a helper function to compute the EM log-likelihood up to an
# additive constant
#
#' @importFrom matrixStats rowMaxs
lLikFull_func <- function(parList) {
    M <- outer(em_y(parList), em_betaVec(parList))
    M <- M - outer(em_depth(parList), exp(em_betaVec(parList)))
    M <- M +
        outer(rep(1, length(em_depth(parList))), em_gamList(parList)$gamVec)
    S <- rowMaxs(M)
    return(sum(log(rowSums(exp(M - S))) + S) -
            length(em_y(parList)) * log(em_gamList(parList)$gamDot) +
            sum(
                em_y(parList) * log(em_depth(parList)) -
                    lfactorial(em_y(parList))
            ))
}


# lLikUpdate_func is a helper function to compute the EM log-likelihood up to
# an addative constant at update pars
#
#' @importFrom matrixStats rowMaxs
lLikUpdate_func <- function(parList) {
    M <- outer(em_y(parList), em_parUpdate(parList)$betaVec)
    M <- M - outer(em_depth(parList), exp(em_parUpdate(parList)$betaVec))
    M <- M + outer(rep(1, length(em_depth(parList))),
                    em_parUpdate(parList)$gamList$gamVec)
    S <- rowMaxs(M)
    return(sum(log(rowSums(exp(M - S))) + S) -
            length(em_y(parList)) * log(em_parUpdate(parList)$gamList$gamDot))
}


# lpMat_func is a helper function to compute the EM log membership probabilities
#
#' @importFrom matrixStats rowMaxs
#' @importFrom matrixStats colMaxs
lpMat_func <- function(parList) {
    lpMat <- outer(em_betaVec(parList), em_y(parList))
    lpMat <- lpMat - outer(exp(em_betaVec(parList)), em_depth(parList))
    lpMat <- lpMat + em_gamList(parList)$gamVec
    lpMat <- lpMat - outer(rep(1, em_J(parList)), colMaxs(lpMat))
    lpMat <- lpMat - outer(rep(1, em_J(parList)), log(colSums(exp(lpMat))))
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
    lpMat <- outer(em_parUpdate(parList)$betaVec, em_y(parList))
    lpMat <- lpMat -
        outer(exp(em_parUpdate(parList)$betaVec), em_depth(parList))
    lpMat <- lpMat + em_parUpdate(parList)$gamList$gamVec
    lpMat <- lpMat - outer(rep(1, em_J(parList)), colMaxs(lpMat))
    lpMat <- lpMat - outer(rep(1, em_J(parList)), log(colSums(exp(lpMat))))
    lpMax <- rowMaxs(lpMat)
    return(list(
        spMat = exp(lpMat - lpMax),
        lpMax = lpMax
    ))
}


# gamUpdate_func is a helper function to compute the EM update to the gamma
# parameter
gamUpdate_func <- function(parList) {
    gamUpdate <- log(rowSums(em_pList(parList)$spMat)) +
        em_pList(parList)$lpMax
    return(gamUpdate - gamUpdate[em_zInd(parList)])
}


# gamUpdate2_func is a helper function to compute the EM update to the gamma
# parameter at update pars
gamUpdate2_func <- function(parList) {
    gamUpdate <- log(rowSums(em_pListUpdate(parList)$spMat)) +
        em_pListUpdate(parList)$lpMax
    return(gamUpdate - gamUpdate[em_zInd(parList)])
}


# betaUpdate_func is a helper function to compute the EM update to the beta
# parameter
betaUpdate_func <- function(parList) {
    betaUpdate <- log(pmax(em_pList(parList)$spMat %*% em_y(parList),
                            sqrt(.Machine$double.xmin))) -
        log(pmax(em_pList(parList)$spMat %*% em_depth(parList),
                    sqrt(.Machine$double.xmin)))
    return(c(betaUpdate))
}


# betaUpdate2_func is a helper function to compute the EM update to the beta
# parameter at update pars
betaUpdate2_func <- function(parList) {
    betaUpdate <- log(pmax(em_pListUpdate(parList)$spMat %*% em_y(parList),
                            sqrt(.Machine$double.xmin))) -
        log(pmax(em_pListUpdate(parList)$spMat %*% em_depth(parList),
                    sqrt(.Machine$double.xmin)))
    return(c(betaUpdate))
}


# gTilde_func is a helper function to compute the EM parameter update deltas
gTilde_func <- function(parList) {
    gamUpdate <- gamUpdate_func(parList)
    betaUpdate <- betaUpdate_func(parList)
    return(list(
        gamTilde = gamUpdate - em_gamList(parList)$gamVec,
        betaTilde = betaUpdate - em_betaVec(parList),
        gTilde = c(gamUpdate - em_gamList(parList)$gamVec,
                    betaUpdate - em_betaVec(parList))
    ))
}


# gTildeUpdate_func is a helper function to compute the EM parameter update
# deltas at update pars
gTildeUpdate_func <- function(parList) {
    gamUpdate <- gamUpdate2_func(parList)
    betaUpdate <- betaUpdate2_func(parList)
    return(list(
        gamTilde = gamUpdate - em_parUpdate(parList)$gamList$gamVec,
        betaTilde = betaUpdate - em_parUpdate(parList)$betaVec,
        gTilde = c(gamUpdate - em_parUpdate(parList)$gamList$gamVec,
                    betaUpdate - em_parUpdate(parList)$betaVec)
    ))
}


# gGamma_func is a helper function to compute the likelihood gradient in gamma
gGamma_func <- function(parList) {
    M <- rowSums(em_pList(parList)$spMat) * exp(em_pList(parList)$lpMax)
    M <- M - sum(M) *
        exp(em_gamList(parList)$gamVec - log(em_gamList(parList)$gamDot))
    M[em_zInd(parList)] <- 0
    return(M)
}


# gGamma2_func is a helper function to compute the likelihood gradient in
# gamma at update pars
gGamma2_func <- function(parList) {
    M <- rowSums(em_pListUpdate(parList)$spMat) *
        exp(em_pListUpdate(parList)$lpMax)
    M <- M - sum(M) * exp(em_parUpdate(parList)$gamList$gamVec -
                            log(em_parUpdate(parList)$gamList$gamDot))
    M[em_zInd(parList)] <- 0
    return(M)
}


# gBeta_func is a helper function to compute the likelihood gradient in beta
gBeta_func <- function(parList) {
    M <- em_pList(parList)$spMat %*% em_y(parList) -
        rowSums(em_pList(parList)$spMat *
                outer(exp(em_betaVec(parList)), em_depth(parList)))
    return(c(exp(em_pList(parList)$lpMax) * M))
}


# gBeta2_func is a helper function to compute the likelihood gradient in beta
# at update pars
gBeta2_func <- function(parList) {
    M <- em_pListUpdate(parList)$spMat %*% em_y(parList) -
        rowSums(em_pListUpdate(parList)$spMat *
                outer(exp(em_parUpdate(parList)$betaVec), em_depth(parList)))
    return(c(exp(em_pListUpdate(parList)$lpMax) * M))
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
    return(c(sign(t(em_d(parList)) %*% em_g(parList)$g)))
}


# parUpdate_func is a helper function to compute the EM parameter updates
parUpdate_func <- function(parList) {
    parVec <- c(em_gamList(parList)$gamVec, em_betaVec(parList)) +
        em_alpha(parList) * em_d(parList)
    return(list(
        gamList = list(
            gamVec = parVec[seq_len(em_J(parList))],
            gamDot = sum(exp(parVec[seq_len(em_J(parList))]))
            ),
        betaVec = parVec[(em_J(parList) + 1):(2 * em_J(parList))],
        parVec = parVec
    ))
}


# emParInit_func is a helper function to initialize the EM parameters from
# raw cout data
emParInit_func <- function(y, depth, depthRep, slope, subInd, emPar) {
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

    # Bound slope to exp(-3) times minimum depth
    minD <- min(depth)
    if(minD < 0) {
        fitSlope <- min(slope, (exp(minD - 3) - 1) / (exp(minD) - 1))
    } else {
        fitSlope <- max(slope, (exp(minD - 3) - 1) / (exp(minD) - 1))
    }

    # Initialize EM parameters with 5 strict EM steps
    parList <- initEM_func(
        y = y[subInd],
        slope = fitSlope,
        subCDF = cdfList,
        depth = depth[subInd],
        stepInit = 5,
        emPar
    )

    return(list(
        parList = parList,
        fitSlope = fitSlope
    ))
}


# initEM_func is a helper function to compute the EM parameter updates
#' @importFrom methods new
initEM_func <- function(y, slope, subCDF, depth, stepInit, emPar) {
    # Bound slope to exp(-3) times minimum depth
    minD <- min(depth)
    if(minD < 0) {
        slope <- min(slope, (exp(minD - 3) - 1) / (exp(minD) - 1))
    } else {
        slope <- max(slope, (exp(minD - 3) - 1) / (exp(minD) - 1))
    }
    parList <- new(
        "emFitPars", y = y, depth = exp(depth) * slope - slope + 1,
        J = min(emPar$maxK, max(2, ceiling(sqrt(sum(y > 0)))))
    )
    parList <- em_gamListSet(parList, list(
        gamVec = rep(0, em_J(parList)),
        gamDot = sum(exp(rep(0, em_J(parList))))
    ))
    parList <- em_zIndSet(parList, which.max(em_gamList(parList)$gamVec))
    # bound initial beta values within eCDF minus 1e-4 on either side
    lBound <- max(1e-4, subCDF$eCDF_fun(-100))
    parList <- em_betaVecSet(parList, subCDF$eQuant_fun(
        seq(lBound, 1 - 1e-4, length = em_J(parList))
    ))

    for(i in seq_len(stepInit)) {
        parList <- em_pListSet(parList, lpMat_func(parList))
        parList <- em_gTildeSet(parList, gTilde_func(parList))
        gamList <- em_gamList(parList)
        gamList$gamVec <- em_gamList(parList)$gamVec +
            em_gTilde(parList)$gamTilde
        parrList <- em_gamListSet(parList, gamList)
        parList <- em_zIndSet(parList, which.max(em_gamList(parList)$gamVec))
        gamList <- em_gamList(parList)
        gamList$gamDot <- sum(exp(em_gamList(parList)$gamVec))
        parList  <- em_gamListSet(parList, gamList)
        parList <- em_betaVecSet(
            parList, em_betaVec(parList) + em_gTilde(parList)$betaTilde
        )
    }

    parList <- em_zIndSet(parList, which.max(em_gamList(parList)$gamVec))

    parList <- em_pListSet(parList, lpMat_func(parList))
    parList <- em_gSet(parList, parGrad_func(parList))
    parList <- em_gTildeSet(parList, gTilde_func(parList))
    parList <- em_SSet(parList, matrix(0, 2 * em_J(parList), 2 * em_J(parList)))

    parList <- em_lLikSet(parList, lLik_func(parList))

    return(parList)
}


# fitEM_func is a helper function to fit the EM parameter estimates
fitEM_func <- function(parList, emPar, subInd) {
    nUpdate <- emPar$maxIter
    deltaStop <- emPar$tol
    iter <- 0

    while(iter < nUpdate) {
        iter <- iter + 1

        parList <- em_dSet(
            parList, as.numeric(
                -em_gTilde(parList)$gTilde + em_S(parList) %*% em_g(parList)$g
            )
        )
        parList <- em_alphaSet(parList, -1)

        # Sanity check on bounds of k+1 parameter estimates
        parList <- emParCheck_func(parList, subInd)
        if(parList$endLoop) {
            parList <- parList$parList
            break
        } else {
            parList <- parList$parList
        }

        # Accelerated EM step
        parList <- takeEMStep_func(parList)

        # Check end condition
        if(em_lLikUpdate(parList) - em_lLik(parList) <= deltaStop) {
            break
        }
        parList <- em_lLikSet(parList, em_lLikUpdate(parList))
    }

    return(parList)
}


# emParCheck_func is a helper function to perform sanity checks on  projected
# EM update values
emParCheck_func <- function(parList, subInd) {
    # Sanity check on bounds of k+1 parameter estimates
    endLoop <- FALSE

    gamInd <- em_gamList(parList)$gamVec +
        em_alpha(parList) * em_d(parList)[seq_len(em_J(parList))]
    if(!all(gamInd >= -75)) {
        if(em_J(parList) > 2) {
            parList <- parListUpdate_func(parList)
            parList <- em_lLikSet(parList, lLik_func(parList))
            parList <- em_SSet(parList, em_S(parList) * 0)
        } else {endLoop <- TRUE}
    }

    betaInd <- em_betaVec(parList) +
        em_alpha(parList) *
        em_d(parList)[(em_J(parList) + 1):(2 * em_J(parList))]
    if(!all(betaInd > -75)) {
        if(em_J(parList) > 2) {
            parList <- parListUpdate_func(parList, doBeta = TRUE)
            parList <- em_lLikSet(parList,  lLik_func(parList))
            parList <- em_SSet(parList, em_S(parList) * 0)
        } else {endLoop <- TRUE}
    }
    if(any(gamInd > 75) | any(betaInd > 75) |
       max(gamInd) - min(gamInd) > log(2 * length(subInd))) {
        parList <- em_SSet(parList, em_S(parList) * 0)
    }

    return(list(
        parList = parList,
        endLoop = endLoop
    ))
}


# parListUpdate_func is a helper function to update parList when a parameter
# is to be removed
parListUpdate_func <- function(parList, doBeta = FALSE) {
    if(doBeta) {
        remInd <- which.min(em_betaVec(parList))
    } else {
        remInd <- which.min(em_gamList(parList)$gamVec)
    }
    parList <- em_JSet(parList, em_J(parList) - 1)
    gamList <- em_gamList(parList)
    gamList$gamVec <- gamList$gamVec[-remInd]
    parList <- em_zIndSet(parList, which.max(gamList$gamVec))
    gamList$gamVec <- gamList$gamVec -
        gamList$gamVec[em_zInd(parList)]
    gamList$gamDot <- sum(exp(gamList$gamVec))
    parList <- em_gamListSet(parList, gamList)
    parList <- em_betaVecSet(parList, em_betaVec(parList)[-remInd])
    parList <- em_pListSet(parList, lpMat_func(parList))
    parList <- em_gSet(parList, parGrad_func(parList))
    parList <- em_gTildeSet(parList, gTilde_func(parList))
    parList <- em_SSet(parList, matrix(0, 2 * em_J(parList), 2 * em_J(parList)))
    return(parList)
}


# takeEMStep_func is a helper function to update parList to the next iterate
takeEMStep_func <- function(parList) {
    parList <- calcStep_func(parList)

    dTheta <- em_parUpdate(parList)$parVec -
        c(em_gamList(parList)$gamVec, em_betaVec(parList))
    dg <- em_gUpdate(parList)$g - em_g(parList)$g
    parList <- em_dgTildeSet(parList, gTildeUpdate_func(parList))
    dgTilde <- em_dgTilde(parList)$gTilde - em_gTilde(parList)$gTilde

    dThetaStar <- -dgTilde + em_S(parList) %*% dg
    parList <- em_SSet(
        parList,
        em_S(parList) + c(1 + (t(dg) %*% dThetaStar) / (t(dg) %*% dTheta)) *
            (dTheta %*% t(dTheta)) / c(t(dg) %*% dTheta) -
            (dThetaStar %*% t(dTheta) + t(dThetaStar %*% t(dTheta))) /
            c(t(dg) %*% dTheta)
    )

    parList <- em_pListSet(parList, em_pListUpdate(parList))
    parList <- em_gamListSet(parList, em_parUpdate(parList)$gamList)
    parList <- em_betaVecSet(parList, em_parUpdate(parList)$betaVec)
    parList <- em_gSet(parList, em_gUpdate(parList))
    parList <- em_gTildeSet(parList, em_dgTilde(parList))

    return(parList)
}


# calcStep_func is a helper function to calculate the step direction/length
# for the modified EM iteration
calcStep_func <- function(parList) {
    # Control parameters for strong Wolfe line search
    # c1: small, strictly positive
    # c2: large, strictly less than 1, strictly greater than c1
    c1 <- 1e-3; c2 <- 0.8

    # step direction and initial step length
    parList <- em_dSet(
        parList, as.numeric(
            -em_gTilde(parList)$gTilde + em_S(parList) %*% em_g(parList)$g
        )
    )
    # Sanity check on step length
    parList <- em_alphaSet(parList, -1)
    if(any(abs(em_d(parList)) > 250)) {
        parList <- em_alphaSet(parList, -250 / max(abs(em_d(parList))))
    }

    # check initial Wolfe conditions
    parList <- em_parUpdateSet(parList, parUpdate_func(parList))
    parList <- em_lLikUpdateSet(parList, lLikUpdate_func(parList))
    w1 <- em_lLikUpdate(parList) >= em_lLik(parList) +
        c1 * em_alpha(parList) * t(em_d(parList)) %*% em_g(parList)$g
    parList <- em_pListUpdateSet(parList, lpMatUpdate_func(parList))
    parList <- em_gUpdateSet(parList, parGradUpdate_func(parList))
    w2 <- abs(t(em_d(parList)) %*% em_gUpdate(parList)$g) <=
        c2 * abs(t(em_d(parList)) %*% em_g(parList)$g)

    if(w1 & w2) {
        # Return after conditions met
        return(parList)
    } else if(t(em_g(parList)$g) %*% em_S(parList) %*% em_g(parList)$g > 0) {
        # Reset S matrix (indefinite) and take standard EM step
        parList <- resetS_func(parList)
    } else {
        parList <- alphaLine_func(parList, w1, w2, c1, c2)
    }

    if(is.null(em_alpha(parList))) {
        # Reset S matrix (indefinite) and take standard EM step
        parList <- resetS_func(parList)
    }

    parList <- em_parUpdateSet(parList, parUpdate_func(parList))
    parList <- em_lLikUpdateSet(parList, lLikUpdate_func(parList))
    parList <- em_pListUpdateSet(parList, lpMatUpdate_func(parList))
    parList <- em_gUpdateSet(parList, parGradUpdate_func(parList))

    return(parList)
}


# zoom_func is a helper function to calculate the acceptable step length
# from bounds
zoom_func <- function(parList, alphaBounds) {
    # Control parameters for strong Wolfe line search
    # c1: small, strictly positive
    # c2: large, strictly less than 1, strictly greater than c1
    c1 <- 1e-3; c2 <- 0.8

    if(alphaBounds[1] != 0) {
        parList <- em_alphaSet(parList, alphaBounds[1])
        parList <- em_parUpdateSet(parList, parUpdate_func(parList))
        lLikLow <- lLikUpdate_func(parList)
    } else {
        lLikLow <- em_lLik(parList)
    }
    iter <- 0
    while(TRUE) {
        iter <- iter + 1
        if(iter == 10) {
            return(NULL)
        }
        parList <- em_alphaSet(parList, mean(alphaBounds))
        parList <- em_parUpdateSet(parList, parUpdate_func(parList))
        parList <- em_lLikUpdateSet(parList, lLikUpdate_func(parList))
        w1 <- em_lLikUpdate(parList) >= em_lLik(parList) +
            c1 * em_alpha(parList) * t(em_d(parList)) %*% em_g(parList)$g
        parList <- em_pListUpdateSet(parList, lpMatUpdate_func(parList))
        parList <- em_gUpdateSet(parList, parGradUpdate_func(parList))
        w2 <- abs(t(em_d(parList)) %*% em_gUpdate(parList)$g) <=
            c2 * abs(t(em_d(parList)) %*% em_g(parList)$g)
        if(!w1 | em_lLikUpdate(parList) <= lLikLow) {
            alphaBounds[2] <- em_alpha(parList)
        } else {
            if(w2) {
                return(em_alpha(parList))
            } else if(
                t(em_d(parList)) %*% em_gUpdate(parList)$g *
                diff(alphaBounds) <= 0
            ) {
                alphaBounds[2] <- alphaBounds[1]
            }
            alphaBounds[1] <- em_alpha(parList)
            lLikLow <- em_lLikUpdate(parList)
        }
    }
}


# resetS_func is a helper function to reset the EM acceleration with a strict
# EM step
resetS_func <- function(parList) {
    # Reset S matrix (indefinite) and take standard EM step
    parList <- em_SSet(parList, em_S(parList) * 0)
    parList <- em_dSet(parList, as.numeric(-em_gTilde(parList)$gTilde))
    parList <- em_parUpdateSet(parList, parUpdate_func(parList))
    parList <- em_lLikUpdateSet(parList, lLikUpdate_func(parList))
    parList <- em_pListUpdateSet(parList, lpMatUpdate_func(parList))
    parList <- em_gUpdateSet(parList, parGradUpdate_func(parList))
    parList <- em_alphaSet(parList, -1)

    return(parList)
}


# alphaLine_func is a helper function to perform line-search for the EM update
alphaLine_func <- function(parList, w1, w2, c1, c2) {
    # Line search for step length
    alphaVec <- c(0, em_alpha(parList))
    # Sanity bound on step size
    alphaMin <- -100
    iter <- 1
    while(TRUE) {
        if(!w1 | (em_lLikUpdate(parList) <= em_lLik(parList) & iter > 1)) {
            parList <- em_alphaSet(parList, zoom_func(
                parList, alphaVec[iter:(iter + 1)]
            ))
            break
        } else if (w2) {
            parList <- em_alphaSet(parList, alphaVec[(iter + 1)])
            break
        } else if (t(em_d(parList)) %*% em_gUpdate(parList)$g >= 0) {
            parList <- em_alphaSet(parList, zoom_func(
                parList, alphaVec[(iter + 1):iter]
            ))
            break
        } else {
            if(em_alpha(parList) == alphaMin) {break}
            iter <- iter + 1
            alphaVec[iter + 1] <- max(alphaVec[iter] * 2, alphaMin)
            parList <- em_alphaSet(parList, alphaVec[iter + 1])
            parList <- em_parUpdateSet(parList, parUpdate_func(parList))
            parList <- em_lLikUpdateSet(parList, lLikUpdate_func(parList))
            w1 <- em_lLikUpdate(parList) >= em_lLik(parList) +
                c1 * em_alpha(parList) *
                t(em_d(parList)) %*% em_g(parList)$g
            parList <- em_pListUpdateSet(parList, lpMatUpdate_func(parList))
            parList <- em_gUpdateSet(parList, parGradUpdate_func(parList))
            w2 <- abs(t(em_d(parList)) %*% em_gUpdate(parList)$g) <=
                c2 * abs(t(em_d(parList)) %*% em_g(parList)$g)
        }
    }

    return(parList)
}


# rqsSample_func is a helper function to perform restricted quantile sampling
# given fitted EM parameters
rqsSample_func <- function(parList, y, depth, fitSlope, prec) {
    lamPar <- pmax(exp(em_betaVec(parList)), .Machine$double.xmin * 1e1)
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
    tau <- exp(em_gamList(parList)$gamVec - max(em_gamList(parList)$gamVec))
    tau <- tau / sum(tau)

    # Update model for all points
    parList <- em_ySet(parList, y)
    parList <- em_depthSet(parList, exp(depth) * fitSlope - fitSlope + 1)
    lowP <- rep(0, length(y))
    highP <- rep(0, length(y))
    nzInd <- y > 0
    for(i in seq_len(em_J(parList))) {
        lowP[nzInd] <- lowP[nzInd] + tau[i] *
            ppois(y[nzInd] - 1, lamPar[i] * em_depth(parList)[nzInd])
        highP <- highP + tau[i] * ppois(y, lamPar[i] * em_depth(parList))
    }

    # Estimate mixture CDF
    # Lower bound (minimum from rounding precision)
    lBound <- max(
        10 ^ (-prec) / 2,
        qgamma(min(lowP), min(lamPar) / phi, scale = phi)
    )
    # Upper bound (minimum to guarantee 2 steps)
    hBound <- max(
        qgamma(max(highP), max(lamPar) / phi, scale = phi), lBound + 1e-2
    )
    qSeq <- seq(lBound, hBound, min(1e-1, (hBound - lBound) / 100))
    pSeq <- rep(0, length(qSeq))
    for(i in seq_len(em_J(parList))) {
        pSeq <- pSeq + tau[i] * pgamma(qSeq, lamPar[i] / phi, scale = phi)
    }
    cdfFun <- approxfun(pSeq, qSeq, yleft = 0, yright = max(qSeq))
    normVec <- cdfFun(runif(length(y), lowP, highP))
    normVec <- round(normVec, prec)

    return(normVec)
}


# posteriorSample_func is a helper function to perform sampling from the
# estimated posterior distribution given fitted EM parameters
posteriorSample_func <- function(parList, y, depth, fitSlope, prec, emPar) {
    lamPar <- pmax(exp(em_betaVec(parList)), .Machine$double.xmin * 1e1)
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
    parList <- em_ySet(parList, y)
    parList <- em_depthSet(parList, exp(depth) * fitSlope - fitSlope + 1)
    parList <- em_pListSet(parList, lpMat_func(parList))
    pList <- em_pList(parList)
    pList$spMat <- em_pList(parList)$spMat *
        exp(em_pList(parList)$lpMax)
    parList <- em_pListSet(parList, pList)
    lamPar <- pmax(exp(em_betaVec(parList)), .Machine$double.xmin * 1e1)

    # Resample counts
    distVec <- apply(
        em_pList(parList)$spMat[-em_J(parList), , drop = FALSE], 2,
        function(pVec) {
            which.min(runif(1) > c(cumsum(pVec), 1))
        }
    )

    shapeVec <- lamPar / phi
    concVec <- emPar$conPar

    normVec <- rgamma(n = length(em_y(parList)),
                      shape = em_y(parList) * concVec + shapeVec[distVec],
                      scale = 1 / (em_depth(parList) * concVec + 1 / phi))
    normVec <- round(normVec, prec)

    return(normVec)
}
