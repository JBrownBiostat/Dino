# S4 Class for organization of internal EM-Algorithm parameters and intermediate
# results
setClass(
    "emFitPars",
    representation(
        y = "numeric", # observed counts

        depth = "numeric", # sequencing depth / library size

        J =  "numeric", # number of mixing components

        gamList = "list", # list of current mixture weights
        # gamList$gamVec: uncentered log-transformed weights
        # gamList$gamDot: the sum for transformation to mixture weights

        zInd = "numeric", # index of largest gamList$gamVec for
        # stability centering

        betaVec = "numeric", # current estimate of (log) Poisson means

        pList = "list", # list of coefficients defining current conditional
        # mixture component membership probabilities
        # pList$spMat: row-scaled (for stability) matrix (J x length(y))
        # of conditional probabilities
        # pList$lpMax: log-scale row maxs for transformation of pList$spMat
        # back to probability space

        gTilde = "list", # list of differences between next "strict" EM iterate
        # and current parameter estimates for "accelerated" EM modification
        # gTilde$gamTilde: vector of new mixture wights minus gamList$gamVec
        # gTilde$betaTilde: vector of new Poisson means minus betaVec
        # gTilde$gTilde: c(gTilde$gamTilde, gTilde$betaTilde)

        g = "list", # list of parameter gradients for "accelerated" EM
        # g$gGam: gradient of log-likelihood in mixture weights
        # g$gBeta: gradient of log-likelihood in  Poisson means
        # g$g: c(g$gGam, g$gBeta)

        S = "matrix", # approximate inverse Hessian (2*J x 2*J) matrix for
        # "accelerated" EM updates

        lLik = "numeric", # current log-likelihood

        alpha = "numeric", # optimization step-length

        d = "numeric", # optimization step-direction

        dgTilde = "list", # list of differences between k+1 parameters and
        # k+1 "strict" EM update
        # dgTilde$gamTilde: gTilde$gamTilde at iterate k+1
        # dgTilde$betaTilde: gTilde$betaTilde at iterate k+1
        # dgTilde$gTilde: c(dgTilde$gamTilde, dgTilde$betaTilde)

        parUpdate = "list", # list of optimization parameters at next iterate
        # parUpdate$gamList: gamList at iterate k+1
        # parUpdate$betaVec: betaVec at iterate k+1
        # parUpdate$parVec: c(parUpdate$gamList$gamVec, parUpdate$betaVec)

        gUpdate = "list", # list of log-likelihood gradients at next iterate
        # gUpdate$gGam: g$gGam at iterate k + 1
        # gUpdate$gBeta: g$gBeta at iterate k + 1
        # gUpdate$g: c(gUpdate$gGam, gUpdate$gBeta)

        pListUpdate = "list", # list of conditional mixture component membership
        # probabilities at next iterate
        # pListUpdate$spMat: pList$spMat at iterate k + 1
        # pListUpdate$lpMax: pList$lpMax at iterate k + 1

        lLikUpdate = "numeric" # log-likelihood at next iterate
    )
)


# emFitPars methods
setGeneric("em_y", function(object) {
    standardGeneric("em_y")
})
setMethod(f = "em_y", signature(object = "emFitPars"), function(object) {
    return(object@y)
})
setGeneric("em_ySet", function(object, value) {
    standardGeneric("em_ySet")
})
setMethod(
    f = "em_ySet", signature(object = "emFitPars"),
    function(object, value) {
        object@y <- value
        return(object)
    }
)

setGeneric("em_depth", function(object) {
    standardGeneric("em_depth")
})
setMethod(f = "em_depth", signature(object = "emFitPars"), function(object) {
    object@depth
})
setGeneric("em_depthSet", function(object, value) {
    standardGeneric("em_depthSet")
})
setMethod(
    f = "em_depthSet", signature(object = "emFitPars"),
    function(object, value) {
        object@depth <- value
        return(object)
    }
)

setGeneric("em_J", function(object) {
    standardGeneric("em_J")
})
setMethod(f = "em_J", signature(object = "emFitPars"), function(object) {
    object@J
})
setGeneric("em_JSet", function(object, value) {
    standardGeneric("em_JSet")
})
setMethod(
    f = "em_JSet", signature(object = "emFitPars"),
    function(object, value) {
        object@J <- value
        return(object)
    }
)

setGeneric("em_gamList", function(object) {
    standardGeneric("em_gamList")
})
setMethod(f = "em_gamList", signature(object = "emFitPars"), function(object) {
    object@gamList
})
setGeneric("em_gamListSet", function(object, value) {
    standardGeneric("em_gamListSet")
})
setMethod(
    f = "em_gamListSet", signature(object = "emFitPars"),
    function(object, value) {
        object@gamList <- value
        return(object)
    }
)

setGeneric("em_zInd", function(object) {
    standardGeneric("em_zInd")
})
setMethod(f = "em_zInd", signature(object = "emFitPars"), function(object) {
    object@zInd
})
setGeneric("em_zIndSet", function(object, value) {
    standardGeneric("em_zIndSet")
})
setMethod(
    f = "em_zIndSet", signature(object = "emFitPars"), function(object, value) {
        object@zInd <- value
        return(object)
    }
)

setGeneric("em_betaVec", function(object) {
    standardGeneric("em_betaVec")
})
setMethod(f = "em_betaVec", signature(object = "emFitPars"), function(object) {
    object@betaVec
})
setGeneric("em_betaVecSet", function(object, value) {
    standardGeneric("em_betaVecSet")
})
setMethod(
    f = "em_betaVecSet", signature(object = "emFitPars"),
    function(object, value) {
        object@betaVec <- value
        return(object)
    }
)

setGeneric("em_pList", function(object) {
    standardGeneric("em_pList")
})
setMethod(f = "em_pList", signature(object = "emFitPars"), function(object) {
    object@pList
})
setGeneric("em_pListSet", function(object, value) {
    standardGeneric("em_pListSet")
})
setMethod(
    f = "em_pListSet", signature(object = "emFitPars"),
    function(object, value) {
        object@pList <- value
        return(object)
    }
)

setGeneric("em_gTilde", function(object) {
    standardGeneric("em_gTilde")
})
setMethod(f = "em_gTilde", signature(object = "emFitPars"), function(object) {
    object@gTilde
})
setGeneric("em_gTildeSet", function(object, value) {
    standardGeneric("em_gTildeSet")
})
setMethod(
    f = "em_gTildeSet", signature(object = "emFitPars"),
    function(object, value) {
        object@gTilde <- value
        return(object)
    }
)

setGeneric("em_g", function(object) {
    standardGeneric("em_g")
})
setMethod(f = "em_g", signature(object = "emFitPars"), function(object) {
    object@g
})
setGeneric("em_gSet", function(object, value) {
    standardGeneric("em_gSet")
})
setMethod(
    f = "em_gSet", signature(object = "emFitPars"),
    function(object, value) {
        object@g <- value
        return(object)
    }
)

setGeneric("em_S", function(object) {
    standardGeneric("em_S")
})
setMethod(f = "em_S", signature(object = "emFitPars"), function(object) {
    object@S
})
setGeneric("em_SSet", function(object, value) {
    standardGeneric("em_SSet")
})
setMethod(
    f = "em_SSet", signature(object = "emFitPars"),
    function(object, value) {
        object@S <- value
        return(object)
    }
)

setGeneric("em_lLik", function(object) {
    standardGeneric("em_lLik")
})
setMethod(f = "em_lLik", signature(object = "emFitPars"), function(object) {
    object@lLik
})
setGeneric("em_lLikSet", function(object, value) {
    standardGeneric("em_lLikSet")
})
setMethod(
    f = "em_lLikSet", signature(object = "emFitPars"),
    function(object, value) {
        object@lLik <- value
        return(object)
    }
)

setGeneric("em_alpha", function(object) {
    standardGeneric("em_alpha")
})
setMethod(f = "em_alpha", signature(object = "emFitPars"), function(object) {
    object@alpha
})
setGeneric("em_alphaSet", function(object, value) {
    standardGeneric("em_alphaSet")
})
setMethod(
    f = "em_alphaSet", signature(object = "emFitPars"),
    function(object, value) {
        object@alpha <- value
        return(object)
    }
)

setGeneric("em_d", function(object) {
    standardGeneric("em_d")
})
setMethod(f = "em_d", signature(object = "emFitPars"), function(object) {
    object@d
})
setGeneric("em_dSet", function(object, value) {
    standardGeneric("em_dSet")
})
setMethod(
    f = "em_dSet", signature(object = "emFitPars"),
    function(object, value) {
        object@d <- value
        return(object)
    }
)

setGeneric("em_dgTilde", function(object) {
    standardGeneric("em_dgTilde")
})
setMethod(f = "em_dgTilde", signature(object = "emFitPars"), function(object) {
    object@dgTilde
})
setGeneric("em_dgTildeSet", function(object, value) {
    standardGeneric("em_dgTildeSet")
})
setMethod(
    f = "em_dgTildeSet", signature(object = "emFitPars"),
    function(object, value) {
        object@dgTilde <- value
        return(object)
    }
)

setGeneric("em_parUpdate", function(object) {
    standardGeneric("em_parUpdate")
})
setMethod(
    f = "em_parUpdate", signature(object = "emFitPars"), function(object) {
        object@parUpdate
    }
)
setGeneric("em_parUpdateSet", function(object, value) {
    standardGeneric("em_parUpdateSet")
})
setMethod(
    f = "em_parUpdateSet", signature(object = "emFitPars"),
    function(object, value) {
        object@parUpdate <- value
        return(object)
    }
)

setGeneric("em_gUpdate", function(object) {
    standardGeneric("em_gUpdate")
})
setMethod(f = "em_gUpdate", signature(object = "emFitPars"), function(object) {
    object@gUpdate
})
setGeneric("em_gUpdateSet", function(object, value) {
    standardGeneric("em_gUpdateSet")
})
setMethod(
    f = "em_gUpdateSet", signature(object = "emFitPars"),
    function(object, value) {
        object@gUpdate <- value
        return(object)
    }
)

setGeneric("em_pListUpdate", function(object) {
    standardGeneric("em_pListUpdate")
})
setMethod(
    f = "em_pListUpdate", signature(object = "emFitPars"), function(object) {
        object@pListUpdate
    }
)
setGeneric("em_pListUpdateSet", function(object, value) {
    standardGeneric("em_pListUpdateSet")
})
setMethod(
    f = "em_pListUpdateSet", signature(object = "emFitPars"),
    function(object, value) {
        object@pListUpdate <- value
        return(object)
    }
)

setGeneric("em_lLikUpdate", function(object) {
    standardGeneric("em_lLikUpdate")
})
setMethod(
    f = "em_lLikUpdate", signature(object = "emFitPars"), function(object) {
        object@lLikUpdate
    }
)
setGeneric("em_lLikUpdateSet", function(object, value) {
    standardGeneric("em_lLikUpdateSet")
})
setMethod(
    f = "em_lLikUpdateSet", signature(object = "emFitPars"),
    function(object, value) {
        object@lLikUpdate <- value
        return(object)
    }
)
