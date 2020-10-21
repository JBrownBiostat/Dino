#' Calculate upper and lower deviations at intercept points
#'
#' @description \code{calcDev} is a helper function to calculate upper and
#'   lower sums of absolute deviations for quantile regression.
#'
#' @usage \code{calcDev(y, depth, depthRep, slope = 1)}
calcDev <- function(y, depth, depthRep, slope = 1) {
  # Initialize return
  logY <- pmax(log(y), log(0.999))

  ret <- data.frame(
    intercept = sort(depth - (logY - log(0.999)) / slope)
  )
  intDup <- duplicated(ret$intercept)
  freq <- c(which(!intDup), length(depth) + 1)
  ret <- ret[!intDup, , drop = FALSE]
  ret$freq <- freq[-1] - freq[-length(freq)]
  rm(list = c("intDup", "freq"))

  # Calculate lower deviations
  ret$lowDev <- sum((log(0.999) + slope * (depth - ret$intercept[1])) - logY) -
    cumsum(c(0, diff(ret$intercept) * slope * rev(cumsum(rev(ret$freq[-1])))))

  # Calculate upper deviations
  intDepthOrd <- order(c(ret$intercept, depthRep$depth))
  ret$highDev <- c(
    0,
    cumsum(
      (cumsum(c(ret$freq, -depthRep$rep)[intDepthOrd])
       [-(nrow(ret) + nrow(depthRep))] *
         diff(c(ret$intercept, depthRep$depth)[intDepthOrd])) *
        slope)[intDepthOrd <= nrow(ret)][-nrow(ret)]
  )

  # Calculate numbers of zeros (censored points)
  zInd <- match(depth[y == 0], ret$intercept)
  ret$zInt <- 0
  ret$zInt[sort(unique(zInd))] <- table(zInd)
  remDepth <- depthRep$rep

  # Remove duplicate quantiles
  if(ret$freq[1] == ret$zInt[1]) {
    ret <- ret[which.max(ret$freq > ret$zInt):nrow(ret), ]
  }

  # Convexify percentile estimates #
  nY <- sum(ret$freq)
  chullInt <- chull(ret$highDev, ret$lowDev)
  chullInt <- chullInt[which.max(chullInt == nrow(ret)):which.max(chullInt == 1)]
  hRet <- ret[rev(chullInt), ]
  hRet$m <- c(1e3 * diff(hRet$lowDev[1:2]) / diff(hRet$highDev[1:2]), diff(hRet$lowDev) / diff(hRet$highDev))
  mFunc <- approxfun(x = hRet$highDev, y = hRet$m)
  m <- mFunc(ret$highDev)
  ret$pct <- -m / (1 - m)
  pctDiff <- pmax(-diff(ret$pct), (ret$freq - ret$zInt)[-nrow(ret)] / nY)
  lZ <- ret$zInt[1]
  for(i in 2:length(pctDiff)) {
    lZ <- lZ + ret$zInt[i] - (pctDiff[i-1] - (ret$freq[i - 1] - ret$zInt[i - 1]) / nY) * nY
    pctDiff[i] <- min(pctDiff[i], (ret$freq[i] - ret$zInt[i] + lZ) / nY)
  }
  ret$pct <- c(1, 1 - cumsum(pctDiff))
  ret$quant <- -ret$intercept * slope + log(0.999)
  ret <- ret[!duplicated(ret$pct), ]
  ret <- ret[!duplicated(ret$quant), ]

  if(nrow(ret) > 100) {
    subLen <- floor(log(nrow(ret)))
    subInd <- (1:(ceiling(nrow(ret) / subLen))) * subLen - (subLen - 1)
    ret <- ret[subInd[1:(length(subInd) - subLen)], ]
  }
  ret <- ret[, 6:7]

  quantFun <- approxfun(ret$pct, ret$quant)
  mu <- sum(
    exp(quantFun(seq(min(ret$pct), 1, length = 1e3)))  *
      ((1 - min(ret$pct)) / 1e3)
  )
  minQ <- exp(min(ret$quant))
  minP <- min(ret$pct)
  out <- optimize(f = function(beta, mu, minQ, minP) {
    (pgamma(minQ, mu / beta, beta) - minP)^2
  }, interval = c(0, 1e2), mu = mu, minQ = minQ, minP = minP,
  tol = 1e-6)
  beta <- out$minimum

  minP <- min(minP, pgamma(minQ, mu / beta, beta))

  ret <- rbind(
    ret,
    data.frame(
      pct = seq(minP, 0, length = ceiling(100 * minP)),
      quant = log(qgamma(seq(minP, 0, length = ceiling(100 * minP)),
                         mu / beta, beta)) -
        (min(ret$quant) - log(qgamma(minP, mu / beta, beta)))
    )
  )
  ret <- ret[!duplicated(ret$pct), ]
  ret <- ret[!is.infinite(ret$quant), ]
  ret <- ret[!duplicated(ret$quant), ]
  if(min(ret$pct) != 0) {
    ret <- rbind(
      ret, c(0, min(ret$quant) - 100)
    )
  }

  return(ret)
}
