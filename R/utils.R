#' @name pfunc
#' @importFrom stats pnorm
#' @noRd

pfunc <- function(z, alternative, draws = NULL) {
  z <- as.numeric(z)
  # analytical variance estimate
  if (is.null(draws)) {
    if (alternative == "greater") {
      p <- pnorm(z, lower.tail = FALSE)
    } else if (alternative == "lower") {
      p <- pnorm(z, lower.tail = TRUE)
    } else {
      p <- 2 * pnorm(abs(z), lower.tail = FALSE)
    }
  } else {
    # simulation-based variance estimate
    # see e.g., North/ Curtis/ Sham (2002) [Am J Hum Genet]
    # or Dary (2011) [Geogr. An.] for the '+1'
    if (alternative == "greater") {
      p <- (sum(draws >= z) + 1) / (length(draws) + 1)
    } else if (alternative == "lower") {
      p <- (sum(draws <= z) + 1) / (length(draws) + 1)
    } else {
      p <- (sum(abs(draws) >= abs(z)) + 1) / (length(draws) + 1)
      # see e.g., Hartwig (2013) [J Clin Trials]
    }
  }
  return(p)
}


#' @name star
#' @noRd

star <- function(p) {
  out <- NULL
  out[p <= .001] <- "***"
  out[p <= .01 & p > .001] <- "**"
  out[p <= .05 & p > .01] <- "*"
  out[p <= .1 & p > .05] <- "."
  out[p > .1] <- " "
  return(out)
}


#' @name candsetsize
#' @noRd

candsetsize <- function(npos, zMI) {
  denominator <- 1 + exp(2.1480 - (6.1808 * (zMI+.6)^.1742) / npos^.1298 + 3.3534 / (zMI + .6)^.1742)
  nc <- npos / denominator
  return(round(nc, 0))
}


#' @name residfun
#' @noRd

residfun <- function(y, fitvals, size = NULL, model) {
  if (!(model %in% c("linear", "probit", "logit", "poisson", "nb"))) {
    stop("'model' must be either 'linear', 'probit', 'logit', 'poisson', or 'nb'")
  }
  # raw residuals
  raw <- y - fitvals
  # pearson & deviance residuals
  if (model %in% c("probit", "logit")) {
    pearson <- (y - fitvals) / sqrt(fitvals * (1 - fitvals))
    sign <- ifelse(y == 1, 1, -1)
    deviance <- sign * sqrt(-2 * (y * log(fitvals) + (1 - y) * log(1 - fitvals)))
  } else if (model == "poisson") {
    pearson <- (y - fitvals) / sqrt(fitvals)
    sign <- ifelse(y > fitvals, 1, -1)
    ratio <- ifelse(y == 0, 1, y / fitvals)
    deviance <- sign * sqrt(2 * (y * log(ratio) - (y - fitvals)))
  } else if (model == "nb") {
    pearson <- (y - fitvals) / sqrt((fitvals + (fitvals^2 / size)))
    deviance <- Inf
  } else if (model == "linear") {
    pearson <- deviance <- raw
  }
  # output
  out <- data.frame(raw, pearson, deviance)
  return(out)
}


#' @name fittedval
#' #' @importFrom stats pnorm
#' @noRd

fittedval <- function(x, params, model) {
  if(model != "nb") mu <- x %*% params
  if (model == "linear") {
    yhat <- mu
  } else if (model == "probit") {
    yhat <- pnorm(mu)
  } else if (model == "logit") {
    yhat <- exp(mu) / (1 + exp(mu))
  } else if (model == "poisson") {
    yhat <- exp(mu)
  } else if (model == "nb"){
    mu <- x %*% params[-length(params)]
    yhat <- exp(mu)
  }
  return(yhat)
}


#' @name conditionNumber
#' @importFrom stats cor
#' @noRd

# function to calculate the condition number - degree of EV multicollinearity
# (see Griffith 2004b, p. 1797 and Griffith/ Amrhein 1997, p. 98)
conditionNumber <- function(evecs = NULL, round = 8) {
  if (!is.null(evecs)) {
    cormat <- cor(evecs)
    corevals <- eigen(cormat)$values
    res <- round(sqrt(corevals[1] / corevals[length(corevals)]), round)
  } else {
    res <- NULL
  }
  return(res)
}


#' @name getICs
#' @noRd

getICs <- function(negloglik, n, df) {
  AIC <- 2 * negloglik + 2 * df
  BIC <- 2 * negloglik + log(n) * df
  out <- data.frame(AIC, BIC)
  return(out)
}


#' @name pseudoR2
#' @noRd

# McFadden's (adjusted) pseudo-R2 (filtered vs unfiltered model)
pseudoR2 <- function(negloglik_n, negloglik_f, nev) {
  R2 <- 1 - (-negloglik_n / -negloglik_f)
  adjR2 <- 1 - ((-negloglik_n - nev) / -negloglik_f)
  out <- data.frame(R2, adjR2)
  return(out)
}


#' @name Zscore
#' @importFrom stats sd
#' @noRd

Zscore <- function(x, na.rm = TRUE) {
  x <- as.matrix(x)
  Z <- apply(x, 2, function(v) (v - mean(v, na.rm = na.rm)) / sd(v, na.rm = na.rm))
  return(Z)
}
