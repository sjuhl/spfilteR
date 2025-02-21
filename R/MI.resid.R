#' @name MI.resid
#'
#' @title Moran Test for Residual Spatial Autocorrelation
#'
#' @description This function assesses the degree of spatial
#' autocorrelation present in regression residuals by means of the Moran
#' coefficient.
#'
#' @param resid residual vector
#' @param x vector/ matrix of regressors (default = NULL)
#' @param W spatial connectivity matrix
#' @param alternative specification of alternative hypothesis as 'greater' (default),
#' 'lower', or 'two.sided'
#' @param boot optional integer specifying the number of simulation iterations to
#' compute the variance. If NULL (default), variance calculated under assumed normality
#' @param na.rm listwise deletion of observations with missing values (TRUE/ FALSE)
#'
#' @return A \code{data.frame} object with the following elements:
#' \describe{
#' \item{\code{I}}{observed value of the Moran coefficient}
#' \item{\code{EI}}{expected value of Moran's I}
#' \item{\code{VarI}}{variance of Moran's I}
#' \item{\code{zI}}{standardized Moran coefficient}
#' \item{\code{pI}}{\emph{p}-value of the test statistic}
#' }
#'
#' @details The function assumes an intercept-only model if \code{x = NULL}.
#' Furthermore, \code{MI.resid} automatically symmetrizes the matrix
#' \emph{\strong{W}} by: 1/2 * (\emph{\strong{W}} + \emph{\strong{W}}').
#'
#' @note Calculations are based on the procedure proposed by Cliff and Ord
#' (1981). See also Cliff and Ord (1972).
#'
#' @author Sebastian Juhl
#'
#' @references Cliff, Andrew D. and John K. Ord (1981): Spatial Processes:
#' Models & Applications. Pion, London.
#'
#' Cliff, Andrew D. and John K. Ord (1972): Testing for Spatial Autocorrelation
#' Among Regression Residuals. Geographical Analysis, 4 (3): pp. 267 - 284
#'
#' @importFrom stats var
#'
#' @examples
#' data(fakedata)
#' y <- fakedataset$x1
#' x <- fakedataset$x2
#'
#' resid <- y - x %*% solve(crossprod(x)) %*% crossprod(x,y)
#' (Moran <- MI.resid(resid = resid, x = x, W = W, alternative = "greater"))
#'
#' # intercept-only model
#' x <- rep(1, length(y))
#' resid2 <- y - x %*% solve(crossprod(x)) %*% crossprod(x,y)
#' intercept <- MI.resid(resid = resid2, W = W, alternative = "greater")
#' # same result with MI.vec for the intercept-only model
#' vec <- MI.vec(x = resid2, W = W, alternative = "greater")
#' rbind(intercept, vec)
#'
#' @seealso \code{\link{lmFilter}}, \code{\link{glmFilter}}, \code{\link{MI.vec}},
#' \code{\link{MI.local}}
#'
#' @export

MI.resid <- function(resid, x = NULL, W, alternative = "greater", boot = NULL, na.rm = TRUE) {

  #####
  # Input
  # Checks
  #####
  if (!(alternative %in% c("greater", "lower", "two.sided"))) {
    stop("Invalid input: 'alternative' must be either 'greater', 'lower', or 'two.sided'")
  }
  if (!any(class(W) %in% c("matrix", "Matrix", "data.frame"))) {
    stop("W must be of class 'matrix' or 'data.frame'")
  }
  if (any(class(W) != "matrix")) {
    W <- as.matrix(W)
  }
  if (!is.null(x) & anyNA(x) & na.rm == FALSE) {
    stop("Missing values detected in x")
  }

  n <- length(resid)
  if (is.null(x)) {
    x <- data.matrix(rep(1, n))
  } else {
    x <- data.matrix(x)
    # remove missing values
    if (na.rm) {
      miss <- apply(x, 1, anyNA)
      x <- data.matrix(x[!miss,])
      W <- W[!miss, !miss]
    }
  }
  
  # add intercept term if required
  #if (!isTRUE(all.equal(x[, 1], rep(1, n)))) {
  if (!isTRUE(any(apply(x, 2, function(c) all(1 == c))))) {
    x <- cbind(1, x)
  }

  df <- n - qr(x)$rank

  # step 1
  W <- .5 * (W + t(W))
  S0 <- crossprod(rep(1, n), W %*% rep(1, n))
  S1 <- .5 * sum((W + t(W))^2)

  # step 2
  Z <- crossprod(W, x)

  # step 3
  C1 <- crossprod(x, Z)
  C2 <- crossprod(Z)

  # steps 4 & 5
  G <- qr.solve(crossprod(x))
  traceA <- sum(diag(crossprod(G, C1)))
  traceB <- sum(diag(crossprod(4 * G, C2)))
  traceA2 <- sum(diag(crossprod(G %*% C1)))
  I <- n / S0 * crossprod(resid, W %*% resid) / crossprod(resid)
  EI <- -(n * traceA) / (df * S0)
  if (!is.null(boot)) {
    boot <- round(boot)
    if (boot < 100) {
      warning(paste0("Number of bootstrap iterations (",boot,") too small. Set to 100"))
      boot <- 100
    }

    boot.I <- NULL
    for (i in 1:boot) {
      ind <- sample(1:n, replace = TRUE)
      boot.I[i] <- n / S0 * crossprod(resid[ind], W %*% resid[ind]) / crossprod(resid[ind])
    }
    VarI <- var(boot.I)
    zI <- (I - EI) / sqrt(VarI)
    pI <- pfunc(z = I, alternative = alternative, draws = boot.I)
  } else {
    VarI <- (n^2 / (S0^2 * df * (df + 2))) * (S1 + 2 * traceA2 - traceB - ((2 * traceA^2) / df))
    zI <- (I - EI) / sqrt(VarI)
    pI <- pfunc(z = zI, alternative = alternative)
  }

  #####
  # Output
  #####
  out <- data.frame(I, EI, VarI, zI, pI, NA)
  colnames(out) <- c("I", "EI", "VarI", "zI", "pI", "")
  out[1, 6] <- star(p = out[1, "pI"])
  return(out)
}
