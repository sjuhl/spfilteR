#' @name MI.decomp
#'
#' @title Decomposition of the Moran Coefficient
#'
#' @description A decomposition of the Moran coefficient in order to separately
#' test for the simultaneous presence of positive and negative autocorrelation
#' in a variable.
#'
#' @param x a vector or matrix
#' @param W spatial connectivity matrix
#' @param nsim number of iterations to simulate the null distribution
#' @param na.rm listwise deletion of observations with missing values (TRUE/ FALSE)
#'
#' @return Returns a \code{data.frame} that contains the following information
#' for each variable:
#' \describe{
#' \item{\code{I+}}{observed value of Moran's I (positive part)}
#' \item{\code{VarI+}}{variance of Moran's I (positive part)}
#' \item{\code{pI+}}{simulated \emph{p}-value of Moran's I (positive part)}
#' \item{\code{I-}}{observed value of Moran's I (negative part)}
#' \item{\code{VarI-}}{variance of Moran's I (negative part)}
#' \item{\code{pI-}}{simulated \emph{p}-value of Moran's I (negative part)}
#' \item{\code{pItwo.sided}}{simulated \emph{p}-value of the two-sided test}
#' }
#'
#' @details If \code{x} is a matrix, this function computes the Moran
#' test for spatial autocorrelation for each column.
#'
#' The \emph{p}-values calculated for \code{I+} and \code{I-} assume
#' a directed alternative hypothesis. Statistical significance is assessed
#' using a permutation procedure to generate a simulated null distribution.
#'
#' @author Sebastian Juhl
#'
#' @references Dary, Stéphane (2011): A New Perspective about Moran’s
#' Coefficient: Spatial Autocorrelation as a Linear Regression Problem.
#' Geographical Analysis, 43 (2): pp. 127 - 141.
#'
#' @seealso \code{\link{MI.vec}}, \code{\link{MI.ev}}, \code{\link{MI.sf}},
#' \code{\link{MI.resid}}, \code{\link{MI.local}}, \code{\link{getEVs}}
#'
#' @importFrom stats var cor sd
#'
#' @examples
#' data(fakedata)
#' X <- cbind(fakedataset$x1, fakedataset$x2,
#' fakedataset$x3, fakedataset$negative)
#'
#' (MI.dec <- MI.decomp(x = X, W = W, nsim = 100))
#'
#' # the sum of I+ and I- equals the observed Moran coefficient:
#' I <- MI.vec(x = X, W = W)[, "I"]
#' cbind(MI.dec[, "I+"] + MI.dec[, "I-"], I)
#'
#' @export

MI.decomp <- function(x, W, nsim = 100, na.rm = TRUE) {

  # convert x to a matrix and save names (if provided)
  x <- data.matrix(x)
  if (!is.null(colnames(x))) {
    nams <- colnames(x)
  }
  x <- unname(x)

  # missing values
  miss <- is.na(x)

  #####
  # Input
  # Checks
  #####
  if (0 %in% apply(x, 2, sd, na.rm = TRUE)) {
    warning("Constant term removed from x")
    x <- data.matrix(x[, apply(x, 2, sd, na.rm = TRUE) != 0])
  }
  if (!any(class(W) %in% c("matrix", "Matrix", "data.frame"))) {
    stop("W must be of class 'matrix' or 'data.frame'")
  }
  if (any(class(W) != "matrix")) {
    W <- as.matrix(W)
  }
  if (anyNA(W)) {
    stop("Missing values in W detected")
  }
  if (!na.rm & anyNA(x)) {
    stop("Missing values in x detected")
  }
  if (nsim < 100) {
    warning(paste0("Number of permutations (",nsim,") too small. Set to 100"))
    nsim <- 100
  }

  # nr of input variables
  nx <- ncol(x)

  # symmetrize W
  W <- .5 * (W + t(W))

  # z-score transformation
  Z <- Zscore(x, na.rm = na.rm)

  # eigendecomposition
  eigen <- getEVs(W = W)

  #####
  # Output
  #####
  out <- data.frame(matrix(NA, nrow = nx, ncol = 10))
  colnames(out) <- c("I+", "VarI+", "pI+", "",
                     "I-", "VarI-", "pI-", "",
                     "pItwo.sided","")

  for (i in seq_len(nx)) {
    # expected value
    n <- sum(!miss[, i])
    EI <- -1 / (n - 1)

    # null distribution
    Ip <- In <- NULL
    for (k in seq_len(nsim)) {
      ind <- sample(which(Z[, i] != miss[, i]), replace = TRUE)
      c2 <- cor(Z[ind, i], eigen$vectors[!miss[, i], !miss[, i]], method = "pearson")^2
      Ip[k] <- c2[, eigen$values[!miss[, i]] > EI] %*% eigen$moran[eigen$values > EI & !miss[, i]]
      In[k] <- c2[, eigen$values[!miss[, i]] < EI] %*% eigen$moran[eigen$values < EI & !miss[, i]]
    }

    # observed
    cor2 <- cor(Z[!miss[, i], i], eigen$vectors[!miss[, i], !miss[, i]], method = "pearson")^2
    out[i, "I+"] <- cor2[, eigen$values[!miss[, i]] > EI] %*% eigen$moran[eigen$values > EI & !miss[, i]]
    out[i, "I-"] <- cor2[, eigen$values[!miss[, i]] < EI] %*% eigen$moran[eigen$values < EI & !miss[, i]]

    # variance
    out[i, "VarI+"] <- var(Ip)
    out[i, "VarI-"] <- var(In)

    # significance
    out[i, "pI+"] <- pfunc(z = out[i, "I+"], alternative = "greater", draws = Ip)
    out[i, "pI-"] <- pfunc(z = out[i, "I-"], alternative = "lower", draws = In)
    out[i, "pItwo.sided"] <- 2 * min(out[i, c("pI+", "pI-")])
    out[i, 4] <- star(p = out[i, "pI+"])
    out[i, 8] <- star(p = out[i, "pI-"])
    out[i, 10] <- star(p = out[i, "pItwo.sided"])
  }

  # add variable names
  if (exists('nams', inherits = FALSE)) {
    rownames(out) <- nams
  }

  return(out)
}
