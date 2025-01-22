#' @name MI.local
#'
#' @title Local Moran Coefficient
#'
#' @description Reports the local Moran Coefficient for each unit.
#'
#' @param x a vector
#' @param W spatial connectivity matrix
#' @param alternative specification of alternative hypothesis as 'greater' (default),
#' 'lower', or 'two.sided'
#' @param na.rm listwise deletion of observations with missing values (TRUE/ FALSE)
#'
#' @return Returns an object of class \code{data.frame} that contains the
#' following information for each variable:
#' \describe{
#' \item{\code{Ii}}{observed value of local Moran's I}
#' \item{\code{EIi}}{expected value of local Moran coefficients}
#' \item{\code{VarIi}}{variance of local Moran's I}
#' \item{\code{zIi}}{standardized local Moran coefficient}
#' \item{\code{pIi}}{\emph{p}-value of the test statistic}
#' }
#'
#' @note The calculation of the statistic and its moments follows
#' Anselin (1995) and Sokal et al. (1998).
#'
#' @author Sebastian Juhl
#'
#' @references Anselin, Luc (1991): Local Indicators of Spatial
#' Association-LISA. Geographical Analysis, 27 (2): pp. 93 - 115.
#'
#' Bivand, Roger S. and David W. S. Wong (2018): Comparing Implementations
#' of Global and Local Indicators of Spatial Association. TEST, 27:
#' pp. 716 - 748.
#'
#' Sokal, Robert R., Neal L. Oden, Barbara A. Thomson (1998): Local
#' Spatial Autocorrelation in a Biological Model. Geographical Analysis,
#' 30 (4): pp. 331 - 354.
#'
#' @seealso \code{\link{MI.vec}}, \code{\link{MI.ev}}, \code{\link{MI.sf}},
#' \code{\link{MI.resid}}, \code{\link{MI.decomp}}
#'
#' @examples
#' data(fakedata)
#' x <- fakedataset$x2
#'
#' (MIi <- MI.local(x = x, W = W, alternative = "greater"))
#'
#' @export

MI.local <- function(x, W, alternative = "greater", na.rm = TRUE) {

  # missing values
  miss <- is.na(x)
  if (na.rm) {
    x <- x[!miss]
    W <- W[!miss, !miss]
  }

  #####
  # Input
  # Checks
  #####
  if (!any(class(W) %in% c("matrix", "Matrix", "data.frame"))) {
    stop("W must be of class 'matrix' or 'data.frame'")
  }
  if (any(class(W) != "matrix")) {
    W <- as.matrix(W)
  }
  if (anyNA(x) | anyNA(W)) {
    stop("Missing values detected")
  }
  if (!(alternative %in% c("greater", "lower", "two.sided"))) {
    stop("Invalid input: 'alternative' must be either 'greater',
         'lower', or 'two.sided'")
  }

  # define variables
  n <- length(x)
  z <- x - mean(x)
  m2 <- sum(z^2) / n
  Wi <- apply(W, 1, sum)
  Wi2 <- apply(W, 1, function(x) sum(x^2))
  b2 <- n * sum(z^4) / sum(z^2)^2

  #####
  # Output
  #####
  out <- data.frame(matrix(NA, nrow = length(miss), ncol = 6))
  colnames(out) <- c("Ii", "EIi", "VarIi", "zIi", "pIi", "")
  # observed local Is
  out[!miss, "Ii"] <- (z / m2) * W %*% z
  # expected
  out[!miss, "EIi"] <- -Wi / (n - 1)
  # variance
  out[!miss, "VarIi"] <- Wi2 * (n - b2) / (n - 1) + (Wi^2 - Wi2) * (2 * b2 - n) / ((n - 1) * (n - 2)) - out[!miss, "EIi"]^2
  # test statistic
  out[!miss, "zIi"] <- apply(out[!miss,], 1, function(x) (x[1] - x[2]) / sqrt(x[3]))
  # pI
  out[!miss, "pIi"] <- vapply(out[!miss, "zIi"], pfunc ,alternative = alternative
                        ,FUN.VALUE = numeric(1))
  out[!miss, 6] <- vapply(out[!miss, "pIi"], star, FUN.VALUE = character(1))

  # return
  return(out)
}
