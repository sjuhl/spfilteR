#' @name vp
#'
#' @title Variance Partitioning with Moran Spectral Randomization
#'
#' @description This function decomposes the variation in an outcome variable
#' into four fractions: a) the influence of covariates, b) joint influence of
#' covariates and space, c) the influence of space, and d) unexplained residual
#' variation. Moran spectral randomization is applied to obtain the expected
#' value of the coefficient of determination adjusted for spurious correlations.
#'
#' @param y outcome vector
#' @param x vector/ matrix of covariates
#' @param evecs selected eigenvectors
#' @param msr number of permutations to compute the expected value under H0
#'
#' @return Returns an object of class \code{vpart} which provides the following
#' information:
#' \describe{
#' \item{\code{R2}}{unadjusted fractions of explained variation}
#' \item{\code{adjR2}}{adjusted fractions (based on Moran spectral randomization)}
#' \item{\code{msr}}{number of permutations to obtain the expected value under H0}
#' }
#'
#' @note The adjusted R-squared values are obtained by: 1 - (1 - R2) / (1 - E(R2|H0)).
#' For fractions [ab] and [a], Moran spectral randomization is used to derive
#' E(R2|H0). To this end, the rows in matrix (or column vector) x are randomly
#' permuted in order to preserve the correlation structure (see e.g., Clappe et
#' al. 2018).
#'
#' @examples
#' data(fakedata)
#' E <- getEVs(W = W, covars = NULL)$vectors
#'
#' (partition <- vp(y = fakedataset$x1, evecs = E[, 1:10], msr = 100))
#'
#' @author Sebastian Juhl
#'
#' @seealso \code{\link{getEVs}}
#'
#' @references
#' Clappe, Sylvie, Dray Stéphane. and Pedro R. Peres-Neto (2018):
#' Beyond neutrality: disentangling the effects of species sorting
#' and spurious correlations in community analysis. Ecology 99 (8):
#' pp. 1737 - 1747.
#'
#' Wagner, Helene H., and Stéphane Dray (2015): Generating spatially
#' constrained null models for irregularly spaced data using Moran spectral
#' randomization methods. Methods in Ecology and Evolution 6 (10):
#' pp. 1169 - 1178.
#'
#' @export

vp <- function(y, x = NULL, evecs = NULL, msr = 100) {
  n <- length(y)
  if (is.null(x)) {
    x <- rep(1, n)
  }
  x <- as.matrix(x)
  if (!isTRUE(all.equal(x[, 1], rep(1, n)))) {
    x <- cbind(1, x)
  }

  ### checks
  if (anyNA(y) | anyNA(x)) {
    stop("Missing values detected")
  }
  if (ncol(x) > 1 & qr(x)$rank != ncol(x)) {
    stop("Perfect multicollinearity in covariates detected")
  }
  if (msr < 100) {
    warning(paste0("Number of permutations (",msr,") too small. Set to 100"))
    msr <- 100
  }

  if (!is.null(evecs)) {
    evecs <- as.matrix(evecs)
    # drop eigenvector with zero eigenvalue. If W is not of full-rank,
    # multiple eigenvectors with zero eigenvalues exist
    # (Wagner/ Dray 2015: 1170)
    null <- which(!vapply(apply(evecs, 2, sum),
                          function(x) isTRUE(all.equal(x, 0, tolerance = 1e-7)),
                          FUN.VALUE = TRUE))
    if (length(null) == 1) {
      evecs <- evecs[, -null]
    } else if (length(null) > 1) {
      sub <- cbind(1, evecs[, null])
      sub <- qr.Q(qr(sub))
      evecs[,null] <- sub[,-ncol(sub)]
      evecs <- evecs[, -null[1]]
    }
  }

  # test nr of supplied covars
  if (ncol(cbind(x, evecs)) >= n) {
    stop("Nr of covariates equals or exceeds n")
  }

  ### unadjusted R-squared
  TSS <- sum((y - mean(y))^2)
  # ab
  resid.x <- y - x %*% solve(crossprod(x)) %*% crossprod(x, y)
  R2.ab <- 1 - (sum(resid.x^2) / TSS)
  # abc
  xe <- cbind(x, evecs)
  resid.xe <- y - xe %*% solve(crossprod(xe)) %*% crossprod(xe, y)
  R2.abc <- 1 - (sum(resid.xe^2) / TSS)
  if (!is.null(evecs)) {
    # bc
    resid.e <- y - cbind(1, evecs) %*% solve(crossprod(cbind(1, evecs))) %*% crossprod(cbind(1, evecs), y)
    R2.bc <- 1 - (sum(resid.e^2) / TSS)
  } else {
    R2.bc <- 0
  }
  # individual fractions
  R2.a <- R2.abc - R2.bc
  R2.c <- R2.abc - R2.ab
  R2.b <- R2.abc - R2.a - R2.c
  R2.d <- 1 - R2.abc

  # E(H0) - Moran Spectral Randomization
  msr <- round(msr)
  r2.ab <- r2.a <- NULL
  for (i in seq_len(msr)) {
    ind <- sample(1:n, replace = TRUE)
    x.msr <- x[ind,]
    resid.msr <- y - x.msr %*% solve(crossprod(x.msr)) %*% crossprod(x.msr, y)
    r2.ab[i] <- 1 - (sum(resid.msr^2) / TSS)
    xe.msr <- cbind(x.msr, evecs)
    resid.xe.msr <- y - xe.msr %*% solve(crossprod(xe.msr)) %*% crossprod(xe.msr, y)
    r2.abc <- 1 - (sum(resid.xe.msr^2) / TSS)
    r2.a[i] <- r2.abc - R2.bc
  }
  Eh0.ab <- mean(r2.ab)
  Eh0.a <- mean(r2.a)

  ### adjusted R-squared
  adjR2.ab <- 1 - (1 / (1 - Eh0.ab)) * (1 - R2.ab)
  adjR2.bc <- ifelse(is.null(evecs), 0, 1 - (1 / (1 - ncol(evecs) / (n - 1))) * (1 - R2.bc))
  adjR2.a <- 1 - (1 / (1 - Eh0.a)) * (1 - R2.a)
  adjR2.b <- adjR2.ab - adjR2.a
  adjR2.c <- adjR2.bc - adjR2.b
  adjR2.abc <- adjR2.ab + adjR2.c
  adjR2.d <- 1 - adjR2.abc

  ### Output
  R2 <- c(ab = R2.ab,
          bc = R2.bc,
          abc = R2.abc,
          a = R2.a,
          b = R2.b,
          c = R2.c,
          d = R2.d)
  adjR2 <- c(ab = adjR2.ab,
             bc = adjR2.bc,
             abc = adjR2.abc,
             a = adjR2.a,
             b = adjR2.b,
             c = adjR2.c,
             d = adjR2.d)
  out <- list(R2 = R2, adjR2 = adjR2, msr =msr)
  class(out) <- "vpart"
  return(out)
}

#' @export
print.vpart <- function(x, ...) {
  res <- data.frame(cbind(format(round(x$adjR2, 7), nsmall = 7),
                          format(round(x$R2, 7), nsmall = 7)),
                          row.names = c("ab", "bc", "abc", "a", "b", "c", "d"))
  colnames(res) <- c("Adj. R2", "R2")
  cat("\n\t - Variation Partitioning -\n\n")
  print(res)
  cat(paste("---\n", "Permutations:", x$msr, "\n"))
}
