#' @name MI.vec
#'
#' @title Local Moran Coefficient
#'
#' @description Calculates the local Moran Coefficient for each unit.
#'
#' @param x a vector
#' @param W spatial connectivity matrix
#' @param alternative specification of alternative hypothesis as 'greater' (default),
#' 'lower', or 'two.sided'
#'
#' @return Returns an object of class \code{data.frame} that contains the
#' following information for each variable:
#' \describe{
#' \item{\code{Ii}}{observed value of local Moran coefficients}
#' \item{\code{EIi}}{expected value of local Moran coefficients}
#' \item{\code{VarIi}}{variance of local Moran coefficients}
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
#' @importFrom stats
#'
#' @examples
#' data(fakedata)
#' x <- fakedataset$x2
#'
#' (MI <- MI.local(x=X,W=W,alternative="greater"))
#'
#' @export

MI.local <- function(x,W,alternative="greater"){

  #####
  # Input
  # Checks
  #####
  if(!any(class(W) %in% c("matrix","Matrix","data.frame"))){
    stop("W must be of class 'matrix' or 'data.frame'")
  }
  if(any(class(W)!="matrix")) W <- as.matrix(W)
  if(anyNA(x) | anyNA(W)) stop("Missing values detected")
  if(!(alternative %in% c("greater","lower", "two.sided"))){
    stop("Invalid input: 'alternative' must be either 'greater',
         'lower', or 'two.sided'")
  }

  # save names (if provided)
  if(!is.null(names(x))) nams <- names(x)

  n <- length(x)
  z <- x - mean(x)
  m2 <- sum(z^2)/n
  Wi <- apply(W,1,sum)

  #####
  # Output
  #####
  out <- data.frame(matrix(NA,nrow=n,ncol=6))
  colnames(out) <- c("Ii","EIi","VarIi","zIi","pIi","")
  # observed local Is
  out[,"Ii"] <- (z/m2) * W%*%z
  # expected
  out[,"EIi"] <- -Wi/(n-1)
  # variance
  out[,"VarIi"] <- NA
  # test statistic
  out[i,"zI"] <- NA #(out[i,"I"]-out[i,"EI"])/sqrt(out[i,"VarI"])
  # pI
  out[i,"pI"] <- NA #pfunc(z=out[i,"zI"],alternative=alternative)
  out[i,6] <- NA #star(p=out[i,"pI"])

  if(!is.null(colnames(x))) rownames(out) <- nams

  return(out)
}
