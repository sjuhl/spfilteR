#' @title Moran Test for Spatial Autocorrelation in a Vector
#'
#' @description Tests for the presence of spatial autocorrelation
#' in a vector as indicated by the Moran coefficient. The variance
#' is calculated under the normality assumtion.
#'
#' @param x a vector
#' @param W spatial connectivity matrix
#' @param alternative specification of alternative hypothesis as 'greater' (default),
#' 'lower', or 'two.sided'
#'
#' @return Returns an object of class \code{data.frame} that contains the
#' following information:
#' \tabular{lcl}{
#' \code{I}\tab \tab observed value of the Moran coefficient\cr
#' \code{EI}\tab\tab expected value of Moran's I\cr
#' \code{VarI}\tab \tab variance of Moran's I (under normality)\cr
#' \code{zI}\tab\tab standardized Moran coefficient\cr
#' \code{pI}\tab\tab \emph{p}-value of the test statistic
#' }
#'
#' @author Sebastian Juhl
#'
#' @references Cliff, Andrew D. and John K. Ord (1981): Spatial Processes:
#' Models & Applications. Pion.
#'
#' Bivand, Roger S. and David W. S. Wong (2018): Comparing
#' Implementations of Global and Local Indicators of Spatial Association.
#' TEST 27: pp. 716 - 748.
#'
#' Tiefelsdorf, Michael and Barry Boots (1995): The Exact Distribution
#' of Moran's I. Environment and Planning A: Economy and Space, 27 (6):
#' pp. 985 - 999.
#'
#' @seealso \code{\link{getMoran}}
#'
#' @export

MI.vec <- function(x,W,alternative="greater"){
  # convert x to a matrix and save names (if provided)
  x <- as.matrix(x)
  if(!is.null(colnames(x))) nams <- colnames(x)

  #####
  # Input
  # Checks
  #####
  if(0 %in% apply(x,2,sd)) warning("Constant term detected in x.")
  if(!any(class(W) %in% c("matrix","Matrix","data.frame"))){
    stop("W must be of class 'matrix' or 'data.frame'")
  }
  if(any(class(W)!="matrix")) W <- as.matrix(W)
  if(anyNA(x) | anyNA(W)) stop("Missing values detected")
  if(!(alternative %in% c("greater","lower", "two.sided"))){
    stop("Invalid input: 'alternative' must be either 'greater',
         'lower', or 'two.sided'")
  }

  #####
  # Additional
  # Variables
  #####
  nx <- ncol(x)
  n <- nrow(W)
  df <- n-1 # only one variable considered
  S0 <- t(rep(1,n)) %*% W %*% rep(1,n)
  S1 <- sum((W*W)+(W*t(W)))
  S2 <- sum((rowSums(W)+colSums(W))^2)

  # projection matrix M
  M <- diag(n)-rep(1,n)%*%t(rep(1,n))/n
  MWM <- M %*% W %*% M

  #####
  # Output
  #####
  out <- data.frame(matrix(NA,nrow=nx,ncol=6))
  colnames(out) <- c("I","EI","VarI","zI","pI","")
  for(i in 1:nx){
    # observed
    out[i,"I"] <- (n/S0) * t(x[,i]) %*% MWM %*% x[,i] / crossprod(x[,i],M) %*% x[,i]
    # expected
    out[i,"EI"] <- -1/(n-1)
    # variance (normality assumption)
    out[i,"VarI"] <- ((n^2*S1 - n*S2 + 3*S0^2) / (S0^2 * (n^2-1))) - out[i,"EI"]^2
    # test statistic
    out[i,"zI"] <- (out[i,"I"]-out[i,"EI"])/sqrt(out[i,"VarI"])
    # pI
    out[i,"pI"] <- pfunc(z=out[i,"zI"],alternative=alternative)
    out[i,6] <- star(p=out[i,"pI"])
  }
  if(!is.null(colnames(x))) rownames(out) <- nams

  return(out)
}
