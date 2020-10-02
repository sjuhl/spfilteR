#' @name MI.vec
#'
#' @title Moran Test for Spatial Autocorrelation in a Vector
#'
#' @description Tests for the presence of spatial autocorrelation
#' in a variable as indicated by the Moran coefficient. The variance
#' is calculated under the normality assumtion.
#'
#' @param x a vector or matrix
#' @param W spatial connectivity matrix
#' @param alternative specification of alternative hypothesis as 'greater' (default),
#' 'lower', or 'two.sided'
#' @param symmetrize symmetrizes the connectivity matrix \emph{\strong{W}}
#' by: 1/2 * (\emph{\strong{W}} + \emph{\strong{W}}') (TRUE/ FALSE).
#'
#' @return Returns an object of class \code{data.frame} that contains the
#' following information for each variable:
#' \describe{
#' \item{\code{I}}{observed value of the Moran coefficient}
#' \item{\code{EI}}{expected value of Moran's I}
#' \item{\code{VarI}}{variance of Moran's I (under normality)}
#' \item{\code{zI}}{standardized Moran coefficient}
#' \item{\code{pI}}{\emph{p}-value of the test statistic}
#' }
#'
#' @note Estimation of the variance (under the normality assumption)
#' follows Cliff and Ord (1981), see also Upton and Fingleton (1985).
#' It assumes the connectivity matrix \emph{\strong{W}} to be symmetric.
#' For inherently non-symmetric matrices, it is recommended to specify
#' \code{symmetrize=TRUE}.
#'
#' @author Sebastian Juhl
#'
#' @references Cliff, Andrew D. and John K. Ord (1981): Spatial Processes:
#' Models & Applications. Pion, London.
#'
#' Upton, Graham J. G. and Bernard Fingleton (1985): Spatial Data Analysis
#' by Example, Volume 1. New York, Wiley.
#'
#' Bivand, Roger S. and David W. S. Wong (2018): Comparing Implementations
#' of Global and Local Indicators of Spatial Association. TEST 27:
#' pp. 716 - 748.
#'
#' @seealso \code{\link{getMoran}}
#'
#' @importFrom stats sd
#'
#' @examples
#' data(fakedata)
#' X <- cbind(fakedataset$x1,fakedataset$x2
#' ,fakedataset$x3,fakedataset$x4)
#'
#' MI <- MI.vec(x=X,W=W,alternative="greater",symmetrize=TRUE)
#' MI
#'
#' @export

MI.vec <- function(x,W,alternative="greater",symmetrize=TRUE){
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
  if(symmetrize) W <- .5 * (W + t(W))
  nx <- ncol(x)
  n <- nrow(W)
  df <- n-1 # only one variable considered at a time
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
