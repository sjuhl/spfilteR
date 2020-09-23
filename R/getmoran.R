#' @title Moran Coefficient for Residual Spatial Autocorrelation
#'
#' @description This function assesses the degree of spatial
#' autocorrelation in regression residuals by means of the Moran
#' coefficient.
#'
#' @param y vector of regressands
#' @param x vector/ matrix of regressors (default = NULL)
#' @param fitted.values a vector of fitted values from a regression model (default = NULL)
#' @param W spatial connectivity matrix
#' @param alternative specification of alternative hypothesis as 'greater' (default),
#' 'lower', or 'two.sided'
#' @param boot optional integer specifying the number of bootstrap iterations
#' to compute the variance. If NULL (default), variance calculated under normality
#'
#' @return A \code{data.frame} object with the following elements:
#' \tabular{lcl}{
#' \code{I}\tab \tab observed value of the Moran coefficient\cr
#' \code{EI}\tab\tab expected value of Moran's I\cr
#' \code{VarI}\tab \tab variance of Moran's I\cr
#' \code{zI}\tab\tab standardized Moran coefficient\cr
#' \code{pI}\tab\tab \emph{p}-value of the test statistic
#' }
#'
#' @author Sebastian Juhl
#'
#' @references Tiefelsdorf, Michael and Barry Boots (1995): The Exact Distribution
#' of Moran's I. Environment and Planning A: Economy and Space, 27 (6):
#' pp. 985 - 999.
#'
#' Tiefelsdorf, Michael and Daniel A. Griffith (2007):
#' Semiparametric filtering of spatial autocorrelation: the eigenvector
#' approach. Environment and Planning A: Economy and Space, 39 (5):
#' pp. 1193 - 1221.
#'
#' @seealso \code{\link{lmFilter}}, \code{\link{MI.vec}}
#'
#' @export

getMoran <- function(y,x=NULL,fitted.values=NULL,W,alternative="greater",boot=NULL){
  if(!(alternative %in% c("greater","lower", "two.sided"))){
    stop("Invalid input: 'alternative' must be either 'greater', 'lower', or 'two.sided'")
  }
  if(!any(class(W) %in% c("matrix","Matrix","data.frame"))){
    stop("W must be of class 'matrix' or 'data.frame'")
  }
  if(any(class(W)!="matrix")) W <- as.matrix(W)
  n <- nrow(W)
  if(is.null(x)) x <- rep(1,n)
  x <- as.matrix(x)
  if (!all(x[,1]==1)) x <- cbind(1,x) # add intercept term
  if(is.null(fitted.values)) fitted.values <- x %*% solve(crossprod(x), crossprod(x, y))
  resid <- y-fitted.values
  if(!isSymmetric(W)) W <- .5 * (W + t(W)) # symmetric connectivity matrix
  I <- n/crossprod(rep(1,n),W%*%rep(1,n)) * crossprod(resid,W%*%resid) / crossprod(resid)
  M <- diag(n)-x%*%qr.solve(crossprod(x),t(x))
  EI <- sum(diag(M%*%W))/(n-ncol(x))
  if(!is.null(boot)){
    if(boot<100){
      warning(paste0("Number of bootstrap iterations (",boot,") too small. Set to 100"))
      boot <- 100
    }
    boot.I <- NULL
    for(i in 1:boot){
      ind <- sample(1:n,replace=T)
      boot.I[i] <- n/crossprod(rep(1,n),W%*%rep(1,n)) * crossprod(resid[ind],W%*%resid[ind]) / crossprod(resid[ind])
    }
    VarI <- var(boot.I)
    zI <- (I-EI)/sqrt(VarI)
    pI <- emp.pfunc(draws=boot.I,obs=I,alternative=alternative)
  } else {
    num <- sum(diag(M%*%W%*%M%*%t(W))) + sum(diag(M%*%W%*%M%*%W)) + sum(diag(M%*%W))^2
    denom <- (n-ncol(x))*(n-ncol(x)+2)
    VarI <- num/denom - EI^2
    if(VarI<=0){
      zI <- 0
    } else zI <- (I-EI)/sqrt(VarI)
    pI <- pfunc(z=zI,alternative=alternative)
  }

  #####
  # Output
  #####
  out <- data.frame(I,EI,VarI,zI,pI,NA)
  colnames(out) <- c("I","EI","VarI","zI","pI","")
  out[1,6] <- star(p=out[1,"pI"])
  return(out)
}
