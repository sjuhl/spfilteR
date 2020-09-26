#' @name vif.ev
#'
#' @title Variance Inflation Factor for Eigenvectors
#'
#' @description Calculate the variance inflation factor (VIF) for
#' the eigenvectors in the spatial filter.
#'
#' @param x vector/ matrix of regressors (default=NULL)
#' @param evecs (selected) eigenvectors
#'
#' @return Returns a vector containing the VIF for each eigenvector.
#'
#' @examples
#' data("fakedata")
#' evecs <- getEVs(W=W)$vectors
#' evecs <- evecs[,1:10]
#' vif.ev(x=fakedataset$x1,evecs=evecs)
#'
#' @author Sebastian Juhl
#'
#' @seealso \code{\link{lmFilter}}, \code{\link{getEVs}}
#'
#' @export

vif.ev <- function(x=NULL,evecs){
  evecs <- as.matrix(evecs)
  if(is.null(x)) x <- rep(1,nrow(evecs))
  x <- as.matrix(x)
  if (!all(x[,1]==1)) x <- cbind(1,x)
  inflate <- NULL
  for(i in 1:ncol(evecs)){
    xi <- evecs[,i]
    tots <- sum((xi - mean(xi))^2)
    ests <- solve(crossprod(x)) %*% crossprod(x,xi)
    re <- xi - x %*% ests
    Rsq <- 1-(sum(crossprod(re))/tots)
    inflate[i] <- 1/(1-Rsq)
  }
  return(inflate)
}
