#' @name vif.ev
#'
#' @title Variance Inflation Factor for Eigenvectors
#'
#' @description Calculate the variance inflation factor (VIF) for
#' the eigenvectors in the spatial filter.
#'
#' @param x vector/ matrix of regressors (default=NULL)
#' @param evecs (selected) eigenvectors
#' @param na.rm remove missing values in supplied covariates (TRUE/ FALSE)
#'
#' @return Returns a vector containing the VIF for each eigenvector.
#'
#' @note This function assumes a linear model.
#'
#' @examples
#' data(fakedata)
#' E <- getEVs(W=W,covars=NULL)$vectors
#' (VIF <- vif.ev(x=fakedataset$x1,evecs=E[,1:10]))
#'
#' @importFrom stats complete.cases
#'
#' @author Sebastian Juhl
#'
#' @seealso \code{\link{lmFilter}}, \code{\link{getEVs}}
#'
#' @export

vif.ev <- function(x=NULL,evecs,na.rm=TRUE){
  evecs <- as.matrix(evecs)
  if(is.null(x)) x <- rep(1,nrow(evecs))
  x <- as.matrix(x)
  if (!all(x[,1]==1)) x <- cbind(1,x)
  if(na.rm){
    evecs <- as.matrix(evecs[complete.cases(x),])
    x <- as.matrix(x[complete.cases(x),])
  }
  if(anyNA(x)) stop("Missing values detected")
  inflate <- NULL
  for(i in seq_len(ncol(evecs))){
    xi <- evecs[,i]
    tots <- sum((xi - mean(xi))^2)
    ests <- solve(crossprod(x)) %*% crossprod(x,xi)
    re <- xi - x %*% ests
    Rsq <- 1-(sum(crossprod(re))/tots)
    inflate[i] <- 1/(1-Rsq)
  }
  return(inflate)
}
