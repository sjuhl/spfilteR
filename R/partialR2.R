#' @name partialR2
#'
#' @title Coefficient of Partial Determination
#'
#' @description This function computes the partial R-squared of all
#' selected eigenvectors in a spatially filtered linear regression model.
#'
#' @param y vector of regressands
#' @param x vector/ matrix of regressors
#' @param evecs (selected) eigenvectors
#'
#' @return Vector of partial R-squared values for each eigenvector.
#'
#' @note The function assumes a linear regression model. Since the
#' eigenvectors are mutually uncorrelated, \code{partialR2} evaluates
#' them sequentially. In generalized linear models, the presence of a link
#' function can corrupt the uncorrelatedness of the eigenvectors.
#'
#' @author Sebastian Juhl
#'
#' @seealso \code{\link{lmFilter}}, \code{\link{getEVs}}
#'
#' @examples
#' data(fakedata)
#' y <- fakedataset$x1
#' x <- fakedataset$x2
#'
#' # get eigenvectors
#' E <-getEVs(W=W,covars=NULL)$vectors
#'
#' (out <- partialR2(y=y,x=x,evecs=E[,1:5]))
#'
#'
#' @export

partialR2 <- function(y,x=NULL,evecs){
  # if no evecs are supplied
  if(is.null(evecs)){ warning("No eigenvectors supplied"); return(pR2 <- NULL)}

  #####
  # Extract Information
  # & Formatting
  #####
  if(is.null(x)) x <- rep(1,length(y))
  x <- as.matrix(x)
  if (!all(x[,1]==1)) x <- cbind(1,x)
  evecs <- as.matrix(evecs)
  nev <- ncol(evecs)
  # store names
  names_evecs <- colnames(evecs)

  #####
  # Input Checks
  #####
  if(anyNA(y) | anyNA(x)) stop("Missing values detected")
  if(qr(x)$rank!=ncol(x)) stop("Perfect multicollinearity in covariates detected")

  # R2 reduced model (without eigenvectors)
  fitvals <- x %*% solve(crossprod(x)) %*% crossprod(x,y)
  resid_wo <- residfun(y=y,fitvals=fitvals,model="linear")$raw
  pR2 <- rep(NA,nev)
  names(pR2) <- names_evecs
  # loop over all eigenvectors
  for(i in 1:nev){
    xev <- cbind(x,evecs[,i])
    fitvals_ev <- xev %*% solve(crossprod(xev)) %*% crossprod(xev,y)
    resid_ev <- residfun(y=y,fitvals=fitvals_ev,model="linear")$raw
    pR2[i] <- (sum(crossprod(resid_wo)) - sum(crossprod(resid_ev)))/sum(crossprod(resid_wo))
  }

  #####
  # Output
  #####
  return(pR2)
}
