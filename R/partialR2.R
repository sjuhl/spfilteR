#' @title Coefficient of Partial Determination
#'
#' @description This function computes the partial R-squared of all
#' selected eigenvectors in a (filtered) linear regression model.
#'
#' @param y vector of regressands
#' @param x matrix of regressors
#' @param evecs (selected) eigenvectors
#'
#' @return Vector of partial R-squared values for each eigenvector.
#'
#' @author Sebastian Juhl
#'
#' @seealso \code{\link{lmFilter}}, \code{\link{getEVs}}
#'
#' @export

partialR2 <- function(y,x,evecs){
  # if no evecs are supplied
  if(is.null(evecs)){ warning("No eigenvectors supplied"); return(pR2 <- NULL)}

  #####
  # Extract Information
  # & Formatting
  #####
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
  resid_wo <- y - x %*% solve(crossprod(x)) %*% crossprod(x,y)
  pR2 <- rep(NA,nev)
  names(pR2) <- names_evecs
  # loop over all eigenvectors
  for(i in 1:nev){
    xev <- cbind(x,evecs[,i])
    resid_ev <- y - xev %*% solve(crossprod(xev)) %*% crossprod(xev,y)
    pR2[i] <- (sum(crossprod(resid_wo)) - sum(crossprod(resid_ev)))/sum(crossprod(resid_wo))
  }

  #####
  # Output
  #####
  return(pR2)
}