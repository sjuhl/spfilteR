#' @name MI.ev
#'
#' @title Moran Coefficients of Eigenvectors
#'
#' @description Calculates the Moran coefficient for each eigenvector
#' included in the spatial filter.
#'
#' @param W spatial connectivity matrix
#' @param evals eigenvalues of selected eigenvectors
#'
#' @return Returns a vector containing the Moran coefficients of each
#' selected eigenvector.
#'
#' @author Sebastian Juhl
#'
#' @references Le Gallo, Julie and Antonio Páez (2013): Using synthetic
#' variables in instrumental variable estimation of spatial series models.
#' Environment and Planning A, 45 (9): pp. 2227 - 2242.
#'
#' Tiefelsdorf, Michael and Barry Boots (1995): The Exact Distribution
#' of Moran's I. Environment and Planning A: Economy and Space, 27 (6):
#' pp. 985 - 999.
#'
#' @seealso \code{\link{lmFilter}}, \code{\link{getEVs}}, \code{\link{MI.sf}}
#'
#' @export

MI.ev <- function(W,evals){
  n <- nrow(W)
  evMI <- sapply(evals,function(evals) n/crossprod(rep(1,n),W%*%rep(1,n)) * evals)
  return(evMI)
}
