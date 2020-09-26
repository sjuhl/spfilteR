#' @name MI.sf
#'
#' @title Moran Coefficient of the Spatial Filter
#'
#' @description Computes the Moran coefficient for the spatial filter.
#'
#' @param gamma regression coefficient associated with the eigenvectors
#' @param evMI Moran coefficient of selected eigenvectors
#'
#' @return Moran coefficient of the spatial filter.
#'
#' @author Sebastian Juhl
#'
#' @references Le Gallo, Julie and Antonio Páez (2013): Using synthetic
#' variables in instrumental variable estimation of spatial series models.
#' Environment and Planning A: Economy and Space, 45 (9): pp. 2227 - 2242.
#'
#' @seealso \code{\link{lmFilter}}, \code{\link{getEVs}}, \code{\link{MI.ev}}
#'
#' @export

MI.sf <- function(gamma,evMI){
  sfMI <- 1/sum(gamma^2) * sum(gamma^2 * evMI)
  return(sfMI)
}
