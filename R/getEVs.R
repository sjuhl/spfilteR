#' @name getEVs
#'
#' @title Eigenfunction Decomposition of a (Transformed) Spatial Connectivity
#' Matrix
#'
#' @description Extract eigenvectors and corresponding eigenvalues from
#' the matrix \emph{\strong{MWM}}, where \emph{\strong{M}} denotes a symmetric
#' and idempotent projection matrix and \emph{\strong{W}} is the spatial
#' connectivity matrix.
#'
#' @param W spatial connectivity matrix
#' @param covars vector/ matrix of regressors included in the construction
#' of the projection matrix \emph{\strong{M}} - see Details below
#'
#' @return A list containing the following objects:
#' \describe{
#' \item{\code{vectors}}{a matrix of \emph{n} eigenvectors}
#' \item{\code{values}}{a vector of the corresponding eigenvalues}
#' \item{\code{moran}}{a vector of the Moran coefficients associated with
#' the eigenvectors}
#' }
#'
#' @details The eigenfunctions obtained by \code{getEVs}
#' can be used to perform supervised eigenvector selection and to
#' manually create a spatial filter. To this end, a candidate set
#' might be determined by 1) the sign of the spatial autocorrelation
#' in the outcome variable and 2) the strength of spatial association
#' found in each eigenvector as indicated by \code{moran}.
#'
#' If \emph{\strong{W}} is not symmetric, \code{getEVs} symmetrizes the
#' matrix by: 0.5 * (\emph{\strong{W}} + \emph{\strong{W}}').
#'
#' If \code{covars} are supplied, the function uses the covariates to construct
#' projection matrix: \emph{\strong{M} = \strong{I} - \strong{X} (\strong{X}'
#' \strong{X})^-1\strong{X}'}. Using this matrix results in a set of
#' eigenvectors that are uncorrelated to each other as well as to the
#' covariates. If \code{covars=NULL}, only the intercept term is used
#' to construct \emph{\strong{M}}. See e.g., Griffith and Tiefelsdorf (2007)
#' for more details on the appropriate choice of \emph{\strong{M}}.
#'
#' @author Sebastian Juhl
#'
#' @references Tiefelsdorf, Michael and Daniel A. Griffith (2007):
#' Semiparametric filtering of spatial autocorrelation: the eigenvector
#' approach. Environment and Planning A: Economy and Space, 39 (5):
#' pp. 1193 - 1221.
#'
#' @seealso \code{\link{lmFilter}}, \code{\link{glmFilter}}, \code{\link{MI.ev}},
#' \code{\link{MI.sf}}, \code{\link{vif.ev}}, \code{\link{partialR2}}
#'
#' @examples
#' data(fakedata)
#'
#' E <-getEVs(W=W,covars=NULL)
#'
#' @export

getEVs <- function(W,covars=NULL){
  n <- nrow(W)
  # symmetric connectivity matrix V
  # Note: if W is symmetric, V == W
  V <- .5 * (W + t(W))

  # projection matrix M (orthogonal eigenvectors)
  if(is.null(covars)) covars <- rep(1,n)
  covars <- as.matrix(covars)
  if (!all(covars[,1]==1)) covars <- cbind(1,covars)

  M <- diag(n)-covars%*%qr.solve(crossprod(covars),t(covars))
  # if no covars, M equals diag(n)-rep(1,n)%*%t(rep(1,n))/n

  # MWM
  MVM <- M %*% V %*% M

  # eigenvectors and eigenvalues
  eigs <- eigen(MVM,symmetric=T)

  # Moran coefficient for eigenvectors
  moran <- MI.ev(W=W,evals=eigs$values)

  # return
  return(list(vectors=eigs$vectors
              ,values=eigs$values
              ,moran=moran))
}
