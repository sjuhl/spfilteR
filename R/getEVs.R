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
#' of the projection matrix \emph{\strong{M}}
#'
#' @return A list containing the following objects:
#' \describe{
#' \item{\code{vectors}}{A matrix containing \emph{n} eigenvectors}
#' \item{\code{values}}{A vector of the corresponding eigenvalues}
#' }
#'
#' @details The eigenfunctions obtained by \code{getEVs}
#' can be used to perform supervised eigenvector selection and to
#' manually create a spatial filter.
#'
#' If \emph{\strong{W}} is not symmetric, \code{getEVs} symmetrizes the
#' matrix by: \eqn{0.5 * (W + t(W))}.
#'
#' @author Sebastian Juhl
#'
#' @seealso \code{\link{lmFilter}}, \code{\link{MI.ev}}, \code{\link{MI.sf}},
#' \code{\link{vif.ev}}, \code{\link{partialR2}}
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

  # return
  return(list(vectors=eigs$vectors
              ,values=eigs$values))
}
