#' @title Extract Eigenvectors and Eigenvalues
#'
#' @description Extract eigenvectors and the corresponding eigenvalues from
#' the matrix \emph{\strong{MWM}}, where \emph{\strong{M}} denotes a symmetric
#' and idempotent projection matrix and \emph{\strong{W}} is the spatial
#' connectivity matrix.
#'
#' @param W Spatial connectivity matrix
#' @param covars Design matrix of regressors included in the construction
#' of the projection matrix \emph{\strong{M}}
#'
#' @return A list containing the following objects:
#' \describe{
#' \item{\code{vectors}}{A matrix containing \emph{n} eigenvectors}
#' \item{\code{values}}{A vector of the corresponding eigenvalues}
#' }
#'
#' @author Sebastian Juhl
#'
#' @seealso \code{\link{lmfilter}}, \code{\link{MI.ev}}, \code{\link{MI.sf}},
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
