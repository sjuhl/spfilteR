#' Fake Dataset
#'
#' An artificially generated cross-sectional dataset together with
#' an accompanying binary connectivity matrix \emph{\strong{W}} to
#' illustrate the use of the package. The \eqn{10 x 10} matrix connects
#' the \eqn{n=100} units based on rook's adjacency definition of contiguity.
#'
#' @docType data
#'
#' @usage data(toydata)
#'
#' @return The file contains two objects:
#' \tabular{lcl}{
#' \code{fakedataset}\tab \tab a synthetic dataset of dimensionality \eqn{100 x 8}\cr
#' \code{W}\tab\tab a binary connectivity matrix
#' }
#'
#' @keywords dataset
#' @examples
#'
#' data(fakedata)
#' head(fakedataset)
#' dim(W)

'fakedataset'

