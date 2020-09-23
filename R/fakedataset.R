#' Fake Dataset
#'
#' An artificially generated cross-sectional dataset together with
#' an accompanying binary connectivity matrix \emph{\strong{W}}. The
#' \eqn{n=100} units are located on a regular grid and \emph{\strong{W}}
#' is defined according to rook's adjacency definition of contiguity.
#' The synthetic data can be used to illustrate the functionality
#' of this package.
#'
#' @usage data(fakedata)
#'
#' @docType data
#'
#' @return The file contains two objects:
#' \tabular{lcl}{
#' \code{fakedataset}\tab \tab a synthetic dataset\cr
#' \code{W}\tab\tab an artificial \eqn{100 x 100} connectivity matrix
#' }
#'
#' @keywords dataset
#' @examples
#'
#' data(fakedata)
#' head(fakedataset)
#' dim(W)

'fakedataset'

