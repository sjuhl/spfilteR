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
#' @seealso \code{\link{lmfilter}}, \code{\link{getEVs}}, \code{\link{MI.ev}},
#' \code{\link{MI.sf}}, \code{\link{vif.ev}}, \code{\link{partialR2}},
#' \code{\link{candsetsize}}
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



#' @title Variance Inflation Factor for Eigenvectors
#'
#' @description Calculate the variance inflation factor (VIF) for
#' the eigenvectors in the spatial filter
#'
#' @param x Matrix of regressors (default = NULL)
#' @param evecs Selected eigenvectors for inclusion in the spatial filter
#'
#' @return Returns a vector containing the VIF for each selected eigenvector
#'
#' @author Sebastian Juhl
#'
#' @seealso \code{\link{lmfilter}}, \code{\link{getEVs}}
#'
#' @export

vif.ev <- function(x=NULL,evecs){
  evecs <- as.matrix(evecs)
  if(is.null(x)) x <- as.matrix(rep(1,nrow(evecs)))
  if (!all(x[,1]==1)) x <- cbind(1,x)
  inflate <- NULL
  for(i in 1:ncol(evecs)){
    xi <- evecs[,i]
    tots <- sum((xi - mean(xi))^2)
    ests <- solve(crossprod(x)) %*% crossprod(x,xi)
    re <- xi - x %*% ests
    Rsq <- 1-(sum(crossprod(re))/tots)
    inflate[i] <- 1/(1-Rsq)
  }
  return(inflate)
}


#' @title Calculate \emph{p}-Value for Standardized Moran Coefficient
#'
#' @description Derive \emph{p}-value for the standardized Moran coefficient
#'
#' @param z Standardized Moran coefficient
#' @param alternative Optional: alternative hypothesis to be assessed. Default
#' is 'two-sided'
#'
#' @return Returns the \emph{p}-value associated with the supplied
#' (standardized) Moran coefficient
#'
#' @author Sebastian Juhl
#'
#' @seealso \code{\link{lmfilter}}, \code{\link{emp.pfunc}}
#'
#' @export

pfunc <- function(z,alternative){
  if(alternative=="greater"){
    p <- pnorm(z, lower.tail=F)
  } else if(alternative=="lower"){
    p <- pnorm(z, lower.tail=T)
  } else p <- 2*pnorm(abs(z), lower.tail=F)
  return(p)
}


#' @title Calculate empirical \emph{p}-Value for Standardized
#' Moran Coefficient
#'
#' @description Derive empirical \emph{p}-value for the standardized
#' Moran coefficient
#'
#' @param draws Number of simulation iterations
#' @param z Observed value for the standardized Moran coefficient
#' @param alternative Optional: alternative hypothesis to be assessed. Default
#' is 'two-sided'
#'
#' @return Returns the empirical \emph{p}-value associated with the supplied
#' (standardized) Moran coefficient based on a simulated null distribution
#'
#' @author Sebastian Juhl
#'
#' @seealso \code{\link{lmfilter}}, \code{\link{pfunc}}
#'
#' @export

emp.pfunc <- function(draws,z,alternative){
  z <- as.numeric(z)
  # see e.g., North/ Curtis/ Sham (2002) [Am J Hum Genet] for the '+1'
  if(alternative=="greater"){
    p <- (sum(draws>=z)+1)/(length(draws)+1)
  } else if(alternative=="lower"){
    p <- (sum(draws<=z)+1)/(length(draws)+1)
  } else {
    p <- (sum(abs(draws)>=abs(z))+1)/(length(draws)+1) # see e.g., Hartwig (2013) [J Clin Trials]
  }
  return(p)
}


#' @title Graphical Illustration of Significance Levels
#'
#' @description Illustrate significance levels for output
#'
#' @param p \emph{p}-value
#'
#' @return Returns the empirical \emph{p}-value associated with the supplied
#' (standardized) Moran coefficient based on a simulated null distribution
#'
#' @author Sebastian Juhl
#'
#' @seealso \code{\link{lmfilter}}, \code{\link{pfunc}}, \code{\link{emp.pfunc}}
#'
#' @export

star <- function(p){
  out <- NULL
  out[p<=.001] <- "***"
  out[p<=.01 & p>.001] <- "**"
  out[p<=.05 & p>.01] <- "*"
  out[p<=.1 & p>.05] <- "."
  out[p>.1] <- " "
  return(out)
}

#' @title Moran Coefficients for Eigenvectors
#'
#' @description Moran coefficient for each eigenvector in the spatial filter
#'
#' @param W Spatial connectivity matrix
#' @param evals Eigenvalues of selected eigenvectors
#'
#' @return Returns a vector containing the Moran coefficients of each
#' selected eigenvector
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
#' @seealso \code{\link{lmfilter}}, \code{\link{getEVs}}, \code{\link{MI.sf}}
#'
#' @export

MI.ev <- function(W,evals){ # evals = eigenvalues
  n <- nrow(W)
  evMI <- sapply(evals,function(evals) n/crossprod(rep(1,n),W%*%rep(1,n)) * evals)
  return(evMI)
}


#' @title Moran Coefficient of the Spatial Filter
#'
#' @description Moran coefficient for the spatial filter
#'
#' @param gamma Regression coefficient associated with the eigenvectors
#' @param evMI Moran coefficient of selected eigenvectors
#'
#' @return Moran coefficient of the spatial filter
#'
#' @author Sebastian Juhl
#'
#' @references Le Gallo, Julie and Antonio Páez (2013): Using synthetic
#' variables in instrumental variable estimation of spatial series models.
#' Environment and Planning A: Economy and Space, 45 (9): pp. 2227 - 2242.
#'
#' @seealso \code{\link{lmfilter}}, \code{\link{getEVs}}, \code{\link{MI.ev}}
#'
#' @export

MI.sf <- function(gamma,evMI){
  sfMI <- 1/sum(gamma^2) * sum(gamma^2 * evMI)
  return(sfMI)
}


#' @title Partial R-squared
#'
#' @description Partial R-squared of selected eigenvectors
#'
#' @param y Vector of regressands
#' @param x Matrix of regressors
#' @param evecs Selected eigenvectors for inclusion in the spatial filter
#'
#' @return Vector of partial R-squared values for each eigenvector
#'
#' @author Sebastian Juhl
#'
#' @seealso \code{\link{lmfilter}}, \code{\link{getEVs}}
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

#' @title Ideal Candidate Set Size
#'
#' @description Determine the ideal size of the candidate set
#'
#' @param npos Regression coefficient associated with the eigenvectors
#' @param zMI Standardized Moran coefficient of regression residuals
#'
#' @return An integer indicating the ideal size of the candidate set
#'
#' @author Sebastian Juhl
#'
#' @references Chun, Yongwan, Daniel A. Griffith, Monghyeon Lee, Parmanand
#' Sinha (2016): Eigenvector selection with stepwise regression techniques
#' to construct eigenvector spatial filters. Journal of Geographical
#' Systems, 18, pp. 67 – 85.
#'
#' @seealso \code{\link{lmfilter}}, \code{\link{getEVs}}
#'
#' @export

candsetsize <- function(npos, zMI){
  denominator <- 1+exp(2.1480-(6.1808*(zMI+.6)^.1742)/npos^.1298 + 3.3534/(zMI+.6)^.1742)
  nc <- npos/denominator
  return(round(nc,0))
}
