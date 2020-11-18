#' @name MI.decomp
#'
#' @title Decomposition of the Moran Coefficient
#'
#' @description A decomposition of the Moran coefficient in order to separately
#' test for the presence of positive and negative autocorrelation in a variable.
#'
#' @param x a vector or matrix
#' @param W spatial connectivity matrix
#' @param nsim number of iterations to simulate the null distribution
#'
#' @return Returns a \code{data.frame} that contains the following information
#' for each variable:
#' \describe{
#' \item{\code{I+}}{observed value of Moran's I (positive part)}
#' \item{\code{VarI+}}{variance of Moran's I (positive part)}
#' \item{\code{pI+}}{simulated \emph{p}-value of Moran's I (positive part)}
#' \item{\code{I-}}{observed value of Moran's I (negative part)}
#' \item{\code{VarI-}}{variance of Moran's I (negative part)}
#' \item{\code{pI-}}{simulated \emph{p}-value of Moran's I (negative part)}
#' \item{\code{pItwo.sided}}{simulated \emph{p}-value of the two-sided test}
#' }
#'
#' @details If \code{x} is a matrix, this function computes the Moran
#' test for spatial autocorrelation for each column.
#'
#' The \emph{p}-values calculated for \code{I+} and \code{I-} assume
#' a directed alternative hypothesis.
#'
#' @author Sebastian Juhl
#'
#' @references Dary, Stéphane (2011): A New Perspective about Moran’s
#' Coefficient: Spatial Autocorrelation as a Linear Regression Problem.
#' Geographical Analysis, 43 (2): pp. 127 - 141.
#'
#' @seealso \code{\link{MI.vec}}, \code{\link{MI.ev}}, \code{\link{MI.sf}},
#' \code{\link{MI.resid}}, \code{\link{getEVs}}
#'
#' @importFrom stats var cor sd
#'
#' @examples
#' data(fakedata)
#' X <- cbind(fakedataset$x1,fakedataset$x2
#' ,fakedataset$x3,fakedataset$negative)
#'
#' (MI.dec <- MI.decomp(x=X,W=W,nsim=100))
#'
#' # the sum of I+ and I- equals the observed Moran coefficient:
#' I <- MI.vec(x=X,W=W)[,"I"]
#' cbind(MI.dec[,"I+"] + MI.dec[,"I-"], I)
#'
#' @export

MI.decomp <- function(x,W,nsim=100){
  # convert x to a matrix and save names (if provided)
  x <- as.matrix(x)
  if(!is.null(colnames(x))) nams <- colnames(x)

  #####
  # Input
  # Checks
  #####
  if(0 %in% apply(x,2,sd)) warning("Constant term detected in x")
  if(!any(class(W) %in% c("matrix","Matrix","data.frame"))){
    stop("W must be of class 'matrix' or 'data.frame'")
  }
  if(any(class(W)!="matrix")) W <- as.matrix(W)
  if(anyNA(x) | anyNA(W)) stop("Missing values detected")
  if(nsim<100){
    warning(paste0("Number of permutations (",nsim,") too small. Set to 100"))
    nsim <- 100
  }

  # symmetrize W
  W <- .5 * (W + t(W))

  # z-score transformation
  nx <- ncol(x)
  Z <- Zscore(x)

  # expected value
  n <- nrow(W)
  EI <- -1/(n-1)

  # eigendecomposition
  eigen <- getEVs(W=W)

  # null distribution
  Ip <- In <- matrix(NA,nrow=nx,ncol=nsim)
  for(i in seq_len(nsim)){
    ind <- sample(1:n,replace=TRUE)
    c2 <- cor(Z[ind,],eigen$vectors, method="pearson")^2
    Ip[,i] <- c2[,eigen$values>EI] %*% eigen$moran[eigen$values>EI]
    In[,i] <- c2[,eigen$values<EI] %*% eigen$moran[eigen$values<EI]
  }

  #####
  # Output
  #####
  out <- data.frame(matrix(NA,nrow=nx,ncol=10))
  colnames(out) <- c("I+","VarI+","pI+",""
                     ,"I-","VarI-","pI-",""
                     ,"pItwo.sided","")
  # observed
  cor2 <- cor(Z,eigen$vectors, method="pearson")^2
  out[,"I+"] <- cor2[,eigen$values>EI] %*% eigen$moran[eigen$values>EI]
  out[,"I-"] <- cor2[,eigen$values<EI] %*% eigen$moran[eigen$values<EI]
  # variance
  out[,"VarI+"] <- apply(Ip,1,var)
  out[,"VarI-"] <- apply(In,1,var)
  # significance
  for(i in seq_len(nx)){
    out[i,"pI+"] <- pfunc(z=out[i,"I+"],alternative="greater",draws=Ip[i,])
    out[i,"pI-"] <- pfunc(z=out[i,"I-"],alternative="lower",draws=In[i,])
    out[i,"pItwo.sided"] <- 2*min(out[i,c("pI+","pI-")])
    out[i,4] <- star(p=out[i,"pI+"])
    out[i,8] <- star(p=out[i,"pI-"])
    out[i,10] <- star(p=out[i,"pItwo.sided"])
  }
  if(!is.null(colnames(x))) rownames(out) <- nams

  return(out)
}
