#' @name lmFilter
#'
#' @title Unsupervised Spatial Filtering with Eigenvectors in Linear Regression Models
#'
#' @description This function implements the eigenvector-based semiparametric
#' spatial filtering approach in a linear regression framework using OLS.
#' Eigenvectors are selected by an unsupervised stepwise regression
#' technique. Supported selection criteria are the minimization of residual
#' autocorrelation, maximization of model fit, significance of residual autocorrelation,
#' and the statistical significance of eigenvectors. Alternatively, all eigenvectors in
#' the candidate set can be included as well.
#'
#' @param y vector of regressands
#' @param x vector/ matrix of regressors (default=NULL)
#' @param W spatial connectivity matrix
#' @param objfn the objective function to be used for eigenvector
#' selection. Possible criteria are: the maximization of the
#' adjusted R-squared ('R2'), minimization of residual autocorrelation ('MI'),
#' significance level of candidate eigenvectors ('p'), significance of residual spatial
#' autocorrelation ('pMI') or all eigenvectors in the candidate set ('all')
#' @param MX covariates used to construct the projection matrix (default=NULL) - see
#' Details
#' @param sig significance level to be used for eigenvector selection
#' if \code{objfn='p'} or \code{objfn='pMI'}
#' @param bonferroni Bonferroni adjustment for the significance level
#' (TRUE/ FALSE) if \code{objfn='p'}. Set to FALSE if \code{objfn='pMI'} -
#' see Details
#' @param positive restrict search to eigenvectors associated with positive
#' levels of spatial autocorrelation (TRUE/ FALSE)
#' @param ideal.setsize if \code{positive=TRUE}, uses the formula proposed in
#' Chun et al. (2016) to determine the ideal size of the candidate set
#' (TRUE/ FALSE)
#' @param alpha a value in (0,1] indicating the range of candidate eigenvectors
#' according to their associated level of spatial autocorrelation, see e.g.,
#' Griffith (2003)
#' @param tol if \code{objfn='MI'}, determines the amount of remaining residual
#' autocorrelation at which the eigenvector selection terminates
#' @param boot.MI number of iterations used to estimate the variance of Moran's I.
#' If \code{boot=NULL} (default), analytical results will be used
#' @param na.rm remove observations with missing values (TRUE/ FALSE)
#' @param object an object of class \code{spfilter}
#' @param EV display summary statistics for selected eigenvectors (TRUE/ FALSE)
#' @param ... additional arguments
#'
#' @return An object of class \code{spfilter} containing the following
#' information:
#' \describe{
#' \item{\code{Estimates}}{summary statistics of the parameter estimates}
#' \item{\code{varcovar}}{estimated variance-covariance matrix}
#' \item{\code{EV}}{a matrix with summary statistics of selected eigenvectors}
#' \item{\code{selvecs}}{vector/ matrix of selected eigenvectors}
#' \item{\code{evMI}}{Moran coefficient of all eigenvectors}
#' \item{\code{moran}}{residual autocorrelation for the initial and the
#' filtered model}
#' \item{\code{fit}}{adjusted R-squared of the initial and the filtered model}
#' \item{\code{residuals}}{initial and filtered model residuals}
#' \item{\code{other}}{a list providing supplementary information:
#' \describe{
#' \item{\code{ncandidates}}{number of candidate eigenvectors considered}
#' \item{\code{nev}}{number of selected eigenvectors}
#' \item{\code{sel_id}}{ID of selected eigenvectors}
#' \item{\code{sf}}{vector representing the spatial filter}
#' \item{\code{sfMI}}{Moran coefficient of the spatial filter}
#' \item{\code{model}}{type of the regression model}
#' \item{\code{dependence}}{filtered for positive or negative spatial dependence}
#' \item{\code{objfn}}{selection criteria specified in the objective function of
#' the stepwise regression procedure}
#' \item{\code{bonferroni}}{TRUE/ FALSE: Bonferroni-adjusted significance level
#' (if \code{objfn='p'})}
#' \item{\code{siglevel}}{if \code{objfn='p'} or \code{objfn='pMI'}: actual
#' (unadjusted/ adjusted) significance level}
#' }
#' }
#' }
#'
#' @details If \emph{\strong{W}} is not symmetric, it gets symmetrized by
#' 1/2 * (\emph{\strong{W}} + \emph{\strong{W}}') the eigenfunction decomposition.
#'
#' If covariates are supplied to \code{MX}, the function uses these regressors
#' to construct the following projection matrix:
#'
#' \emph{\strong{M} = \strong{I} - \strong{X} (\strong{X}'\strong{X})^-1\strong{X}'}
#'
#' Eigenvectors from \emph{\strong{MWM}} using this specification of
#' \emph{\strong{M}} are not only mutually uncorrelated but also orthogonal
#' to the regressors specified in \code{MX}. Alternatively, if \code{MX=NULL}, the
#' projection matrix becomes \emph{\strong{M} = \strong{I} - \strong{11}'/\emph{n}},
#' where \emph{\strong{1}} is a vector of ones and \emph{n} represents the number of
#' observations. Griffith and Tiefelsdorf (2007) show how the choice of the appropriate
#' \emph{\strong{M}} depends on the underlying process that generates the spatial
#' dependence.
#'
#' The Bonferroni correction is only possible if eigenvector selection is based on
#' the significance level of the eigenvectors (\code{objfn='p'}). It is set to
#' FALSE if eigenvectors are added to the model until the residuals exhibit no
#' significant level of spatial autocorrelation (\code{objfn='pMI'}).
#'
#' @examples
#' data(fakedata)
#' y <- fakedataset$x1
#' X <- cbind(fakedataset$x2,fakedataset$x3,fakedataset$x4)
#'
#' res <- lmFilter(y=y,x=X,W=W,objfn='MI',positive=FALSE)
#' print(res)
#' summary(res,EV=TRUE)
#'
#' E <- res$selvecs
#' (ols <- coef(lm(y~X+E)))
#' coef(res)
#'
#' @references Tiefelsdorf, Michael and Daniel A. Griffith (2007):
#' Semiparametric filtering of spatial autocorrelation: the eigenvector
#' approach. Environment and Planning A: Economy and Space, 39 (5):
#' pp. 1193 - 1221.
#'
#' Griffith, Daniel A. (2003): Spatial Autocorrelation and Spatial Filtering:
#' Gaining Understanding Through Theory and Scientific Visualization.
#' Berlin/ Heidelberg, Springer.
#'
#' Chun, Yongwan, Daniel A. Griffith, Monghyeon Lee, Parmanand
#' Sinha (2016): Eigenvector selection with stepwise regression techniques
#' to construct eigenvector spatial filters. Journal of Geographical
#' Systems, 18, pp. 67 – 85.
#'
#' Le Gallo, Julie and Antonio Páez (2013): Using synthetic
#' variables in instrumental variable estimation of spatial series models.
#' Environment and Planning A: Economy and Space, 45 (9): pp. 2227 - 2242.
#'
#' Tiefelsdorf, Michael and Barry Boots (1995): The Exact Distribution
#' of Moran's I. Environment and Planning A: Economy and Space, 27 (6):
#' pp. 985 - 999.
#'
#' @importFrom stats pt sd
#'
#' @seealso \code{\link{glmFilter}}, \code{\link{getEVs}}, \code{\link{MI.resid}}
#'
#' @export

lmFilter <- function(y,x=NULL,W,objfn="MI",MX=NULL,sig=.05
                     ,bonferroni=TRUE,positive=TRUE,ideal.setsize=FALSE
                     ,alpha=.25,tol=.1,boot.MI=NULL,na.rm=TRUE){

  if(!is.null(MX)) MX <- as.matrix(MX)
  if(!is.null(x)) x <- as.matrix(x)
  if(!is.null(colnames(x))) nams <- colnames(x) else nams <- NULL

  # missing values
  if(na.rm){
    if(!is.null(x)){
      miss <- apply(cbind(y,x),1,anyNA)
      x <- as.matrix(x[!miss,])
    } else miss <- is.na(y)
    y <- y[!miss]
    W <- W[!miss,!miss]
    if(!is.null(MX)) MX[!miss,]
  }

  # number of observations
  n <- length(y)

  # add intercept if not included in x
  if(is.null(x)) x <- as.matrix(rep(1,n))
  if(!all(x[,1]==1)) x <- cbind(1,x)
  if(!is.null(MX) && any(apply(MX,2,sd)==0)) MX <- as.matrix(MX[,apply(MX,2,sd)!=0])
  nx <- ncol(x)

  #####
  # Input Checks
  #####
  if(anyNA(y) | anyNA(x) | anyNA(W)) stop("Missing values detected")
  if(alpha==0) alpha <- 1e-07
  if(alpha<1e-07 | alpha>1){
    stop("Invalid argument: 'alpha' must be in the interval (0,1]")
  }
  if(qr(x)$rank!=ncol(x)) stop("Perfect multicollinearity in covariates detected")
  if(!any(class(W) %in% c("matrix","Matrix","data.frame"))){
    stop("W must be of class 'matrix' or 'data.frame'")
  }
  if(any(class(W)!="matrix")) W <- as.matrix(W)
  if(!(objfn %in% c("R2","p","MI","pMI","all"))){
    stop("Invalid argument: objfn must be one of 'R2', 'p', 'MI', 'pMI', or'all'")
  }
  if(positive==FALSE & ideal.setsize==TRUE){
    stop("Estimating the ideal set size is only valid for positive spatial autocorrelation")
  }

  # no bonferroni adjustment for 'pMI'
  if(objfn=="pMI" & bonferroni) bonferroni <- FALSE

  #####
  # Objective Function
  #####
  objfunc <- function(y,xe,n,W,objfn,boot.MI,alternative){
    resid <- y-xe%*%solve(crossprod(xe)) %*% crossprod(xe,y)
    if(objfn=="R2"){
      TSS <- sum((y - mean(y))^2)
      R2 <- 1-(sum(resid^2)/TSS)
      test <- -( 1-(1-R2)*(n-1)/(n-(ncol(xe))) ) # negative adjusted R-squared
    }
    if(objfn=="p"){
      est <- (solve(crossprod(xe)) %*% crossprod(xe,y))[ncol(xe),]
      se <- sqrt( (solve(t(xe)%*%xe)[ncol(xe),ncol(xe)]*sum(resid^2))/(nrow(W)-ncol(xe)) )
      test <- 2*pt(abs(est/se),df=(n-ncol(xe)),lower.tail=FALSE) # p-value
    }
    if(objfn=="MI"){
      test <- abs(MI.resid(resid=resid,x=xe,W=W,boot=boot.MI)$zI) # (absolute) standardized Moran's I
    }
    if(objfn=="pMI"){
      test <- -(MI.resid(resid=resid,x=xe,W=W,boot=boot.MI,alternative=alternative)$pI)
    }
    return(test)
  }

  #####
  # Eigenvectors and
  # Eigenvalues
  #####
  eigs <- getEVs(W,covars=MX)
  evecs <- eigs$vectors
  evals <- eigs$values

  # MI of eigenvectors
  evMI <- MI.ev(W=W,evals=evals)

  #####
  # Nonspatial
  # OLS Regression
  #####
  TSS <- sum((y - mean(y))^2)
  coefs_init <- solve(crossprod(x)) %*% crossprod(x,y)
  fitvals <- fittedval(x=x,params=coefs_init,model="linear")
  resid_init <- residfun(y=y,fitvals=fitvals,model="linear")$raw
  R2 <- 1-(sum(crossprod(resid_init))/TSS)
  adjR2_init <- 1-(1-R2)*(n-1)/(n-nx)
  zMI_init <- MI.resid(resid=resid_init,x=x,W=W,boot=boot.MI)$zI
  if(objfn=="MI") oldZMI <- abs(zMI_init)
  if(objfn=="R2") adjR2 <- -adjR2_init

  #####
  # Eigenvector Selection:
  # Candidate Set
  #####
  if(positive | zMI_init >= 0){
    if(ideal.setsize){
      # avoids problems of NaN if positive=TRUE but zMI < 0:
      csize <- candsetsize(npos=length(evals[evals > 1e-07])
                             ,zMI=ifelse(zMI_init<0,0,zMI_init))
      sel <- evals %in% evals[1:csize]
      dep <- "positive"
    } else{
      sel <- evMI/evMI[1] >= alpha
      dep <- "positive"
    }
  } else {
    sel <- evMI/evMI[n] >= alpha
    dep <- "negative"
  }

  # number of feasible eigenvectors in the candidate set
  ncandidates <- sum(sel)

  # Bonferroni adjustment (Griffith/ Chun 2014: 1490)
  if(objfn=="p" | objfn=="pMI"){
    if(bonferroni & ncandidates>0) sig <- sig/ncandidates
  } else {
    sig <- bonferroni <- NULL
  }

  if(objfn=="pMI") oldpMI <- -(MI.resid(resid=resid_init,x=x,W=W,boot=boot.MI
                                        ,alternative=ifelse(dep=="positive"
                                                            ,"greater","lower")
                                        )$pI)

  #####
  # Search Algorithm:
  # Stepwise Regression
  #####
  if(objfn=="all"){
    sel_id <- which(sel)
  } else {
    sel_id <- NULL
    selset <- which(sel)

    # start forward search
    for(i in which(sel)){
      if(objfn=="pMI") if(abs(oldpMI) > sig) break
      ref <- Inf
      sid <- NULL

      # select candidate eigenvector
      for(j in selset){
        xe <- cbind(x,evecs[,sel_id],evecs[,j])
        test <- objfunc(y=y,xe=xe,n=n,W=W,objfn=objfn,boot.MI=boot.MI
                        ,alternative=ifelse(dep=="positive","greater","lower"))
        if(test<ref){
          sid <- j
          ref <- test
        }
      }

      # stopping rules
      if(objfn=="R2"){
        if(ref < adjR2){
          adjR2 <- ref
          sel_id <- c(sel_id,sid)
        } else break
      }
      if(objfn=="p"){
        if(ref < sig){
          sel_id <- c(sel_id,sid)
        } else break
      }
      if(objfn=="MI"){
        if(ref < oldZMI){
          oldZMI <- ref
          sel_id <- c(sel_id,sid)
        } else break
        if(oldZMI < tol) break
      }
      if(objfn=="pMI"){
        if(abs(ref) > oldpMI){
          oldpMI <- abs(ref)
          sel_id <- c(sel_id,sid)
        }
        if(abs(ref) > sig)  break
      }

      # remove selected eigenvectors from search set
      selset <- selset[!(selset %in% sel_id)]
    } # end search
  }

  # covariates & eigenvectors
  xev <- cbind(x,evecs[,sel_id])

  # number of selected EVs
  count <- length(sel_id)

  #####
  # Filtered
  # OLS Results
  #####
  coefs <- solve(crossprod(xev)) %*% crossprod(xev,y)
  vcov <- (solve(crossprod(xev)) * sum((y-xev %*% coefs)^2))/(n-ncol(xev))
  se <- sqrt(diag(vcov))
  p.val <- 2*pt(abs(coefs/se),df=(n-ncol(xev)),lower.tail=FALSE)

  # fit & spatial autocorrelation
  fitvals <- fittedval(x=xev,params=coefs,model="linear")
  resid <- residfun(y=y,fitvals=fitvals,model="linear")$raw
  R2 <- 1-(sum(crossprod(resid))/TSS)
  adjR2 <- 1-(1-R2)*(n-1)/(n-ncol(xev))
  MI_filtered <- MI.resid(resid=resid,x=xev,W=W,boot=boot.MI
                          ,alternative=ifelse(dep=="positive","greater","lower"))
  MI_init <- MI.resid(resid=resid_init,x=x,W=W,boot=boot.MI
                      ,alternative=ifelse(dep=="positive","greater","lower"))

  #####
  # Output
  #####
  # OLS estimates (filtered)
  est <- cbind(coefs[1:nx],se[1:nx],p.val[1:nx])
  colnames(est) <- c("Estimate", "SE", "p-value")
  varcovar <- vcov[1:nx,1:nx]
  if(nx==1){
    rownames(est) <- names(varcovar) <- "(Intercept)"
  } else {
    if(!is.null(nams)){
      rownames(est) <- rownames(varcovar) <- colnames(varcovar) <- c("(Intercept)",nams)
    } else {
      rownames(est) <- c("(Intercept)",paste0("beta_",1:(nx-1)))
      rownames(varcovar) <- colnames(varcovar) <- rownames(est)
    }
  }

  # selected eigenvectors & eigenvalues
  if(count!=0){
    # selected eigenvectors
    selvecs <- as.matrix(evecs[,sel_id])
    colnames(selvecs) <- paste0("evec_",sel_id)
    # EV matrix
    gammas <- coefs[(nx+1):(nx+count)]
    gse <- se[(nx+1):(nx+count)]
    gp <- 2*pt(abs(gammas/gse),df=(n-ncol(xev)),lower.tail=FALSE)
    pR2 <- partialR2(y=y,x=xev[,1:nx],evecs=xev[,(nx+1):(nx+count)])
    vif <- vif.ev(x=xev[,1:nx],evecs=xev[,(nx+1):(nx+count)])
    EV <- cbind(gammas,gse,gp,pR2,vif,evMI[sel_id])
    colnames(EV) <- c("Estimate","SE","p-value","partialR2","VIF","MI")
    rownames(EV) <- paste0("ev_", sel_id)
    # construct spatial filter
    sf <- selvecs %*% gammas
    sfMI <- MI.sf(gamma=gammas, evMI=evMI[sel_id])
  } else selvecs <- EV <- sf <- sfMI <- NULL

  # Moran's I
  moran <- rbind(MI_init[,colnames(MI_init)!=""],MI_filtered[,colnames(MI_filtered)!=""])
  rownames(moran) <- c("Initial", "Filtered")
  colnames(moran) <- c("Observed","Expected","Variance","z","p-value")

  # model fit
  fit <- c(adjR2_init,adjR2)
  names(fit) <- c("Initial","Filtered")

  # residuals
  residuals <- cbind(resid_init,resid)
  colnames(residuals) <- c("Initial","Filtered")

  # output list
  out_list <- list(estimates=est
                   ,varcovar=varcovar
                   ,EV=EV
                   ,selvecs=selvecs
                   ,evMI=evMI
                   ,moran=moran
                   ,fit=fit
                   ,residuals=residuals
                   ,other=list(ncandidates=ncandidates
                               ,nev=count
                               ,sel_id=sel_id
                               ,sf=sf
                               ,sfMI=sfMI
                               ,model="linear"
                               ,dependence=dep
                               ,objfn=objfn
                               ,bonferroni=bonferroni
                               ,siglevel=sig
                   )
  )

  # define class
  class(out_list) <- "spfilter"

  # return
  return(out_list)
}
