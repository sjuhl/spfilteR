#' @name lmFilter
#' @title Eigenvector-Based Spatial Filtering in Linear Regression Models
#'
#' @description This function implements the eigenvector-based semiparametric
#' spatial filtering approach in a linear regression framework using OLS.
#' Eigenvectors are selected by an unsupervised steppwise regression
#' technique. Supported selection criteria are the minimization of residual
#' autocorrelation, maximization of model fit, and the statistical significance
#' of eigenvectors. Alternatively, all eigenvectors in the candidate set
#' can be included as well.
#'
#' @param y vector of regressands
#' @param x vector/ matrix of regressors (default = NULL)
#' @param W spatial connectivity matrix
#' @param objfn specifies the objective function to be used for eigenvector
#' selection. Possible criteria are: the maximization of the
#' adjusted R-squared ('R2'), minimization of residual autocorrelation ('MI'),
#' significance level of candidate eigenvectors ('p'), or
#' all eigenvectors in the candidate set ('all')
#' @param MX covariates used to construct the projection matrix (TRUE/ FALSE)
#' @param sig significance level to be used for eigenvector selection
#' if \code{objfn='p'}
#' @param bonferroni Bonferroni adjustment for the significance level
#' (TRUE/ FALSE)
#' @param positive restrict search to eigenvectors associated with positive
#' levels of spatial autocorrelation (TRUE/ FALSE)
#' @param ideal.setsize if \code{positive=TRUE}, uses the formula proposed in
#' Chun et al. (2016) to determine the ideal size of the candidate set
#' (TRUE/ FALSE)
#' @param tol if \code{objfn='MI'}, determines the amount of remaining residual
#' autocorrelation at which the eigenvector selection terminates
#' @param na.rm removes missing values in covariates (TRUE/ FALSE)
#'
#' @return An object of class \code{spfilter} containing the following
#' information:
#' \tabular{lcl}{
#' \code{Estimates}\tab \tab summary statistics of the parameter estimates\cr
#' \code{varcovar}\tab\tab estimated variance-covariance matrix\cr
#' \code{EV}\tab\tab a matrix with summary statistics of selected eigenvectors\cr
#' \code{selvecs}\tab\tab matrix of selected eigenvectors\cr
#' \code{evMI}\tab\tab Moran coefficient of all eigenvectors\cr
#' \code{moran}\tab\tab residual autocorrelation for the initial and the
#' filtered model\cr
#' \code{fit}\tab\tab adjusted R-squared of the initial and the filtered model\cr
#' \code{residuals}\tab\tab initial and filtered model residuals\cr
#' \code{other}\tab\tab a list providing supplementary information:\cr
#' \tab\tab \describe{
#' \item{\code{ncandidates}}{number of candidate eigenvectors considered}
#' \item{\code{nev}}{number of selected eigenvectors}
#' \item{\code{sel_id}}{ID of selected eigenvectors}
#' \item{\code{sf}}{vector representing the spatial filter}
#' \item{\code{sfMI}}{Moran coefficient of the spatial filter}
#' \item{\code{model}}{class of the regression model}
#' \item{\code{MX}}{TRUE/ FALSE: include covariates in the projection matrix
#' \emph{\strong{M}}}
#' \item{\code{dependence}}{filtered for positive or negative spatial dependence}
#' \item{\code{objfn}}{selection criteria specified in the objective function of
#' the stepwise regression procedure}
#' \item{\code{bonferroni}}{TRUE/ FALSE: Bonferroni-adjusted significance level
#' (if \code{objfn='p'})}
#' \item{\code{siglevel}}{if \code{objfn='p'}: actual (unadjusted/ adjusted)
#' significance level}
#' }
#' }
#'
#' @details If \emph{\strong{W}} is not symmetric, it gets symmetrized by
#' \eqn{0.5 * (W + t(W))} the eigenfunction decomposition.
#'
#' @examples
#' data(fakedata)
#' y <- fakedataset$x1
#' X <- cbind(fakedataset$count,fakedataset$x2,fakedataset$x5)
#'
#' res <- lmFilter(y=y,x=X,W=W,objfn='R2')
#' print(res)
#' summary(res,EV=F)
#' plot(res)
#'
#' E <- res$selvecs
#' (lm(y~X+E))
#'
#' @references Tiefelsdorf, Michael and Daniel A. Griffith (2007):
#' Semiparametric filtering of spatial autocorrelation: the eigenvector
#' approach. Environment and Planning A: Economy and Space, 39 (5):
#' pp. 1193 - 1221.
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
#' @seealso \code{\link{getEVs}}, \code{\link{getMoran}}
#'
#' @export

lmFilter <- function(y,x=NULL,W,objfn="MI",MX=FALSE,sig=.05
                     ,bonferroni=TRUE,positive=TRUE,ideal.setsize=FALSE
                     ,alpha=.25,tol=.1,na.rm=TRUE){

  if(!is.null(x)) x <- as.matrix(x)

  # missing values
  if(na.rm){
    if(!is.null(x)){
      miss <- apply(cbind(y,x),1,anyNA)
      x <- as.matrix(x[!miss,])
    } else miss <- is.na(y)
    y <- y[!miss]
    W <- W[!miss,!miss]
  }

  # number of observations
  n <- length(y)
  if(is.null(x)) x <- as.matrix(rep(1,n))

  # add intercept if not included in x
  if(!all(x[,1]==1)) x <- cbind(1,x)
  nx <- ncol(x)

  #####
  # Input Checks
  #####
  if(anyNA(y) | anyNA(x) | anyNA(W)) stop("Missing values detected")
  if(alpha==0) alpha <- 1e-07
  if(alpha<1e-07 | alpha>1){
    stop("Invalid argument: 'alpha' must satisfy the restriction: 0 < alpha <= 1.")
  }
  if(qr(x)$rank!=ncol(x)) stop("Perfect multicollinearity in covariates detected")
  if(!any(class(W) %in% c("matrix","Matrix","data.frame"))){
    stop("W must be of class 'matrix' or 'data.frame'")
  }
  if(any(class(W)!="matrix")) W <- as.matrix(W)
  if(!(objfn %in% c("R2","p","MI","all"))){
    stop("Invalid argument: objfn must be one of 'R2', 'p', 'MI', or'all'")
  }
  if(positive==F & ideal.setsize==T){
    stop("Estimating the ideal set size is only valid for positive spatial autocorrelation")
  }

  #####
  # Objective Function
  #####
  objfunc <- function(y,xe,n,W,objfn){
    resid <- y-xe%*%solve(crossprod(xe)) %*% crossprod(xe,y)
    if(objfn=="R2"){
      R2 <- 1-(sum(resid^2)/TSS)
      test <- -( 1-(1-R2)*(n-1)/(n-(ncol(xe))) ) # (negative) adjusted R-squared
    }
    if(objfn=="p"){
      est <- (solve(crossprod(xe)) %*% crossprod(xe,y))[ncol(xe),]
      se <- sqrt( (solve(t(xe)%*%xe)[ncol(xe),ncol(xe)]*sum(resid^2))/(nrow(W)-ncol(xe)) )
      test <- 2*pt(abs(est/se),df=(n-ncol(xe)),lower.tail=F) # p-value
    }
    if(objfn=="MI"){
      test <- abs(getMoran(y=y,x=xe,W=W)$zI) # (absolute) standardized Moran's I
    }
    return(test)
  }

  #####
  # Eigenvectors and
  # Eigenvalues
  #####
  if(MX){
    eigs <- getEVs(W,covars=x)
  } else eigs <- getEVs(W)

  evecs <- eigs$vectors
  evals <- eigs$values

  # MI of eigenvectors
  evMI <- MI.ev(W=W,evals=evals)

  #####
  # Nonspatial
  # OLS Regression
  #####
  TSS <- sum((y - mean(y))^2)
  resid_init <- y - x %*% solve(crossprod(x)) %*% crossprod(x,y)
  R2 <- 1-(sum(crossprod(resid_init))/TSS)
  adjR2 <- adjR2_init <- 1-(1-R2)*(n-1)/(n-nx)
  MI_init <- getMoran(y=y,x=x,W=W)
  if(objfn=="MI") oldZMI <- abs(MI_init$zI)

  #####
  # Eigenvector Selection:
  # Candidate Set
  #####
  if(positive | MI_init$zI >= 0){
    if(ideal.setsize){
      # avoids problems of NaN if positive=T but zMI < 0:
      if(MI_init$zI < 0){
        csize <- candsetsize(npos=length(evals[evals > 1e-07]),zMI=0)
      } else {
        csize <- candsetsize(npos=length(evals[evals > 1e-07]),zMI=MI_init$zI)
      }
      if(csize==0) sel <- rep(F,length(evals))
      else sel <- evals %in% evals[1:csize]
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
  if(objfn=="p"){
    if(bonferroni & ncandidates>0) sig <- sig/ncandidates
  } else {
    sig <- bonferroni <- NULL
  }

  #####
  # Search Algorithm:
  # Stepwise Regression
  #####
  if(objfn=="all"){
    sel_id <- which(sel)
  } else {
    sel_id <- NULL # list of selected eigenvectors
    selset <- which(sel)

    # start outer loop
    for(i in which(sel)){
      oldtest <- Inf
      sid <- NULL

      # start inner loop
      for(j in selset){
        xe <- cbind(x,evecs[,sel_id],evecs[,j])
        test <- objfunc(y=y,xe=xe,n=n,W=W,objfn=objfn)
        if(test<oldtest){
          sid <- j
          oldtest <- test
        }
      } # end inner loop

      # stopping rules
      if(objfn=="R2"){
        if(oldtest < adjR2){
          adjR2 <- oldtest
          sel_id <- c(sel_id,sid)
        } else break
      }
      if(objfn=="p"){
        if(oldtest < sig){
          sel_id <- c(sel_id,sid)
        } else break
      }
      if(objfn=="MI"){
        if(oldtest < oldZMI){
          oldZMI <- oldtest
          sel_id <- c(sel_id,sid)
        } else break
        if(oldZMI < tol) break
      }

      # remove selected eigenvectors from search set
      selset <- selset[!(selset %in% sel_id)]

    } # end outer loop
  }

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
  p.val <- 2*pt(abs(coefs/se),df=(n-ncol(xev)),lower.tail=F)

  # fit & spatial autocorrelation
  resid <- y - xev %*% coefs
  R2 <- 1-(sum(crossprod(resid))/TSS)
  adjR2 <- 1-(1-R2)*(n-1)/(n-ncol(xev))
  MI_filtered <- getMoran(y=y,x=xev,W=W)

  #####
  # Output
  #####
  # OLS estimates (filtered)
  est <- cbind(coefs[1:nx],se[1:nx],p.val[1:nx])
  colnames(est) <- c("Estimate", "SE", "p-value")
  if(nx==1){
    rownames(est) <- "(Intercept)"
  } else rownames(est) <- c("(Intercept)",paste0("beta_",1:(nx-1)))
  varcovar <- vcov[1:nx,1:nx]
  rownames(varcovar) <- colnames(varcovar) <- rownames(est)

  # selected eigenvectors & eigenvalues
  if(count!=0){
    # selected eigenvectors
    selvecs <- as.matrix(evecs[,sel_id])
    colnames(selvecs) <- paste0("evec_",sel_id)
    # EV matrix
    gammas <- coefs[(nx+1):(nx+count)]
    gse <- se[(nx+1):(nx+count)]
    gp <- 2*pt(abs(gammas/gse),df=(n-ncol(xev)),lower.tail=F)
    pR2 <- partialR2(y=y,x=xev[,1:nx],evec=xev[,(nx+1):(nx+count)])
    vif <- vif.ev(x=xev[,1:nx],evec=xev[,(nx+1):(nx+count)])
    EV <- cbind(gammas,gse,gp,pR2,vif,evMI[sel_id])
    colnames(EV) <- c("Estimate","SE","p-value","partialR2","VIF","MI")
    rownames(EV) <- paste0("ev_", sel_id)
    # construct spatial filter
    sf <- selvecs %*% gammas
    sfMI <- MI.sf(gamma=gammas, evMI=evMI[sel_id])
  } else selvecs <- EV <- sf <- sfMI <- NULL

  # Moran's I
  if(dep=="positive"){
    MI_init$pI <- pnorm(MI_init$zI,lower.tail=F)
    MI_filtered$pI <- pnorm(MI_filtered$zI,lower.tail=F)
  } else{
    MI_init$pI <- pnorm(MI_init$zI,lower.tail=T)
    MI_filtered$pI <- pnorm(MI_filtered$zI,lower.tail=T)
  }
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
                               ,MX=MX
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
