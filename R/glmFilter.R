#' @name glmFilter
#'
#' @title Unsupervised Spatial Filtering with Eigenvectors in Generalized
#' Linear Regression Models
#'
#' @description This function implements the eigenvector-based semiparametric
#' spatial filtering approach in a generalized linear regression framework using MLE.
#' Eigenvectors are selected by an unsupervised steppwise regression
#' technique. Supported selection criteria are the minimization of residual
#' autocorrelation, maximization of model fit, and the statistical significance
#' of eigenvectors. Alternatively, all eigenvectors in the candidate set
#' can be included as well.
#'
#' @param y vector of regressands
#' @param x vector/ matrix of regressors (default=NULL)
#' @param W spatial connectivity matrix
#' @param objfn specifies the objective function to be used for eigenvector
#' selection. Possible criteria are: the maximization of the
#' adjusted R-squared ('AIC' or 'BIC'), minimization of residual autocorrelation ('MI'),
#' significance level of candidate eigenvectors ('p'), or
#' all eigenvectors in the candidate set ('all')
#' @param MX covariates used to construct the projection matrix (TRUE/ FALSE)
#' @param model a character string indicating the model to be estimated
#' @param optim.method a character specifying the optimization method
#' @param sig significance level to be used for eigenvector selection
#' if \code{objfn='p'}
#' @param bonferroni Bonferroni adjustment for the significance level
#' (TRUE/ FALSE)
#' @param positive restrict search to eigenvectors associated with positive
#' levels of spatial autocorrelation (TRUE/ FALSE)
#' @param min.reduction if \code{objfn} is either 'AIC' or 'BIC'. A value in the
#' interval [0,1) that determines the minimum reduction in AIC/ BIC a candidate
#' eigenvector need to achieve in order to be selected
#' @param boot.MI number of iterations used to estimate the variance of Moran's I
#' (default=100). Alternatively, if \code{boot=NULL}, analytical results will
#' be used
#' @param resid.type character string specifying the residual type to be used.
#' Options are 'raw', 'deviance', and 'pearson' (default)
#' @param alpha a value in (0,1] indicating the range of candidate eigenvectors
#' according to their associated level of spatial autocorrelation, see e.g.,
#' Griffith (2003)
#' @param tol if \code{objfn='MI'}, determines the amount of remaining residual
#' autocorrelation at which the eigenvector selection terminates
#' @param na.rm remove missing values in covariates (TRUE/ FALSE)
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
#' \item{\code{condnum}}{condition number to assess the degree of multicollinearity
#' among the eigenvectors induced by the link function, see e.g., Griffith/ Amrhein
#' (1997)}
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
#' \item{\code{resid.type}}{residual type used ('raw', 'deviance', or 'pearson')}
#' \item{\code{pseudoR2}}{McFadden's pseudo R-squared (filtered vs. unfiltered model)}
#' }
#' }
#'
#' @details If \emph{\strong{W}} is not symmetric, it gets symmetrized by
#' 0.5 * (\emph{\strong{W}} + \emph{\strong{W}}') the eigenfunction decomposition.
#'
#' If \code{MX=TRUE}, the function uses the covariates specified in the argument
#' \code{x} to construct the following projection matrix:
#'
#' \emph{\strong{M} = \strong{I} - \strong{X} (\strong{X}'\strong{X})^-1\strong{X}'}
#'
#' Eigenvectors from \emph{\strong{MWM}} using this specification of
#' \emph{\strong{M}} are not only mutually uncorrelated but also orthogonal
#' to the included regressors. Alternatively, if \code{MX=FALSE}, the projection
#' matrix becomes \emph{\strong{M} = \strong{I} - \strong{1}
#' (\strong{1}'\strong{1})^-1\strong{1}'}, where \emph{\strong{1}} is a vector of
#' ones. Griffith and Tiefelsdorf (2007) show how the choice of the appropriate
#' \emph{\strong{M}} depends on the underlying process that generates the spatial
#' dependence.
#'
#' @note If the condition number (\code{condnum}) suggests high levels of multicollinearity,
#' problematic eigenvectors can be manually removed from \code{selvecs} and
#' the model can be reestimated using the \code{glm} function. Moreover, if other models that
#' are currently not implemented here need to be estimated (e.g., quasi-binomials),
#' users can extract eigenvectors using the function \code{getEVs} and perform a
#' supervised eigenvector search using the \code{glm} function.
#'
#' @examples
#' data(fakedata)
#'
#' # poisson model
#' y_pois <- fakedataset$count
#' poisson <- glmFilter(y=y_pois,x=NULL,W=W,objfn="MI",positive=FALSE
#' ,model="poisson",boot.MI=100)
#' print(poisson)
#' summary(poisson,EV=FALSE)
#'
#' # probit model - summarize EVs
#' y_prob <- fakedataset$indicator
#' probit <- glmFilter(y=y_prob,x=NULL,W=W,objfn="MI",positive=FALSE
#' ,model="probit",boot.MI=100)
#' print(probit)
#' summary(probit,EV=TRUE)
#'
#' # logit model - AIC objective function
#' y_logit <- fakedataset$indicator
#' logit <- glmFilter(y=y_logit,x=NULL,W=W,objfn="AIC",positive=FALSE
#' ,model="logit",min.reduction=.05)
#' print(logit)
#' summary(logit,EV=FALSE)
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
#' Griffith, Daniel A. and Carl G. Amrhein (1997): Multivariate Statistical
#' Analysis for Geographers. Englewood Cliffs, Prentice Hall.
#'
#' @importFrom stats pnorm dpois optim pt
#'
#' @seealso \code{\link{lmFilter}}, \code{\link{getEVs}}, \code{\link{getMoran}}
#'
#' @export

glmFilter <- function(y,x=NULL,W,objfn="MI",MX=F,model,optim.method="BFGS"
                      ,sig=.05,bonferroni=TRUE,positive=TRUE,min.reduction=.05
                      ,boot.MI=100,resid.type="pearson",alpha=.25,tol=.1
                      ,na.rm=T){

  if(!is.null(x)) x <- as.matrix(x)
  if(!is.null(colnames(x))) nams <- colnames(x)

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

  # add intercept if not included in x
  if(is.null(x)) x <- as.matrix(rep(1,n))
  if(!all(x[,1]==1)) x <- cbind(1,x)
  nx <- ncol(x)

  #####
  # Input Checks
  #####
  if(anyNA(y) | anyNA(x) | anyNA(W)) stop("Missing values detected")
  if(alpha==0) alpha <- 1e-07
  if(alpha<1e-07 | alpha>1){
    stop("Invalid argument: 'alpha' must be in the interval (0,1]")
  }
  if(min.reduction<0 | min.reduction>=1){
    stop("Invalid argument: 'min.reduction' must be in the interval [0,1)")
  }
  if(qr(x)$rank!=ncol(x)) stop("Perfect multicollinearity in covariates detected")
  if(!any(class(W) %in% c("matrix","Matrix","data.frame"))){
    stop("W must be of class 'matrix' or 'data.frame'")
  }
  if(any(class(W)!="matrix")) W <- as.matrix(W)
  if(!(model %in% c("probit","logit","poisson"))){
    stop("'model' must be either 'probit', 'logit', or 'poisson'")
  }
  if(!(objfn %in% c("p","MI","AIC","BIC","all"))){
    stop("Invalid argument: objfn must be one of 'p', 'MI', 'AIC', 'BIC', or 'all'")
  }
  if(!(resid.type %in% c("raw","pearson","deviance"))){
    stop("Invalid argument: resid.type must be one of 'raw', 'pearson', or 'deviance'")
  }

  #####
  # Log-Likelihood Functions
  #####
  # loglik function
  loglik <- function(theta,y,x,model){
    if(model=="probit"){
      p <- pnorm(x%*%theta)
      ll <- sum(y*log(p) + (1-y)*log(1-p))
    }
    if(model=="logit"){
      p <- exp(x%*%theta)/(1+exp(x%*%theta))
      ll <- sum(y*log(p) + (1-y)*log(1-p))
    }
    if(model=="poisson"){
      mu <- exp(x%*%theta)
      #ll <- sum(y*log(mu) - mu)
      ll <- sum(dpois(y, lambda=mu,log=TRUE))
    }
    # return negative log-likelihood
    return(-ll)
  }

  # objective function to evaluate
  objfunc <- function(y,xe,n,W,objfn,model,optim.method,boot.MI,resid.type){
    inits <- rep(0,ncol(xe))
    o <- optim(par=inits,fn=loglik,x=xe,y=y,model=model
               ,method=optim.method,hessian=T)
    if(objfn=="p"){
      est <- o$par[ncol(xe)]
      se <- sqrt(diag(solve(o$hessian)))[ncol(xe)]
      test <- 2*pt(abs(est/se),df=(n-ncol(xe)),lower.tail=F) # p-value
    }
    if(objfn=="AIC"){
      test <- getICs(negloglik=o$value,n=n,df=ncol(xe))$AIC
    }
    if(objfn=="BIC"){
      test <- getICs(negloglik=o$value,n=n,df=ncol(xe))$BIC
    }
    if(objfn=="MI"){
      fitvals <- fittedval(x=xe,params=o$par,model=model)
      resid <- residfun(y=y,fitvals=fitvals,model=model)[,resid.type]
      test <- abs(getMoran(resid=resid,x=xe,W=W,boot=boot.MI)$zI)
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

  # MI of each eigenvector
  evMI <- MI.ev(W=W,evals=evals)

  #####
  # Nonspatial Regression
  #####
  inits <- rep(0,nx)
  opt <- optim(par=inits,fn=loglik,x=x,y=y,model=model
               ,method=optim.method,hessian=T)
  coefs_init <- opt$par
  se_init <- sqrt(diag(solve(opt$hessian)))
  ll_init <- opt$value
  ICs_init <- getICs(negloglik=ll_init, n=n, df=length(coefs_init))
  yhat_init <- fittedval(x=x,params=coefs_init,model=model)
  resid_init <- residfun(y=y,fitvals=yhat_init,model=model)[,resid.type]
  MI_init <- getMoran(resid=resid_init,x=x,W=W,boot=boot.MI)
  if(objfn=="MI") oldZMI <- abs(MI_init$zI)
  if(objfn %in% c("AIC","BIC")){
    IC <- ICs_init[,objfn]
    mindiff <- abs(IC*min.reduction)
  } else IC <- mindiff <- NULL

  #####
  # Eigenvector Selection:
  # Candidate Set
  #####
  if(positive | MI_init$zI>=0){
    sel <- evMI/evMI[1] >= alpha
    dep <- "positive"
  } else {
    sel <- evMI/evMI[n] >= alpha
    dep <- "negative"
  }

  # number of feasible eigenvectors
  ncandidates <- sum(sel)

  # Bonferroni adjustment
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

    # start forward selection
    for(i in which(sel)){
      ref <- Inf
      sid <- NULL

      # identify next test eigenvector
      for(j in selset){
        xe <- cbind(x,evecs[,sel_id],evecs[,j])
        test <- objfunc(y=y,xe=xe,n=n,W=W,objfn=objfn,model=model
                        ,optim.method=optim.method,boot.MI=boot.MI
                        ,resid.type=resid.type)
        if(test<ref){
          sid <- j
          ref <- test
        }
      }

      # stopping rules
      if(objfn %in% c("AIC","BIC")){
        if(ref < IC & abs(IC-ref) >= mindiff){
          IC <- ref
          sel_id <- c(sel_id,sid)
          mindiff <- abs(IC*min.reduction)
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

      # remove selected eigenvectors from candidate set
      #selset <- which(sel)[!(which(sel) %in% sel_id)] # re-evaluate excluded eigenvectors
      selset <- selset[!(selset %in% sel_id)] # evaluate just forward selected eigenvectors
    } # end selection
  }

  # number of selected EVs
  count <- length(sel_id)

  # filtered regression
  xev <- cbind(x,evecs[,sel_id])
  inits <- rep(0,ncol(xev))
  opt <- optim(par=inits,fn=loglik,x=xev,y=y,model=model
               ,method=optim.method,hessian=T)
  coefs_out <- opt$par
  se_out <- sqrt(diag(solve(opt$hessian)))
  p.val <- 2*pt(abs(coefs_out/se_out),df=(n-ncol(xev)),lower.tail=F)
  ll_out <- opt$value
  ICs_out <- getICs(negloglik=ll_out, n=n, df=length(coefs_out))
  yhat_out <- fittedval(x=xev,params=coefs_out,model=model)
  resid_out <- residfun(y=y,fitvals=yhat_out,model=model)[,resid.type]
  MI_out <- getMoran(resid=resid_out,x=xev,W=W,boot=boot.MI)

  #####
  # Output
  #####
  # Estimates
  est <- cbind(coefs_out[1:nx],se_out[1:nx],p.val[1:nx])
  colnames(est) <- c("Estimate", "SE", "p-value")
  varcovar <- solve(opt$hessian)[1:nx,1:nx]
  if(nx==1){
    rownames(est) <- names(varcovar) <- "(Intercept)"
  } else {
    if(!is.null(colnames(x))){
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
    condnum <- conditionNumber(evecs=selvecs,round=8) # multicollinearity: condnum > 1
    # EV matrix
    gammas <- coefs_out[(nx+1):(nx+count)]
    gse <- se_out[(nx+1):(nx+count)]
    gp <- 2*pt(abs(gammas/gse),df=(n-length(coefs_out)),lower.tail=F)
    EV <- cbind(gammas,gse,gp,evMI[sel_id])
    colnames(EV) <- c("Estimate","SE","p-value","MI")
    rownames(EV) <- paste0("ev_", sel_id)
    # generate spatial filter
    sf <- selvecs %*% gammas
    sfMI <- MI.sf(gamma=gammas, evMI=evMI[sel_id])
  } else selvecs <- EV <- sf <- sfMI <- condnum <- NULL

  # Moran's I
  moran <- rbind(MI_init[,colnames(MI_init)!=""],MI_out[,colnames(MI_out)!=""])
  rownames(moran) <- c("Initial", "Filtered")
  colnames(moran) <- c("Observed","Expected","Variance","z","p-value")

  # model fit
  fit <- rbind(c(-ll_init,ICs_init),c(-ll_out,ICs_out))
  rownames(fit) <- c("Initial", "Filtered")
  colnames(fit) <- c("logL","AIC","BIC")

  # McFadden's pseudo-R2
  pRsqr <- pseudoR2(negloglik_n=ll_init,negloglik_f=ll_out,nev=count)

  # residuals
  residuals <- cbind(resid_init,resid_out)
  colnames(residuals) <- c("Initial","Filtered")

  # define output
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
                               ,condnum=condnum
                               ,sel_id=sel_id
                               ,sf=sf
                               ,sfMI=sfMI
                               ,model=model
                               ,MX=MX
                               ,dependence=dep
                               ,objfn=objfn
                               ,bonferroni=bonferroni
                               ,siglevel=sig
                               ,resid.type=resid.type
                               ,pseudoR2=pRsqr
                   )
  )

  # define class
  class(out_list) <- "spfilter"

  # return
  return(out_list)
}
