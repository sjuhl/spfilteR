#' @name glmFilter
#'
#' @title Unsupervised Spatial Filtering with Eigenvectors in Generalized
#' Linear Regression Models
#'
#' @description This function implements the eigenvector-based semiparametric
#' spatial filtering approach in a generalized linear regression framework using
#' maximum likelihood estimation (MLE). Eigenvectors are selected by an unsupervised
#' stepwise regression technique. Supported selection criteria are the minimization of
#' residual autocorrelation, maximization of model fit, significance of residual
#' autocorrelation, and the statistical significance of eigenvectors. Alternatively,
#' all eigenvectors in the candidate set can be included as well.
#'
#' @param y response variable
#' @param x vector/ matrix of regressors (default = NULL)
#' @param W spatial connectivity matrix
#' @param objfn the objective function to be used for eigenvector
#' selection. Possible criteria are: the maximization of model fit
#' ('AIC' or 'BIC'), minimization of residual autocorrelation ('MI'),
#' significance level of candidate eigenvectors ('p'), significance of residual spatial
#' autocorrelation ('pMI'), or all eigenvectors in the candidate set ('all')
#' @param MX covariates used to construct the projection matrix (default = NULL) - see
#' Details
#' @param model a character string indicating the type of model to be estimated.
#' Currently, 'probit', 'logit', and 'poisson' are valid inputs
#' @param optim.method a character specifying the optimization method used by
#' the \code{optim} function
#' @param sig significance level to be used for eigenvector selection
#' if \code{objfn = 'p'} or \code{objfn = 'pMI'}
#' @param bonferroni Bonferroni adjustment for the significance level
#' (TRUE/ FALSE) if \code{objfn = 'p'}. Set to FALSE if \code{objfn = 'pMI'} -
#' see Details
#' @param positive restrict search to eigenvectors associated with positive
#' levels of spatial autocorrelation (TRUE/ FALSE)
#' @param ideal.setsize if \code{positive = TRUE}, uses the formula proposed by
#' Chun et al. (2016) to determine the ideal size of the candidate set
#' (TRUE/ FALSE)
#' @param min.reduction if \code{objfn} is either 'AIC' or 'BIC'. A value in the
#' interval [0,1) that determines the minimum reduction in AIC/ BIC (relative to the
#' current AIC/ BIC) a candidate eigenvector need to achieve in order to be selected
#' @param boot.MI number of iterations used to estimate the variance of Moran's I
#' (default = 100). Alternatively, if \code{boot.MI = NULL}, analytical results will
#' be used
#' @param resid.type character string specifying the residual type to be used.
#' Options are 'raw', 'deviance', and 'pearson' (default)
#' @param alpha a value in (0,1] indicating the range of candidate eigenvectors
#' according to their associated level of spatial autocorrelation, see e.g.,
#' Griffith (2003)
#' @param tol if \code{objfn = 'MI'}, determines the amount of remaining residual
#' autocorrelation at which the eigenvector selection terminates
#' @param na.rm remove observations with missing values (TRUE/ FALSE)
#'
#' @return An object of class \code{spfilter} containing the following
#' information:
#' \describe{
#' \item{\code{estimates}}{summary statistics of the parameter estimates}
#' \item{\code{varcovar}}{estimated variance-covariance matrix}
#' \item{\code{EV}}{a matrix containing the summary statistics of selected eigenvectors}
#' \item{\code{selvecs}}{vector/ matrix of selected eigenvectors}
#' \item{\code{evMI}}{Moran coefficient of all eigenvectors}
#' \item{\code{moran}}{residual autocorrelation in the initial and the
#' filtered model}
#' \item{\code{fit}}{adjusted R-squared of the initial and the filtered model}
#' \item{\code{residuals}}{initial and filtered model residuals}
#' \item{\code{other}}{a list providing supplementary information:
#' \describe{
#' \item{\code{ncandidates}}{number of candidate eigenvectors considered}
#' \item{\code{nev}}{number of selected eigenvectors}
#' \item{\code{condnum}}{condition number to assess the degree of multicollinearity
#' among the eigenvectors induced by the link function, see e.g., Griffith/ Amrhein
#' (1997)}
#' \item{\code{sel_id}}{ID of selected eigenvectors}
#' \item{\code{sf}}{vector representing the spatial filter}
#' \item{\code{sfMI}}{Moran coefficient of the spatial filter}
#' \item{\code{model}}{type of the regression model}
#' \item{\code{dependence}}{filtered for positive or negative spatial dependence}
#' \item{\code{objfn}}{selection criterion specified in the objective function of
#' the stepwise regression procedure}
#' \item{\code{bonferroni}}{TRUE/ FALSE: Bonferroni-adjusted significance level
#' (if \code{objfn='p'})}
#' \item{\code{siglevel}}{if \code{objfn = 'p'} or \code{objfn = 'pMI'}: actual
#' (unadjusted/ adjusted) significance level}
#' \item{\code{resid.type}}{residual type ('raw', 'deviance', or 'pearson')}
#' \item{\code{pseudoR2}}{McFadden's pseudo R-squared (filtered vs. unfiltered model)}
#' }
#' }
#' }
#'
#' @details If \emph{\strong{W}} is not symmetric, it gets symmetrized by
#' 1/2 * (\emph{\strong{W}} + \emph{\strong{W}}') before the decomposition.
#'
#' If covariates are supplied to \code{MX}, the function uses these regressors
#' to construct the following projection matrix:
#'
#' \emph{\strong{M} = \strong{I} - \strong{X} (\strong{X}'\strong{X})^-1\strong{X}'}
#'
#' Eigenvectors from \emph{\strong{MWM}} using this specification of
#' \emph{\strong{M}} are not only mutually uncorrelated but also orthogonal
#' to the regressors specified in \code{MX}. Alternatively, if \code{MX = NULL}, the
#' projection matrix becomes \emph{\strong{M} = \strong{I} - \strong{11}'/\emph{n}},
#' where \emph{\strong{1}} is a vector of ones and \emph{n} represents the number of
#' observations. Griffith and Tiefelsdorf (2007) show how the choice of the appropriate
#' \emph{\strong{M}} depends on the underlying process that generates the spatial
#' dependence.
#'
#' The Bonferroni correction is only possible if eigenvector selection is based on
#' the significance level of the eigenvectors (\code{objfn = 'p'}). It is set to
#' FALSE if eigenvectors are added to the model until the residuals exhibit no
#' significant level of spatial autocorrelation (\code{objfn = 'pMI'}).
#'
#' @note If the condition number (\code{condnum}) suggests high levels of
#' multicollinearity, eigenvectors can be sequentially removed from \code{selvecs}
#' and the model can be re-estimated using the \code{glm} function in order to
#' identify and manually remove the problematic eigenvectors. Moreover, if other
#' models that are currently not implemented here need to be estimated
#' (e.g., quasi-binomial models), users can extract eigenvectors using the function
#' \code{getEVs} and perform a supervised eigenvector search using the \code{glm}
#' function.
#'
#' In contrast to eigenvector-based spatial filtering in linear regression models,
#' Chun (2014) notes that only a limited number of studies address the problem
#' of measuring spatial autocorrelation in generalized linear model residuals.
#' Consequently, eigenvector selection may be based on an objective function that
#' maximizes model fit rather than minimizes residual spatial autocorrelation.
#'
#' @examples
#' data(fakedata)
#'
#' # poisson model
#' y_pois <- fakedataset$count
#' poisson <- glmFilter(y = y_pois, x = NULL, W = W, objfn = "MI", positive = FALSE,
#' model = "poisson", boot.MI = 100)
#' print(poisson)
#' summary(poisson, EV = FALSE)
#'
#' # probit model - summarize EVs
#' y_prob <- fakedataset$indicator
#' probit <- glmFilter(y = y_prob, x = NULL, W = W, objfn = "p", positive = FALSE,
#' model = "probit", boot.MI = 100)
#' print(probit)
#' summary(probit, EV = TRUE)
#'
#' # logit model - AIC objective function
#' y_logit <- fakedataset$indicator
#' logit <- glmFilter(y = y_logit, x = NULL, W = W, objfn = "AIC", positive = FALSE,
#' model = "logit", min.reduction = .05)
#' print(logit)
#' summary(logit, EV = FALSE)
#'
#' @references Chun, Yongwan (2014): Analyzing Space-Time Crime Incidents Using
#' Eigenvector Spatial Filtering: An Application to Vehicle Burglary.
#' Geographical Analysis 46 (2): pp. 165 - 184.
#'
#' Tiefelsdorf, Michael and Daniel A. Griffith (2007):
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
#' @importFrom stats pnorm dpois optim pt sd
#'
#' @seealso \code{\link{lmFilter}}, \code{\link{getEVs}}, \code{\link{MI.resid}},
#' \code{\link[stats]{optim}}
#'
#' @export

glmFilter <- function(y, x = NULL, W, objfn = "AIC", MX = NULL, model, optim.method = "BFGS",
                      sig = .05, bonferroni = TRUE, positive = TRUE, ideal.setsize = FALSE,
                      min.reduction = .05, boot.MI = 100, resid.type = "pearson",
                      alpha = .25, tol = .1, na.rm = TRUE) {

  if (!is.null(MX)) {
    MX <- as.matrix(MX)
  }
  if (!is.null(x)) {
    x <- as.matrix(x)
  }
  if (!is.null(colnames(x))) {
    nams <- colnames(x)
  } else {
    nams <- NULL
  }

  # missing values
  if (na.rm) {
    if (!is.null(x)) {
      miss <- apply(cbind(y, x), 1, anyNA)
      x <- as.matrix(x[!miss,])
    } else {
      miss <- is.na(y)
    }
    y <- y[!miss]
    W <- W[!miss, !miss]
    if (!is.null(MX)) {
      MX[!miss,]
    }
  }

  # number of observations
  n <- length(y)

  # add intercept if not included in x
  if (is.null(x)) {
    x <- as.matrix(rep(1, n))
  }
  if (!isTRUE(all.equal(x[, 1], rep(1, n)))) {
    x <- cbind(1, x)
  }
  if (!is.null(MX) && any(apply(MX, 2, sd) == 0)) {
    MX <- as.matrix(MX[, apply(MX, 2, sd) != 0])
  }
  nx <- ncol(x)

  #####
  # Input Checks
  #####
  if (anyNA(y) | anyNA(x) | anyNA(W)) {
    stop("Missing values detected")
  }
  if (alpha == 0) {
    alpha <- 1e-07
  }
  if (alpha < 1e-07 | alpha > 1) {
    stop("Invalid argument: 'alpha' must be in the interval (0,1]")
  }
  if (min.reduction < 0 | min.reduction >= 1) {
    stop("Invalid argument: 'min.reduction' must be in the interval [0,1)")
  }
  if (qr(x)$rank != ncol(x)) {
    stop("Perfect multicollinearity in covariates detected")
  }
  if (!any(class(W) %in% c("matrix", "Matrix", "data.frame"))) {
    stop("W must be of class 'matrix' or 'data.frame'")
  }
  if (any(class(W) != "matrix")) {
    W <- as.matrix(W)
  }
  if (!(model %in% c("probit", "logit", "poisson"))) {
    stop("'model' must be either 'probit', 'logit', or 'poisson'")
  }
  if (!(objfn %in% c("p", "MI", "pMI", "AIC", "BIC", "all"))) {
    stop("Invalid argument: objfn must be one of 'p', 'MI','pMI, 'AIC', 'BIC', or 'all'")
  }
  if (!(resid.type %in% c("raw", "pearson", "deviance"))) {
    stop("Invalid argument: resid.type must be one of 'raw', 'pearson', or 'deviance'")
  }
  if (positive == FALSE & ideal.setsize == TRUE) {
    stop("Estimating the ideal set size is only valid for positive spatial autocorrelation")
  }

  # no bonferroni adjustment for 'pMI'
  if (objfn == "pMI" & bonferroni) {
    bonferroni <- FALSE
  }

  #####
  # Log-Likelihood Functions
  #####
  # loglik function
  loglik <- function(theta, y, x, model) {
    if (model == "probit") {
      p <- pnorm(x %*% theta)
      ll <- sum(y * log(p) + (1 - y) * log(1 - p))
    } else if (model == "logit") {
      p <- exp(x %*% theta) / (1 + exp(x %*% theta))
      ll <- sum(y * log(p) + (1 - y) * log(1 - p))
    } else if (model == "poisson") {
      mu <- exp(x %*% theta)
      #ll <- sum(y*log(mu) - mu)
      ll <- sum(dpois(y, lambda = mu, log = TRUE))
    }
    # return negative log-likelihood
    return(-ll)
  }

  # objective function to evaluate
  objfunc <- function(y, xe, n, W, objfn, model, optim.method, boot.MI,
                      resid.type, alternative) {
    inits <- rep(0, ncol(xe))
    o <- optim(par = inits, fn = loglik, x = xe, y = y, model = model,
               method = optim.method, hessian = TRUE)
    if (objfn == "p") {
      est <- o$par[ncol(xe)]
      se <- sqrt(diag(solve(o$hessian)))[ncol(xe)]
      test <- 2 * pt(abs(est / se), df = (n - ncol(xe)), lower.tail = FALSE) # p-value
    } else if (objfn == "AIC") {
      test <- getICs(negloglik = o$value, n = n, df = ncol(xe))$AIC
    } else if (objfn == "BIC") {
      test <- getICs(negloglik = o$value, n = n, df = ncol(xe))$BIC
    } else if (objfn == "MI") {
      fitvals <- fittedval(x = xe, params = o$par, model = model)
      resid <- residfun(y = y, fitvals = fitvals, model = model)[, resid.type]
      test <- abs(MI.resid(resid = resid, x = xe, W = W, boot = boot.MI)$zI)
    } else if (objfn == "pMI") {
      fitvals <- fittedval(x = xe, params = o$par, model = model)
      resid <- residfun(y = y, fitvals = fitvals, model = model)[, resid.type]
      test <- -(MI.resid(resid = resid, x = xe, W = W, boot = boot.MI,
                         alternative = alternative)$pI)
    }
    return(test)
  }

  #####
  # Eigenvectors and
  # Eigenvalues
  #####
  eigs <- getEVs(W, covars = MX)
  evecs <- eigs$vectors
  evals <- eigs$values

  # MI of each eigenvector
  evMI <- MI.ev(W = W, evals = evals)

  #####
  # Nonspatial Regression
  #####
  inits <- rep(0, nx)
  opt <- optim(par = inits, fn = loglik, x = x, y = y, model = model,
               method = optim.method, hessian = TRUE)
  coefs_init <- opt$par
  se_init <- sqrt(diag(solve(opt$hessian)))
  ll_init <- opt$value
  ICs_init <- getICs(negloglik = ll_init, n = n, df = length(coefs_init))
  yhat_init <- fittedval(x = x, params = coefs_init, model = model)
  resid_init <- residfun(y = y, fitvals = yhat_init, model = model)[, resid.type]
  zMI_init <- MI.resid(resid = resid_init, x = x, W = W, boot = boot.MI)$zI
  if (objfn == "MI") {
    oldZMI <- abs(zMI_init)
  } else if (objfn %in% c("AIC", "BIC")) {
    IC <- ICs_init[, objfn]
    mindiff <- abs(IC * min.reduction)
  } else {
    IC <- mindiff <- NULL
  }

  #####
  # Eigenvector Selection:
  # Candidate Set
  #####
  if (positive | zMI_init >= 0) {
    if (ideal.setsize) {
      # avoids problems of NaN if positive=TRUE but zMI < 0:
      csize <- candsetsize(npos = length(evals[evals > 1e-07]),
                           zMI=ifelse(zMI_init < 0, 0, zMI_init))
      sel <- evals %in% evals[1:csize]
      dep <- "positive"
    } else {
      sel <- evMI / evMI[1] >= alpha
      dep <- "positive"
    }
  } else {
    sel <- evMI / evMI[n] >= alpha
    dep <- "negative"
  }

  # number of feasible eigenvectors
  ncandidates <- sum(sel)

  # Bonferroni adjustment
  if (objfn == "p" | objfn == "pMI") {
    if (bonferroni & ncandidates > 0) {
      sig <- sig / ncandidates
    }
  } else {
    sig <- bonferroni <- NULL
  }

  if (objfn == "pMI") {
    oldpMI <- -(MI.resid(resid = resid_init, x = x, W = W, boot = boot.MI,
                         alternative = ifelse(dep == "positive", "greater", "lower")
                         )$pI)
  }

  #####
  # Search Algorithm:
  # Stepwise Regression
  #####
  if (objfn == "all") {
    sel_id <- which(sel)
  } else {
    sel_id <- NULL
    selset <- which(sel)

    # start forward selection
    for (i in which(sel)) {
      if (objfn == "pMI") {
        if (abs(oldpMI) > sig) {
          break
        }
      }
      ref <- Inf
      sid <- NULL

      # identify next test eigenvector
      for (j in selset) {
        xe <- cbind(x, evecs[, sel_id], evecs[, j])
        test <- objfunc(y = y, xe = xe, n = n, W = W, objfn = objfn, model = model,
                        optim.method = optim.method, boot.MI = boot.MI,
                        resid.type = resid.type, alternative = ifelse(dep == "positive",
                                                                      "greater", "lower"))
        if (test < ref) {
          sid <- j
          ref <- test
        }
      }

      # stopping rules
      if (objfn %in% c("AIC", "BIC")) {
        if (ref < IC & abs(IC - ref) >= mindiff) {
          IC <- ref
          sel_id <- c(sel_id, sid)
          mindiff <- abs(IC * min.reduction)
        } else {
          break
        }
      } else if (objfn == "p") {
        if (ref < sig) {
          sel_id <- c(sel_id, sid)
        } else {
          break
        }
      } else if (objfn == "MI") {
        if(ref < oldZMI){
          oldZMI <- ref
          sel_id <- c(sel_id, sid)
        } else {
          break
        }
        if (oldZMI < tol) {
          break
        }
      } else if (objfn == "pMI") {
        if (abs(ref) > oldpMI) {
          oldpMI <- abs(ref)
          sel_id <- c(sel_id, sid)
        }
        if (abs(ref) > sig) {
          break
        }
      }

      # remove selected eigenvectors from candidate set
      selset <- selset[!(selset %in% sel_id)]
    } # end selection
  }

  # number of selected EVs
  count <- length(sel_id)

  # filtered regression
  xev <- cbind(x, evecs[, sel_id])
  inits <- rep(0, ncol(xev))
  opt <- optim(par = inits, fn = loglik, x = xev, y = y, model = model,
               method = optim.method, hessian = TRUE)
  coefs_out <- opt$par
  se_out <- sqrt(diag(solve(opt$hessian)))
  p.val <- 2 * pt(abs(coefs_out / se_out), df = (n - ncol(xev)), lower.tail = FALSE)
  ll_out <- opt$value
  ICs_out <- getICs(negloglik = ll_out, n = n, df = length(coefs_out))
  yhat_out <- fittedval(x = xev, params = coefs_out, model = model)
  resid_out <- residfun(y = y, fitvals = yhat_out, model = model)[, resid.type]
  MI_out <- MI.resid(resid = resid_out, x = xev, W = W, boot = boot.MI,
                     alternative = ifelse(dep == "positive", "greater", "lower"))
  MI_init <- MI.resid(resid = resid_init, x = x, W = W, boot = boot.MI,
                      alternative = ifelse(dep == "positive", "greater", "lower"))

  #####
  # Output
  #####
  # Estimates
  est <- cbind(coefs_out[1:nx], se_out[1:nx], p.val[1:nx])
  colnames(est) <- c("Estimate", "SE", "p-value")
  varcovar <- solve(opt$hessian)[1:nx, 1:nx]
  if (nx == 1) {
    rownames(est) <- names(varcovar) <- "(Intercept)"
  } else {
    if (!is.null(nams)) {
      rownames(est) <- rownames(varcovar) <- colnames(varcovar) <- c("(Intercept)", nams)
    } else {
      rownames(est) <- c("(Intercept)", paste0("beta_", 1:(nx - 1)))
      rownames(varcovar) <- colnames(varcovar) <- rownames(est)
    }
  }

  # selected eigenvectors & eigenvalues
  if (count != 0) {
    # selected eigenvectors
    selvecs <- as.matrix(evecs[, sel_id])
    colnames(selvecs) <- paste0("evec_", sel_id)
    condnum <- conditionNumber(evecs = selvecs, round = 8) # multicollinearity: condnum > 1
    # EV matrix
    gammas <- coefs_out[(nx + 1):(nx + count)]
    gse <- se_out[(nx + 1):(nx + count)]
    gp <- 2 * pt(abs(gammas / gse), df = (n - length(coefs_out)), lower.tail = FALSE)
    EV <- cbind(gammas, gse, gp, evMI[sel_id])
    colnames(EV) <- c("Estimate", "SE", "p-value", "MI")
    rownames(EV) <- paste0("ev_", sel_id)
    # generate spatial filter
    sf <- selvecs %*% gammas
    sfMI <- MI.sf(gamma = gammas, evMI = evMI[sel_id])
  } else {
    selvecs <- EV <- sf <- sfMI <- condnum <- NULL
  }

  # Moran's I
  moran <- rbind(MI_init[, colnames(MI_init) != ""], MI_out[, colnames(MI_out) != ""])
  rownames(moran) <- c("Initial", "Filtered")
  colnames(moran) <- c("Observed", "Expected", "Variance", "z", "p-value")

  # model fit
  fit <- rbind(c(-ll_init, ICs_init), c(-ll_out, ICs_out))
  rownames(fit) <- c("Initial", "Filtered")
  colnames(fit) <- c("logL", "AIC", "BIC")

  # McFadden's pseudo-R2
  pRsqr <- pseudoR2(negloglik_n = ll_init, negloglik_f = ll_out, nev = count)

  # residuals
  residuals <- cbind(resid_init, resid_out)
  colnames(residuals) <- c("Initial", "Filtered")

  # define output
  out_list <- list(estimates = est,
                   varcovar = varcovar,
                   EV = EV,
                   selvecs = selvecs,
                   evMI = evMI,
                   moran = moran,
                   fit = fit,
                   residuals = residuals,
                   other = list(ncandidates = ncandidates,
                                nev = count,
                                condnum = condnum,
                                sel_id = sel_id,
                                sf = sf,
                                sfMI = sfMI,
                                model = model,
                                dependence = dep,
                                objfn = objfn,
                                bonferroni = bonferroni,
                                siglevel = sig,
                                resid.type = resid.type,
                                pseudoR2 = pRsqr
                   )
  )

  # define class
  class(out_list) <- "spfilter"

  # return
  return(out_list)
}
