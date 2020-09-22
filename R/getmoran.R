
# Moran Coefficient for residuals
getmoran <- function(y,x=NULL,fitted.values=NULL,W,alternative="greater",boot=NULL){
  if(!(alternative %in% c("greater","lower", "two.sided"))){
    stop("Invalid input: 'alternative' must be either 'greater', 'lower', or 'two.sided'")
  }
  n <- nrow(W)
  if(is.null(x)) x <- rep(1,n)
  x <- as.matrix(x)
  if (!all(x[,1]==1)) x <- cbind(1,x) # add intercept term
  if(is.null(fitted.values)) fitted.values <- x %*% solve(crossprod(x), crossprod(x, y))
  resid <- y-fitted.values
  if(!isSymmetric(W)) W <- .5 * (W + t(W)) # symmetric connectivity matrix
  I <- n/crossprod(rep(1,n),W%*%rep(1,n)) * crossprod(resid,W%*%resid) / crossprod(resid)
  M <- diag(n)-x%*%qr.solve(crossprod(x),t(x))
  EI <- sum(diag(M%*%W))/(n-ncol(x))
  if(!is.null(boot)){
    if(boot<100){
      warning(paste0("Number of bootstrap iterations (",boot,") too small. Set to 100"))
      boot <- 100
    }
    boot.I <- NULL
    for(i in 1:boot){
      ind <- sample(1:n,replace=T)
      boot.I[i] <- n/crossprod(rep(1,n),W%*%rep(1,n)) * crossprod(resid[ind],W%*%resid[ind]) / crossprod(resid[ind])
    }
    VarI <- var(boot.I)
    ZI <- (I-EI)/sqrt(VarI)
    pI <- emp.pfunc(draws=boot.I,obs=I,alternative=alternative)
  } else {
    num <- sum(diag(M%*%W%*%M%*%t(W))) + sum(diag(M%*%W%*%M%*%W)) + sum(diag(M%*%W))^2
    denom <- (n-ncol(x))*(n-ncol(x)+2)
    VarI <- num/denom - EI^2
    if(VarI<=0){
      ZI <- 0
    } else ZI <- (I-EI)/sqrt(VarI)
    pI <- pfunc(z=ZI,alternative=alternative)
  }
  return(list(I=I,EI=EI,VarI=VarI,ZI=ZI,pI=pI))
}
