library(testthat)
library(spfilteR)

# load data
data(fakedata)

#####
# getEVs()
#####
test_that("getEVs() returns a list", {
  out <- getEVs(W=W,covars=NULL)
  expect_is(out, "list")
})

test_that("getEVs() adds an intercept term", {
  covars <- fakedataset$x3
  EVs1 <- getEVs(W=W,covars=covars)$vectors[,1]
  covars2 <- cbind(1,fakedataset$x3)
  EVs2 <- getEVs(W=W,covars=covars2)$vectors[,1]
  expect_equal(EVs1, EVs2)
})


#####
# MI.vec()
#####
test_that("MI.vec() works with vector input and returns a data.frame", {
  x <- fakedataset$x1
  out <- MI.vec(x=x,W=W)
  expect_is(out, "data.frame")
})

test_that("MI.vec() works with matrix input and returns a data frame with
          number of columns equal to number of columns in x", {
  X <- cbind(fakedataset$x1,fakedataset$x2,fakedataset$x3)
  out <- MI.vec(x=X,W=W)
  expect_equal(nrow(out), ncol(X))
})

test_that("MI.vec() detects NAs and returns an error message", {
  x <- fakedataset$x1
  x[1] <- NA
  expect_error(MI.vec(x=x,W=W), "Missing values detected")
})

test_that("MI.vec() works with named input", {
  x <- as.matrix(fakedataset$x1)
  colnames(x) <- "name"
  out <- MI.vec(x=x,W=W)
  expect_equal(rownames(out), colnames(x))
})

test_that("MI.vec() warns if a constant is supplied", {
  X <- cbind(1,fakedataset$x1)
  expect_warning(MI.vec(x=X,W=W), "Constant term detected in x")
})

test_that("check the permissible attributes for 'alternative'", {
  x <- fakedataset$x1
  alternative <- c("greater","lower","two.sided","else")
  expect <- c(TRUE,TRUE,TRUE,FALSE)
  out <- NULL
  for(i in 1:4){
    res <- try(MI.vec(x=x,W=W,alternative=alternative[i])
               ,silent=TRUE)
    out[i] <- class(res)!="try-error"
  }
  expect_equal(out, expect)
})

test_that("W must be of class 'matrix', 'Matrix', or 'data.frame'", {
  W2 <- as.vector(W)
  expect_error(MI.vec(x=fakedataset$x1,W=W2), "W must be of class 'matrix' or 'data.frame'")
})


#####
# MI.local()
#####
test_that("MI.local() works with vector input and returns a data.frame", {
  x <- fakedataset$x1
  out <- MI.local(x=x,W=W,alternative="greater")
  expect_is(out, "data.frame")
})

test_that("MI.local() stores vector names", {
  x <- fakedataset$x1
  names(x) <- paste0("name_",1:length(x))
  out <- MI.local(x=x,W=W)
  expect_equal(rownames(out), names(x))
})

test_that("MI.local() detects NAs and returns an error message", {
  x <- fakedataset$x1
  x[1] <- NA
  expect_error(MI.local(x=x,W=W), "Missing values detected")
})

test_that("check the permissible attributes for 'alternative'", {
  x <- fakedataset$x1
  alternative <- c("greater","lower","two.sided","else")
  expect <- c(TRUE,TRUE,TRUE,FALSE)
  out <- NULL
  for(i in 1:4){
    res <- try(MI.local(x=x,W=W,alternative=alternative[i])
               ,silent=TRUE)
    out[i] <- class(res)!="try-error"
  }
  expect_equal(out, expect)
})

test_that("W must be of class 'matrix', 'Matrix', or 'data.frame'", {
  W2 <- as.vector(W)
  expect_error(MI.local(x=fakedataset$x1,W=W2), "W must be of class 'matrix' or 'data.frame'")
})


#####
# MI.decomp()
#####
test_that("MI.decomp() works with matrix input and returns a data frame with
          number of columns equal to number of columns in x", {
            X <- cbind(fakedataset$x1,fakedataset$x2,fakedataset$x3)
            out <- MI.decomp(x=X,W=W,nsim=100)
            expect_equal(nrow(out), ncol(X))
          })

test_that("MI.decomp() detects NAs and returns an error message", {
  x <- fakedataset$x1
  x[1] <- NA
  expect_error(MI.decomp(x=x,W=W), "Missing values detected")
})

test_that("MI.decomp() works with named input", {
  x <- as.matrix(fakedataset$x1)
  colnames(x) <- "name"
  out <- MI.decomp(x=x,W=W)
  expect_equal(rownames(out), colnames(x))
})

test_that("MI.decomp() warns if a constant is supplied", {
  X <- cbind(1,fakedataset$x1)
  expect_warning(MI.decomp(x=X,W=W), "Constant term detected in x")
})

test_that("W must be of class 'matrix', 'Matrix', or 'data.frame'", {
  W2 <- as.vector(W)
  expect_error(MI.decomp(x=fakedataset$x1,W=W2), "W must be of class 'matrix' or 'data.frame'")
})

test_that("MI.decomp() warns if nsim < 100", {
  nsim <- 99
  test <- tryCatch(MI.decomp(x=fakedataset$x1,W=W,nsim=nsim),warning=function(t) TRUE)
  expect_true(test)
})

test_that("If nsim < 100, MI.decomp() sets nsim=100", {
  x <- fakedataset$x1
  nsim1 <- 90
  nsim2 <- 100
  set.seed(123)
  MI1 <- suppressWarnings(MI.decomp(x=x,W=W,nsim=nsim1))
  set.seed(123)
  MI2 <- MI.decomp(x=x,W=W,nsim=nsim2)
  expect_equal(MI1[,"VarI+"],MI2[,"VarI+"])
})


#####
# MI.resid()
#####
test_that("The argument 'alternative' in MI.resid() needs to be either
          'greater', 'lower', or 'two.sided'", {
  y <- fakedataset$x1
  x <- fakedataset$x2
  resid <- y - x %*% solve(crossprod(x)) %*% crossprod(x,y)
  alternative <- c("greater","lower","two.sided","nothing")
  out <- NULL
  for(i in 1:4){
    res <- try(MI.resid(resid=resid,x=x,W=W,alternative=alternative[i])
               ,silent=TRUE)
    out[i] <- class(res)!="try-error"
  }
  expect_equal(out, c(TRUE,TRUE,TRUE,FALSE))
})

test_that("MI.resid() works without supplying x (intercept-only model)", {
  y <- fakedataset$x1
  x <- as.matrix(rep(1,length(y)))
  resid <- y - x %*% solve(crossprod(x)) %*% crossprod(x,y)
  out <- MI.resid(resid=resid,W=W)
  expect_is(out, "data.frame")
})

test_that("W must be of class 'matrix', 'Matrix', or 'data.frame'", {
  y <- fakedataset$x1
  x <- as.matrix(rep(1,length(y)))
  resid <- y - x %*% solve(crossprod(x)) %*% crossprod(x,y)
  W2 <- as.vector(W)
  expect_error(MI.resid(resid=resid,W=W2), "W must be of class 'matrix' or 'data.frame'")
})

test_that("MI.resid() warns if boot<100", {
  y <- fakedataset$x1
  x <- as.matrix(rep(1,length(y)))
  resid <- y - x %*% solve(crossprod(x)) %*% crossprod(x,y)
  boot <- 95
  test <- tryCatch(MI.resid(resid=resid,W=W,boot=boot),warning=function(t) TRUE)
  expect_true(test)
})

test_that("If boot < 100, MI.resid() sets boot=100", {
  y <- fakedataset$x1
  x <- as.matrix(rep(1,length(y)))
  resid <- y - x %*% solve(crossprod(x)) %*% crossprod(x,y)
  boot1 <- 90
  boot2 <- 100
  set.seed(123)
  moran1 <- suppressWarnings(MI.resid(resid=resid,W=W,boot=boot1))
  set.seed(123)
  moran2 <- MI.resid(resid=resid,W=W,boot=boot2)
  expect_equal(moran1$VarI,moran2$VarI)
})


#####
# partialR2()
#####
test_that("partialR2() produces a warning if no eigenvectors were
          supplied/ selected", {
  y <- fakedataset$x1
  x <- fakedataset$x2
  expect_warning(partialR2(y=y,x=x,evecs=NULL), "No eigenvectors supplied")
})

test_that("partialR2() works without covariates (intercept-only model)", {
  y <- fakedataset$x1
  E <- getEVs(W=W)$vectors[,1:5]
  res <- try(partialR2(y=y,x=NULL,evecs=E),silent=TRUE)
  out <- class(res)!="try-error"
  expect_true(out)
})

test_that("partialR2() preserves number/names of eigenvectors (if supplied)", {
  y <- fakedataset$x1
  x <- fakedataset$x2
  E <- getEVs(W=W)$vectors[,c(2,4,6)]
  colnames(E) <- c("2","4","6")
  out <- partialR2(y=y,x=x,evecs=E)
  expect_equal(names(out), c("2","4","6"))
})

test_that("partialR2() detects perfect multicollinearity", {
  y <- fakedataset$x1
  X <- cbind(1,fakedataset$x2,fakedataset$x2)
  E <- getEVs(W=W)$vectors[,c(2,4,6)]
  expect_error(partialR2(y=y,x=X,evecs=E)
               ,"Perfect multicollinearity in covariates detected")
})

test_that("partialR2() detects missings", {
  y <- fakedataset$x1
  X <- cbind(1,fakedataset$x2)
  X[1,2] <- NA
  E <- getEVs(W=W)$vectors[,1:5]
  expect_error(partialR2(y=y,x=X,evecs=E)
               ,"Missing values detected")
})


#####
# vif.ev()
#####
test_that("vif.ev() returns an object with length equal to the number of evecs", {
  E <- getEVs(W=W)$vectors[,1:2]
  expect_equal(length(vif.ev(evecs=E)), ncol(E))
})

test_that("vif.ev() removes missings in supplied covariates", {
  X <- cbind(fakedataset$x1,fakedataset$x2)
  X[3,1] <- NA
  E <- getEVs(W=W)$vectors[,1:2]
  vif.ev(x=X,evecs=E,na.rm=TRUE)
  expect_equal(length(vif.ev(evecs=E)), ncol(E))
})

test_that("vif.ev() returns error if missings in covariates but na.rm=FALSE", {
  X <- cbind(fakedataset$x1,fakedataset$x2)
  X[3,1] <- NA
  E <- getEVs(W=W)$vectors[,1:2]
  expect_error(vif.ev(x=X,evecs=E,na.rm=FALSE), "Missing values detected")
})


#####
# lmFilter()
#####
test_that("lmFilter() detects perfect multicollinearity in covariates", {
  y <- fakedataset$x1
  X <- cbind(fakedataset$x2,fakedataset$x3,fakedataset$x2)
  expect_error(lmFilter(y=y,x=X,W=W)
               ,"Perfect multicollinearity in covariates detected")
})

test_that("The argument 'ideal.setsize' is only valid if positive=TRUE", {
  y <- fakedataset$x1
  expect_error(lmFilter(y=y,W=W,positive=FALSE,ideal.setsize=TRUE)
               ,"Estimating the ideal set size is only valid for positive spatial autocorrelation")
})

test_that("lmFilter() works if missings are present in y", {
  y <- fakedataset$x1
  y[6] <- NA
  res <- try(lmFilter(y=y,W=W,na.rm=TRUE),silent=TRUE)
  out <- class(res)!="try-error"
  expect_true(out)
})

test_that("lmFilter() returns an error message if 'alpha' is not contained
          in (0,1]", {
  y <- fakedataset$x1
  alpha <- c(-1,.001,.5,1,1.0001)
  expected <- c(TRUE,FALSE,FALSE,FALSE,TRUE)
  out <- NULL
  for(i in seq_along(alpha)){
    res <- try(lmFilter(y=y,W=W,alpha=alpha[i]),silent=TRUE)
    out[i] <- class(res)=="try-error"
  }
  expect_equal(out,expected)
})

test_that("'objfn=all' selects all eigenvectors in the candidate set", {
  eigen <- getEVs(W=W)
  alpha <- .25
  inselset <- eigen$moran/eigen$moran[1]>=alpha
  E <- eigen$vectors[,inselset]
  sf <- lmFilter(y=fakedataset$x1,W=W,objfn="all",alpha=alpha)
  expect_equal(sf$other$nev,sum(inselset))
})

test_that("The argument 'objfn' works with 'p', 'MI', 'R2', 'pMI', and 'all' ", {
  y <- fakedataset$x1
  objfn <- c("MI","p","R2","all","else","pMI")
  expected <- c(TRUE,TRUE,TRUE,TRUE,FALSE,TRUE)
  out <- NULL
  for(i in seq_along(objfn)){
    res <- try(lmFilter(y=y,W=W,objfn=objfn[i]),silent=TRUE)
    out[i] <- class(res)!="try-error"
  }
  expect_equal(out,expected)
})

test_that("If 'objfn=p': selects fewer EVs if 'bonferroni=T'", {
  bonferroni <- lmFilter(y=fakedataset$x1,W=W,objfn="p",bonferron=TRUE)
  nobonferroni <- lmFilter(y=fakedataset$x1,W=W,objfn="p",bonferron=FALSE)
  expect_true(bonferroni$other$nev < nobonferroni$other$nev)
})

test_that("W must be of class 'matrix', 'Matrix', or 'data.frame'", {
  W2 <- as.vector(W)
  expect_error(lmFilter(y=fakedataset$x1,W=W2,objfn="all",na.rm=F)
               ,"W must be of class 'matrix' or 'data.frame'")
})

test_that("lmFilter() works with named input", {
  y <- fakedataset$x1
  x <- as.matrix(fakedataset$x2)
  colnames(x) <- "name"
  out <- lmFilter(y=y,x=x,W=W,objfn="p",bonferroni=FALSE,sig=.05)
  expect_equal(rownames(out$estimates)[2], colnames(x))
})

test_that("lmFilter() returns an error message if 'alpha' is not contained
          in (0,1]", {
            y <- fakedataset$x1
            alpha <- c(-1,.001,.5,1,1.0001)
            expected <- c(TRUE,FALSE,FALSE,FALSE,TRUE)
            out <- NULL
            for(i in seq_along(alpha)){
              res <- try(lmFilter(y=y,W=W,objfn="all",alpha=alpha[i])
                         ,silent=TRUE)
              out[i] <- class(res)=="try-error"
            }
            expect_equal(out,expected)
          })

test_that("If 'alpha=0', lmFilter() sets 'alpha' to 1e-07", {
  y <- fakedataset$x1
  res <- try(lmFilter(y=y,W=W,objfn="all",alpha=0),silent=TRUE)
  out <- class(res)!="try-error"
  expect_true(out)
})

test_that("W must be of class 'matrix', 'Matrix', or 'data.frame'", {
  W2 <- as.vector(W)
  expect_error(lmFilter(y=fakedataset$x1,W=W2,objfn="all",na.rm=F)
               ,"W must be of class 'matrix' or 'data.frame'")
})

test_that("lmFilter() accepts covariates in 'MX'", {
  y <- fakedataset$x1
  X <- cbind(1,fakedataset$x3)
  MX <- X
  res <- try(lmFilter(y=y,x=X,W=W,MX=X,objfn="all",positive=F),silent=TRUE)
  out <- class(res)!="try-error"
  expect_true(out)
})

test_that("lmFilter() also works with negative autocorrelation", {
  y <- fakedataset$negative
  res <- try(lmFilter(y=y,W=W,objfn="all",positive=F),silent=TRUE)
  out <- class(res)!="try-error"
  expect_true(out)
})

test_that("lmFilter() breaks if missings are present but na.rm=FALSE", {
  y <- fakedataset$x2
  y[6] <- NA
  expect_error(lmFilter(y=y,W=W,na.rm=FALSE)
               ,"Missing values detected")
})

test_that("check ideal candidate set size", {
  sf <- lmFilter(y=fakedataset$x3,W=W,positive=TRUE
                 ,ideal.setsize=TRUE)
  expect_is(sf, "spfilter")
})

test_that("check 'tol' in lmFilter()", {
  tol <- c(.5,.01)
  nev <- NULL
  for(i in seq_along(tol)){
    sf <- lmFilter(y=fakedataset$x1,objfn="MI",alpha=.05,W=W,tol=tol[i])
    nev[i] <- sf$other$nev
  }
  expect_true(nev[1]<nev[2])
})

test_that("check 'pMI' with positive autocorrelation in lmFilter()", {
  sf <- lmFilter(y=fakedataset$x1,W=W,objfn="pMI",positive=TRUE
                 ,bonferroni=FALSE)
  expect_is(sf, "spfilter")
})

test_that("lmFilter sets 'bonferroni=FALSE' for 'pMI'", {
  sf <- lmFilter(y=fakedataset$x1,W=W,objfn="pMI",positive=TRUE
                 ,bonferroni=TRUE)
  expect_false(sf$other$bonferroni)
})

test_that("check 'pMI' with negative autocorrelation in lmFilter()", {
  sf <- lmFilter(y=fakedataset$negative,W=W,objfn="pMI",positive=FALSE)
  expect_is(sf, "spfilter")
})

test_that("selects no EVs if objfn=='pMI' and initial residuals are insignificant", {
  sf <- lmFilter(y=fakedataset$x3,W=W,objfn="pMI",positive=TRUE
                 ,bonferroni=TRUE)
  expect_equal(sf$other$nev,0)
})



#####
# glmFilter()
#####
test_that("glmFilter() estimates poisson models", {
  y <- fakedataset$count
  X <- cbind(fakedataset$x1,fakedataset$x2,fakedataset$x3)
  out <- glmFilter(y=y,x=X,W=W,objfn="MI",model="poisson")
  expect_is(out, "spfilter")
})

test_that("glmFilter() estimates probit models", {
  y <- fakedataset$indicator
  X <- cbind(fakedataset$x1,fakedataset$x2,fakedataset$x3,fakedataset$x3)
  out <- glmFilter(y=y,x=NULL,W=W,objfn="MI",model="probit")
  expect_is(out, "spfilter")
})

test_that("glmFilter() estimates logit models", {
  y <- fakedataset$indicator
  out <- glmFilter(y=y,W=W,objfn="MI",model="logit")
  expect_is(out, "spfilter")
})

test_that("check ideal candidate set size", {
  y <- fakedataset$indicator
  sf <- glmFilter(y=y,W=W,objfn="p",model="logit",positive=TRUE
                  ,ideal.setsize=TRUE)
  expect_is(sf, "spfilter")
})

test_that("The argument 'ideal.setsize' is only valid if positive=TRUE", {
  y <- fakedataset$indicator
  expect_error(glmFilter(y=y,W=W,objfn="p",model="logit",positive=FALSE,ideal.setsize=TRUE)
               ,"Estimating the ideal set size is only valid for positive spatial autocorrelation")
})

test_that("The argument 'min.reduction' works (for AIC & BIC) - fewer EVs are
          selected for higher values", {
  y <- fakedataset$count
  lowerAIC <- glmFilter(y=y,W=W,objfn="AIC",model="poisson",min.reduction=0)
  higherAIC <- glmFilter(y=y,W=W,objfn="AIC",model="poisson",min.reduction=.1)
  outAIC <- lowerAIC$other$nev > higherAIC$other$nev
  lowerBIC <- glmFilter(y=y,W=W,objfn="BIC",model="poisson",min.reduction=0)
  higherBIC <- glmFilter(y=y,W=W,objfn="BIC",model="poisson",min.reduction=.1)
  outBIC <- lowerBIC$other$nev > higherBIC$other$nev
  expect_true(all(outAIC,outBIC))
})

test_that("glmFilter() works with named input", {
  y <- fakedataset$count
  x <- as.matrix(fakedataset$x1)
  colnames(x) <- "name"
  out <- glmFilter(y=y,x=x,W=W,objfn="p",bonferroni=FALSE
                   ,model="poisson",sig=.05)
  expect_equal(rownames(out$estimates)[2], colnames(x))
})

test_that("glmFilter() breaks if missings are present but na.rm=FALSE", {
  y <- fakedataset$count
  y[6] <- NA
  expect_error(glmFilter(y=y,W=W,model="poisson",na.rm=FALSE)
               ,"Missing values detected")
})

test_that("glmFilter() works even if no EVs are selected", {
  y <- fakedataset$count
  res <- glmFilter(y=y,W=W,model="poisson",objfn="AIC",min.reduction=.9)
  expect_equal(res$other$nev,0)
})

test_that("glmFilter() returns an error message if 'alpha' is not contained
          in (0,1]", {
            y <- fakedataset$indicator
            alpha <- c(-1,.001,.5,1,1.0001)
            expected <- c(TRUE,FALSE,FALSE,FALSE,TRUE)
            out <- NULL
            for(i in seq_along(alpha)){
              res <- try(glmFilter(y=y,W=W,model="logit",objfn="p",alpha=alpha[i])
                         ,silent=TRUE)
              out[i] <- class(res)=="try-error"
            }
            expect_equal(out,expected)
})

test_that("If 'alpha=0', glmFilter() sets 'alpha' to 1e-07", {
  y <- fakedataset$indicator
  x <- fakedataset$x3
  res <- try(glmFilter(y=y,W=W,model="probit",objfn="BIC",alpha=0),silent=TRUE)
  out <- class(res)!="try-error"
  expect_true(out)
})

test_that("glmFilter() detects perfect multicollinearity", {
  y <- fakedataset$indicator
  X <- cbind(1,fakedataset$x2,fakedataset$x2)
  expect_error(glmFilter(y=y,x=X,W=W,model="probit",objfn="p")
               ,"Perfect multicollinearity in covariates detected")
})

test_that("glmFilter() works if all EVs should be selected", {
  y <- fakedataset$count
  expect_is(glmFilter(y=y,W=W,model="poisson",objfn="all")
               ,"spfilter")
})

test_that("if 'objfn=AIC' or 'objfn=BIC', 'min.reduction' needs to be in [0,1)", {
  y <- fakedataset$indicator
  reduce <- c(0,.3,1)
  expect <- c(TRUE,TRUE,FALSE)
  out <- NULL
  for(i in 1:3){
    res <- try(glmFilter(y=y,W=W,model="logit",objfn="AIC",min.reduction=reduce[i])
               ,silent=TRUE)
    out[i] <- class(res)!="try-error"
  }
  expect_equal(out,expect)
})

test_that("W must be of class 'matrix', 'Matrix', or 'data.frame'", {
  W2 <- as.vector(W)
  expect_error(glmFilter(y=fakedataset$count,W=W2,model="poisson",objfn="all",na.rm=F)
               ,"W must be of class 'matrix' or 'data.frame'")
})

test_that("'model' must be one of 'poisson', 'probit', or 'logit'", {
  model <- c("probit","logit","poisson","else")
  expect <- c(TRUE,TRUE,TRUE,FALSE)
  out <- NULL
  for(i in 1:4){
    if(model[i] %in% c("probit","logit")) y <- fakedataset$indicator
    else y <- fakedataset$count
    res <- try(glmFilter(y=y,W=W,model=model[i],objfn="BIC")
               ,silent=TRUE)
    out[i] <- class(res)!="try-error"
  }
  expect_equal(out,expect)
})

test_that("The argument 'objfn' works with 'p', 'MI', 'pMI', 'AIC', 'BIC', and 'all'", {
  y <- fakedataset$count
  objfn <- c("MI","p","AIC","BIC","all","else","pMI")
  expected <- c(TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE)
  out <- NULL
  for(i in seq_along(objfn)){
    res <- try(glmFilter(y=y,W=W,model="poisson",objfn=objfn[i]),silent=TRUE)
    out[i] <- class(res)!="try-error"
  }
  expect_equal(out,expected)
})

test_that("The argument 'resid.type' must be 'raw','pearson', or 'deviance'", {
  y <- fakedataset$count
  resid.type <- c("raw","pearson","deviance","else")
  expected <- c(TRUE,TRUE,TRUE,FALSE)
  out <- NULL
  for(i in 1:4){
    res <- try(glmFilter(y=y,W=W,model="poisson",objfn="BIC",resid.type=resid.type[i])
               ,silent=TRUE)
    out[i] <- class(res)!="try-error"
  }
  expect_equal(out,expected)
})

test_that("eigenvectors are orthogonal to covariates specified in MX", {
  y <- fakedataset$indicator
  X <- cbind(1,fakedataset$x2)
  res <- glmFilter(y=y,x=X,W=W,MX=X,model="probit",objfn="MI")
  correlation <- cor(cbind(X[,2],res$selvecs),method="pearson")
  offdiag <- correlation - diag(1,nrow(correlation))
  max <- max(offdiag)
  expect_true(max<1e-07)
})

test_that("check that glmFilter() works with negative autocorrelation", {
  y <- fakedataset$negcount
  res <- glmFilter(y=y,W=W,model="poisson",objfn="MI",positive=FALSE)
  expect_equal(res$other$dependence,"negative")
})

test_that("check 'pMI' with positive autocorrelation in glmFilter()", {
  sf <- glmFilter(y=fakedataset$count,W=W,objfn="pMI",model="poisson"
                  ,positive=TRUE,bonferroni=FALSE)
  expect_is(sf, "spfilter")
})

test_that("lmFilter sets 'bonferroni=FALSE' for 'pMI'", {
  sf <- glmFilter(y=fakedataset$count,W=W,objfn="pMI",model="poisson"
                  ,positive=TRUE,bonferroni=TRUE)
  expect_false(sf$other$bonferroni)
})

test_that("check 'pMI' with negative autocorrelation in lmFilter()", {
  sf <- glmFilter(y=fakedataset$negcount,W=W,objfn="pMI",model="poisson"
                 ,positive=FALSE)
  expect_equal(sf$other$dependence, "negative")
})

test_that("selects EVs if objfn=='pMI' and initial residuals are significant", {
  y <- fakedataset$indicator
  X <- cbind(1,fakedataset$x3)
  sf <- glmFilter(y=y,x=X,W=W,objfn="pMI",model="poisson",sig=.15,positive=TRUE
                 ,bonferroni=FALSE,boot.MI=NULL,resid.type="deviance")
  expect_true(sf$other$nev>0)
})


#####
# vp()
#####
test_that("If msr<100, vp() gives a warning", {
  y <- fakedataset$x1
  msr <- c(100,99,101,1)
  expect <- c(FALSE,TRUE,FALSE,TRUE)
  out <- NULL
  for(i in 1:length(expect)){
    res <- tryCatch(vp(y=y,evecs=NULL,msr=msr[i]),warning=function(x) TRUE)
    out[i] <- isTRUE(res)
  }
  expect_equal(out,expect)
})

test_that("If msr<100, vp() sets it to 100", {
  y <- fakedataset$x1
  out <- suppressWarnings(vp(y=y,x=NULL,evecs=NULL,msr=1))
  expect_equal(out$msr,100)
})

test_that("vp() detects perfect multicollinearity", {
  X <- cbind(fakedataset$x2,fakedataset$x2*2)
  expect_error(vp(y=fakedataset$x1,x=X,evecs=NULL,msr=100)
               ,"Perfect multicollinearity in covariates detected")
})

test_that("vp() detects NAs", {
  y <- fakedataset$x2
  y[1] <- NA
  expect_error(vp(y=y,x=NULL,evecs=NULL,msr=100)
               ,"Missing values detected")
})

test_that("vp() works if a single eigenvalue equals zero", {
  evecs <- getEVs(W=W)$vectors
  y <- fakedataset$x3
  zeroes <- which(!vapply(apply(evecs,2,sum)
                          ,function(x) isTRUE(all.equal(x,0,tolerance=1e-7))
                          ,FUN.VALUE = TRUE))
  part <- vp(y=y,x=NULL,evecs=cbind(evecs[,1:20],evecs[,zeroes[1]]),msr=100)
  expect_equal(class(part),"vpart")
})

test_that("vp() works if multiple eigenvalues equal zero", {
  evecs <- getEVs(W=W)$vectors
  y <- fakedataset$x3
  part <- vp(y=y,x=NULL,evecs=evecs[,20:80],msr=100)
  expect_output(print(part))
})

test_that("Error if nr of covars >= n", {
  evecs <- getEVs(W=W)$vectors
  y <- fakedataset$x3
  expect_error(vp(y=y,x=NULL,evecs=evecs,msr=100)
               ,"Nr of covariates equals or exceeds n")
})


#####
# methods
#####
test_that("check the summary function", {
  sf <- lmFilter(y=fakedataset$x1,W=W,objfn="R2")
  expect_output(summary(sf))
})

test_that("summary function for glmFilter()", {
  y <- fakedataset$count
  X <- cbind(1,fakedataset$x1)
  sf <- glmFilter(y=fakedataset$count,x=X,W=W,model="poisson",objfn="MI")
  expect_output(summary(sf))
})

test_that("summary function shows significance level without adjustment", {
  sf <- lmFilter(y=fakedataset$x3,W=W,objfn="p",bonferroni=FALSE)
  expect_output(summary(sf))
})

test_that("summary function shows significance level with adjustment", {
  sf <- lmFilter(y=fakedataset$x3,W=W,objfn="p",bonferroni=TRUE)
  expect_output(summary(sf,EV=TRUE))
})

test_that("summary function summarizes eigenvectors", {
  sf <- lmFilter(y=fakedataset$x3,W=W,objfn="all")
  expect_output(summary(sf,EV=TRUE))
})

test_that("plot method", {
  sf <- lmFilter(y=fakedataset$x1,W=W,objfn="all")
  out <- plot(sf)
  expect_is(out, "NULL")
})

test_that("check the print method", {
  sf <- lmFilter(y=fakedataset$x1,W=W,objfn="R2")
  expect_output(print(sf))
})

test_that("coef() gives the correct number of coefs", {
  X <- cbind(fakedataset$x2,fakedataset$x3)
  out <- lmFilter(y=fakedataset$x1,x=X,W=W,objfn="R2")
  expect_equal(length(coef(out)),(ncol(X)+1))
})

test_that("correct dimensionality of vcov()", {
  X <- cbind(fakedataset$x2,fakedataset$x3)
  out <- lmFilter(y=fakedataset$x1,x=X,W=W,objfn="R2")
  dim1 <- dim(vcov(out))[1]
  dim2 <- dim(vcov(out))[2]
  out <- dim1==dim2 & dim1==(ncol(X)+1)
  expect_true(out)
})


#####
# utils
#####
test_that("pfunc() works with empirical p-values and 'alternative=lower'", {
  z <- .5
  draws <- rnorm(100,.5,.2)
  alternative <- c("else","lower")
  out <-NULL
  for(i in 1:2){
    out[i] <- is.numeric(pfunc(z=z,alternative=alternative[i],draws=draws))
  }
  expect_true(all(out))
})

test_that("candsetsize() returns an integer", {
  zMI <- .5
  npos <- sum(getEVs(W)$moran>0)
  out <- candsetsize(npos=npos,zMI=zMI)
  expect_true(out==round(out))
})

test_that("residfun() takes as model types 'linear', 'probit'
          ,'logit', or 'poisson' ", {
            model <- c("linear","probit","logit","poisson","else")
            expect <- c(TRUE,TRUE,TRUE,TRUE,FALSE)
            out <- NULL
            for(i in 1:5){
              if(model[i]=="linear"){
                y <- fakedataset$x1
                fit <- fakedataset$x1+rnorm(length(y),0,.3)
              } else if(model[i] %in% c("probit","logit")){
                y <- fakedataset$indicator
                fit <- runif(length(y),.001,.999)
              } else {
                y <- fakedataset$count
                fit <- round(fakedataset$x1+rnorm(length(y),0,.3))
              }
              res <- try(residfun(y=y,fitvals=fit,model=model[i])
                         ,silent=TRUE)
              out[i] <- class(res)!="try-error"
            }
            expect_equal(out,expect)
})

test_that("conditionNumber() returns NULL if no EVs supplied", {
  out <- conditionNumber(evecs=NULL)
  expect_null(out)
})

test_that("Zscore() returns standardized values", {
  x <- fakedataset$x1
  z <- Zscore(x)
  expect_true(all(round(mean(z),10)==0 & sd(z)==1))
})

