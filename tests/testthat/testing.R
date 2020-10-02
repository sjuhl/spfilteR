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
  covars <- fakedataset$x4
  EVs1 <- getEVs(W=W,covars=covars)$vectors
  covars2 <- cbind(1,fakedataset$x4)
  EVs2 <- getEVs(W=W,covars=covars2)$vectors
  out <- all(EVs1==EVs2)
  expect_true(out)
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
# getMoran()
#####
test_that("The argument 'alternative' in getMoran() needs to be either
          'greater', 'lower', or 'two.sided'", {
  y <- fakedataset$x1
  x <- fakedataset$x2
  resid <- y - x %*% solve(crossprod(x)) %*% crossprod(x,y)
  alternative <- c("greater","lower","two.sided","nothing")
  out <- NULL
  for(i in 1:4){
    res <- try(getMoran(resid=resid,x=x,W=W,alternative=alternative[i])
               ,silent=TRUE)
    out[i] <- class(res)!="try-error"
  }
  expect_equal(out, c(TRUE,TRUE,TRUE,FALSE))
})

test_that("getMoran() works without supplying x (intercept-only model)", {
  y <- fakedataset$x1
  x <- as.matrix(rep(1,length(y)))
  resid <- y - x %*% solve(crossprod(x)) %*% crossprod(x,y)
  out <- getMoran(resid=resid,W=W)
  expect_is(out, "data.frame")
})

test_that("W must be of class 'matrix', 'Matrix', or 'data.frame'", {
  y <- fakedataset$x1
  x <- as.matrix(rep(1,length(y)))
  resid <- y - x %*% solve(crossprod(x)) %*% crossprod(x,y)
  W2 <- as.vector(W)
  expect_error(getMoran(resid=resid,W=W2), "W must be of class 'matrix' or 'data.frame'")
})

test_that("getMoran() warns if boot<100", {
  y <- fakedataset$x1
  x <- as.matrix(rep(1,length(y)))
  resid <- y - x %*% solve(crossprod(x)) %*% crossprod(x,y)
  boot <- 95
  test <- tryCatch(getMoran(resid=resid,W=W,boot=boot),warning=function(t) TRUE)
  expect_true(test)
})

test_that("If boot < 100, getMoran() sets boot=100", {
  y <- fakedataset$x1
  x <- as.matrix(rep(1,length(y)))
  resid <- y - x %*% solve(crossprod(x)) %*% crossprod(x,y)
  boot1 <- 90
  boot2 <- 100
  set.seed(123)
  moran1 <- suppressWarnings(getMoran(resid=resid,W=W,boot=boot1))
  set.seed(123)
  moran2 <- getMoran(resid=resid,W=W,boot=boot2)
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
  for(i in 1:length(alpha)){
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
  filter <- lmFilter(y=fakedataset$x1,W=W,objfn="all",alpha=alpha)
  expect_equal(filter$other$nev,sum(inselset))
})

test_that("The argument 'objfn' works with 'p', 'MI', 'R2', and 'all' ", {
  y <- fakedataset$x1
  objfn <- c("MI","p","R2","all","else")
  expected <- c(TRUE,TRUE,TRUE,TRUE,FALSE)
  out <- NULL
  for(i in 1:length(objfn)){
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


#####
# glmFilter()
#####
test_that("glmFilter() estimates poisson models", {
  y <- fakedataset$count
  X <- cbind(fakedataset$x1,fakedataset$x2,fakedataset$x3,fakedataset$x4)
  out <- glmFilter(y=y,x=X,W=W,objfn="MI",model="poisson")
  expect_is(out, "spfilter")
})

test_that("glmFilter() estimates probit models", {
  y <- fakedataset$indicator
  X <- cbind(fakedataset$x1,fakedataset$x2,fakedataset$x3,fakedataset$x4)
  out <- glmFilter(y=y,x=NULL,W=W,objfn="MI",model="probit")
  expect_is(out, "spfilter")
})

test_that("glmFilter() estimates logit models", {
  y <- fakedataset$indicator
  X <- cbind(fakedataset$x1,fakedataset$x2,fakedataset$x3,fakedataset$x4)
  out <- glmFilter(y=y,x=X,W=W,objfn="MI",model="logit")
  expect_is(out, "spfilter")
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
            for(i in 1:length(alpha)){
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


#####
# methods
#####
test_that("Check the summary function", {
  filter <- lmFilter(y=fakedataset$x1,W=W,objfn="R2")
  expect_output(summary(filter))
})

test_that("Check the print method", {
  filter <- lmFilter(y=fakedataset$x1,W=W,objfn="R2")
  expect_output(print(filter))
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
                y <- fakedataset$x3
                fit <- fakedataset$x3+rnorm(length(y),0,.3)
              } else if(model[i] %in% c("probit","logit")){
                y <- fakedataset$indicator
                fit <- runif(length(y),.001,.999)
              } else {
                y <- fakedataset$count
                fit <- round(fakedataset$x3+rnorm(length(y),0,.3))
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

