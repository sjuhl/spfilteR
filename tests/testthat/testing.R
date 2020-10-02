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

