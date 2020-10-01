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
