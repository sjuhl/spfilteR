library(testthat)
library(spfilteR)

data(fakedata)

test_that("getEVs() returns a list", {
  out <- getEVs(W=W,covars=NULL)
  expect_is(out, "list")
})
