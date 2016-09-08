context("regression")

test_that("regression predicts correctly", {
  x.train <- seq(-5, 5, .1)
  y.train <- rnorm(length(x.train))
  x.new <- x.train
  pars <- list(var=1, inv.scale=2, gamma=2, noise=.1, kernel="gamma.exp")
  pred <- lvgpr(x.train, y.train, x.new, pars)
  expect_equal(pred$y.predict, y.train, tolerance = 0.05)
})
