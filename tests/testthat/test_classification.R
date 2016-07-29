library(lvgpp)

context("classification")

test_that("classification predicts correctly", {
  x.train <- seq(-5, 5, .1)
  c.train <- c(rep(0,25), rep(1, length(x.train)-50), rep(0,25))
  x.new <- x.train
  pars <- list(var=1, inv.scale=2, gamma=2, noise=.1, kernel="gamma.exp")
  pred <- lvgpc(x.train, c.train, x.new, pars)$mean.c.predict
  pred[pred < .5] <- 0
  pred[pred >= .5] <- 1
  expect_equal(pred, c.train, tolerance = 0.005 )
})
