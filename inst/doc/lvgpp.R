## ------------------------------------------------------------------------
library(lvgpp)

## ---- eval=F-------------------------------------------------------------
#  lvgpp::demo.regression()
#  lvgpp::demo.bin.classification()

## ------------------------------------------------------------------------
 x.train <- seq(-5, 5, .1)
 y.train <- rnorm(length(x.train))

## ------------------------------------------------------------------------
 x.new <- rnorm(100)

## ------------------------------------------------------------------------
 pars <- list(var=1, inv.scale=2, gamma=2, noise=.1, kernel="gamma.exp")

## ------------------------------------------------------------------------
 pred <- lvgpr(x.train, y.train, x.new, pars)

## ---- fig.width=6--------------------------------------------------------
 plot(x.train, y.train, type="l", xlab="X", ylab="Y")
 points(x.new, pred$y.predict, col="blue")

## ------------------------------------------------------------------------
 x.train <- seq(-5, 5, .1)
 c.train <- c(rep(0,25), rep(1, length(x.train)-50), rep(0,25))
 x.new <- rnorm(100, 0, 2)

## ------------------------------------------------------------------------
 pars <- list(var=1, inv.scale=2, gamma=2, noise=.1, kernel="gamma.exp")

## ---- , fig.width=6------------------------------------------------------
 pred <- lvgpp::lvgpc(x.train, c.train, x.new, pars)
 plot(x.train, c.train, type="l", xlab="X", ylab="C")
 points(x.new, pred$c.label, col="blue")

## ------------------------------------------------------------------------
 x <- rnorm(10)
 K <- covariance.function(x, x, list(var=1, inv.scale=2, gamma=2, noise=0, kernel="gamma.exp"))
 K[1:5, 1:5]

## ------------------------------------------------------------------------
 x = seq(-10, 10, .25)
 m = rep(0, length(x))
 K = covariance.function(x, x, list(var=1, inv.scale=2, gamma=2, noise=0, kernel="gamma.exp"))
 sample.from.gp(x, m, K)[1:5]

