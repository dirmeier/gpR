## ------------------------------------------------------------------------
library(gpR)

## ---- eval=F-------------------------------------------------------------
#  gpR::demo.regression()
#  gpR::demo.bin.classification()

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
 pred <- lvgpc(x.train, c.train, x.new, pars)
 plot(x.train, c.train, type="l", xlab="X", ylab="C")
 points(x.new, pred$c.label, col="blue")

