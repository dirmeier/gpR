# gpR: Gaussian processes in R
#
# Copyright (C) 2015 - 2016 Simon Dirmeier
#
# This file is part of gpR.
#
# gpR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gpR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with gpR. If not, see <http://www.gnu.org/licenses/>.

context("binomial approximation")

predict.r <- function(K, c.train, sig, x.new, train, test)
{
  cov.train.train <- K[train, train]
  c.mean <- c.train - sig
  Dinv <- diag(1 / (sig * (1 - sig)))
  pred.means <- pred.samples <- c()
  Cu <- solve(cov.train.train + Dinv)
  for (i in seq(length(x.new)))
  {
    cov.new.train <- t(K[train, test[i]])
    cov.new.new <- K[test[i], test[i]]
    m <- cov.new.train %*% c.mean
    var <- cov.new.new - cov.new.train %*%
      solve(cov.train.train + Dinv) %*% t(cov.new.train)
    me <- .logistic.normal.integral(m, var)
    rn <- .sample.gaussian(m, var)
    pre <- .sigmoid(rn)
    pred.samples <- c(pred.samples , pre)
    pred.means   <- c(pred.means, me)
  }
  m.new <- (K[test, train] %*% c.mean)[,1]
  cov.new.new <- K[test, test] - K[test, train] %*% Cu %*% K[train, test]
  list(m.new=m.new,
       cov.new.new=cov.new.new,
       pred.samples=pred.samples,
       pred.means=pred.means)
}

predict.f <- function(K, c.train, sig, x.new, train, test)
{
  ntrain <- length(train)
  nnew   <- length(test)
  n <- ntrain + nnew
  me <- pred.sample <- pred.means <- rep(0, nnew)
  K.new.new <- matrix(0, nnew, nnew)
  res <- .Fortran("predict_binomial",
                  na=as.integer(n),
                  ntrain=as.integer(ntrain),
                  nnew=as.integer(nnew),
                  K=K,
                  ctrain=c.train,
                  sigtrain=sig,
                  xnew=x.new,
                  trainidx=as.integer(train),
                  newidx=as.integer(test),
                  mnew=me,
                  ps=pred.sample,
                  pm=pred.means,
                  Kn=K.new.new,
                  PACKAGE="gpR")
  list(m.new=res$mnew,
       cov.new.new=res$Kn,
       pred.samples=res$ps,
       pred.means=res$pm)
}

test_that("fortran binomial prediction", {
  ntrain <- 100
  nnew <- 10
  x.train <- rnorm(ntrain)
  x.new  <- rnorm(nnew)
  x.norm <- c(x.train, x.new)
  K <- covariance.function(x.norm, x.norm)
  c.train <- c(rep(0, ntrain/2), rep(1, ntrain/2))
  sig <- c(runif(ntrain/2, 0, 0.5), runif(ntrain/2, 0.5, 1))
  train <- 1:ntrain
  test <- 1:nnew + ntrain
  res.r <- predict.r(K, c.train, sig, x.new, train, test)
  res.f <- predict.f(K, c.train, sig, x.new, train, test)
  expect_equal(res.r$m.new, res.f$m.new, tolerance=.01)
  expect_equal(res.r$pred.means, res.f$pred.means, tolerance=.01)
})
