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

context("classification")

test_that("classification predicts correctly", {
  x.train <- seq(-5, 5, .1)
  c.train <- c(rep(0,25), rep(1, length(x.train)-50), rep(0,25))
  x.new <- x.train
  pars <- list(var=1, inv.scale=2, gamma=2, noise=.1, kernel="gamma.exp")
  pred <- lvgpc(x.train, c.train, x.new, pars)$mean.c.predict
  pred[pred < .5] <- 0
  pred[pred >= .5] <- 1
  expect_equal(pred, c.train, tolerance = 0.005)
})
