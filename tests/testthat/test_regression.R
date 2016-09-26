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

context("regression")

test_that("regression predicts correctly", {
  x.train <- seq(-5, 5, .1)
  y.train <- rnorm(length(x.train))
  x.new <- x.train
  pars <- list(var=1, inv.scale=2, gamma=2, noise=.1, kernel="gamma.exp")
  pred <- lvgpr(x.train, y.train, x.new, pars)
  expect_equal(pred$y.predict, y.train, tolerance = 0.05)
})
