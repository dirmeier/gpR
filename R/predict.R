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

#' @noRd
setGeneric("predict", function(obj, ...) standardGeneric("predict"))

#' @noRd
setMethod("predict", signature(obj = "lvgpr.data"), function(obj, ...) {
  predict.gaussian(obj, ...)
})

#' @noRd
setMethod("predict", signature(obj = "lvgpc.data"), function(obj, ...) {
  predict.binomial(obj, ...)
})
