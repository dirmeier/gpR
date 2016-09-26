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

#' Data wrapper for Gausian process regression
#'
#' @noRd
#' @slot x.train  the training data for the explanatory variables
#' @slot y.train  the training data for the dependant variables
#' @slot x.new  the data for which the response should be predicted
#' @slot pars  parameter list
setClass(
  "lvgpr.data",
  representation(x.train="numeric", y.train="numeric",
                 x.new="numeric", pars="list"),
  validity=function(object)length(object@x.train) == length(object@y.train)
)

#' Data wrapper for Gausian process classification
#'
#' @noRd
#' @slot x.train  the training data for the explanatory variable
#' @slot c.train  the training data for the labels
#' @slot x.new  the data for which labels should be predicted
#' @slot pars  parameter list
setClass(
  "lvgpc.data",
  representation(x.train="numeric", c.train="numeric",
                 x.new="numeric", pars="list"),
  validity=function(object) all(object@c.train %in% c(0,1)) &
    (length(object@x.train) == length(object@c.train))
)

