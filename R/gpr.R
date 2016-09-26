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

#' @title Prediction of normal-distributed random variables using GP regression
#'
#' @description Prediction of normal-distributed random variables using GP regression
#' @export
#' @param x.train  vector of independent variables used for training
#' @param y.train  vector of dependent variables used for training
#' @param x.new  vector of variables for which a response should be predicted
#' @param pars  a list containing the hyper-parameters and kernel specifications
#' @param ... additional parameters (not specified)
#' @return an object of class \code{lvgpr.pred}
#' \item{y.predict }{ the predicted y* values given the \code{lvgpr} object}
#' \item{mean }{ the posterior mean values}
#' \item{cov  }{ the posterior covariance/kernel}
#' \item{call  }{ the function call}
#' @references
#'  Rasmussen C.E. and Williams C.K.I. (2006), \emph{Gaussian Processes for Machine Learning}, MIT Press \cr
#'  \url{http://www.gaussianprocess.org/gpml/} \cr \cr
#'  Barber D. (2013), \emph{Bayesian Reasoning and Machine Learning}, Cambridge University Press \cr
#'  \url{http://web4.cs.ucl.ac.uk/staff/D.Barber/pmwiki/pmwiki.php?n=Brml.HomePage}
#' @examples
#'  x.train <- seq(-5, 5, .1)
#'  y.train <- rnorm(length(x.train))
#'  x.new <- rnorm(100)
#'  pars <- list(var=1, inv.scale=2, gamma=2, noise=.1, kernel="gamma.exp")
#'  pred <- lvgpr(x.train, y.train, x.new, pars)
lvgpr <- function(x.train=NULL, y.train=NULL, x.new=NULL,
                  pars=list(var=1, inv.scale=2, gamma=2, noise=.1,
                            kernel="gamma.exp"),
                  ...)
{
  UseMethod("lvgpr")
}

#' @export
#' @importFrom methods new
lvgpr.default <- function(x.train=NULL, y.train=NULL, x.new=NULL,
                          pars=list(var=1, inv.scale=2, gamma=2, noise=.1,
                                    kernel="gamma.exp"),
                          ...)
{
  # create an regression object
  lvgpr.obj <- methods::new("lvgpr.data",
                            x.train=x.train, y.train=y.train,
                            x.new=x.new, pars=pars)
  # predict the posterior means, covariance,
  # and predicted value for y*
  pred.posterior <- predict(lvgpr.obj, ...)
  # create a list containing the predicted values
  obj <- base::list(y.predict=pred.posterior$y.predict,
              cov=pred.posterior$cov,
              mean=pred.posterior$mean)
  # add function call to object
  obj$call <- match.call()
  # cast to class "lvgpr.pred"
  class(obj) <- "lvgpr.pred"
  obj
}
