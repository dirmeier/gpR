#' Data wrapper for Gausian process regression
#'
#' @slot x.train  the training data for the explanatory variables
#' @slot y.train  the training data for the dependant variables
#' @slot x.new  the data for which the response should be predicted
#' @slot pars  parameter list
setClass(
  "lvgpr.data",
  representation(x.train="numeric", y.train="numeric",x.new="numeric", pars="list"),
  validity=function(object)length(object@x.train) == length(object@y.train)
)

#' Data wrapper for Gausian process classification
#'
#' @slot x.train  the training data for the explanatory variable
#' @slot c.train  the training data for the labels
#' @slot x.new  the data for which labels should be predicted
#' @slot pars  parameter list
setClass(
  "lvgpc.data",
  representation(x.train="numeric", c.train="numeric",x.new="numeric", pars="list"),
  validity=function(object) all(object@c.train %in% c(0,1)) &
    (length(object@x.train) == length(object@c.train))
)

