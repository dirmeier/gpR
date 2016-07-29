#' Data wrapper for Gausian process regression
#'
#' @noRd
setClass(
  "lvgpr.data",
  representation(x.train="numeric", y.train="numeric",x.new="numeric", pars="list"),
  validity=function(object)length(object@x.train) == length(object@y.train)
)

#' Data wrapper for Gausian process classification
#'
#' @noRd
setClass(
  "lvgpc.data",
  representation(x.train="numeric", c.train="numeric",x.new="numeric", pars="list"),
  validity=function(object) all(object@c.train %in% c(0,1)) &
    (length(object@x.train) == length(object@c.train))
)

