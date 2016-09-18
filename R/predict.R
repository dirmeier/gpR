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
