#' @title Logistic transfer-function
#'
#' @noRd
#' @examples
#'  .sigmoid(1)
.sigmoid <- function(x) (1 / (1 + exp(-x)))
