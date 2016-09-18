#' @title Covariance function
#'
#' @noRd
#' @description Calculates the covariance/kernel of a set of input points
#' @param x1  a vector of random real numbers
#' @param x2  a vector of random real numbers
#' @param pars  a list containing hyperparameters and kernel specification
#' @return the calculated kernel/covariance matrix of size (\code{length(x1)} x \code{length(x2)})
covariance.function <- function(x1, x2, pars)
{
  fun <-  switch(pars$kernel, gamma.exp = .gamma.exp, NULL)
  if (is.null(fun))
    stop("Wrong kernel provided! Pls check method documentation for details.")
  K <- matrix(outer(x1, x2, FUN = function(x, y) fun(x, y, pars)),
              length(x1), length(x2))
  if (length(pars) == 5 && length(x1) == length(x2))
    K <- K + pars$noise * diag(length(x1))
  invisible(K)
}

#' @title Gamma-exponential kernel function
#'
#' @noRd
#' @examples
#' .gamma.exp(1:5, 1:5, pars=list(inv.scale=1, gamma=1))
.gamma.exp <-
function(x1, x2, pars)
{
  d <- x1 - x2
  pars$var * exp(- (1/(pars$inv.scale)) * (d**pars$gamma))
}
