#' @title Covariance function
#'
#' @description Calculates the covariance/kernel of a set of input points
#' @export
#' @param x1   a vector of random real numbers
#' @param x2   a vector of random real numbers
#' @param pars  a list containing hyperparameters and kernel specification
#' @return the calculated kernel/covariance matrix of size (\code{length(x1)} x \code{length(x2)})
#' @examples
#'  x = rnorm(10)
#'  K = covariance.function(x, x, list(var=1, inv.scale=2, gamma=2, noise=0, kernel="gamma.exp"))
covariance.function <-
function(x1, x2, pars)
{
  fun <-  switch(pars$kernel, gamma.exp = gamma.exp, NULL)
  if (is.null(fun)) stop("Wrong kernel provided! Pls check method documentation for details.")
  K <- matrix(outer(x1, x2, FUN = function(x, y) fun(x, y, pars)), length(x1), length(x2))
  if (length(pars) == 5 && length(x1) == length(x2)) K <- K + pars$noise * diag(length(x1))
  invisible(K)
}

#' @title Gamma-exponential kernel function
#'
#' @description Implementation of the gamma-exponential kernel function, which is the squared-exponential kernel for pars$gamma=2
#' @noRd
gamma.exp <-
function(x1, x2, pars)
{
  d <- x1 - x2
  pars$var * exp(- (1/(pars$inv.scale)) * (d**pars$gamma))
}
