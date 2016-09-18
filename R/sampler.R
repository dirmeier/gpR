#' Sample function values from a Gaussian process
#'
#' @noRd
#' @importFrom stats rnorm
#' @param x  a vector of real values
#' @param m  a vector of mean values, e.g. a vector of zeroes
#' @param K  a symmetric positive semi-definite matrix of size \code{length(x)} x \code{length(x)}, e.g. calculated using \code{convariance.function}
#' @return a vector of normal-distributed values
#' @examples
#'  sample.from.gp()
sample.from.gp <- function(x=seq(-10, 10, .25), m=NULL, K=NULL)
{
  xlen <- length(x)
  if (is.null(m)) m <- rep(0)
  if (is.null(K))
    K <- covariance.function(x, x,
                             list(var=1, inv.scale=2, gamma=2,
                                  noise=0, kernel="gamma.exp"))
  # small pseudo-count for numerical stability
  C <- 1e-5 * diag(xlen)
  # standard transformation of Gaussian random variables: X <- \mu + \sigma Z
  # the cholesky decomposition of a matrix K is its square root
  # the square root of the covariance yields the standard deviation
  # K has to be symmetric and semi-definite.
  # Since this is a requirement for a kernel-function,
  # we are safe to use chol().
  # thus the 'sampling' from a Gaussian process is just generating
  # xlen random variables and transforming them accordingly
  y <- m + t(chol(K + C)) %*% stats::rnorm(xlen)
  as.vector(y)
}

#' Sample a Gaussian random variable
#'
#' @noRd
#' @param m  the mean of the Gausian
#' @param var  the variance of the Gaussian
#' @return a sample from a Gaussian with mean \code{m} and variance \code{var}
.sample.gaussian <- function(m, var) m + var * rnorm(1)

#' Sample class label probabilities from a logistic mapping of a Gaussian process
#'
#' @noRd
#' @param x  a vector of real values
#' @param m  a vector of mean values, e.g. a vector of zeroes
#' @param K  a symmetric positive semi-definite matrix of size \code{length(x)} x \code{length(x)}, e.g. calculated using \code{convariance.function}
#' @return a vector of class mapping probabilities
#' @examples
#'  x=seq(-10, 10, .25)
#'  m=rep(0, length(x))
#'  K=covariance.function(x, x, list(var=1, inv.scale=2, gamma=2, noise=0, kernel="gamma.exp"))
#'  .sample.from.sigmoid.gp(x, m, K)
.sample.from.sigmoid.gp <-
function(x=seq(-10, 10, .25),m=NULL, K=NULL)
{
  xlen <- length(x)
  if (is.null(m)) m <- rep(0, xlen)
  if (is.null(K))
    K <- covariance.function(
      x, x, list(var=1, inv.scale=2, gamma=2, noise=0, kernel="gamma.exp"))
  y.latent <- sample.from.gp(x, m, K)
  y <- .sigmoid(y.latent)
  as.vector(y)
}
