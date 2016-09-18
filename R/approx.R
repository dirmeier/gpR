#' @include classes.R
NULL

#' @title Approximate a non-Gaussian distribution to a Gaussian normal
#'  distribution using Laplace's method using Newton updates
#'
#' @noRd
#' @param object  an \code{lvgpc.data} object
#' @return a list containing the MAP estimates for \code{y}, predicted labels \code{c},
#'  the covariance \code{cov} for train and test inputs and the hessian \code{D} of the likelihood
#' \item{map }{ the approximated MAP-estimated of the \code{y} values (also the mean since it is Gaussian)}
#' \item{log.transform }{ the logistic transformation of the \code{map} values}
#' \item{cov  }{ the joint covariance matrix of the train and test inputs}
#' \item{D }{ the Hessian matrix of the log-likelihood p(y|f)}
#' @references
#'  Rasmussen C.E. and Williams C.K.I. (2006), \emph{Gaussian Processes for Machine Learning}, MIT Press \cr
#'  \url{http://www.gaussianprocess.org/gpml/} \cr \cr
#'  Barber D. (2015), \emph{Bayesian Reasoning and Machine Learning}, Cambridge University Press \cr
#'  \url{http://web4.cs.ucl.ac.uk/staff/D.Barber/pmwiki/pmwiki.php?n=Brml.Online}
setGeneric("approx.posterior", function(obj, ...)
  standardGeneric("approx.posterior"))

#' @noRd
setMethod("approx.posterior", signature(obj = "lvgpc.data"),
          function(obj, ...) {
            .approx.posterior(obj, ...)
          }
)

#' @noRd
.approx.posterior <- function(obj, ...)
{
    x.train <- obj@x.train
    c.train <- obj@c.train
    x.new   <- obj@x.new
    x <- c(x.train, x.new)
    # kernel hyperparameters and kernel function (i.e. the covariance)
    pars <- obj@pars
    n <- length(x.train)
    train <- 1:n
    # create the covariance
    K <- covariance.function(x1=x, x2=x, pars=pars)
    # posterior and logistic transform values (estimated in next lines)
    y <- sig <- rep(0, n)
    # Hessian of log-likelihood (estimated in the next lines)
    D <- diag(sig)
    # covariance of training points
    K.train <- K[train, train]
    # estimate the mean of the posterior
    # since in a Gaussian distribution the mean equals the mode,
    # we can apply a Laplace approximation.
    # (see references: Barber p577, Section~2.82)
    # the update is done using Newton's method,
    #  where y.new = y.old - (\nabla \nabla \Psi)^{-1} \nabla \Psi
    # (see references: Barber p405, Equation~19.5.14;
    #                  Rasmussen p43, Equation~3.18)
    # \Psi is the unnormalized log posterior GP = log p(y|f) + log(f)
    res <- .laplace.approximation(n, c.train, K.train, y, sig, D)
    list(map=res$y,             # mode/mean/map-estimator of posterior
         log.transform=res$sig, # map class-label prediction
         cov=K,             # joint covariance matrix
         D=res$D)               # Hessian of log-likelihood
}

#' @noRd
.laplace.approximation <- function(n, c.train, K.train, y, sig, D)
{
  .Fortran("laplace_approximation",
            n=as.integer(n),
            c=c.train,
            K=K.train,
            y=y,
            sig=sig,
            D=D,
            PACKAGE="gpR")
}

#' Approximation of logistic-normal integral using the error function
#'
#' @noRd
#' @import pracma
#' @references
#'  \url{http://threeplusone.com/Crooks-LogisticNormal.pdf}
#'  @examples
#'  .logistic.normal.integral(0, 1)
.logistic.normal.integral <- function(m, sig)
{
  lam <- sqrt(pi)/4
  0.5 + 0.5 * pracma::erf(lam * m / sqrt(1 + 2 * lam**2 * sig))
}
