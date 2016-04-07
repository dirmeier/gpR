#' @title Approximate a non-Gaussian distribution to a Gaussian normal distribution using Laplace's method using Newton updates
#'
#' @noRd
#' @param object  an \code{lvgpc} object
#' @return a list containing the MAP estimates for \code{y}, predicted labels \code{c}, the covariance \code{cov} for train and test inputs and the hessian \code{D} of the likelihood
#' \item{map }{ the approximated MAP-estimated of the \code{y} values (also the mean since it is Gaussian)}
#' \item{log.transform }{ the logistic transformation of the \code{map} values}
#' \item{cov  }{ the joint covariance matrix of the train and test inputs}
#' \item{D }{ the Hessian matrix of the log-likelihood p(y|f)}
#' @references
#'  Rasmussen C.E. and Williams C.K.I. (2006), \emph{Gaussian Processes for Machine Learning}, MIT Press \cr
#'  \url{http://www.gaussianprocess.org/gpml/} \cr \cr
#'  Barber D. (2015), \emph{Bayesian Reasoning and Machine Learning}, Cambridge University Press \cr
#'  \url{http://web4.cs.ucl.ac.uk/staff/D.Barber/pmwiki/pmwiki.php?n=Brml.Online}
approx.posterior <-
function(object)
{
  if (class(object) != "lvgpc") stop("Please provide an lvgpc object!")
  x.train <- object$x.train
  c.train <- object$c.train
  x.new   <- object$x.new
  x <- c(x.train, x.new)
  # kernel hyperparameters and kernel function (i.e. the covariance)
  pars <- object$pars
  n <- length(x.train)
  train <- 1:n
  test <- 1:length(x.new) + n
  # create the covariance
  K <- covariance.function(x1=x, x2=x, pars=pars)
  diagn <- diag(n)
  # posterior and logistic transform values (estimated in next lines)
  y <- sig <- rep(0, n)
  # Hessian of log-likelihood (estimated in the next lines)
  D <- diag(sig)
  # covariance of training points
  K.train <- K[train, train]
  y.old <- Inf
  niter <- 1000
  iter <- 1
  # estimate the mean of the posterior
  # since in a Gaussian distribution the mean equals the mode, we can apply a Laplace approximation :) (see references: Barber p577, Section~2.82)
  # the update is done using Newton's method, where y.new = y.old - (\nabla \nabla \Psi)^{-1} \nabla \Psi (see references: Barber p405, Equation~19.5.14; Rasmussen p43, Equation~3.18)
  # \Psi is the unnormalized log posterior GP = log p(y|f) + log(f)
  while (mean(abs(y - y.old)) > 10e-5 && iter <= niter)
  {
    # compute class mapping for y
    sig <- as.vector(sigmoid(y))
    # compute Hessian of log-likelihood (ONLY IN THIS CASE WHERE THE LOG-TRANSFORM IS USED!)
    D <- diag(sig*(1-sig))
    y.old <- y
    # compute Newton update (see references: Barber p406, Equation~19.5.19; Rasmussen p43, Equation~3.18)
    y <- as.vector(K.train %*% solve(diagn + D %*% K.train) %*% (D%*%y + c.train  - sig))
    iter <- iter + 1
  }
  list(map=y,             # mode/mean/map-estimator of posterior
       log.transform=sig, # map class-label prediction
       cov=K,             # joint covariance matrix
       D=D)               # Hessian of log-likelihood
}

#' Approximation of logistic-normal integral using the error function
#'
#' @import pracma
#' @noRd
#' @references
#'  \url{http://threeplusone.com/Crooks-LogisticNormal.pdf}
logistic.normal.integral <-
function(m, sig)
{
  lam <- sqrt(pi)/4
  0.5 + 0.5 * pracma::erf(lam * m / sqrt(1 + 2*lam**2*sig))
}
