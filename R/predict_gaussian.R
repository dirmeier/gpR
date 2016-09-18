#' @title Predict a new set of data points x.new using Gaussian process regression
#'
#' @noRd
#' @param object  an \code{lvgpr.data} object
#' @return an \code{lvgpr.pred} object
#' \item{y.predict }{ the predicted y* values given the \code{lvgpr} object}
#' \item{mean }{ the posterior mean values}
#' \item{cov  }{ the posterior covariance/kernel}
#' @references
#'  Rasmussen C.E. and Williams C.K.I. (2006), \emph{Gaussian Processes for Machine Learning}, MIT Press \cr
#'  \url{http://www.gaussianprocess.org/gpml/} \cr \cr
#'  Barber D. (2015), \emph{Bayesian Reasoning and Machine Learning}, Cambridge University Press \cr
#'  \url{http://web4.cs.ucl.ac.uk/staff/D.Barber/pmwiki/pmwiki.php?n=Brml.Online}
predict.gaussian <- function(obj, ...)
{
  x.train <- obj@x.train
  y.train <- obj@y.train
  x.new   <- obj@x.new
  pars    <- obj@pars
  # small pseudo-count for numerical stability
  C <- 1e-5 * diag(length(x.train))
  # covariance of training inputs and new inputs
  cov.new.train <- covariance.function(x1=x.new, x2=x.train, pars=pars)
  # inverse covariance of training points
  cov.inv.train.train <- solve(C + covariance.function(x1=x.train,
                                                       x2=x.train,
                                                       pars=pars))
  # posterior mean value (see references: Barber p396, Equation~19.2.5;
  #                                       Rasmussen p16, Equation~2.23)
  m <- cov.new.train %*% cov.inv.train.train %*% y.train
  # posterior covariance of new inputs
  #  (see references: Barber p396, Equation~19.2.5;
  #                   Rasmussen p16, Equation~2.24)
  cov.new.new <- covariance.function(x1=x.new, x2=x.new, pars=pars) -
    cov.new.train %*% cov.inv.train.train %*% t(cov.new.train)
  # sample from the posterior Gaussian process to get y*
  y   <- sample.from.gp(x=x.new, m=m, K=cov.new.new)
  ret <- list(y.predict=y, cov=cov.new.new, mean=m)
  ret$call <- match.call()
  class(ret) <- "lvgpr.pred"
  ret
}
