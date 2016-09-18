#' @title Predict a new set of data points x.new using Gaussian process classification
#'
#' @noRd
#' @param object  an \code{lvgpc.data} object
#' @return an \code{lvgpc.pred} object
#' \item{c.predict }{ the predicted c* values given the \code{lvgpc.data} object}
#' \item{c.labels }{ the predicted c* labale given the \code{lvgpc.data} object}
#' \item{mean.c.predict }{ the predicted mean c* values given the \code{lvgpc.data} object}
#' \item{mean }{ the (approximated) posterior mean values}
#' \item{cov  }{ the (approximated) posterior covariance/kernel}
#' @references
#'  Rasmussen C.E. and Williams C.K.I. (2006), \emph{Gaussian Processes for Machine Learning}, MIT Press \cr
#'  \url{http://www.gaussianprocess.org/gpml/} \cr \cr
#'  Barber D. (2015), \emph{Bayesian Reasoning and Machine Learning}, Cambridge University Press \cr
#'  \url{http://web4.cs.ucl.ac.uk/staff/D.Barber/pmwiki/pmwiki.php?n=Brml.Online}
predict.binomial <- function(obj, ...)
{
  x.train <- obj@x.train
  c.train <- obj@c.train
  x.new   <- obj@x.new
  n <- length(x.train)
  train <- 1:n
  test <- 1:length(x.new) + n
  # approximate the mode/mean of the posterior, its covariance,
  #  the estimated class mappings and the Hessian of the log-likelihood
  approx <- approx.posterior(obj, ...)
  # logistic transform of posterior mean values
  sig <- approx$log.transform
  # joint covariance matrix of x.train and x.new
  K <- approx$cov
  cov.train.train <- K[train, train]
  # compute log probability of the class mapping of c.train given
  # the posterior mean
  # (see references: Barber p406, Equation~19.5.24;
  #                  Rasmussen p44, Equation~3.15 & Equation~3.21)
  c.mean <- c.train - sig
  # compute the inverse of the Hessian of the log-likelihood
  Dinv <- diag(1 / (sig * (1 - sig)))
  pred.means <- pred.samples <- c()
  # make a prediction for all new input points
  for (i in seq(length(x.new)))
  {
    # get the covariance of the training points and a new input
    cov.new.train <- t(K[train, test[i]])
    # get the variance of the new input point
    cov.new.new <- K[test[i], test[i]]
    # compute the posterior mean of the predictive distribution
    # (see references: Barber p406, Equation~19.5.24;
    #                  Rasmussen p44, Equation~3.21)
    m <- cov.new.train %*% c.mean
    # compute the posterior variance of the predictive distribution
    # (see references: Barber p406, Equation~19.5.26;
    #                  Rasmussen p44, Equation~3.24)
    var <- cov.new.new - cov.new.train %*%
      solve(cov.train.train + Dinv) %*% t(cov.new.train)
    # compute the mean class probability mapping
    # (see references: Barber p406, Equation~19.5.27;
    #                  Rasmussen p44, Equation~3.25)
    me <- .logistic.normal.integral(m, var)
    # compute a sample from the posterior class mapping
    pre <- .sigmoid(.sample.gaussian(m, var))
    # store the predicted class mapping in a vector
    pred.samples <- base::c(pred.samples , pre)
    pred.means   <- base::c(pred.means, me)
  }
  m.new <- K[test, train] %*% c.mean
  cov.new.new <- K[test, test] - K[test, train] %*%
    solve(cov.train.train + Dinv) %*% K[train, test]
  c.labels <- rep(0, length(pred.samples))
  c.labels[pred.samples >= 0.5] <- 1
  ret <- list(mean.c.predict=pred.means,
              c.predict=pred.samples,
              c.labels = c.labels,
              cov=cov.new.new,
              mean=m.new)
  ret$call <- match.call()
  class(ret) <- "lvgpc.pred"
  ret
}
