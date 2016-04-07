#' @title Predict a new set of data points x.new using Gaussian process regression
#'
#' @noRd
#' @param object  an \code{lvgpr} object
#' @return an \code{lvgpr.pred} object
#' \item{y.predict }{ the predicted y* values given the \code{lvgpr} object}
#' \item{mean }{ the posterior mean values}
#' \item{cov  }{ the posterior covariance/kernel}
#' @references
#'  Rasmussen C.E. and Williams C.K.I. (2006), \emph{Gaussian Processes for Machine Learning}, MIT Press \cr
#'  \url{http://www.gaussianprocess.org/gpml/} \cr \cr
#'  Barber D. (2015), \emph{Bayesian Reasoning and Machine Learning}, Cambridge University Press \cr
#'  \url{http://web4.cs.ucl.ac.uk/staff/D.Barber/pmwiki/pmwiki.php?n=Brml.Online}
predict.lvgpr <-
function(object, ...)
{
  if (class(object) != "lvgpr") stop("Please provide an lvgpr object!")
  x.train <- object$x.train
  y.train <- object$y.train
  x.new   <- object$x.new
  pars    <- object$pars
  # small pseudo-count for numerical stability
  C <- 1e-5 * diag(length(x.train))
  # covariance of training inputs and new inputs
  cov.new.train <- covariance.function(x1=x.new, x2=x.train, pars=pars)
  # inverse covariance of training points
  cov.inv.train.train <- solve(C + covariance.function(x1=x.train, x2=x.train, pars=pars))
  # posterior mean value (see references: Barber p396, Equation~19.2.5; Rasmussen p16, Equation~2.23)
  m <- cov.new.train %*% cov.inv.train.train %*% y.train
  # posterior covariance of new inputs (see references: Barber p396, Equation~19.2.5; Rasmussen p16, Equation~2.24)
  cov.new.new <- covariance.function(x1=x.new, x2=x.new, pars=pars) - cov.new.train %*% cov.inv.train.train %*% t(cov.new.train)
  # sample from the posterior Gaussian process to get y*
  y <- sample.from.gp(x=x.new, m=m, K=cov.new.new)
  obj <- list(y.predict=y,
              cov=cov.new.new,
              mean=m)
  class(obj) <- "lvgpr.pred"
  obj
}

#' @title Predict a new set of data points x.new using Gaussian process classification
#'
#' @noRd
#' @param object  an \code{lvgpc} object
#' @return an \code{lvgpc.pred} object
#' \item{c.predict }{ the predicted c* values given the \code{lvgpc} object}
#' \item{mean.c.predict }{ the predicted mean c* values given the \code{lvgpc} object}
#' \item{mean }{ the (approximated) posterior mean values}
#' \item{cov  }{ the (approximated) posterior covariance/kernel}
#' @references
#'  Rasmussen C.E. and Williams C.K.I. (2006), \emph{Gaussian Processes for Machine Learning}, MIT Press \cr
#'  \url{http://www.gaussianprocess.org/gpml/} \cr \cr
#'  Barber D. (2015), \emph{Bayesian Reasoning and Machine Learning}, Cambridge University Press \cr
#'  \url{http://web4.cs.ucl.ac.uk/staff/D.Barber/pmwiki/pmwiki.php?n=Brml.Online}
predict.lvgpc <-
function(object, ...)
{
  if (class(object) != "lvgpc") stop("Please provide an lvgpc object!")
  x.train <- object$x.train
  c.train <- object$c.train
  x.new   <- object$x.new
  n <- length(x.train)
  train <- 1:n
  test <- 1:length(x.new) + n
  # approximate the mode/mean of the posterior, its covariance, the estimated class mappings and the Hessian of the log-likelihood
  approx <- approx.posterior(object)
  # logistic transform of posterior mean values
  sig <- approx$log.transform
  # joint covariance matrix of x.train and x.new
  K <- approx$cov
  cov.train.train <- K[train, train]
  # compute log probability of the class mapping of c.train given the posterior mean (see references: Barber p406, Equation~19.5.24; Rasmussen p44, Equation~3.15 & Equation~3.21)
  c.mean <- c.train - sig
  # compute the inverse of the Hessian of the log-likelihood
  Dinv <- diag(1/(sig*(1-sig)))
  pred.means <- pred.samples <- c()
  # make a prediction for all new input points
  for (i in seq(length(x.new)))
  {
    # get the covariance of the training points and a new input
    cov.new.train <- t(K[train, test[i]])
    # get the variance of the new input point
    cov.new.new <- K[test[i], test[i]]
    # compute the posterior mean of the predictive distribution (see references: Barber p406, Equation~19.5.24; Rasmussen p44, Equation~3.21)
    m <- cov.new.train %*% c.mean
    # compute the posterior variance of the predictive distribution (see references: Barber p406, Equation~19.5.26; Rasmussen p44, Equation~3.24)
    var <- cov.new.new - cov.new.train %*% solve(cov.train.train + Dinv) %*% t(cov.new.train)
    # compute the mean class probability mapping (see references: Barber p406, Equation~19.5.27; Rasmussen p44, Equation~3.25)
    me <- logistic.normal.integral(m, var)
    # compute a sample from the posterior class mapping
    pre <- sigmoid(sample.gaussian(m, var))
    # store the predicted class mapping in a vector
    pred.samples <- c(pred.samples , pre)
    pred.means <- c(pred.means, me)
  }
  m.new <- K[test, train] %*% c.mean
  cov.new.new <- K[test, test] - K[test, train] %*% solve(cov.train.train + Dinv) %*% K[train, test]
  c.labels <- rep(0, length(pred.samples))
  c.labels[pred.samples >= 0.5] <- 1
  obj <- list(mean.c.predict=pred.means,
              c.predict=pred.samples,
              c.labels = c.labels,
              cov=cov.new.new,
              mean=m.new)
  obj$call <- match.call()
  obj
}
