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
  res <- .predict.binomial(K, c.train, sig, x.new, train, test)
  c.labels <- rep(0, length(res$pred.samples))
  c.labels[res$pred.samples >= 0.5] <- 1
  ret <- base::list(mean.c.predict=res$pred.means,
                    c.predict=res$pred.samples,
                    c.labels = c.labels,
                    cov=res$cov.new.new,
                    mean=res$m.new)
  ret$call <- match.call()
  ret
}

#' @noRd
.predict.binomial <- function(K, c.train, sig, x.new, train, test)
{
  ntrain <- length(train)
  nnew   <- length(test)
  n      <- ntrain + nnew
  me <- pred.sample <- pred.means <- rep(0, nnew)
  K.new.new <- matrix(0, nnew, nnew)
  res <- .Fortran("predict_binomial",
                  na=as.integer(n),
                  ntrain=as.integer(ntrain),
                  nnew=as.integer(nnew),
                  K=K,
                  ctrain=c.train,
                  sigtrain=sig,
                  xnew=x.new,
                  trainidx=as.integer(train),
                  newidx=as.integer(test),
                  mnew=me,
                  ps=pred.sample,
                  pm=pred.means,
                  Kn=K.new.new,
                  PACKAGE="gpR")
  list(m.new=res$mnew,
       cov.new.new=res$Kn,
       pred.samples=res$ps,
       pred.means=res$pm)
}
