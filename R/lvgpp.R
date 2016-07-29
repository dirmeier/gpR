#' @title Prediction of normal-distributed random variables using GP regression
#'
#' @description Prediction of normal-distributed random variables using GP regression
#' @export
#' @param x.train  vector of independent variables used for training
#' @param y.train  vector of dependent variables used for training
#' @param x.new  vector of variables for which a response should be predicted
#' @param pars  a list containing the hyper-parameters and kernel specifications
#' @param ... additional parameters (not specified)
#' @return an object of class \code{lvgpr.pred}
#' \item{y.predict }{ the predicted y* values given the \code{lvgpr} object}
#' \item{mean }{ the posterior mean values}
#' \item{cov  }{ the posterior covariance/kernel}
#' \item{call  }{ the function call}
#' @references
#'  Rasmussen C.E. and Williams C.K.I. (2006), \emph{Gaussian Processes for Machine Learning}, MIT Press \cr
#'  \url{http://www.gaussianprocess.org/gpml/} \cr \cr
#'  Barber D. (2013), \emph{Bayesian Reasoning and Machine Learning}, Cambridge University Press \cr
#'  \url{http://web4.cs.ucl.ac.uk/staff/D.Barber/pmwiki/pmwiki.php?n=Brml.HomePage}
#' @examples
#'  x.train <- seq(-5, 5, .1)
#'  y.train <- rnorm(length(x.train))
#'  x.new <- rnorm(100)
#'  pars <- list(var=1, inv.scale=2, gamma=2, noise=.1, kernel="gamma.exp")
#'  pred <- lvgpr(x.train, y.train, x.new, pars)
lvgpr <-
function(x.train=NULL, y.train=NULL, x.new=NULL,
         pars=list(var=1, inv.scale=2, gamma=2, noise=.1, kernel="gamma.exp"), ...)
  UseMethod("lvgpr")

#' @title Prediction of binomial-distributed random variables using GP classification
#'
#' @description Prediction of binomial-distributed random variables using GP classification
#' @export
#' @param x.train  vector of independent variables used for training
#' @param c.train  vector of dependent variables used for training
#' @param x.new  vector of variables for which a response should be predicted
#' @param pars  a list containing the hyper-parameters and kernel specifications
#' @param ... additional parameters (not specified)
#' @return An object of class \code{lvgpc.pred}
#' @return an \code{lvgpc.pred} object
#' \item{c.predict }{ the predicted c* values given the \code{lvgpc} object}
#' \item{mean.c.predict }{ the predicted mean c* values given the \code{lvgpc} object}
#' \item{mean }{ the (approximated) posterior mean values}
#' \item{cov  }{ the (approximated) posterior covariance/kernel}
#' \item{call  }{ the function call}
#' @references
#'  Rasmussen C.E. and Williams C.K.I. (2006), \emph{Gaussian Processes for Machine Learning}, MIT Press \cr
#'  \url{http://www.gaussianprocess.org/gpml/} \cr \cr
#'  Barber D. (2013), \emph{Bayesian Reasoning and Machine Learning}, Cambridge University Press \cr
#'  \url{http://web4.cs.ucl.ac.uk/staff/D.Barber/pmwiki/pmwiki.php?n=Brml.HomePage}
#' @examples
#'  x.train <- seq(-5, 5, .1)
#'  c.train <- c(rep(0,25), rep(1, length(x.train)-50), rep(0,25))
#'  x.new <- rnorm(100)
#'  pars <- list(var=1, inv.scale=2, gamma=2, noise=.1, kernel="gamma.exp")
#'  pred <- lvgpc(x.train, c.train, x.new, pars)
lvgpc <-
function(x.train=NULL, c.train=NULL, x.new=NULL,
         pars=list(var=1, inv.scale=2, gamma=2, noise=.1, kernel="gamma.exp"),
         ...)
  UseMethod("lvgpc")

#' @export
#' @importFrom methods new
lvgpr.default <-
function(x.train=NULL, y.train=NULL, x.new=NULL,
         pars=list(var=1, inv.scale=2, gamma=2, noise=.1, kernel="gamma.exp"),
         ...)
{
  # create an regression object
  lvgpr.obj <- methods::new("lvgpr.data",
                            x.train=x.train, y.train=y.train,
                            x.new=x.new, pars=pars)
  # predict the posterior means, covariance,
  # and predicted value for y*
  pred.posterior <- predict(lvgpr.obj, ...)
  # create a list containing the predicted values
  obj <- base::list(y.predict=pred.posterior$y.predict,
              cov=pred.posterior$cov,
              mean=pred.posterior$mean)
  # add function call to object
  obj$call <- match.call()
  # cast to class "lvgpr.pred"
  class(obj) <- "lvgpr.pred"
  obj
}

#' @export
#' @importFrom methods new
lvgpc.default <-
function(x.train=NULL, c.train=NULL, x.new=NULL,
         pars=list(var=1, inv.scale=2, gamma=2, noise=.1, kernel="gamma.exp"), ...)
{
    # create an classification object
    lvgpc.obj <-  methods::new("lvgpc.data",
                               x.train=x.train, c.train=c.train,
                               x.new=x.new, pars=pars)
    # predict the posterior means, covariance,
    # and class mapping probabilities (i.e. the probability that c.new=1)
    pred.posterior <- predict(lvgpc.obj, ...)
    obj <- base::list(mean.c.predict=pred.posterior$mean.c.predict,
                c.predict=pred.posterior$c.predict,
                cov=pred.posterior$cov,
                mean=pred.posterior$mean,
                c.label=pred.posterior$c.labels)
    # add function call to object
    obj$call <- match.call()
    # cast to class "lvgpc.pred"
    class(obj) <- "lvgpc.pred"
    obj
}
