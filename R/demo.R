#' @title Demonstration of GPR
#'
#' @description Run a GPR example. Generate samples using the prior GP, then learn the posterior GP and predict a new set of data points.
#   Plot both some prior samples and posterior samples and the confidence interval of the posterior.
#' @export
#' @importFrom graphics par plot lines polygon points abline
#' @importFrom methods new
#' @importFrom stats rnorm
#' @examples
#'  demo.regression()
demo.regression <-
function()
{
  # hyperparameters and kernel specification
  #  (these should be calculated by maximizing the marg. likelihood)
  params <- base::list(var=1, inv.scale=2, gamma=2, noise=0, kernel='gamma.exp')
  x <- y <- base::seq(-5, 5, .2)
  graphics::par(mfrow=c(2,1))
  # create empty plot
  graphics::plot(type="n", x, y, col=1, main = "Prior GP",
                 ylim=c(-3,3), xlim=c(-5,5))
  # sample three times from prior GP
  for (i in 1:3)
    graphics::lines(x,sample.from.gp(x, base::rep(0, base::length(x)),
                                     covariance.function(x, x, params)),
                    col=(i+1))
  # create training set (just use the sampled y values as true values.
  # This is no issue, since we set the noise to zero)
  x.train <- base::c(-4:-1, 1)
  y.train <- sample.from.gp(x.train,
                            base::rep(0, base::length(x.train)),
                            covariance.function(x.train, x.train, params))
  # create a new set of covariables
  # btw.: this can also be higher dimensional.
  # the one-dimensionality is only chosen for convencience.
  # for higher dimensions the mean and cov functions have to be
  #  altered to produce a mean vector and a gram matrix
  x.new <- seq(-5, 5, .2)
  # sample from the posterior, which is basically an update to the
  #  previous mean function and covariance function
  # plot training points
  graphics::plot(type="p", x.train, y.train, col=1,
                 main = "Posterior GP",
                 ylim=base::c(-3,3), xlim=base::c(-5,5),
                 xlab="x", ylab="y")
  for (i in 1:3)
  {
    # make the prediction
    pred <- lvgpr(x.train=x.train, y.train=y.train, x.new=x.new, pars=params)
    graphics::lines(x.new, pred$y.predict, col=(i+1))
    # plot confidence intervals
    if (i==1)
    {
      conf.lower <- pred$mean - 1.96 * sqrt(diag(pred$cov))
      conf.upper <- pred$mean + 1.96 * sqrt(diag(pred$cov))
      graphics::polygon(base::c((x.new), base::rev(x.new)),
                        base::c(conf.lower, base::rev(conf.upper)))
    }
  }
}

#' @title Demonstration of binary GPC
#'
#' @description Run a GPC example.
#'  Generate class labels using the sigmoid prior GP,
#'  then learn the posterior GP and predict a new set of data points.
#'  Plot both samples from the sigmoid prior and the sigmoid posterior
#'  and the mean prediction of the posterior.
#' @export
#' @importFrom graphics par plot lines polygon points abline
#' @importFrom methods new
#' @importFrom stats rnorm
#' @examples
#'  demo.bin.classification()
demo.bin.classification <- function()
{
    # hyperparameters and kernel specification
    #  (these should be calculated by maximizing the marg. likelihood)
    params <- list(var=1, inv.scale=2, gamma=2, noise=0, kernel='gamma.exp')
    graphics::par(mfrow=c(2,1))
    x <- y <- seq(-5, 5, .2)
    graphics::par(mfrow=c(2,1))
    # create empty plot
    graphics::plot(type="n", x, y,
                   col=1, main = "Prior sigmoid-GP",
                   ylim=c(0, 1),
                   xlim=c(-5,5))
    graphics::abline(h = 0.5, col = "gray60")
    # sample three times from prior sigmoid-GP
    for (i in 1:3)
    {
      c <-  .sample.from.sigmoid.gp(x, rep(0, length(x)),
                                   covariance.function(x, x, params))
      # plot all that is greater .5 as red, else blue
      cols <- ifelse(c >= 0.5, 2, 4)
      # use different point styles and plot small points to fit into the plot
      graphics::points(x, c, col=cols, pch=(i+1), cex=.5)
    }
    # create training data
    x.train <- c(seq(-5, 5, length.out=20), seq(-10, -5, length.out=20),
                 seq(5, 10, length.out=20) )
    c.train <- c(rep(1, 20), rep(0, 40))
    # create testing data
    x.new <- sort(c(stats::rnorm(20, -5),
                    stats::rnorm(20, 5),
                    stats::rnorm(20)))
    # specify label colors and point style
    cols <- base::ifelse(c.train == 1, 2, 4)
    pch  <- base::ifelse(c.train == 1, "o", "+")
    # plot training points
    graphics::plot(x.train, c.train, type="p", col=cols, cex=.5, pch=pch,
                   xlab="x" , ylab="c",
                   ylim=base::c(0,1), xlim=base::c(-10,10),
                   main = "Posterior sigmoid-GP")
    graphics:: abline(h = 0.5, col = "gray60")
    # sample from the posterior, which is basically an update to the
    #  previous mean function and covariance function
    for (i in 1:3)
    {
      pred <- lvgpc(x.train=x.train, c.train=c.train,
                    x.new=x.new, pars=params)
      c <- pred$c.predict
      cols <- base::ifelse(c >= 0.5, 2, 4)
      graphics::points(x.new, c, col=cols, pch=(i+1), cex=.5)
      if (i==1)
      {
        c.mean <- pred$mean.c.predict
        graphics::lines(x.new, c.mean)
      }
    }
}
