context("laplace approximation")

laplace.r <- function(n, K.train, c.train)
{
  y <- rep(1, n)
  y.old <- sig <- rep(0, n)
  iter <- 1
  niter <- 10000
  thresh <- 0.00001
  diagn <- diag(n)
  while (sum(abs(y - y.old)) > thresh && iter < niter)
  {
    sig <- as.vector(.sigmoid(y))
    D <- diag(sig*(1-sig))
    y.old <- y
    y <- base::as.vector(K.train %*% solve(diagn + D %*% K.train) %*%
                           (D%*%y + c.train  - sig))
    iter <- iter + 1
  }
  list(y=y, sig=sig)
}

test_that("fortran laplace approximation", {
  n <- 100
  y <- c(rep(0, n/2), rep(1, n/2))
  diagn <- diag(n)
  pars <- list(var=1, inv.scale=2, gamma=2, noise=.1, kernel="gamma.exp")
  x <- c(rnorm(n/2), rnorm(n/2, 1))
  K.train <- covariance.function(x1=x, x2=x, pars=pars)
  c.train <- c(rep(0, n/2) , rep(1, n/2))
  D <- matrix(0, n, n)
  sig <- as.vector(.sigmoid(y))
  res.r <- laplace.r(n, K.train, c.train)
  res.f <- .Fortran("laplace_approximation",
                    n=as.integer(n),
                    c=c.train, K=K.train,
                    y=y, sig=sig, D=D,
                    PACKAGE="gpR")
  expect_equal(res.r$y, res.f$y)
  expect_equal(res.r$sig, res.f$sig)
})

test_that("fortran y posterior", {
  n <- 100
  y <- c(rep(0, n/2), rep(1, n/2))
  diagn <- diag(n)
  pars <- list(var=1, inv.scale=2, gamma=2, noise=.1, kernel="gamma.exp")
  x <- c(rnorm(n/2), rnorm(n/2, 1))
  K.train <- covariance.function(x1=x, x2=x, pars=pars)
  c.train <- c(rep(0, n/2) , rep(1, n/2))
  D <- matrix(0, n, n)
  DIAG <- diag(n)
  sig <- as.vector(.sigmoid(y))
  D <- diag(sig*(1-sig))
  res.f <- .Fortran("y_posterior",
                   n=as.integer(n),
                   c=c.train, K=K.train,
                   y=y, sig=sig, D=D,
                   DIAG=diag(n),
                   PACKAGE="gpR")$y
  res.r <- as.vector(K.train %*% solve(diagn + D %*% K.train) %*%
                       (D%*%y + c.train  - sig))
  expect_equal(res.r, res.f)
})


