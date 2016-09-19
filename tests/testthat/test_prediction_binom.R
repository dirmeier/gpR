.predict.a <- function(K, c.train, sig, x.new, train, test)
{
  cov.train.train <- K[train, train]
  c.mean <- c.train - sig
  Dinv <- diag(1 / (sig * (1 - sig)))
  pred.means <- pred.samples <- c()
  for (i in seq(length(x.new)))
  {
    cov.new.train <- t(K[train, test[i]])
    cov.new.new <- K[test[i], test[i]]
    m <- cov.new.train %*% c.mean
    var <- cov.new.new - cov.new.train %*%
    solve(cov.train.train + Dinv) %*% t(cov.new.train)
    me <- .logistic.normal.integral(m, var)
    pre <- .sigmoid(.sample.gaussian(m, var))
    pred.samples <- c(pred.samples , pre)
    pred.means   <- c(pred.means, me)
  }
  m.new <- K[test, train] %*% c.mean
  cov.new.new <- K[test, test] - K[test, train] %*%
    solve(cov.train.train + Dinv) %*% K[train, test]

}

(na, ntrain, nnew, K, ctrain, sigtrain, xnew, trainidx, newidx, mnew, ps, pm, Kn)

#' @noRd
.predict.binomial <- function(K, c.train, sig, x.new, train, test)
{
  ntrain <- length(train)
  nnew   <- length(test)
  n <- ntrain + nnew
  me <- pred.sample <- pred.means <- rep(0, nnew)
  K.new.new <- matrix(0, nnew, nnew)
  res <- .Fortran("predict_binomial",
                  n=as.integer(n),
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
}
