.predict.a <- function(K, c.train, sig, x.new, train, test)
{
  # cov.train.train <- K[train, train]
  # c.mean <- c.train - sig
  # Dinv <- diag(1 / (sig * (1 - sig)))
  # pred.means <- pred.samples <- c()
  # for (i in seq(length(x.new)))
  # {
  #   cov.new.train <- t(K[train, test[i]])
  #   cov.new.new <- K[test[i], test[i]]
  #   m <- cov.new.train %*% c.mean
  #   var <- cov.new.new - cov.new.train %*%
  #   solve(cov.train.train + Dinv) %*% t(cov.new.train)
  #   me <- .logistic.normal.integral(m, var)
  #   pre <- .sigmoid(.sample.gaussian(m, var))
  #   pred.samples <- c(pred.samples , pre)
  #   pred.means   <- c(pred.means, me)
  # }
  # m.new <- K[test, train] %*% c.mean
  # cov.new.new <- K[test, test] - K[test, train] %*%
  #   solve(cov.train.train + Dinv) %*% K[train, test]
  # list(m.new=m.new, cov.new.new=cov.new.new,
  #      pred.means=pred.means, pred.samples=pred.samples)
}

#' @noRd
.predict.binomial <- function(K, c.train, sig, x.new, train, test)
{
  # me <- pred.sample <- pred.means <- rep(0, n)
  # K.new.new <- matrix(0, n, n)
  # .Fortran("predict_binomial",
  #          n=
  #            K=K,
  #          ct=c.train,
  #          sig=sig,
  #          xn=x.new,
  #          trainidx=train,
  #          testidx=test,
  #          me=me,
  #          ps=pred.sample,
  #          pm=pred.means,
  #          Kn=K.new.new,
  #          PACKAGE="gpR")
}
