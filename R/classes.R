#' @title Description for class \code{lvgpr}
#'
#' @noRd
new.lvgpr.obj <-
function(x.train, y.train, x.new, pars)
{
    obj <- list(x.train=x.train, y.train=y.train, x.new=x.new, pars=pars)
    class(obj) <- "lvgpr"
    obj
}

#' @title Description for class \code{lvgpc}
#'
#' @noRd
new.lvgpc.obj <-
function(x.train, c.train, x.new, pars)
{
  obj <- list(x.train=x.train, c.train=c.train, x.new=x.new, pars=pars)
  class(obj) <- "lvgpc"
  obj
}
