myMLESpatialProcess=function (x, y, lambda.start = 0.5, theta.start = NULL, cov.function = "stationary.cov", 
          cov.args = list(Covariance = "Matern", smoothness = 1), Distance = "rdist", 
          verbose = FALSE, optim.args = NULL, ...) 
{
  if (supportsArg(cov.function, arg = "Distance") && is.null(cov.args$Distance)) 
    cov.args = c(cov.args, list(Distance = Distance))
  if (is.null(theta.start)) {
    pairwiseD <- do.call(Distance, list(x))
    pairwiseD <- pairwiseD[col(pairwiseD) > row(pairwiseD)]
    theta.start <- median(pairwiseD)
  }
  if (is.null(optim.args)) {
    optim.args = list(method = "BFGS", control = list(fnscale = -1, 
                                                      parscale = c(0.5, 0.5), ndeps = c(0.05, 0.05)))
  }
  opt = do.call("myMKrigMLEJoint", c(list(x, y, lambda.guess = lambda.start, 
                                        cov.params.guess = list(theta = theta.start), cov.fun = cov.function, 
                                        cov.args = cov.args, optim.args = optim.args, verbose = verbose), 
                                   list(...)))
  return(opt)
}

myMKrigMLEJoint = function (x, y, weights = rep(1, nrow(x)), lambda.guess = 1, 
                          cov.params.guess = NULL, cov.fun = "stationary.cov", cov.args = NULL, 
                          Z = NULL, optim.args = NULL, find.trA.MLE = FALSE, ..., verbose = FALSE) 
{
  if (is.null(optim.args)) 
    optim.args = list(method = "BFGS", control = list(fnscale = -1, 
                                                      ndeps = rep(log(1.1), length(cov.params.guess) + 
                                                                    1), reltol = 1e-06, maxit = 15))
  supportsDistMat = supportsArg(cov.fun, "distMat")
  if (supportsDistMat) {
    Dist.fun = c(cov.args, list(...))$Distance
    Dist.args = c(cov.args, list(...))$Dist.args
    if (is.null(Dist.fun)) {
      Dist.fun = "rdist"
      if (is.null(Dist.args)) 
        Dist.args = list(compact = TRUE)
    }
    distMat = do.call(Dist.fun, c(list(x), Dist.args))
    cov.args = c(cov.args, list(distMat = distMat, onlyUpper = TRUE))
  }
  mKrig.args <- c(list(x = x, y = y, weights = weights, Z = Z, 
                       cov.fun = cov.fun), list(...))
  mKrig.args$find.trA = find.trA.MLE
  ncolSummary = 8 + length(cov.params.guess)
  summary <- matrix(NA, nrow = 1, ncol = ncolSummary)
  dimnames(summary) <- list(NULL, c("EffDf", "lnProfLike", 
                                    "GCV", "sigma.MLE", "rho.MLE", "llambda.MLE", names(cov.params.guess), 
                                    "counts eval", "counts grad"))
  lnProfileLike.max <- -Inf
  temp.fn <- function(parameters) {
    lambda = exp(parameters[1])
    if (length(parameters) > 1) {
      otherParams = as.list(exp(parameters[2:length(parameters)]))
      names(otherParams) = names(cov.params.guess)
    }
    else otherParams = NULL
    cov.args.temp = c(cov.args, otherParams)
    hold <- do.call("mKrig", c(mKrig.args, list(lambda = lambda), 
                               cov.args.temp))
    if (hold$lnProfileLike.FULL > lnProfileLike.max) {
      out <<- hold
      lnProfileLike.max = hold$lnProfileLike.FULL
    }
    hold = hold[c("rho.MLE.FULL", "sigma.MLE.FULL", "lnProfileLike.FULL")]
    temp.eval <- get("capture.evaluations")
    assign("capture.evaluations", rbind(temp.eval, c(parameters, 
                                                     unlist(hold))), envir = capture.env)
    return(hold$lnProfileLike.FULL)
  }
  capture.evaluations <- matrix(NA, ncol = 4 + length(cov.params.guess), 
                                nrow = 1, dimnames = list(NULL, c("lambda", names(cov.params.guess), 
                                                                  "rho.MLE", "sigma.MLE", "lnProfileLike.FULL")))
  capture.env <- environment()
  init.guess = log(unlist(c(lambda.guess, cov.params.guess)))
  look <- do.call(optim, c(list(par = init.guess), list(temp.fn), 
                           optim.args))
  optim.counts <- look$counts
  llambda.opt <- look$par[1]
  lambda.opt <- exp(llambda.opt)
  if (length(look$par) > 1) {
    params.opt <- exp(look$par[2:length(look$par)])
    params.opt <- as.list(params.opt)
    names(params.opt) <- names(cov.params.guess)
  }
  else params.opt = NULL
  lnLik.eval <- capture.evaluations[-1, ]
  lnLik.eval[, 1:length(look$par)] = exp(lnLik.eval[, 1:length(look$par)])
  find.trA = list(...)$find.trA
  if (is.null(find.trA) || find.trA) {
    iseed = list(...)$iseed
    NtrA = list(...)$NtrA
    if (is.null(iseed)) 
      iseed = 123
    if (is.null(NtrA)) 
      NtrA = 20
    out2 <- mKrig.trace(out, iseed, NtrA)
    out$eff.df <- out2$eff.df
    out$trA.info <- out2$trA.info
    np <- nrow(x)
    out$GCV <- (sum(out$residuals^2)/np)/(1 - out2$eff.df/np)^2
    if (NtrA < np) 
      out$GCV.info <- (sum(out$residuals^2)/np)/(1 - out2$trA.info/np)^2
    else out$GCV.info <- NA
  }
  summary[1, 1:ncolSummary] <- unlist(c(out$eff.df, out$lnProfileLike.FULL, 
                                        out$GCV, out$sigma.MLE.FULL, out$rho.MLE.FULL, llambda.opt, 
                                        params.opt, optim.counts))
  
  hess = optimHess(look$par, temp.fn, control=optim.args$control)
  
  if (verbose) {
    cat("Summary: ", 1, summary[1, 1:ncolSummary], fill = TRUE)
  }
  out = c(out, list(summary = summary, lnLik.eval = lnLik.eval, 
                    optim.obj = look, hess=hess))
  class(out) = "mKrig"
  return(out)
}
