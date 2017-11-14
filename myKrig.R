myKrig = function (x, Y, cov.function = "stationary.cov", lambda = NA, 
                   df = NA, GCV = FALSE, Z = NULL, cost = 1, knots = NA, weights = NULL, 
                   m = 2, nstep.cv = 200, scale.type = "user", x.center = rep(0, 
                                                                              ncol(x)), x.scale = rep(1, ncol(x)), rho = NA, sigma2 = NA, 
                   method = "REML", verbose = FALSE, mean.obj = NA, sd.obj = NA, 
                   null.function = "Krig.null.function", wght.function = NULL, 
                   offset = 0, na.rm = TRUE, cov.args = NULL, chol.args = NULL, 
                   null.args = NULL, wght.args = NULL, W = NULL, give.warnings = TRUE, 
                   lambda.grid=NA, ...) 
{
  out <- list()
  out$call <- match.call()
  if (options()$warn < 0) {
    give.warnings <- FALSE
  }
  if (!is.character(cov.function)) {
    out$cov.function.name <- as.character(substitute(cov.function))
  }
  else {
    out$cov.function.name <- cov.function
  }
  out$null.function.name <- as.character(substitute(null.function))
  if (is.null(wght.function)) {
    out$wght.function.name <- NULL
  }
  else {
    out$wght.function.name <- as.character(substitute(wght.function))
  }
  out$W <- W
  if (verbose) {
    print(out$cov.function.name)
    print(out$null.function.name)
    print(out$wght.function.name)
  }
  C.arg.missing <- all(names(formals(get(out$cov.function.name))) != 
                         "C")
  if (C.arg.missing) 
    stop("Need to have C argument in covariance function\nsee Exp.cov.simple as an example")
  if (!is.null(cov.args)) 
    out$args <- c(cov.args, list(...))
  else out$args <- list(...)
  out$null.args <- null.args
  if (out$null.function.name == "Krig.null.function") {
    out$null.args <- list(m = m)
    out$m <- m
  }
  if (is.null(chol.args)) {
    out$chol.args <- list(pivot = FALSE)
  }
  else {
    out$chol.args <- chol.args
  }
  out$wght.args <- wght.args
  out$offset <- offset
  out$cost <- cost
  out$lambda <- lambda
  out$eff.df <- df
  out$sigma2 <- sigma2
  out$rho <- rho
  out$method <- method
  out$GCV <- GCV
  out$mean.obj <- mean.obj
  out$sd.obj <- sd.obj
  out$correlation.model <- !(is.na(mean.obj[1]) & is.na(sd.obj[1]))
  out$scale.type <- scale.type
  out$x.center <- x.center
  out$x.scale <- x.scale
  if (verbose) {
    cat("  Cov function arguments in call  ", fill = TRUE)
    print(out$args)
    cat(" covariance function used is : ", fill = TRUE)
    print(out$cov.function.name)
  }
  out2 <- Krig.check.xY(x, Y, Z, weights, na.rm, verbose = verbose)
  out <- c(out, out2)
  if (out$correlation.model) {
    out$y <- Krig.cor.Y(out, verbose = verbose)
  }
  out2 <- Krig.transform.xY(out, knots, verbose = verbose)
  out <- c(out, out2)
  out2 <- Krig.which.lambda(out)
  out[names(out2)] <- out2
  out$nondiag.W <- (!is.null(wght.function)) | (!is.null(W))
  if (out$nondiag.W) {
    if (out$knot.model | out$fixed.model) {
      stop("Non diagonal weight matrix for observations\n                    not supported\nwith knots or fixed lambda.")
    }
    if (!is.na(out$shat.pure.error)) {
      stop("Non diagonal weight matrix not implemented\n                    with replicate locations")
    }
  }
  out <- c(out, Krig.make.W(out, verbose = verbose))
  if (out$fixed.model) {
    out$matrices <- Krig.engine.fixed(out, verbose = verbose)
    out$eff.df <- NA
  }
  if (!out$fixed.model) {
    if (out$knot.model) {
      out$matrices <- Krig.engine.knots(out, verbose = verbose)
      out$pure.ss <- out$matrices$pure.ss
    }
    else {
      out$matrices <- Krig.engine.default(out, verbose = verbose)
    }
  }
  out$nt <- out$matrices$nt
  out$np <- out$matrices$np
  out$decomp <- out$matrices$decomp
  if (is.null(out$Z)) {
    out$ind.drift <- rep(TRUE, out$nt)
  }
  else {
    mZ <- ncol(out$ZM)
    out$ind.drift <- c(rep(TRUE, out$nt - mZ), rep(FALSE, 
                                                   mZ))
  }
  if (verbose) {
    cat("null df: ", out$nt, "drift df: ", sum(out$ind.drift), 
        fill = TRUE)
  }
  if (!out$fixed.model | out$GCV) {
    if (verbose) {
      cat("call to gcv.Krig", fill = TRUE)
    }
    gcv.out <- gcv.Krig(out, nstep.cv = nstep.cv, verbose = verbose, lambda.grid = lambda.grid, 
                        cost = out$cost, offset = out$offset, give.warnings = FALSE)
    out$gcv.grid <- gcv.out$gcv.grid
    out$lambda.est <- gcv.out$lambda.est
    out$warningTable <- gcv.out$warningTable
    if (verbose) {
      cat("summaries from grid search/optimization", fill = TRUE)
      print(out$lambda.est)
      print(out$warningTable)
    }
    if (give.warnings) {
      printGCVWarnings(gcv.out$warningTable, method = method)
    }
    if (out$method != "user") {
      out$lambda <- gcv.out$lambda.est[out$method, 1]
      out$eff.df <- out$lambda.est[out$method, 2]
    }
    else {
      if (!is.na(out$eff.df)) {
        out$lambda <- Krig.df.to.lambda(out$eff.df, out$matrices$D)
      }
      else {
        out$eff.df <- Krig.ftrace(out$lambda, out$matrices$D)
      }
    }
  }
  out2 <- Krig.coef(out, yM = out$yM, verbose = verbose)
  out <- c(out, out2)
  out$fitted.values <- predict.Krig(out, x = out$x, Z = out$Z, 
                                    eval.correlation.model = FALSE)
  out$residuals <- out$y - out$fitted.values
  Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
                                               list(x = out$x, Z = out$Z)))
  out$fitted.values.null <- as.matrix(Tmatrix) %*% out$d
  if (verbose) {
    cat("residuals", out$residuals, fill = TRUE)
  }
  out2 <- Krig.parameters(out)
  out <- c(out, out2)
  passed.sigma2 <- (!is.na(out$sigma2))
  if (out$method == "user" & passed.sigma2) {
    out$best.model <- c(out$lambda, out$sigma2, out$rho)
  }
  else {
    out$best.model <- c(out$lambda, out$shat.MLE^2, out$rhohat)
  }
  class(out) <- c("Krig")
  return(out)
}



testFun = function(thetas, gpsDat) {
  for(i in 1:length(thetas)) {
    print(paste0("i: ", i))
    out = myKrig(cbind(gpsDat$lon, gpsDat$lat), gpsDat$slip, theta=thetas[i], GCV=TRUE, method="GCV", 
                 cov.args=list(Distance="rdist.earth", Dist.args=list(miles=FALSE)), m=2, 
                 lambda.grid=seq(0, 5e-2, l=200))
    bestModI = which.min(out$gcv.grid$GCV)
    lnLik = -gcv.grid$`-lnLike Prof`[bestModI]
    print(paste0("lnLik: ", lnLik))
    if(lnLik > maxLnLik) {
      maxPar = out$best.model
      maxTheta = thetas[i]
      maxLnLik = lnLik
      maxBeta = out$d
    }
  }
  5
}