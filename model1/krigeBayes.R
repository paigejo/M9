krige.bayes = function (geodata, coords = geodata$coords, data = geodata$data, 
                        locations = "no", borders, model, prior, output) 
{
  if (missing(geodata)) 
    geodata <- list(coords = coords, data = data)
  if (missing(borders)) 
    borders <- geodata$borders
  call.fc <- match.call()
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    warning(".Random.seed not initialised. Creating it with runif(1)")
    runif(1)
  }
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  locations <- .check.locations(locations)
  do.prediction <- ifelse(all(locations == "no"), FALSE, TRUE)
  if (do.prediction && any(!is.numeric(locations))) 
    stop("krige.bayes: non numeric coordinates passed to the argument \"locations\"")
  base.env <- sys.frame(sys.nframe())
  message.prediction <- character()
  phidist <- list()
  "krige.bayes.counter" <- function(.temp.ap, n.points) {
    if (n.points <= 50) 
      cat(paste(.temp.ap, ", ", sep = ""))
    if (n.points > 50 & n.points <= 500) 
      if (.temp.ap%%10 == 1) 
        cat(paste(.temp.ap, ", ", sep = ""))
    if (n.points > 500) 
      if (.temp.ap%%100 == 1) 
        cat(paste(.temp.ap, ", ", sep = ""))
    if (n.points == .temp.ap) 
      cat("\n")
  }
  kb <- list(posterior = list(beta = list(), sigmasq = list(), 
                              phi = list(), tausq.rel = list()), predictive = list(mean = NULL, 
                                                                                   variance = NULL, distribution = NULL))
  oldClass(kb$posterior) <- c("posterior.krige.bayes", "variomodel")
  oldClass(kb$predictive) <- "predictive.krige.bayes"
  pred.env <- new.env()
  if (missing(model)) 
    model <- model.control()
  else {
    if (length(class(model)) == 0 || class(model) != "model.geoR") {
      if (!is.list(model)) 
        stop("krige.bayes: the argument model only takes a list or an output of the function model.control")
      else {
        model.names <- c("trend.d", "trend.l", "cov.model", 
                         "kappa", "aniso.pars", "lambda")
        model.user <- model
        model <- list()
        if (length(model.user) > 0) {
          for (i in 1:length(model.user)) {
            n.match <- match.arg(names(model.user)[i], 
                                 model.names)
            model[[n.match]] <- model.user[[i]]
          }
        }
        if (is.null(model$trend.d)) 
          model$trend.d <- "cte"
        if (is.null(model$trend.l)) 
          model$trend.l <- "cte"
        if (is.null(model$cov.model)) 
          model$cov.model <- "matern"
        if (is.null(model$kappa)) 
          model$kappa <- 0.5
        if (is.null(model$aniso.pars)) 
          model$aniso.pars <- NULL
        if (is.null(model$lambda)) 
          model$lambda <- 1
        model <- model.control(trend.d = model$trend.d, 
                               trend.l = model$trend.l, cov.model = model$cov.model, 
                               kappa = model$kappa, aniso.pars = model$aniso.pars, 
                               lambda = model$lambda)
      }
    }
  }
  cov.model <- model$cov.model
  cov.model.number <- .cor.number(cov.model)
  kappa <- model$kappa
  lambda <- model$lambda
  if (missing(prior)) 
    prior <- prior.control()
  else {
    if (length(class(prior)) == 0 || class(prior) != "prior.geoR") {
      if (!is.list(prior)) 
        stop("krige.bayes: the argument prior only takes a list or an output of the function prior.control")
      else {
        prior.names <- c("beta.prior", "beta", "beta.var.std", 
                         "sigmasq.prior", "sigmasq", "df.sigmasq", "phi.prior", 
                         "phi", "phi.discrete", "tausq.rel.prior", "tausq.rel", 
                         "tausq.rel.discrete")
        prior.user <- prior
        prior <- list()
        if (length(prior.user) > 0) {
          for (i in 1:length(prior.user)) {
            n.match <- match.arg(names(prior.user)[i], 
                                 prior.names)
            prior[[n.match]] <- prior.user[[i]]
          }
        }
        if (is.null(prior$beta)) 
          prior$beta <- NULL
        if (is.null(prior$beta.prior)) 
          prior$beta.prior <- c("flat", "normal", "fixed")
        if (is.null(prior$beta.var.std)) 
          prior$beta.var.std <- NULL
        if (is.null(prior$sigmasq)) 
          prior$sigmasq <- NULL
        if (is.null(prior$sigmasq.prior)) 
          prior$sigmasq.prior <- c("reciprocal", "uniform", 
                                   "sc.inv.chisq", "fixed")
        if (is.null(prior$df.sigmasq)) 
          prior$df.sigmasq <- NULL
        if (is.null(prior$phi)) 
          prior$phi <- NULL
        if (is.null(prior$phi.prior)) 
          prior$phi.prior <- c("uniform", "exponential", 
                               "fixed", "squared.reciprocal", "reciprocal")
        if (is.null(prior$phi.discrete)) 
          prior$phi.discrete <- NULL
        if (is.null(prior$tausq.rel)) 
          prior$tausq.rel <- 0
        if (is.null(prior$tausq.rel.prior)) 
          prior$tausq.rel.prior <- c("fixed", "uniform", 
                                     "reciprocal")
        if (is.null(prior$tausq.rel.discrete)) 
          prior$tausq.rel.discrete <- NULL
        prior <- prior.control(beta.prior = prior$beta.prior, 
                               beta = prior$beta, beta.var.std = prior$beta.var.std, 
                               sigmasq.prior = prior$sigmasq.prior, sigmasq = prior$sigmasq, 
                               df.sigmasq = prior$df.sigmasq, phi.prior = prior$phi.prior, 
                               phi = prior$phi, phi.discrete = prior$phi.discrete, 
                               tausq.rel.prior = prior$tausq.rel.prior, tausq.rel = prior$tausq.rel, 
                               tausq.rel.discrete = prior$tausq.rel.discrete)
      }
    }
  }
  kb$prior <- prior$priors.info
  kb$model <- model
  oldClass(kb$prior) <- "prior.geoR"
  if (prior$dep.prior) {
    npr <- length(prior$sigmasq)
    nphipr <- nrow(as.matrix(prior$sigmasq))
    ntaupr <- ncol(as.matrix(prior$sigmasq))
  }
  else nphipr <- ntaupr <- npr <- 1
  beta <- prior$beta
  if (prior$beta.prior == "fixed") 
    beta.fixed <- beta
  if (prior$beta.prior == "normal") {
    nbeta <- attr(prior$beta.var.std, "Size")
    betares <- list()
    for (j in 1:ntaupr) {
      for (i in 1:nphipr) {
        beta <- array(prior$beta, dim = c(nphipr, ntaupr, 
                                          nbeta))[i, j, ]
        beta.var.std <- array(prior$beta.var.std, dim = c(nphipr, 
                                                          ntaupr, nbeta^2))[i, j, ]
        beta.var.std <- matrix(beta.var.std, nbeta, nbeta)
        ind.pos <- (j - 1) * nphipr + i
        betares[[ind.pos]] <- list(iv = .solve.geoR(beta.var.std), 
                                   ivm = drop(.solve.geoR(beta.var.std, beta)), 
                                   mivm = drop(crossprod(beta, .solve.geoR(beta.var.std, 
                                                                           beta))))
      }
    }
  }
  if (prior$sigmasq.prior != "fixed") 
    S2.prior <- prior$sigmasq
  else sigmasq.fixed <- S2.prior <- prior$sigmasq
  df.sigmasq.prior <- prior$df.sigmasq
  phi.discrete <- prior$phi.discrete
  exponential.par <- prior$phi
  tausq.rel.fixed <- tausq.rel <- prior$tausq.rel
  exponential.tausq.rel.par <- prior$tausq.rel
  if (tausq.rel.fixed > 2) 
    print("WARNING: relative (NOT absolute) nugget should be specified.")
  tausq.rel.discrete <- prior$tausq.rel.discrete
  n <- length(data)
  if (is.vector(coords)) {
    coords <- cbind(coords, 0)
    warning("krige.bayes: vector of coordinates: assuming one spatial dimension (transect)")
  }
  coords <- as.matrix(coords)
  dists.env <- new.env()
  assign("data.dist", as.vector(dist(coords)), envir = dists.env)
  data.dist.range <- range(get("data.dist", envir = dists.env))
  data.dist.min <- data.dist.range[1]
  data.dist.max <- data.dist.range[2]
  if (round(1e+12 * data.dist.min) == 0) 
    stop("krige.bayes: this function does not allow two data at same location")
  if (missing(output)) 
    output <- output.control()
  else {
    if (length(class(output)) == 0 || class(output) != "output.geoR") {
      if (!is.list(output)) 
        stop("krige.bayes: the argument output only takes a list or an output of the function output.control")
      else {
        output.names <- c("n.posterior", "n.predictive", 
                          "moments", "n.back.moments", "simulations.predictive", 
                          "mean.var", "quantile", "threshold", "signal", 
                          "messages.screen")
        output.user <- output
        output <- list()
        if (length(output.user) > 0) {
          for (i in 1:length(output.user)) {
            n.match <- match.arg(names(output.user)[i], 
                                 output.names)
            output[[n.match]] <- output.user[[i]]
          }
        }
        if (is.null(output$n.posterior)) 
          output$n.posterior <- 1000
        if (is.null(output$n.predictive)) 
          output$n.predictive <- NULL
        if (is.null(output$moments)) 
          output$moments <- TRUE
        if (is.null(output$n.back.moments)) 
          output$n.back.moments <- 1000
        if (is.null(output$simulations.predictive)) {
          if (is.null(output$n.predictive)) 
            output$simulations.predictive <- NULL
          else output$simulations.predictive <- ifelse(output$n.predictive > 
                                                         0, TRUE, FALSE)
        }
        if (is.null(output$mean.var)) 
          output$mean.var <- NULL
        if (is.null(output$quantile)) 
          output$quantile <- NULL
        if (is.null(output$threshold)) 
          output$threshold <- NULL
        if (is.null(output$sim.means)) 
          output$sim.means <- NULL
        if (is.null(output$sim.vars)) 
          output$sim.vars <- NULL
        if (is.null(output$signal)) 
          output$signal <- NULL
        if (is.null(output$messages.screen)) 
          output$messages.screen <- TRUE
        output <- output.control(n.posterior = output$n.posterior, 
                                 n.predictive = output$n.predictive, moments = output$moments, 
                                 n.back.moments = output$n.back.moments, simulations.predictive = output$simulations.predictive, 
                                 mean.var = output$mean.var, quantile = output$quantile, 
                                 threshold = output$threshold, sim.means = output$sim.means, 
                                 sim.vars = output$sim.vars, signal = output$signal, 
                                 messages = output$messages.screen)
      }
    }
  }
  n.posterior <- output$n.posterior
  messages.screen <- output$messages.screen
  if (!do.prediction) {
    if (prior$beta.prior != "fixed" & prior$sigmasq.prior != 
        "fixed" & prior$phi.prior != "fixed" & output$messages.screen) {
      cat("krige.bayes: no prediction locations provided.\n")
      cat("             Only samples of the posterior for the parameters will be returned.\n")
    }
  }
  else {
    n.predictive <- output$n.predictive
    if (is.null(n.predictive)) 
      n.predictive <- ifelse(prior$phi.prior == "fixed", 
                             0, n.posterior)
    simulations.predictive <- output$simulations.predictive
    if (is.null(simulations.predictive)) 
      simulations.predictive <- ifelse(prior$phi.prior == 
                                         "fixed", FALSE, TRUE)
    keep.simulations <- output$keep.simulations
    if (is.null(keep.simulations)) 
      keep.simulations <- simulations.predictive
    mean.estimator <- output$mean.estimator
    if (is.null(mean.estimator)) 
      mean.estimator <- ifelse(simulations.predictive, 
                               TRUE, FALSE)
    moments <- output$moments
    if (is.null(moments) | prior$phi.prior == "fixed") 
      moments <- TRUE
    n.back.moments <- output$n.back.moments
    sim.means <- output$sim.means
    if (is.null(sim.means)) 
      sim.means <- ifelse(simulations.predictive, TRUE, 
                          FALSE)
    sim.vars <- output$sim.vars
    if (is.null(sim.vars)) 
      sim.vars <- FALSE
    signal <- ifelse(is.null(output$signal), TRUE, output$signal)
    quantile.estimator <- output$quantile.estimator
    probability.estimator <- output$probability.estimator
    if (simulations.predictive & n.predictive == 0) 
      n.predictive <- 1000
  }
  if (abs(lambda - 1) > 0.001) {
    if (messages.screen) 
      cat(paste("krige.bayes: Box-Cox's transformation performed for lambda =", 
                round(lambda, digits = 3), "\n"))
    data <- BCtransform(x = data, lambda = lambda)$data
  }
  dimnames(coords) <- list(NULL, NULL)
  if (nrow(coords) != length(data)) 
    stop("krige.bayes: number of data is different of number of data locations (coordinates)")
  if (class(model$trend.d) == "trend.spatial") 
    trend.data <- unclass(model$trend.d)
  else trend.data <- unclass(trend.spatial(trend = model$trend.d, 
                                           geodata = geodata))
  if (nrow(trend.data) != nrow(coords)) 
    stop("trend specification not compatible with the length of the data")
  beta.size <- ncol(trend.data)
  if (beta.size > 1) 
    beta.names <- paste("beta", (0:(beta.size - 1)), sep = "")
  else beta.names <- "beta"
  if (prior$beta.prior == "normal" | prior$beta.prior == "fixed") {
    if (beta.size != length(beta)) 
      stop("size of beta incompatible with the trend model (covariates)")
  }
  if (do.prediction) {
    locations <- .check.locations(locations)
    if (!is.null(borders)) {
      nloc0 <- nrow(locations)
      ind.loc0 <- .geoR_inout(locations, borders)
      locations <- locations[ind.loc0, , drop = TRUE]
      if (nrow(locations) == 0) {
        warning("\nkrige.bayes: no prediction to be performed.\n There are no prediction locations inside the borders")
        do.prediction <- FALSE
      }
      if (messages.screen) 
        cat("krige.bayes: results will be returned only for prediction locations inside the borders\n")
    }
    krige1d <- ifelse((length(unique(locations[, 1])) == 
                         1 | length(unique(locations[, 2])) == 1), TRUE, FALSE)
    if (inherits(model$trend.d, "formula") | inherits(model$trend.l, 
                                                      "formula")) {
      if ((inherits(model$trend.d, "formula") == FALSE) | 
          (inherits(model$trend.l, "formula") == FALSE)) 
        stop("krige.bayes: model$trend.d and model$trend.l must have similar specification\n")
    }
    else {
      if ((length(class(model$trend.d)) > 0 && class(model$trend.d) == 
           "trend.spatial") & (length(class(model$trend.l)) > 
                               0 && class(model$trend.l) == "trend.spatial")) {
        if (ncol(model$trend.d) != ncol(model$trend.l)) 
          stop("krige.bayes: trend.d and trend.l do not have the same number of columns")
      }
      else if (model$trend.d != model$trend.l) 
        stop("krige.bayes: especification of model$trend.l and model$trend.d must be compatible")
    }
    if (messages.screen) {
      cat(switch(as.character(model$trend.d)[1], cte = "krige.bayes: model with constant mean", 
                 `1st` = "krige.bayes: model with mean given by a 1st order polynomial on the coordinates", 
                 `2nd` = "krige.bayes: model with mean given by a 2nd order polynomial on the coordinates", 
                 "krige.bayes: model with mean defined by covariates provided by the user"))
      cat("\n")
    }
    dimnames(locations) <- list(NULL, NULL)
    if (class(model$trend.l) == "trend.spatial") 
      assign("trend.loc", unclass(model$trend.l), envir = pred.env)
    else assign("trend.loc", unclass(trend.spatial(trend = model$trend.l, 
                                                   geodata = list(coords = locations))), envir = pred.env)
    ni <- nrow(get("trend.loc", envir = pred.env))
    if (!is.null(borders)) 
      if (ni == nloc0) 
        assign("trend.loc", get("trend.loc", envir = pred.env)[ind.loc0, 
                                                               , drop = FALSE], envir = pred.env)
    if (nrow(locations) != ni) 
      stop("trend.l is not compatible with number of prediction locations")
    expect.env <- new.env()
    assign("expect", 0, envir = expect.env)
    assign("expect2", 0, envir = expect.env)
  }
  if (!is.null(model$aniso.pars)) {
    if (abs(model$aniso.pars[2] - 1) > 0.001) {
      if (messages.screen) 
        cat("krige.bayes: anisotropy parameters provided and assumed to be constants\n")
      coords <- coords.aniso(coords = coords, aniso.pars = model$aniso.pars)
      if (do.prediction) 
        locations <- coords.aniso(coords = locations, 
                                  aniso.pars = model$aniso.pars)
      remove(dists.env)
      dists.env <- new.env()
      assign("data.dist", as.vector(dist(coords)), envir = dists.env)
    }
  }
  if (do.prediction) {
    assign("d0", loccoords(coords = coords, locations = locations), 
           envir = pred.env)
    loc.coincide <- apply(get("d0", envir = pred.env), 2, 
                          function(x) {
                            any(x < 1e-12)
                          })
    if (any(loc.coincide)) 
      loc.coincide <- (1:ni)[loc.coincide]
    else loc.coincide <- NULL
    if (!is.null(loc.coincide)) {
      temp.f <- function(x, data) {
        return(data[x < 1e-10])
      }
      data.coincide <- apply(get("d0", envir = pred.env)[, 
                                                         loc.coincide, drop = FALSE], 2, temp.f, data = data)
    }
    else data.coincide <- NULL
    n.loc.coincide <- length(loc.coincide)
    assign("loc.coincide", loc.coincide, envir = pred.env)
    assign("data.coincide", data.coincide, envir = pred.env)
    remove(data.coincide, loc.coincide)
  }
  beta.info <- list()
  sigmasq.info <- list()
  for (i in 1:npr) {
    beta.info[[i]] <- switch(prior$beta.prior, fixed = list(mivm = 0, 
                                                            ivm = 0, iv = Inf, beta.fixed = beta.fixed, p = 0), 
                             flat = list(mivm = 0, ivm = 0, iv = 0, p = beta.size), 
                             normal = list(mivm = betares[[i]]$mivm, ivm = betares[[i]]$ivm, 
                                           iv = betares[[i]]$iv, p = 0))
    sigmasq.info[[i]] <- switch(prior$sigmasq.prior, fixed = list(df.sigmasq = Inf, 
                                                                  n0S0 = 0, sigmasq.fixed = sigmasq.fixed), reciprocal = list(df.sigmasq = 0, 
                                                                                                                              n0S0 = 0), uniform = list(df.sigmasq = -2, n0S0 = 0), 
                                sc.inv.chisq = list(df.sigmasq = df.sigmasq.prior, 
                                                    n0S0 = df.sigmasq.prior * S2.prior[i]))
  }
  beta.info$p <- switch(prior$beta.prior, fixed = 0, flat = beta.size, 
                        normal = 0)
  sigmasq.info$df.sigmasq <- switch(prior$sigmasq.prior, fixed = Inf, 
                                    reciprocal = 0, uniform = -2, sc.inv.chisq = df.sigmasq.prior)
  if (prior$phi.prior == "fixed") {
    phi.fixed <- prior$phi
    bsp <- geoR:::beta.sigmasq.post(n = n, beta.info = beta.info[[1]], 
                             sigmasq.info = sigmasq.info[[1]], env.dists = dists.env, 
                             model = list(cov.model = model$cov.model, kappa = model$kappa), 
                             xmat = trend.data, y = data, phi = phi.fixed, tausq.rel = tausq.rel.fixed, 
                             do.prediction.moments = (do.prediction && moments), 
                             do.prediction.simulations = (do.prediction && simulations.predictive), 
                             env.pred = pred.env, signal = (do.prediction && signal))
    if (prior$beta.prior == "fixed") 
      kb$posterior$beta <- list(status = "fixed", fixed.value = beta.fixed)
    else {
      if (prior$sigmasq.prior == "fixed") 
        kb$posterior$beta <- list(distribution = "normal")
      else kb$posterior$beta <- list(distribution = "t", 
                                     conditional = "normal")
      kb$posterior$beta$pars <- list(mean = bsp$beta.post, 
                                     var = bsp$S2.post * bsp$beta.var.std.post)
      attr(kb$posterior$beta$pars$var, "Size") <- beta.size
      class(kb$posterior$beta$pars$var) <- "betavar"
    }
    if (prior$sigmasq.prior == "fixed") 
      kb$posterior$sigmasq <- list(status = "fixed", fixed.value = sigmasq.fixed)
    else kb$posterior$sigmasq <- list(distribution = "sc.inv.chisq", 
                                      pars = list(df = bsp$df.post, S2 = bsp$S2.post))
    kb$posterior$phi <- list(status = "fixed", fixed.value = phi.fixed)
    kb$posterior$tausq.rel <- list(status = "fixed", fixed.value = tausq.rel.fixed)
    kb$predictive$mean <- bsp$pred.mean
    kb$predictive$variance <- bsp$pred.var
    kb$predictive$distribution <- ifelse(prior$sigmasq.prior == 
                                           "fixed", "normal", "t")
    bsp[c("pred.mean", "pred.var")] <- NULL
    if (do.prediction && simulations.predictive && n.predictive > 
        0) {
      phidist$s2 <- as.matrix(bsp$S2.post)
      phidist$probphitausq <- as.matrix(1)
      phidist$beta <- array(bsp$beta.post, c(1, 1, beta.size))
      phidist$varbeta <- array(bsp$beta.var.std.post, c(1, 
                                                        1, beta.size^2))
      phi.unique <- phidist$phitausq <- t(c(phi.fixed, 
                                            tausq.rel.fixed))
      df.model <- bsp$df.post
      ind.length <- 1
      inv.lower <- array(bsp$inv.lower, dim = c(1, 1, (n * 
                                                         (n - 1)/2)))
      inv.diag <- array(bsp$inv.diag, dim = c(1, 1, n))
      ind.table <- n.predictive
      phi.discrete <- phi.fixed
      tausq.rel.discrete <- tausq.rel.fixed
    }
  }
  else {
    if (messages.screen) 
      cat("krige.bayes: computing the discrete posterior of phi/tausq.rel\n")
    if (is.null(phi.discrete)) {
      phi.discrete <- seq(0, 2 * data.dist.max, l = 51)[-1]
      if (messages.screen) 
        cat("krige.bayes: argument `phi.discrete` not provided, using default values\n")
    }
    if (mode(phi.discrete) != "numeric") 
      stop("non-numerical value provided in phi.discrete")
    if (length(phi.discrete) == 1) 
      stop("only one value provided in phi.discrete. Use prior.phi=`fixed`")
    n.phi.discrete <- length(phi.discrete)
    n.tausq.rel.discrete <- length(tausq.rel.discrete)
    phi.names <- paste("phi", phi.discrete, sep = "")
    tausq.rel.names <- paste("tausqrel", tausq.rel.discrete, 
                             sep = "")
    phidist$phitausq <- as.matrix(expand.grid(phi.discrete, 
                                              tausq.rel.discrete))
    if (prior$phi.prior == "user" | prior$tausq.rel.prior == 
        "user") {
      if (prior$tausq.rel.prior == "fixed") 
        phidist$phitausq <- cbind(phidist$phitausq, prior$priors.info$phi$probs, 
                                  1)
      else {
        if (is.null(prior$joint.phi.tausq)) 
          phidist$phitausq <- cbind(phidist$phitausq, 
                                    as.matrix(expand.grid(prior$priors.info$phi$probs, 
                                                          prior$priors.info$tausq.rel$probs)))
        else phidist$phitausq <- cbind(phidist$phitausq, 
                                       as.vector(prior$joint.phi.tausq.rel), 1)
      }
    }
    dimnames(phidist$phitausq) <- list(NULL, NULL)
    df.model <- ifelse(sigmasq.info$df.sigmasq == Inf, Inf, 
                       (n + sigmasq.info$df.sigmasq - beta.info$p))
    phi.tausq.rel.post <- function(phinug) {
      par.set <- get("parset", envir = counter.env)
      if (messages.screen) 
        krige.bayes.counter(.temp.ap = par.set, n.points = ntocount)
      on.exit(assign("parset", get("parset", envir = counter.env) + 
                       1, envir = counter.env))
      phi <- phinug[1]
      tausq.rel <- phinug[2]
      if (prior$beta.prior == "normal" && npr > 1) 
        info.id <- par.set
      else info.id <- 1
      bsp <- geoR:::beta.sigmasq.post(n = n, beta.info = beta.info[[info.id]], 
                               sigmasq.info = sigmasq.info[[info.id]], env.dists = dists.env, 
                               model = list(cov.model = model$cov.model, kappa = model$kappa), 
                               xmat = trend.data, y = data, phi = phi, tausq.rel = tausq.rel, 
                               dets = TRUE, do.prediction.moments = (do.prediction && 
                                                                       moments), do.prediction.simulations = (do.prediction && 
                                                                                                                simulations.predictive), env.pred = pred.env, 
                               signal = signal)
      logprobphitausq <- (-0.5) * log(bsp$det.XiRX) - (bsp$log.det.to.half) - 
        (bsp$df.post/2) * log(bsp$S2.post)
      if (prior$phi.prior == "user") {
        if (phinug[3] > 0) 
          logprobphitausq <- logprobphitausq + log(phinug[3])
        else logprobphitausq <- -Inf
      }
      if (prior$phi.prior == "reciprocal") {
        if (phi > 0) 
          logprobphitausq <- logprobphitausq - log(phi)
        else logprobphitausq <- -Inf
      }
      if (prior$phi.prior == "squared.reciprocal") {
        if (phi > 0) 
          logprobphitausq <- logprobphitausq - 2 * log(phi)
        else logprobphitausq <- -Inf
      }
      if (prior$phi.prior == "exponential") {
        logprobphitausq <- logprobphitausq - log(exponential.par) - 
          (phi/exponential.par)
      }
      if (prior$tausq.rel.prior == "user") {
        if (phinug[4] > 0) 
          logprobphitausq <- logprobphitausq + log(phinug[4])
        else logprobphitausq <- -Inf
      }
      if (prior$tausq.rel.prior == "reciprocal") {
        if (tausq.rel > 0) 
          logprobphitausq <- logprobphitausq - log(tausq.rel)
        else logprobphitausq <- -Inf
      }
      if (get("add.cte", envir = counter.env) && is.finite(logprobphitausq)) {
        assign("cte", logprobphitausq, envir = counter.env)
        assign("add.cte", FALSE, envir = counter.env)
      }
      logprobphitausq <- logprobphitausq - get("cte", envir = counter.env)
      bsp$probphitausq <- drop(exp(logprobphitausq))
      if (do.prediction && moments) {
        assign("expect", (get("expect", envir = expect.env) + 
                            (bsp$pred.mean * bsp$probphitausq)), envir = expect.env)
        assign("expect2", (get("expect2", envir = expect.env) + 
                             ((bsp$pred.var + (bsp$pred.mean^2)) * bsp$probphitausq)), 
               envir = expect.env)
      }
      phi.ind <- which.min(abs(phi.discrete - phi))
      nug.ind <- which.min(abs(tausq.rel.discrete - tausq.rel))
      assign("pn.ind", c(phi.ind, nug.ind), envir = fn.frame)
      assign("bsp", bsp, envir = fn.frame)
      eval(expression(phidist$s2[pn.ind[1], pn.ind[2]] <- bsp$S2.post), 
           envir = fn.frame)
      eval(expression(phidist$probphitausq[pn.ind[1], pn.ind[2]] <- bsp$probphitausq), 
           envir = fn.frame)
      eval(expression(phidist$beta[pn.ind[1], pn.ind[2], 
                                   ] <- drop(bsp$beta.post)), envir = fn.frame)
      eval(expression(phidist$varbeta[pn.ind[1], pn.ind[2], 
                                      ] <- drop(bsp$beta.var.std.post)), envir = fn.frame)
      if (do.prediction && simulations.predictive) {
        eval(expression(inv.lower[pn.ind[1], pn.ind[2], 
                                  ] <- bsp$inv.lower), envir = fn.frame)
        eval(expression(inv.diag[pn.ind[1], pn.ind[2], 
                                 ] <- bsp$inv.diag), envir = fn.frame)
      }
      return(invisible())
    }
    fn.frame <- sys.frame(sys.nframe())
    phidist$s2 <- matrix(NA, n.phi.discrete, n.tausq.rel.discrete)
    dimnames(phidist$s2) <- list(phi.names, tausq.rel.names)
    phidist$probphitausq <- matrix(NA, n.phi.discrete, n.tausq.rel.discrete)
    phidist$beta <- array(NA, dim = c(n.phi.discrete, n.tausq.rel.discrete, 
                                      beta.size))
    dimnames(phidist$beta) <- list(phi.names, tausq.rel.names, 
                                   beta.names)
    phidist$varbeta <- array(NA, dim = c(n.phi.discrete, 
                                         n.tausq.rel.discrete, beta.size^2))
    dimnames(phidist$varbeta) <- list(phi.names, tausq.rel.names, 
                                      NULL)
    if (do.prediction && simulations.predictive) {
      inv.lower <- array(NA, dim = c(n.phi.discrete, n.tausq.rel.discrete, 
                                     (n * (n - 1))/2))
      inv.diag <- array(NA, dim = c(n.phi.discrete, n.tausq.rel.discrete, 
                                    n))
      frame.inv <- sys.frame(sys.nframe())
    }
    counter.env <- new.env()
    assign("parset", 1, envir = counter.env)
    assign("add.cte", TRUE, envir = counter.env)
    assign("cte", 0, envir = counter.env)
    if (messages.screen) {
      ntocount <- nrow(phidist$phitausq)
      cat(paste("krige.bayes: computing the posterior probabilities.\n             Number of parameter sets: ", 
                ntocount, "\n"))
    }
    temp.res <- apply(phidist$phitausq, 1, phi.tausq.rel.post)
    remove("bsp")
    if (messages.screen) {
      remove(counter.env)
      cat("\n")
    }
    phidist$sum.prob <- sum(phidist$probphitausq)
    phidist$probphitausq <- phidist$probphitausq/phidist$sum.prob
    kb$posterior$beta <- list(conditional.distribution = "normal", 
                              pars = list(mean = phidist$beta, var = phidist$varbeta))
    attr(kb$posterior$beta$pars$var, "Size") <- beta.size
    kb$posterior$sigmasq <- list(conditional.distribution = "sc.inv.chisq", 
                                 pars = list(df = df.model, S2 = drop(phidist$s2)))
    kb$posterior$phi$distribution <- drop(rowSums(phidist$probphitausq))
    names(kb$posterior$phi$distribution) <- prior$phi.discrete
    if (prior$tausq.rel.prior != "fixed") {
      kb$posterior$tausq.rel$distribution <- drop(colSums(phidist$probphitausq))
      names(kb$posterior$tausq.rel$distribution) <- tausq.rel.discrete
    }
    else {
      kb$posterior$tausq.rel <- list(status = "fixed", 
                                     fixed.value = tausq.rel.fixed)
    }
    if (prior$phi.prior != "fixed" | prior$tausq.rel.prior != 
        "fixed") {
      kb$posterior$joint.phi.tausq.rel <- phidist$probphitausq
      dimnames(kb$posterior$joint.phi.tausq.rel) <- list(phi.names, 
                                                         tausq.rel.names)
    }
    if (n.posterior > 0) {
      if (messages.screen) 
        cat("krige.bayes: sampling from posterior distribution\n")
      n.points <- length(phidist$probphitausq)
      ind <- sample((1:n.points), n.posterior, replace = TRUE, 
                    prob = as.vector(phidist$probphitausq))
      phi.sam <- phidist$phitausq[ind, ]
      ind.unique <- sort(unique(ind))
      ind.length <- length(ind.unique)
      ind.table <- table(ind)
      names(ind.table) <- NULL
      phi.unique <- phidist$phitausq[ind.unique, , drop = FALSE]
      if (messages.screen) {
        cat("krige.bayes: sample from the (joint) posterior of phi and tausq.rel\n")
        print(rbind(phi = phi.unique[, 1], tausq.rel = phi.unique[, 
                                                                  2], frequency = ind.table))
        cat("\n")
      }
      vecpars.back.order <- order(order(ind))
      sigmasq.sam <- rinvchisq(n = n.posterior, df = df.model, 
                               scale = rep(as.vector(phidist$s2)[ind.unique], 
                                           ind.table))
      if (beta.size == 1) {
        vec.beta <- rep(as.vector(phidist$beta)[ind.unique], 
                        ind.table)
        vec.vbeta <- rep(as.vector(phidist$varbeta)[ind.unique], 
                         ind.table)
        beta.sam <- vec.beta + sqrt(sigmasq.sam * vec.vbeta) * 
          rnorm(n.posterior, mean = 0, sd = 1)
      }
      else {
        ind.beta <- matrix(phidist$beta, ncol = beta.size)[ind.unique, 
                                                           , drop = FALSE]
        ind.beta <- ind.beta[rep(1:ind.length, ind.table), 
                             , drop = FALSE]
        ind.vbeta <- matrix(phidist$varbeta, ncol = beta.size^2)[ind.unique, 
                                                                 , drop = FALSE]
        ind.vbeta <- ind.vbeta[rep(1:ind.length, ind.table), 
                               , drop = FALSE] * sigmasq.sam
        temp.res <- apply(ind.vbeta, 1, rMVnorm, beta.size = beta.size)
        beta.sam <- ind.beta + t(temp.res)
        remove("temp.res")
      }
      if (beta.size == 1) {
        trend.mean <- mean(beta.sam)
        trend.median <- median(beta.sam)
      }
      else {
        trend.mean <- colMeans(beta.sam)
        trend.median <- apply(beta.sam, 2, median)
      }
      S2.mean <- mean(sigmasq.sam)
      S2.median <- median(sigmasq.sam)
      phi.marg <- rowSums(phidist$probphitausq)
      .marg <- phi.marg/(sum(phi.marg))
      phi.mean <- phi.discrete %*% phi.marg
      phi.median <- median(phi.sam[, 1])
      phi.mode <- phi.discrete[which.min(abs(phi.marg - 
                                               max(phi.marg)))]
      tausq.rel.marg <- colSums(phidist$probphitausq)
      tausq.rel.marg <- tausq.rel.marg/(sum(tausq.rel.marg))
      tausq.rel.mean <- tausq.rel.discrete %*% tausq.rel.marg
      tausq.rel.median <- median(phi.sam[, 2])
      tausq.rel.mode <- tausq.rel.discrete[which.min(abs(tausq.rel.marg - 
                                                           max(tausq.rel.marg)))]
      mode.ind <- which(phidist$probphitausq == max(phidist$probphitausq))
      phi.tausq.rel.mode <- phidist$phitausq[mode.ind, 
                                             1:2, drop = FALSE]
      if (nrow(phi.tausq.rel.mode) > 1) {
        cat("krige.bayes: WARNING: multiple modes for phi.tausq.rel. Using the first one to compute modes of beta and sigmasq.\n")
        cat("krige.bayes: modes found are:\n")
        print(phi.tausq.rel.mode)
        phi.tausq.rel.mode <- phi.tausq.rel.mode[1, ]
      }
      if (prior$beta.prior == "normal" && npr > 1) 
        info.id <- mode.ind
      else info.id <- 1
      modes <- geoR:::beta.sigmasq.post(n = n, beta.info = beta.info[[info.id]], 
                                 sigmasq.info = sigmasq.info[[info.id]], env.dists = dists.env, 
                                 model = list(cov.model = model$cov.model, kappa = model$kappa), 
                                 xmat = trend.data, y = data, phi = phi.tausq.rel.mode[1], 
                                 tausq.rel = phi.tausq.rel.mode[2], dets = FALSE, 
                                 do.prediction.moments = FALSE, do.prediction.simulations = FALSE, 
                                 env.pred = pred.env, signal = signal)
      beta.mode.cond <- modes$beta.post
      S2.mode.cond <- modes$S2.post
      rm(modes)
      if (beta.size == 1) 
        kb$posterior$beta$summary <- c(mean = trend.mean, 
                                       median = trend.median, mode.cond = beta.mode.cond)
      else kb$posterior$beta$summary <- cbind(mean = trend.mean, 
                                              median = trend.median, mode.cond = beta.mode.cond)
      kb$posterior$sigmasq$summary <- c(mean = S2.mean, 
                                        median = S2.median, mode.cond = S2.mode.cond)
      kb$posterior$phi$summary <- c(mean = phi.mean, median = phi.median, 
                                    mode = phi.mode)
      if (prior$tausq.rel.prior != "fixed") 
        kb$posterior$tausq.rel$summary <- c(mean = tausq.rel.mean, 
                                            median = tausq.rel.median, mode = tausq.rel.mode)
      kb$posterior$sample <- as.data.frame(cbind(drop(as.matrix(beta.sam)[vecpars.back.order, 
                                                                          ]), sigmasq.sam[vecpars.back.order], phi.sam[, 
                                                                                                                       1]))
      beta.sam <- sigmasq.sam <- NULL
      names(kb$posterior$sample) <- c(beta.names, "sigmasq", 
                                      "phi")
      kb$posterior$sample$tausq.rel <- phi.sam[, 2]
      phi.lev <- unique(phidist$phitausq[, 1])
      kb$posterior$phi$phi.marginal <- data.frame(phi = phi.lev, 
                                                  expected = rowSums(phidist$probphitausq), sampled = as.vector(table(factor(phi.sam[, 
                                                                                                                                     1], levels = phi.lev)))/n.posterior)
      tausq.rel.lev <- unique(phidist$phitausq[, 2])
      if (prior$tausq.rel.prior != "fixed") 
        kb$posterior$tausq.rel$tausq.rel.marginal <- data.frame(tausq.rel = tausq.rel.lev, 
                                                                expected = colSums(phidist$probphitausq), sampled = as.vector(table(factor(phi.sam[, 
                                                                                                                                                   2], levels = tausq.rel.lev)))/n.posterior)
    }
    if (do.prediction) {
      kb$predictive$distribution <- "obtained by numerical approximation"
      if (messages.screen) 
        cat("krige.bayes: starting prediction at the provided locations\n")
      if (n.predictive == n.posterior) {
        include.it <- FALSE
        phi.sam <- phidist$phitausq[ind, ]
        message.prediction <- c(message.prediction, "krige.bayes: phi/tausq.rel samples for the predictive are same as for the posterior")
        if (messages.screen) 
          cat(message.prediction, "\n")
      }
      else {
        include.it <- TRUE
        ind <- sample((1:(dim(phidist$phitausq)[1])), 
                      n.predictive, replace = TRUE, prob = as.vector(phidist$probphitausq))
        ind.unique <- sort(unique(ind))
        ind.length <- length(ind.unique)
        ind.table <- table(ind)
        phi.unique <- phidist$phitausq[ind.unique, , 
                                       drop = FALSE]
        message.prediction <- c(message.prediction, "krige.bayes: phi/tausq.rel samples for the predictive are NOT the same as for the posterior ")
        if (messages.screen) {
          cat(message.prediction, "\n")
          cat("krige.bayes: samples and their frequencies from the distribution of  phi and tau.rel when drawing from the predictive distribution\n")
          print(rbind(phi = phi.unique[, 1], tausq.rel = phi.unique[, 
                                                                    2], frequency = ind.table))
        }
        phi.sam <- phidist$phitausq[ind, ]
        vecpars.back.order <- order(order(ind))
      }
      if (moments) {
        if (messages.screen) 
          cat("krige.bayes: computing moments of the predictive distribution\n")
        kb$predictive$mean <- get("expect", envir = expect.env)/phidist$sum.prob
        remove("expect", envir = expect.env)
        kb$predictive$variance <- (get("expect2", envir = expect.env)/phidist$sum.prob) - 
          ((kb$predictive$mean)^2)
        remove("expect2", envir = expect.env)
      }
    }
  }
  if ((do.prediction && moments) & (abs(lambda - 1) > 0.001)) {
    kb$predictive <- backtransform.moments(lambda = lambda, 
                                           mean = kb$predictive$mean, variance = kb$predictive$variance, 
                                           distribution = kb$predictive$distribution, n.simul = n.back.moments)
  }
  if (do.prediction && simulations.predictive) {
    if (is.R()) {
      if (cov.model.number > 12) 
        stop("simulation in krige.bayes not implemented for the choice of correlation function")
    }
    else if (cov.model.number > 10) 
      stop("simulation in krige.bayes not implemented for the correlation function chosen")
    krige.bayes.aux20 <- function(phinug) {
      iter <- get("counter", envir = counter.env)
      if (messages.screen & prior$phi.prior != "fixed") 
        krige.bayes.counter(.temp.ap = iter, n.points = ind.length)
      phinug <- as.vector(phinug)
      phi <- phinug[1]
      tausq.rel <- phinug[2]
      phi.ind <- which.min(abs(phi.discrete - phi))
      nug.ind <- which.min(abs(tausq.rel.discrete - tausq.rel))
      v0 <- cov.spatial(obj = get("d0", envir = pred.env), 
                        cov.model = cov.model, kappa = kappa, cov.pars = c(1, 
                                                                           phi))
      b <- .bilinearformXAY(X = as.vector(cbind(data, trend.data)), 
                            lowerA = as.vector(inv.lower[phi.ind, nug.ind, 
                                                         , drop = TRUE]), diagA = as.vector(inv.diag[phi.ind, 
                                                                                                     nug.ind, , drop = TRUE]), Y = as.vector(v0))
      tv0ivdata <- drop(b[1, ])
      b <- t(get("trend.loc", envir = pred.env)) - b[-1, 
                                                     , drop = FALSE]
      tmean <- tv0ivdata + drop(crossprod(b, as.vector(phidist$beta[phi.ind, 
                                                                    nug.ind, ])))
      tv0ivdata <- NULL
      Nsims <- ind.table[iter]
      if (signal) 
        Dval <- 1
      else Dval <- 1 + tausq.rel
      iter.env <- sys.frame(sys.nframe())
      coincide.cond <- (((tausq.rel < 1e-12) | !signal) & 
                          !is.null(get("loc.coincide", envir = pred.env)))
      if (coincide.cond) {
        nloc <- ni - n.loc.coincide
        ind.not.coincide <- -(get("loc.coincide", envir = pred.env))
        v0 <- v0[, ind.not.coincide, drop = FALSE]
        tmean <- tmean[ind.not.coincide]
        b <- b[, ind.not.coincide, drop = FALSE]
      }
      else {
        nloc <- ni
        ind.not.coincide <- TRUE
      }
      par.set <- get("parset", envir = counter.env)
      if (prior$beta.prior == "normal" && npr > 1) 
        info.id <- par.set
      else info.id <- 1
      if (any(beta.info[[info.id]]$iv == Inf)) 
        vbetai <- matrix(0, ncol = beta.size, nrow = beta.size)
      else vbetai <- matrix(drop(phidist$varbeta[phi.ind, 
                                                 nug.ind, ]), ncol = beta.size, nrow = beta.size)
      simul <- matrix(NA, nrow = ni, ncol = Nsims)
      if (nloc > 0) 
        simul[ind.not.coincide, ] <- .cond.sim(env.loc = base.env, 
                                               env.iter = iter.env, loc.coincide = get("loc.coincide", 
                                                                                       envir = pred.env), coincide.cond = coincide.cond, 
                                               tmean = tmean, Rinv = list(lower = drop(inv.lower[phi.ind, 
                                                                                                 nug.ind, ]), diag = drop(inv.diag[phi.ind, 
                                                                                                                                   nug.ind, ])), mod = list(beta.size = beta.size, 
                                                                                                                                                            nloc = nloc, Nsims = Nsims, n = n, Dval = Dval, 
                                                                                                                                                            df.model = df.model, s2 = phidist$s2[phi.ind, 
                                                                                                                                                                                                 nug.ind], cov.model.number = cov.model.number, 
                                                                                                                                                            phi = phi, kappa = kappa, nugget = tausq.rel), 
                                               vbetai = vbetai, fixed.sigmasq = (sigmasq.info$df.sigmasq == 
                                                                                   Inf))
      if (coincide.cond) 
        simul[get("loc.coincide", envir = pred.env), 
              ] <- rep(get("data.coincide", envir = pred.env), 
                       Nsims)
      remove("v0", "b", "tmean")
      assign("counter", (iter + 1), envir = counter.env)
      assign("parset", get("parset", envir = counter.env) + 
               1, envir = counter.env)
      if (abs(lambda - 1) > 0.001) {
        return(BCtransform(x = simul, lambda = lambda, 
                           inverse = TRUE)$data)
      }
      else return(simul)
    }
    counter.env <- new.env()
    if (messages.screen) {
      cat("krige.bayes: sampling from the predictive\n")
      if (prior$phi.prior != "fixed") 
        cat(paste("             Number of parameter sets: ", 
                  ind.length, "\n"))
    }
    assign("counter", 1, envir = counter.env)
    assign("parset", 1, envir = counter.env)
    kb$predictive$simulations <- matrix(unlist(apply(phi.unique, 
                                                     1, krige.bayes.aux20)), ncol = n.predictive)
    remove("inv.lower", "inv.diag", "counter.env", "pred.env")
    if (messages.screen) 
      if (abs(lambda - 1) > 0.001) 
        cat("krige.bayes: Box-Cox data transformation performed.\n             Simulations back-transformed to the original scale\n")
    if (messages.screen) 
      cat("krige.bayes: preparing summaries of the predictive distribution\n")
    kb$predictive <- c(kb$predictive, statistics.predictive(simuls = kb$predictive$simulations, 
                                                            mean.var = mean.estimator, quantile = quantile.estimator, 
                                                            threshold = probability.estimator, sim.means = sim.means, 
                                                            sim.vars = sim.vars))
    if (sim.means && exists("vecpars.back.order")) 
      kb$predictive$sim.means[vecpars.back.order]
    if (sim.vars && exists("vecpars.back.order")) 
      kb$predictive$sim.vars[vecpars.back.order]
    if (keep.simulations) {
      if (prior$phi.prior != "fixed") 
        kb$predictive$simulations <- kb$predictive$simulations[, 
                                                               vecpars.back.order]
    }
    else {
      kb$predictive$simulations <- NULL
    }
    if (prior$phi.prior != "fixed") {
      if (include.it) {
        phi.lev <- unique(phidist$phitausq[, 1])
        kb$predictive$phi.marginal <- data.frame(phi = phi.lev, 
                                                 expected = rowSums(phidist$probphitausq), sampled = as.vector(table(factor(phi.sam[, 
                                                                                                                                    1], levels = phi.lev)))/n.predictive)
        tausq.rel.lev <- unique(phidist$phitausq[, 2])
        if (prior$tausq.rel.prior != "fixed") 
          data.frame(tausq.rel = tausq.rel.lev, expected = colSums(phidist$probphitausq), 
                     sampled = as.vector(table(factor(phi.sam[, 
                                                              2], levels = tausq.rel.lev)))/n.predictive)
        else kb$predictive$tausq.rel.marginal <- paste("fixed tausq.rel with value =", 
                                                       tausq.rel)
        kb$predictive$tausq.rel.marginal <- data.frame(tausq.rel = tausq.rel.lev, 
                                                       expected = colSums(phidist$probphitausq), sampled = as.vector(table(factor(phi.sam[, 
                                                                                                                                          2], levels = tausq.rel.lev)))/n.predictive)
      }
    }
  }
  if (!do.prediction) 
    kb$predictive <- "no prediction locations provided"
  kb$.Random.seed <- seed
  kb$max.dist <- data.dist.max
  kb$call <- call.fc
  attr(kb, "prediction.locations") <- call.fc$locations
  attr(kb, "parent.env") <- parent.frame()
  if (!is.null(call.fc$coords)) 
    attr(kb, "data.locations") <- call.fc$coords
  else attr(kb, "data.locations") <- substitute(a$coords, list(a = substitute(geodata)))
  if (do.prediction) 
    attr(kb, "sp.dim") <- ifelse(krige1d, "1d", "2d")
  if (!is.null(call.fc$borders)) 
    attr(kb, "borders") <- call.fc$borders
  oldClass(kb) <- c("krige.bayes", "variomodel")
  return(kb)
}