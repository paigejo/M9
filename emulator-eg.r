#emulator example

g <- function(x, theta) {
  outer(x, theta, function(x, t) exp(-x) * sin(t - pi * x))
}
x <- c(0.1, 0.5, 0.7)
theta <- c(0, 1, 2, 3.5, 5) # uneven spacing more interesting
Y <- g(x, theta)

## little picture

matplot(theta, t(Y), xlim = c(0, 2 * pi), ylim = c(-1, 1),
        type = 'p', pch = 1:3, bty = 'n',
        xlab = "Theta in [0, 2pi]", ylab = "g(x, theta)",
        main = "True function and evaluations")
tfull <- seq(from = 0, to = 2*pi, len = 101)
matplot(tfull, t(g(x, tfull)), type = "l", lty = 2, add = TRUE)
legend('topright', legend = paste('x =', x),
       col = 1:3, lty = 2, pch = 1:3, bty = 'n')

## Set up the regressors and variance functions: polynomials for the
## runs regressors (should be Legendre polynomials really, shifted
## onto [0, 1]); Fourier terms for outputs regressors; power
## exponential for the runs variance function; circular correlation
## for the outputs variance matrix (note that pi cannot be too small
## or this variance is singular)

## put rownames in gr and on Gs, just for clarity

gr <- function(x) {
  robj <- cbind(1, 2*x - 1, x^2)
  rownames(robj) <- paste('x', seq(along = x), sep = '')
  robj
}

kappar <- function(x, xp = x, range = 0.5)
  exp(-abs(outer(x, xp, '-') / range)^(3/2))

Gs <- cbind(1, sin(theta), cos(theta), sin(2 * theta), cos(2 * theta))
rownames(Gs) <- paste('th', seq(along = theta), sep = '')

circular <- function(ang1, ang2 = ang1, range = pi / 1.1)
{
  smallestAngle <- function(a, b) {
    dd <- outer(a, b, '-')
    pmin(abs(dd), abs(2*pi + dd), abs(dd - 2*pi))
  }
  
  angles <- smallestAngle(ang1, ang2)
  ifelse(angles < range, 1 - angles / range, 0)
}

Ws <- circular(theta)

## Set up a minimal prior for the NIG (in general, thought is required
## here!)

local({
  vr <- length(gr(0))
  vs <- ncol(Gs)
  m <- rep(0, vr * vs)
  V <- diag(1^2, vr * vs)
  a <- 1
  d <- 1^2
  NIG <<- list(m = m, V = V, a = a, d = d)
})

## Now we're ready to initialise our OPE

myOPE <- initOPE(gr = gr, kappar = kappar, Gs = Gs, Ws = Ws, NIG = NIG)

xnew <- 0.4
pp0 <- predictOPE(myOPE, Rp = xnew) # prior prediction

## Adjust with the evaluations

myOPE <- adjustOPE(myOPE, R = x, Y = Y)

## Sanity check: predict the points we already have

pp1 <- predictOPE(myOPE, R = x)
stopifnot(
  all.equal.numeric(pp1$mu, Y, check.attributes = FALSE),
  all.equal.numeric(pp1$Sigma, array(0, dim(pp1$Sigma)))
) # phew!

## Make a prediction at some new x values, and add to the plot as
## error bars

pp2 <- predictOPE(myOPE, Rp = xnew)
pp2$mu <- c(pp2$mu) # reshape for convenience
dim(pp2$Sigma) <- rep(length(pp2$mu), 2) # 

mu <- pp2$mu
sig <- sqrt(diag(pp2$Sigma))
arrows(theta, mu + sig * qt(0.025, df = pp2$df),
       theta, mu + sig * qt(0.975, df = pp2$df),
       code = 3, angle = 90, length = 0.1, col = 'blue')
lines(tfull, g(xnew, tfull), col = 'blue')

## Add on some sampled values, interpolated using splines

rsam <- sampleOPE(myOPE, Rp = xnew, N = 10)
if (require(splines)) {
  for (i in 1:nrow(rsam)) {
    pispl <- periodicSpline(theta, rsam[i, ], period = 2*pi)
    lines(predict(pispl, tfull), col = 'darkgrey')
  }
  legend('topleft', legend = c(paste('x =', xnew, '(predicted)'), 'sampled'),
         col = c('blue', 'darkgrey'), lty = 1, pch = NA, bty = 'n')
}

## A more complicated prediction

xnew <- c(xnew, 0.8)
pp3 <- predictOPE(myOPE, Rp = xnew, type = 'EV')

#my code:
mu <- pp3$mu
dim(pp3$Sigma) <- rep(length(pp3$mu), 2)
sig <- sqrt(diag(pp3$Sigma))
arrows(theta, mu + sig * qt(0.025, df = pp3$df),
       theta, mu + sig * qt(0.975, df = pp3$df),
       code = 3, angle = 90, length = 0.1, col = 'green')
lines(tfull, g(xnew, tfull), col = 'green')