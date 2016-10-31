dexpMax = function(x, n=1, lambda=1) {
  n*(1 - exp(-lambda*x))^(n-1) * lambda*exp(-lambda*x)
}
rexpMax = function(x, n=1, lambda=1) {
  rexps = matrix(rexp(x*n, rate=lambda), nrow=n)
  return(c(apply(rexps, 2, max)))
}
pexpMax = function(q, n=1, lambda=1) {
  (1 - exp(-q*lambda))^n
}
qexpMax = function(p, n=1, lambda=1) {
  -log(1 - p^(1/n))*(1/lambda)
}

testDexpMax = function(niter=1000, n=50, lambda=1) {
  rexpMaxes = rexpMax(niter, n, lambda)
  hist(rexpMaxes, breaks=50, freq = FALSE)
  xs = seq(min(rexpMaxes), max(rexpMaxes), l=200)
  lines(xs, dexpMax(xs, n, lambda), col="blue")
}

testDexpMaxMLE = function(n=50, lambda=1, niter=10000) {
  rexpMaxes = rexpMax(niter, n, lambda)
  lambdaHats = log(n)/rexpMaxes
  hist(lambdaHats, breaks=50)
  abline(v=lambda, col="blue")
}