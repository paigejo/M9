
# ideas for tapers
# h = function(x) {ifelse(x <= 0, 0, exp(-1/x))}
# g = function(x, a=1) { a*h(x)/(a*h(x) + h(1-x))}
# 
# gInv = function(x, a=1) {
#   ans = x
#   ans[x == 0] = 0
#   ans[x == 1] = 1
#   inds = (x != 0) & (x != 1)
#   
#   term1 = (log(-(a*x - a)/x) + 2)/(2*log(-(a*x - a)/x))
#   term2 = sqrt(log(-(a*x - a)/x)^2 + 4)/(2*log(-(a*x - a)/x))
#   ans[inds] = term1[inds] - term2[inds]
#   
#   ans
# }

# Randy's exponential taper, renormalized
# taper = function(d, lam=1, dStar=1) {
#   ans = (1 - exp(-lam*(dStar-d)/dStar))/(1 - exp(-lam))
#   ans[x > dStar] = 0
#   return(ans)
# }

# power exponential taper, renormalized
taper = function(d, lambda=1, alpha=2, dStar=21000, normalize=TRUE) {
  scaledD = abs((dStar - d)/dStar - 1)^alpha * lambda^alpha
  if(normalize)
    ans = 1 - (1 - exp(-scaledD))/(1 - exp(-lambda^alpha))
  else
    ans = exp(-scaledD)
  ans[d > dStar] = 0
  return(ans)
}