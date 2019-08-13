function (N) 
{
  S = t(N$Post - N$Pre)
  v = ncol(S)
  return(function(x0, t0, deltat, ...) {
    t = t0
    x = x0
    termt = t0 + deltat
    repeat {
      h = N$h(x, t, ...)
      h0 = sum(h)
      if (h0 < 1e-10) t = 1e+99 else if (h0 > 1e+06) {
        t = 1e+99
        warning("Hazard too big - terminating simulation!")
      } else t = t + rexp(1, h0)
      if (t >= termt) return(x)
      j = sample(v, 1, prob = h)
      x = x + S[, j]
    }
  })
}
