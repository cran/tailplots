s_lnorm = Vectorize(function(sigma) {
  integrand = function(t) tanh( sigma * t) * exp(-t^2)
  val = integrate(integrand, lower=0, upper=Inf)$value
  return( 2 / sqrt(pi) * val )
})


sigma_root <- Vectorize(function(t) {
  f = function(u, t) s_lnorm(u) - t
  uniroot(f, c(1e-8, 100), extendInt="upX", t=t)$root
})


s_vec = function(x,u) sapply(u, s_scalar, x=x)