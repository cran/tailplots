###################################################################
## Paper: A Pareto Tail Plot Without Moment Restrictions by Klar ##
###################################################################

# This function returs the theoretical value of t for the pareto distribution with shape a.  
t_pareto = function(a) {
    ifelse(a<=0, 1, a * (digamma((a+1)/2) - digamma(a/2)) - 1) 
  
}

#This function is used for root search. 
alpha_root <- function(t) {
  f <- function(u) t_pareto(u) - t
  eps <- sqrt(.Machine$double.eps)
  uniroot(f, c(eps, 100), extendInt = "yes")$root
}


# This function applies a vector to t.scalar. 
t_vec = function(x, u) sapply(u, t_scalar, x=x)