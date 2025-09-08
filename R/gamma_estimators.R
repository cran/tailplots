############################################################################
## Paper: A gamma tail statistic and its asymptotics by Iwashita and Klar ##
############################################################################


# This function returs the theoretical value of g for the gamma distribution with shape a. See Remark 1 of the paper. 
g_gamma = function(a) {
    ifelse(a<=0, 1, ( 2^{2*a-1}*a*beta(a,a) )^{-1} )
}

#This function is used for the root search. 
alpha_root_gamma <- function(g) {
  f <- function(a) g_gamma(a) - g
  eps <- sqrt(.Machine$double.eps)
  uniroot(f, c(eps, 100), extendInt = "yes")$root
}

# This function applies a vector to g.scalar. 
g_vec = function(x,d) sapply(d, g_scalar, x=x)
