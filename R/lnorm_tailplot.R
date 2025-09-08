#' Estimate of tail functional s and confidence intervals for s and sigma
#' 
#' This function computes the estimate of \eqn{s} and the associated confidence interval for \eqn{s} as well as the standard deviation sigma on the log-scale of the lognormal distribution. Three methods are implemented to compute the confidence intervals: a method based on the unbiased variance estimators of the underlying U-statistics and two resampling methods (jackknife and bootstrap).
#' 
#' @param x a vector containing the sample data. 
#' @param u the threshold for the computation of s.
#' @param confint a boolean value indicating whether the confidence interval should be computed.
#' @param method the method used for computing the confidence intervals (options include unbiased variance estimator, jackknife, and bootstrap).
#' @param R the number of the bootstrap replicates.
#' @param conf.level the confidence level for the interval. 
#' 
#' 
#' @return A matrix containing:
#' \item{threshold}{The value of the threshold u.}
#' \item{s.estimate}{Estimate of the tail functional s.}
#' \item{s.ci1}{The lower bound of the confidence interval for s (if `confint = TRUE`).}
#' \item{s.ci2}{The upper bound of the confidence interval for s (if `confint = TRUE`).}
#' \item{sigma}{Estimate of the scale parameter under a lognormal model.}
#' \item{sigma.ci1}{The lower bound of the confidence interval for sigma (if `confint = TRUE`).}
#' \item{sigma.ci2}{The upper bound of the confidence interval for sigma (if `confint = TRUE`).}
#'
#' @details 
#' The function \eqn{s}, defined by 
#'  \deqn{
#'   s(u) = \mathbb{E}\!\left[
#'     \frac{|X_{1} - X_{2}|}{X_{1} + X_{2}} \;\middle|\; X_{1}\,X_{2} > u
#'   \right], 
#' }
#' where \eqn{X_1,X_2} are independent and identically distributed (i.i.d.) positive random variables, takes a constant value if and only if \eqn{X_1} follows a lognormal distribution. 
#' Thus, \eqn{s} can be used to detect distributions with lognormal tails. The characterization of the lognormal distribution is based on the work of Mosimann (1970). 
#' This function estimates \eqn{s(u)} using U-statistics, similarly as in Iwashita and Klar (2024).
#' 
#' @references 
#' Mosimann, J. E. (1970). Size allometry: size and shape variables with characterizations of the lognormal and generalized gamma distributions. Journal of the American Statistical Association, 65(330):930–945. 
#' \doi{https://doi.org/10.2307/2284599}
#' 
#' Iwashita, T. & Klar, B. (2024). A gamma tail statistic and its asymptotics. Statistica Neerlandica 78:2, 264-280.  
#' \doi{https://doi.org/10.1111/stan.12316}
#' 
#' @examples 
#' x = rlnorm(1e3, 2, 2)
#' u = round( quantile(x, 0.98) )
#' lnorm_tail(x, u, confint = FALSE)
#' 
#' 
#' @importFrom stats uniroot qnorm var
#' @export 
lnorm_tail = function(x, u, confint = FALSE, method = c("unbiased", "bootstrap", "jackknife"),
                      R = 1000, conf.level = 0.95) {
  method = match.arg(method)
  if (!(method %in% c("unbiased", "bootstrap", "jackknife"))) 
    stop("method can only be unbiased, bootstrap or jackknife")
  
  svec = s_vec(x, u)
  sigma = sapply(svec, sigma_root)
  
  if (confint == FALSE) {
    result = cbind(u, svec, sigma)
    colnames(result) = c("threshold", "s.estimate", "sigma")
    rownames(result) = NULL 
  } else {
    conf.int = ci_lnorm(x, u, method, R, conf.level)
    conf.int.a = apply(conf.int, 1:2, sigma_root)
    result = cbind(u, svec, t(conf.int), sigma, conf.int.a[1,], conf.int.a[2,])
    colnames(result) = c("threshold", "s.estimate", "s.ci1", "s.ci2",
                         "sigma", "sigma.ci1", "sigma.ci2")
    rownames(result) = NULL
  }
  return( round(result, 4) )
}



#' Plot the estimated s and the corresponding confidence intervals
#' 
#' This function produces a tail plot for the estimate \eqn{\hat{s}} over a range of thresholds for a given sample, including confidence intervals computed by one of three methods (unbiased, bootstrap or jackknife). The function also allows a choice between original and log scale.
#' 
#' @param x a vector containing the sample data. 
#' @param method the method used for computing the confidence intervals (options include unbiased variance estimator, jackknife, and bootstrap).
#' @param R the number of the bootstrap replicates.
#' @param conf.level the confidence level for the interval. 
#' @param ci.points the number of thresholds used in the calculation of the confidence intervals.
#' @param xscale the scale of the x-axis (options include "o" = original, "l" = log scale, "b" = both).
#' 
#' 
#' @return 
#' A plot showing the estimated \eqn{s(u)} versus threshold \eqn{u}, optionally on a logarithmic x-axis and including confidence intervals. Note that on the right side of the plot, one can observe the corresponding sigma values, which indicate the the standard deviation on the log-scale of the lognormal distribution associated with the estimated s-values.
#' 
#' @references 
#' Mosimann, J. E. (1970). Size allometry: size and shape variables with characterizations of the lognormal and generalized gamma distributions. Journal of the American Statistical Association, 65(330):930–945. 
#' \doi{https://doi.org/10.2307/2284599}
#' 
#' Iwashita, T. & Klar, B. (2024). A gamma tail statistic and its asymptotics. Statistica Neerlandica 78:2, 264-280.  
#' \doi{https://doi.org/10.1111/stan.12316}
#' 
#' @examples
#' \donttest{
#' x = rlnorm(2e2, 2, 2)
#' lnorm_tailplot(x, method="unbiased", xscale="o")
#' }
#' 
#' @import graphics 
#' @import stats 
#' @export 
lnorm_tailplot = function(x, method = c("unbiased", "bootstrap", "jackknife"), 
                          R = 1000, conf.level = 0.95, ci.points=101, xscale="o") {
  method = match.arg(method)
  if (!(method %in% c("unbiased", "bootstrap", "jackknife"))) 
    stop("method can only be unbiased, bootstrap or jackknife")
  x = sort(x)
  n = length(x)
  svec = s_vec(x, x[-n]) 
  ub = sort(x)[n-4]
  u.vec = seq(min(x), ub, length.out=ci.points)
  conf.int = ci_lnorm(x, u.vec, method, R, conf.level)
  
  if (xscale == "o" | xscale == "b") { #original scale
    plot.stepfun( stepfun( x[1:(n-2)], svec, right=TRUE), xaxs="i", xlim=c(min(x), ub), 
                  ylim=c(0,1), lwd=3, col="black", main="", xlab="", ylab="")
    axis(4, at=s_lnorm(c(20,5,3,2,1.5,1,0.6,0.3,0.1)), 
         labels=c(20,5,3,2,1.5,1,0.6,0.3,0.1))
    mtext(expression(hat(s)[n]), side=2, line=2)
    mtext(expression(sigma), side=4, line=2)
    mtext( "Threshold", side=1, line=2)
    abline(h=0.5,lty=3)
    lines(u.vec, conf.int[1,], lwd=2, lty=2, col="black")
    lines(u.vec, conf.int[2,], lwd=2, lty=2, col="black")
  }
  if (xscale == "l" | xscale == "b") { #log scale 
    plot( stepfun( x[1:(n-2)], svec, right=TRUE), xaxs="i", log="x", xlim=c(min(x), ub), 
          ylim=c(0,1), lwd=3, col="black", main="", xlab="", ylab="")
    axis(4, at=s_lnorm(c(20,5,3,2,1.5,1,0.6,0.3,0.1)), 
         labels=c(20,5,3,2,1.5,1,0.6,0.3,0.1))
    mtext(expression(hat(s)[n]), side=2, line=2)
    mtext(expression(sigma), side=4, line=2)
    mtext( "Threshold", side=1, line=2)
    abline(h=0.5,lty=3)
    lines(u.vec, conf.int[1,], lwd=2, lty=2, col="black")
    lines(u.vec, conf.int[2,], lwd=2, lty=2, col="black")
  }
}