#' Estimate of tail functional g and confidence intervals for g and alpha
#' 
#' This function computes the estimate of \eqn{g} and the associated confidence interval for \eqn{g} as well as \eqn{alpha}, the corresponding shape parameter under the assumption of a gamma model, according to Iwashita and Klar (2024). Three methods are implemented to compute the confidence intervals: a method based on the unbiased variance estimators of the underlying U-statistics, and two resampling methods (jackknife and bootstrap).
#' 
#' @param x a vector containing the sample data. 
#' @param d the threshold for the computation of g.
#' @param confint a boolean value indicating whether a confidence interval should be computed.
#' @param method the method used for computing the confidence intervals (options include unbiased variance estimator, jackknife, and bootstrap).
#' @param R the number of the bootstrap replicates.
#' @param conf.level the confidence level for the interval. 
#' 
#' @details 
#' The function \eqn{g} introduced by Asmussen and Lehtomaa (2017) is used to distinguish between 
#' log-concave and log-convex tail behavior. It is defined as:
#' 
#' \deqn{ g(d) = E\left[ \frac{|X_1 - X_2|}{X_1 + X_2} \bigg| X_1 + X_2 > d \right] }
#' 
#' where \eqn{X_1, X_2} are independent and identically distributed (i.i.d.) positive random variables.
#' For gamma distributions, \eqn{g} takes a constant value, making it a useful tool for detecting gamma-tailed distributions.
#'
#' This function estimates \eqn{g(d)} using U-statistics. The estimator \eqn{\hat{g}(d)} is given by:
#' 
#' \deqn{
#' \hat{g}(d) = \frac{ U^{(1)}_n (d) }{ U^{(2)}_n (d) }, \quad d > 0,  
#' }
#' 
#' where
#' 
#' \deqn{
#' U^{(1)}_n (d) = \frac{2}{n(n-1)} \sum_{1 \leq i < j \leq n} \frac{|X_i - X_j|}{X_i + X_j} 1(X_i + X_j > d),
#' }
#' 
#' \deqn{
#' U^{(2)}_n (d) = \frac{2}{n(n-1)} \sum_{1 \leq i < j \leq n} 1(X_i + X_j > d).
#' }
#'
#' Confidence intervals for \eqn{g(d)}, based on the following variance estimation methods, are also provided:
#' - Unbiased Variance Estimator
#' - Bootstrap Resampling
#' - Jackknife Resampling
#'
#' The \eqn{(1-\gamma)} confidence interval for \eqn{\hat{g}_{n}(d)} 
#' is given by:
#' 
#' \deqn{
#'   \left[
#'     \max\!\Bigl\{
#'       \hat{g}_{n}(d)\;-\;
#'       z_{1 - \gamma/2}
#'       \,\frac{\hat{\sigma}_{d}}{
#'         \sqrt{n\,U^{(2)}_{n}(d)}
#'       },
#'       \;0
#'     \Bigr\},
#'     \;\;
#'     \min\!\Bigl\{
#'       \hat{g}_{n}(d)\;+\;
#'       z_{1 - \gamma/2}
#'       \,\frac{\hat{\sigma}_{d}}{
#'         \sqrt{n\,U^{(2)}_{n}(d)}
#'       },
#'       \;1
#'     \Bigr\}
#'   \right].
#' }{
#' CI formula for g(d).  
#' }
#' Here, 
#' \eqn{z_{1 - \gamma/2} = \Phi^{-1}(1 - \tfrac{\gamma}{2})} is the 
#' appropriate quantile of the standard normal distribution and \eqn{\hat{\sigma}_d} is an estimator of the standard deviation based on one of the methods above.
#'
#' @return A matrix containing:
#' \item{threshold}{The value of the threshold d.}
#' \item{g.estimate}{Estimate of g.}
#' \item{g.ci1}{The lower bound of the confidence interval for g (if `confint = TRUE`).}
#' \item{g.ci2}{The upper bound of the confidence interval for g (if `confint = TRUE`).}
#' \item{alpha}{Estimate of the shape parameter under a gamma model.}
#' \item{alpha.ci1}{The lower bound of the confidence interval for alpha (if `confint = TRUE`).}
#' \item{alpha.ci2}{The upper bound of the confidence interval for alpha (if `confint = TRUE`).}
#'
#' @references 
#' Iwashita, T. & Klar, B. (2024). A gamma tail statistic and its asymptotics. Statistica Neerlandica 78:2, 264-280.  
#' \doi{https://doi.org/10.1111/stan.12316}
#' 
#' Asmussen, S. & Lehtomaa, J. (2017). Distinguishing Log-Concavity from Heavy Tails. *Risks 2017, 5, 10*. \doi{https://doi.org/10.3390/risks5010010}
#'
#'
#' @examples 
#' x <- rgamma(100, shape = 2, scale = 1)
#' gamma_tail(x, d = 2, confint = FALSE, method = "unbiased", R = 1000)
#' 
#' @importFrom stats var qnorm uniroot
#' @export 
gamma_tail = function(x, d, confint = FALSE, method = c("unbiased", "bootstrap", "jackknife"),
                      R = 1000, conf.level = 0.95) {
    method = match.arg(method)
    if (!(method %in% c("unbiased", "bootstrap", "jackknife"))) 
        stop("method can only be unbiased, bootstrap or jackknife")
    if(any(x < 0))
        stop("The sample values must be positive")

    gvec = g_vec(x, d)
    alpha = sapply(gvec, alpha_root_gamma)
    
    if (confint == FALSE) {
        result = cbind(d, gvec, alpha)
        colnames(result) = c("threshold", "g.estimate", "alpha")
        rownames(result) = NULL 
    } else {
        conf.int = ci_gamma(x, d, method, R, conf.level)
        conf.int.a = apply(conf.int, 1:2, alpha_root_gamma)
        result = cbind(d, gvec, t(conf.int), alpha, conf.int.a[2,], conf.int.a[1,])
        colnames(result) = c("threshold", "g.estimate", "g.ci1", "g.ci1",
                            "alpha", "alpha.ci1", "alpha.ci2")
        rownames(result) = NULL
    }
    return( round(result, 4) )
}



#' Plot the estimated g and the corresponding confidence intervals
#' 
#' This function produces a tail plot for the estimate \eqn{\hat{g}} over a range of thresholds for a given sample, including confidence intervals computed by one of three methods (unbiased, bootstrap or jackknife). The function also allows a choice between original and log scale.
#' 
#' @param x a vector containing the sample data. 
#' @param method the method used for computing the confidence intervals (options include unbiased variance estimator, jackknife, and bootstrap).
#' @param R the number of the bootstrap replicates.
#' @param conf.level the confidence level for the interval. 
#' @param ci.points the number of thresholds used in the calculation of the confidence intervals.
#' @param xscale the scale of the x-axis (options include "o" = original, "l" = log scale, "b" = both).
#' 
#' @details
#' For more details about the estimator \eqn{\hat{g}} and the computation of the confidence intervals see \link{gamma_tail}. 
#' 
#' @return 
#' A plot showing the estimated \eqn{g(d)} versus threshold \eqn{d}, optionally on a logarithmic x-axis and including confidence intervals.
#' 
#' @references 
#' Iwashita, T. & Klar, B. (2024). A gamma tail statistic and its asymptotics. Statistica Neerlandica 78:2, 264-280.  
#' \doi{https://doi.org/10.1111/stan.12316}
#' 
#' @examples
#' \donttest{
#' x <- rgamma(2e2, 0.5, 0.2)
#' gamma_tailplot(x, method="unbiased", xscale="o")
#' }
#' 
#' @import graphics
#' @import stats
#' @export 
gamma_tailplot = function(x, method = c("unbiased", "bootstrap", "jackknife"), 
                          R = 1000, conf.level = 0.95, ci.points=101, xscale="o") {
    method = match.arg(method)
    if (!(method %in% c("unbiased", "bootstrap", "jackknife"))) 
        stop("method can only be unbiased, bootstrap or jackknife")
    if(any(x < 0))
        stop("The sample values must be positive")
        
    x = sort(x)
    n = length(x)
    gvec = g_vec(x, x[-n]) 
    ub = sort(x)[n-4]
    u.vec = seq(min(x), ub, length.out=ci.points)
    conf.int = ci_gamma(x, u.vec, method, R, conf.level)
    
    #original scale
    if (xscale == "o" | xscale == "b") { 
        plot.stepfun( stepfun( x[1:(n-2)], gvec, right=TRUE), xaxs="i",
                    xlim=c(min(x), ub), ylim=c(0,1), lwd=3, col="black", main="", xlab="", ylab="")
        axis(4, at=g_gamma(c(50,10,3,1,0.5,0.25,0.1,0.01)), labels=c(50,10,3,1,0.5,0.25,0.1,0.01))
        mtext(expression(hat(g)[n]), side=2, line=2)
        mtext(expression(alpha), side=4, line=2)
        mtext( "Threshold", side=1, line=2)
        abline(h=0.5,lty=3)
        lines(u.vec, conf.int[1,], lwd=2, lty=2, col="black")
        lines(u.vec, conf.int[2,], lwd=2, lty=2, col="black")
    }
    
    #log scale
    if (xscale == "l" | xscale == "b") {  
        plot( stepfun( x[1:(n-2)], gvec, right=TRUE), xaxs="i", log="x",
            xlim=c(min(x), ub), ylim=c(0,1), lwd=3, col="black", main="", xlab="", ylab="")
        axis(4, at=g_gamma(c(50,10,3,1,0.5,0.25,0.1,0.01)), labels=c(50,10,3,1,0.5,0.25,0.1,0.01))
        mtext(expression(hat(t)[n]), side=2, line=2)
        mtext(expression(alpha), side=4, line=2)
        mtext( "Threshold", side=1, line=2)
        abline(h=0.5,lty=3)
        lines(u.vec, conf.int[1,], lwd=2, lty=2, col="black")
        lines(u.vec, conf.int[2,], lwd=2, lty=2, col="black")
    }
}
