#' Estimate of tail functional t and confidence intervals for t and alpha
#' 
#' This function computes the estimate of \eqn{t} and the associated confidence interval for \eqn{t} as well as \eqn{alpha}, the corresponding shape parameter under the assumption of a Pareto model according to Klar (2024). Three methods are implemented to compute the confidence intervals: a method based on the unbiased variance estimators of the underlying U-statistics and two resampling methods (jackknife and bootstrap).
#' 
#' @param x a vector containing the sample data. 
#' @param u the threshold for the computation of t.
#' @param confint a boolean value indicating whether the confidence interval should be computed.
#' @param method the method used for computing the confidence intervals (options include unbiased variance estimator, jackknife, and bootstrap).
#' @param R the number of the bootstrap replicates.
#' @param conf.level the confidence level for the interval. 
#' @param alpha.max  the upper limit of the interval to be searched for the root in an internal routine (the default value of 100 should be increased in case of error).
#'
#' @details
#' In Klar (2024) 
#' the function
#' \deqn{
#'   t_X(u) 
#'   \;=\; 
#'   \mathbb{E}\!\biggl[
#'     \frac{\lvert X_1 - X_2 \rvert}{X_1 + X_2} 
#'     \;\Big|\; 
#'     \min\{X_1, X_2\} \,\ge u
#'   \biggr]
#' }
#' is proposed as a tool for detecting Pareto-type tails, where \eqn{X_1, X_2, X} 
#' are \eqn{i.i.d.} random variables from an absolutely continuous distribution supported on \eqn{[x_m,\infty)}.
#' Theorem 1 in Klar (2024) shows that \eqn{t_X(u)} is constant in 
#' \eqn{u} if and only if \eqn{X} has a Pareto distribution.
#'
#' The estimator \eqn{\hat{t}_n\bigl(X_{(k)}\bigr)} can be computed 
#' recursively.  For \eqn{k = 2,\ldots,n-1},
#'
#' \deqn{
#'   \hat{t}_n\bigl(X_{(k)}\bigr)
#'   \;=\;
#'   \frac{n-k+2}{n-k}\,\hat{t}_n\bigl(X_{(k-1)}\bigr)
#'   \;-\;
#'   \frac{1}{\binom{\,n-k+1\,}{2}}
#'   \sum_{j=k}^{n}
#'   \frac{X_{(j)} - X_{(k-1)}}{X_{(j)} + X_{(k-1)}}\,,
#' }
#'
#' which can be evaluated efficiently starting from 
#' \eqn{\hat{t}_n\bigl(X_{(n-1)}\bigr) = \bigl(X_{(n)} - X_{(n-1)}\bigl)/\bigl(X_{(n)} + X_{(n-1)}\bigl)}, where \eqn{X_{(k)}} denotes the \eqn{k}-th order statistic.  
#' 
#' Confidence intervals for \eqn{t(u)} based on the following methods for variance estimation are also provided:
#' \itemize{
#'   \item Unbiased variance estimator
#'   \item Bootstrap resampling
#'   \item Jackknife resampling
#' }
#' 
#' A two-sided \eqn{(1 - \gamma)} confidence interval 
#' for the estimator \eqn{\hat{t}_n(u)} is :
#' \deqn{
#'   \left[
#'     \max\!\Bigl\{
#'       \hat{t}_n(u)
#'       \;-\;
#'       z_{1 - \frac{\gamma}{2}}
#'       \,\frac{\hat{\sigma}_{u}}{
#'         \sqrt{n\,U_n^{(2)}(u)}
#'       },
#'       \;0
#'     \Bigr\},
#'     \,
#'     \min\!\Bigl\{
#'       \hat{t}_n(u)
#'       \;+\;
#'       z_{1 - \frac{\gamma}{2}}
#'       \,\frac{\hat{\sigma}_{u}}{
#'         \sqrt{n\,U_n^{(2)}(u)}
#'       },
#'       \;1
#'     \Bigr\}
#'   \right],
#' }
#' where \eqn{z_{1 - \frac{\gamma}{2}} = \Phi^{-1}(1 - \tfrac{\gamma}{2})} is the appropriate quantile of the standard normal distribution, \eqn{\hat{\sigma}_u} is an estimator of the standard deviation of \eqn{c\,\hat{t}_n(u)}, for a constant c specified in section 4.1. of Klar (2024), and
#' \eqn{U_n^{(2)}(u)} is a U-statistic given by
#'       \deqn{
#'         U_n^{(2)}(u)
#'         \;=\;
#'         \frac{2}{n\,(n-1)}
#'         \sum_{i = 1}^n
#'         (n - i)
#'         1\{X_{(i)} \,\ge\, u\}.
#'       }
#' 
#' @return A matrix containing:
#' \item{threshold}{The value of the threshold u.}
#' \item{t.estimate}{Estimate of the tail functional t.}
#' \item{t.ci1}{The lower bound of the confidence interval for t (if `confint = TRUE`).}
#' \item{t.ci2}{The upper bound of the confidence interval for t (if `confint = TRUE`).}
#' \item{alpha}{Estimate of the shape parameter under a Pareto model.}
#' \item{alpha.ci1}{The lower bound of the confidence interval for alpha (if `confint = TRUE`).}
#' \item{alpha.ci2}{The upper bound of the confidence interval for alpha (if `confint = TRUE`).}
#'
#' 
#' @references 
#' Klar, B. (2024). A Pareto tail plot without moment restrictions. *The American Statistician*. \doi{https://doi.org/10.1080/00031305.2024.2413081}
#'
#' @examples 
#' x <- actuar::rpareto1(1e3, shape=1, min=1)
#' pareto_tail(x, round( quantile(x, c(0.1, 0.5, 0.75, 0.9, 0.95, 0.99)) ), confint = FALSE) 
#' 
#' @importFrom stats uniroot qnorm var
#' @export 
pareto_tail = function(x, u, confint=FALSE, method = c("unbiased", "bootstrap", "jackknife"), 
                       R = 1000, conf.level = 0.95, alpha.max = 100) {
    method = match.arg(method)
    if (!(method %in% c("unbiased", "bootstrap", "jackknife"))) 
        stop("method can only be unbiased, bootstrap or jackknife")
    if(any(x < 0))
        stop("The sample values must be positive")
        
    tvec = t_vec(x, u)
    alpha = sapply(tvec, alpha_root, alpha.max=alpha.max)

    if (confint == FALSE) {
        result = cbind(u, tvec, alpha)
        colnames(result) = c("threshold", "t.estimate", "alpha")
        rownames(result) = NULL 
    } else {
        conf.int = ci(x, u, method, R, conf.level)
        conf.int.a = apply(conf.int, 1:2, alpha_root, alpha.max=alpha.max)
        result = cbind(u, tvec, t(conf.int), alpha, conf.int.a[2,], conf.int.a[1,])
        colnames(result) = c("threshold", "t.estimate", "t.ci1", "t.ci1",
                            "alpha", "alpha.ci1", "alpha.ci2")
        rownames(result) = NULL
        }
    return( round(result, 4) )
}

#' Plot the estimated t and the corresponding confidence intervals
#' 
#' This function produces a tail plot for the estimate \eqn{\hat{t}} over a range of thresholds for a given sample, including confidence intervals computed by one of three methods (unbiased, bootstrap or jackknife). The function also allows a choice between original and log scale.
#' 
#' @param x a vector containing the sample data. 
#' @param method the method used for computing the confidence intervals (options include unbiased variance estimator, jackknife, and bootstrap).
#' @param R the number of the bootstrap replicates.
#' @param conf.level the confidence level for the interval. 
#' @param ci.points the number of thresholds used in the calculation of the confidence intervals.
#' @param xscale the scale of the x-axis (options include "o" = original, "l" = log scale, "b" = both).
#' 
#' @details
#' For more details about the estimator \eqn{\hat{t}} and the computation of the confidence intervals see \link{pareto_tail}. 
#' 
#' @return 
#' A plot showing the estimated \eqn{t(u)} versus threshold \eqn{u}, optionally on a logarithmic x-axis and including confidence intervals. Note that on the right side of the plot, one can observe the corresponding alpha values, which indicate the shape parameter of the Pareto distribution associated with the estimated t-values.
#' 
#' @references 
#' Klar, B. (2024). A Pareto tail plot without moment restrictions. *The American Statistician*. \doi{https://doi.org/10.1080/00031305.2024.2413081}
#' 
#' @examples
#' \donttest{
#' x <- actuar::rpareto1(1e3, shape=1, min=1)
#' pareto_tailplot(x, method="unbiased", xscale="o")
#' }
#' 
#' @import graphics 
#' @import stats 
#' @export 
pareto_tailplot = function(x, method = c("unbiased", "bootstrap", "jackknife"), 
                      R = 1000, conf.level = 0.95, ci.points=101, xscale="b") {
    method = match.arg(method)
    if (!(method %in% c("unbiased", "bootstrap", "jackknife"))) 
        stop("method can only be unbiased, bootstrap or jackknife")
    if(any(x < 0))
        stop("The sample values must be positive")
        
    n = length(x)
    ub = sort(x)[n-4]
    trec = t_rec(x)
    u.vec = seq(min(x), ub, length.out=ci.points)
    conf.int = ci(x, u.vec, method, R, conf.level)
    
    if (xscale == "o" | xscale == "b") { #original scale
        plot.stepfun( stepfun(trec[1:(n-2),1], trec[,2], right=TRUE), xaxs="i",
                    xlim=c(min(x), ub), ylim=c(0,1), lwd=3, col="black", main="", xlab="", ylab="")
        axis(4, at=t_pareto(c(10,3:1,0.5,0.25,0.1)), labels=c(10,3:1,0.5,0.25,0.1))
        mtext(expression(hat(t)[n]), side=2, line=2)
        mtext(expression(alpha), side=4, line=2)
        mtext( "Threshold", side=1, line=2)
        abline(h=2*log(2)-1,lty=3); abline(h=3-4*log(2),lty=3)
        lines(u.vec, conf.int[1,], lwd=2, lty=2, col="black")
        lines(u.vec, conf.int[2,], lwd=2, lty=2, col="black")
    }
    if (xscale == "l" | xscale == "b") { #log scale 
        plot( stepfun( trec[1:(n-2),1], trec[,2], right=TRUE), xaxs="i", log="x",
            xlim=c(min(x), ub), ylim=c(0,1), lwd=3, col="black", main="", xlab="", ylab="")
        axis(4, at=t_pareto(c(10,3:1,0.5,0.25,0.1)), labels=c(10,3:1,0.5,0.25,0.1))
        mtext(expression(hat(t)[n]), side=2, line=2)
        mtext(expression(alpha), side=4, line=2)
        mtext( "Threshold", side=1, line=2)
        abline(h=2*log(2)-1,lty=3); abline(h=3-4*log(2),lty=3)
        lines(u.vec, conf.int[1,], lwd=2, lty=2, col="black")
        lines(u.vec, conf.int[2,], lwd=2, lty=2, col="black")
    }
}
