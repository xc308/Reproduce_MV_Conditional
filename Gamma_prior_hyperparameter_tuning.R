### ------------------
### Distribution tuning
### ------------------

#' @title Quantile difference when tuning Gamma distribution hyper-parameters
#'
#' @description See details below.
#'
#' @param pars the shape and rate parameters (in that order) of the Gamma distribution
#' @param p5 the 5-th percentile of the desired error distribution
#' @param p95 the 95-th percentile of the desired error distribution
#' @return A measure of discrepancy between the Gamma distribution with parameters \code{par} and the desired 5-th 95-th percentiles of the error signal.
#' @details Assume that in our model there is an error term (for example, observation error) 
#'   \eqn{e \sim \mathcal{N}(0,\sigma^2)} and we are required to estimate the precision, 
#'   \eqn{\sigma^{-2}}. Assume further that we have sufficient prior knowledge to know that
#'   \eqn{F_e(0.05)} = \code{p5} and \eqn{F_e(0.95)} = \code{p95}, 
#'   where \code{F_e} is the cumulative distribution function of \eqn{e} and 
#'   \code{p5} and \code{p95} are known. This function helps configure a prior distribution 
#'   over \eqn{\sigma^{-2}} which reflects this prior belief. Specifically, 
#'   let \eqn{\Lambda_{\alpha,\beta}} be the cumulative distribution function of the Gamma distribution 
#'   with parameters \eqn{\alpha} and \eqn{\beta}. Then this function returns the discrepancy
#'   \deqn{\frac{\Lambda_{\alpha,\beta}(0.05) - 1/(p95)^2}{1/(p95)^2}  + \frac{\Lambda_{\alpha,\beta}(0.95) - 1/(p5)^2}{1/(p5)^2}}
#'
#'
#' For example, we can see how representative a Gamma distribution over 
#'   the precision with shape and rate parameters in \code{pars} is of our prior belief 
#'   that \code{p5} and \code{p95} are 1 and 5, respectively, by calling 
#'   \code{Findalphabeta_gamma(pars,1,5)}. This function can hance be passed to an optimisation routine 
#'   to find better values for \code{pars}.
#' @keywords Gamma distribution, prior elicitation
#' @export
#' @examples
#'
#' require(actuar)
#' # Find the Gamma distribution over the precision corresponding to the
#' # prior belief of the error (1/sqrt(precision)) lying between p5=2, p95=5
#' initpars <- c(5,0.1)
#' hyp_pars <- optim(par=initpars,Findalphabeta_gamma, p5=1, p95=5)
#'
#' # Now simulate from a Gamma with these parameters and verify quantiles
#' X <- rgamma(shape = hyp_pars$par[1], rate = hyp_pars$par[2],n=10000)
#' print( quantile(1/sqrt(X),c(0.05,0.95)))
Findalphabeta_gamma <- function(pars,p5,p95) {
  if(any(pars<0)) {
    return(Inf)
  } else {
    return( sum((qgamma(c(0.05,0.95),shape=pars[1],rate=pars[2]) - c(1/p95^2,1/p5^2))^2/c(1/p95^2,1/p5^2)))
  }
}


