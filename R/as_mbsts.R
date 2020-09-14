#' Definition and estimation of a Multivariate Bayesian Structural Time Series model (MBSTS)
#'
#' The function creates a multivariate Bayesian structural time series model. It then estimates
#' the model, samples from the joint posterior distribution of its parameters, and outputs an
#' object of class \code{mbsts}.
#'
#' @param y t x d data.frame (or matrix) of observations, where d is the number of time series
#' in the multivariate model.
#' @param components Character vector specifying the components of the multivariate structural
#' time series model. Possible values are c("trend", "slope", "seasonal", "cycle").
#' @param seas.period Length of the seasonal pattern, if present.
#' @param cycle.period Length of the cycle pattern, if present.
#' @param X Optional t x N data frame (or matrix) of N predictors.
#' @param H P x P variance-covariance matrix of the regression coefficients. Set by
#' default to H = c(X'X)^(-1) which is akin to the Zellner's g-prior. The value of
#' the scaling factor is set to \code{c = 1}. Alternative priors could be
#' H = c*diag((X'X)^(-1)) or H = c*I. See also Smith & Kohn, 1995 that suggest
#' setting \code{c} in the range [10,1000].
#' @param nu0.r Degrees of freedom of the Inverse-Wishart prior for each element of
#' Sigma.r, a vector of errors for state r.
#' Set by default to d + 2 (must be greater than d - 1).
#' @param s0.r Scale matrix of the Inverse-Wishart prior for each Sigma.r, a vector
#' of errors for state r. Must be a (d x d) positive definite. Default set to the
#' variance-covariance matrix of y multiplied by a scaling factor of 0.01.
#' @param nu0.eps Degrees of freedom of the Inverse-Wishart prior for Sigma.eps,
#' a vector of observation errors for each time series. Set by default to d + 2
#' (must be greater than d - 1).
#' @param s0.eps Scale matrix of the Inverse-Wishart prior for Sigma.eps, a vector
#' of observation errors for each time series. Must be a (d x d) positive definite.
#' Default set to Default set to the variance-covariance matrix of y multiplied by
#' a scaling factor of 0.01..
#' @param niter Number of MCMC iterations.
#' @param burn Desired burn-in, set by default to 0.1 * \code{niter}.
#' @param ping A status message is printed every \code{ping} iteration. Default
#'   set to 0.1 * \code{niter}. Set to 0 to not track the status.
#'
#' @return An object of class 'mbsts' which is a list with the following components:
#' \describe{
#'   \item{eta.samples}{(\code{niter}- \code{burn}) draws from the distribution of eta_r.}
#'   \item{eps.samples}{(\code{niter}- \code{burn}) draws from the distribution of eps.}
#'   \item{states.samples}{(\code{niter}- \code{burn}) draws from p(alpha_t | Y_{1:T}).}
#'   \item{Sigma.r}{(\code{niter}- \code{burn}) draws from the posterior distribution of Sigma.r.}
#'   \item{Sigma.eps}{(\code{niter}- \code{burn}) draws from the posterior distribution of Sigma.eps.}
#'   \item{Z.beta}{(\code{niter}- \code{burn}) x P matrix of the models selected at each
#'   iteration (if a matrix of predictors is provided).}
#'   \item{beta}{ P x d x (\code{niter}- \code{burn}) ) array of the draws from the posterior
#'   distribution of the regression coefficient matrix (if a matrix of predictors is provided).}
#'   \item{X}{Predictor matrix (if provided).}
#'   \item{y}{Matrix of observations.}
#'   \item{Z}{(d x m) selection matrix of the observation equation.}
#'   \item{T}{(m x m) matrix of the state equation.}
#'   \item{R}{(m x r) matrix selecting the state disturbances.}
#'   \item{niter}{Number of mcmc iterations.}
#'   \item{burn}{Burn-in.}
#'   }
#' @export
#'
#' @examples
#' ## Example 1 : local level + seasonal (d = 3)
#' y <- cbind(seq(0.5,200,by=0.5)*0.1 + rnorm(400),
#'            seq(100.25,200,by=0.25)*0.05 + rnorm(400),
#'            rnorm(400, 5,1))
#' mbsts.1 <- as.mbsts(y = y, components = c("trend", "seasonal"), seas.period = 7,
#'                     s0.r = diag(3), s0.eps = diag(3), niter = 100, burn = 10)
#'
#' ## Example 2 : local level + seasonal + covariates (d = 2)
#' y <- cbind(rnorm(100), rnorm(100, 2, 3))
#' X <- cbind(rnorm(100, 0.5, 1) + 5, rnorm(100, 0.2, 2) - 2)
#' mbsts.2 <- as.mbsts(y = y, components = c("trend", "seasonal"), , seas.period = 7,
#'                     X = X, s0.r = diag(2), s0.eps = diag(2), niter = 100, burn = 10)

as.mbsts <- function(y, components, seas.period = NULL, cycle.period = NULL, X = NULL,
                     H = NULL, nu0.r = NULL, s0.r = 0.01 * var(y, na.rm = T), nu0.eps = NULL,
                     s0.eps = 0.01 * var(y, na.rm = T), niter, burn, ping = NULL){

  ## Parameter checks
  if(!is.matrix(y) && !is.data.frame(y)) stop("`y` must be a matrix or data.frame")
  if(all(!components %in% c("trend", "slope", "seasonal", "cycle")))
    stop("`components` must be one of 'trend', 'slop', 'seasonal', or 'cycle'")
  if(!missing(seas.period) && (!is.numeric(seas.period) ||length(seas.period) != 1))
    stop("`seas.period` must be a numeric vector of length one")
  if(!missing(cycle.period) && (!is.numeric(cycle.period) ||length(cycle.period) != 1))
    stop("`cycle.period` must be a numeric vector of length one")
  if(!missing(X)) {
    if(!is.matrix(X) && !is.data.frame(X)) stop("`X` must be a matrix or data.frame")
    if(nrow(X) != nrow (y)) stop("nrow(X) != nrow(X)")
  }
  if(!missing(s0.r) && (!is.matrix(s0.r) || !all(dim(s0.r) == ncol(y))))
    stop("`s0.r` must be a d x d matrix")
  if(!missing(nu0.eps)  && (!is.numeric(nu0.eps) || length(nu0.eps) != 1 || nu0.eps <= nrow(y)))
    stop("`nu0.eps` must be a length-one numeric vector value >= to d")
  if(!missing(s0.eps) && (!is.matrix(s0.eps) || !all(dim(s0.r) == ncol(y))))
    stop("`s0.eps` must be a d x d matrix")
  if(!missing(niter) && !is.numeric(niter)) stop("`niter` must be a length-one numeric vector")
  if(!missing(ping) && !is.numeric(ping)) stop("`ping` must be a length-one numeric vector")
  # Model definition
  Smodel <- model(y = y, components = components, seas.period = seas.period, cycle.period = cycle.period)

  # Model estimation (MCMC)
  est <- mcmc(Smodel = Smodel, X = X, H = H, nu0.r = nu0.r, s0.r = s0.r, nu0.eps = nu0.eps,
              s0.eps = s0.eps, niter = niter, burn = burn, ping = ping)

  class(est) <- "mbsts"
  return(est)
}
