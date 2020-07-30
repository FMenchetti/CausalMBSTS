
#' Multivariate structural time series model definition and estimation
#'
#' @param y t x d data.frame (or matrix) of observations, where d is the number of time series in the multivariate model.
#' @param components Character vector specifying the components of the multivariate structural time series model. Possible values in c("trend", "slope", "seasonal", "cycle").
#' @param seas.period Length of the seasonal pattern.
#' @param cycle.period Length of the cycle pattern.
#' @param X Optional t x N data frame (or matrix) of predictors.
#' @param H P x P variance-covariance matrix of the regression coefficients. Set by default to Zellner's g-prior, H = (X'X)^(-1).
#' @param nu0.r Degrees of freedom of the Inverse-Wishart prior for each element of Sigma.r, a vector of errors for state r.
#' Set by default to d + 2 (must be greater than d - 1).
#' @param s0.r Scale matrix of the Inverse-Wishart prior for each Sigma.r, a vector of errors for state r. Must be a (d x d)
#' positive definite. Default set to ???.
#' @param nu0.eps Degrees of freedom of the Inverse-Wishart prior for Sigma.eps, a vector of observation errors for each time
#' series. Set by default to d + 2 (must be greater than d - 1).
#' @param s0.eps Scale matrix of the Inverse-Wishart prior for Sigma.eps, a vector of observation errors for each time series.
#' Must be a (d x d) positive definite. Default set to ???.
#' @param niter Number of MCMC iterations.
#' @param burn Desired burn-in, set by default to 0.1 * niter.
#' @param ping A status message is printed every 'ping' iteration. Default set to 0.1 * 'niter'.
#'
#' @return An object of class 'mbsts' which is a list with the following components
#' \describe{
#'   \item{eta.samples}{'niter' draws from the distribution of eta_r.}
#'   \item{eps.samples}{'niter' draws from the distribution of eps.}
#'   \item{states.samples}{draws from p(alpha_t | Y_{1:T}).}
#'   \item{Sigma.r}{'niter' draws from the posterior distribution of Sigma.r.}
#'   \item{sigma.eps}{'niter' draws from the posterior distribution of Sigma.eps.}
#'   \item{Z.beta}{('niter'- 'burn') x P matrix of the models selected at each iteration.}
#'   \item{beta}{ P x d x ('niter' - 'burn') ) array of the draws from the posterior distribution of the regression coefficient matrix.}
#'   \item{X}{Predictor matrix.}
#'   \item{y}{Matrix of observations.}
#'   \item{Z}{(1 x m) selection matrix of the observation equation.}
#'   \item{T}{(m x m) matrix of the state equation.}
#'   \item{R}{(1 x r) matrix selecting the state disturbances.}
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
#' mbsts.1 <- as.mbsts(y = y, components = c("trend", "seasonal"), seas.period = 7, s0.r = diag(3), s0.eps = diag(3), niter = 100, burn = 10)
#'
#' ## Example 2 : local level + seasonal + covariates (d = 2)
#' y <- cbind(rnorm(100), rnorm(100, 2, 3))
#' X <- cbind(rnorm(100, 0.5, 1) + 5, rnorm(100, 0.2, 2) - 2)
#' mbsts.2 <- as.mbsts(y = y, components = c("trend", "seasonal"), , seas.period = 7, X = X, s0.r = diag(2), s0.eps = diag(2), niter = 100, burn = 10)

as.mbsts <- function(y, components, seas.period = NULL, cycle.period = NULL, X = NULL, H = NULL, nu0.r = NULL, s0.r, nu0.eps = NULL, s0.eps = NULL, niter, burn, ping = NULL){

  # Model definition
  Smodel <- model(y = y, components = components, seas.period = seas.period, cycle.period = cycle.period)

  # Model estimation (MCMC)
  est <- mcmc(Smodel = Smodel, X = X, H = H, nu0.r = nu0.r, s0.r = s0.r, nu0.eps = nu0.eps, s0.eps = s0.eps, niter = niter, burn = burn, ping = ping)

  class(est) <- "mbsts"
  return(est)
}
