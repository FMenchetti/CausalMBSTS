######################################################################################
######################################################################################
####  Author:           Fiammetta Menchetti                                       ####
####                                                                              ####
####  Date last update: 2020-06-23                                                ####
####                                                                              ####
####  Content:          MCMC for a given MBSTS model                              ####
####                                                                              ####
####  Main function :   mbsts.mcmc                                                ####
####  Dependencies:     lpy.X , block.m                                           ####
####                                                                              ####
######################################################################################
######################################################################################


#' MCMC for a given MBSTS model
#'
#' MCMC to sample from the joint posterior of model parameters in an mbsts model.
#'
#' @import KFAS
#' @importFrom CholWishart rInvWishart
#' @importFrom CholWishart lmvgamma
#' @importFrom Matrix bdiag
#' @importFrom MixMatrix rmatrixnorm
#' @param Smodel A multivariate state space model of class 'SSModel'.
#' @param X t x P optional matrix of predictors.
#' @param H P x P variance-covariance matrix of the regression coefficients. Set by default to H = (X'X)^(-1).
#' @param nu0.k Degrees of freedom of the Inverse-Wishart prior for each Sigma_k. Set by default to n0.k = d + 2 where d is the number of time series in the multivariate model.
#' @param s0.k Scale matrix of the Inverse-Wishart prior for each Sigma_k.
#' @param nu0.eps Degrees of freedom of the Inverse-Wishart prior for Sigma_eps. Set by default to d + 2.
#' @param s0.eps Scale matrix of the Inverse-Wishart prior for Sigma.eps.
#' @param niter Number of MCMC iteration.
#' @param burn Desidered burn-in, set by default to 0.1 * niter.
#' @param ping A status message it's printed every 'ping' iteration, defaults to 0.1 * 'niter'.
#'
#' @return An object of class 'mbsts' which is a list with the following components
#' \describe{
#'   \item{eta.samples}{'niter' draws from the distribution of eta_k.}
#'   \item{eps.samples}{'niter' draws from the distribution of eps.}
#'   \item{states.samples}{draws from p(alpha_t | Y_{1:T}).}
#'   \item{sigma.eta}{'niter' draws from the posterior distribution of Sigma_eta.}
#'   \item{sigma.eps}{'niter' draws from the posterior distribution of Sigma_eps.}
#'   \item{Z.beta}{('niter'- 'burn') x P matrix of the models selected at each iteration.}
#'   \item{beta}{ P x d x ('niter' - 'burn') ) array of the draws from the posterior distribution of the regression coefficient matrix.}
#'   \item{X}{?}
#'   \item{y}{?}
#'   \item{Z}{?}
#'   \item{T}{?}
#'   \item{R}{?}
#'   \item{niter}{?}
#'   \item{burn}{?}
#' }
#' @export
#'
#' @examples
#' ## Example 1 : local level + seasonal
#' y <- cbind(seq(0.5,200,by=0.5)*0.1 + rnorm(400),
#'            seq(100.25,200,by=0.25)*0.05 + rnorm(400),
#'            rnorm(400, 5,1))
#' model.1 <- SSModel(y ~ SSMtrend(degree = 1, Q = matrix(NA)) + SSMseasonal(period=7, Q = matrix(NA)))
#' mcmc.1 <- mbsts.mcmc(model.1, s0.k = diag(3), s0.eps = diag(3), niter = 100, burn = 10)
#'
#' ## Example 2 : local level + seasonal + covariates
#' y <- cbind(rnorm(100), rnorm(100, 2, 3))
#' X <- cbind(rnorm(100, 0.5, 1) + 5, rnorm(100, 0.2, 2) - 2)
#' model.2 <- SSModel(y ~ SSMtrend(degree = 1, Q = matrix(NA,2,2)) + SSMseasonal(period=7, Q = matrix(NA,2,2)))
#' mcmc.2 <- mbsts.mcmc(model.2, X = X, s0.k = diag(2), s0.eps = diag(2), niter = 100, burn = 10)

mbsts.mcmc <- function(Smodel, X = NULL, H = NULL, nu0.k = NULL, s0.k, nu0.eps = NULL, s0.eps, niter, burn,
                       ping = NULL) {

  # MCMC to sample from the joint posterior of model parameters in an mbsts model.
  # For now, the case without contemporaneous predictors variables is not allowed
  #
  # Args:
  #   Smodel  : a multivariate state space model of class 'SSModel'
  #   X       : an optional T x N matrix of predictors
  #   H       : N x N variance-covariance matrix between regression coefficients.
  #             The default is Zellner's g-prior, H = (X'X)^(-1)
  #   nu0.k   : degrees of freedom of the Inverse-Wishart prior for each Sigma_k.
  #             The default is the smallest integer such that the expectation of eta_k exists,
  #             that is, nu0.k = p + 2 where p is the number of time series in the
  #             multivariate model
  #   nu0.eps : degrees of freedom of the Inverse-Wishart prior for Sigma_eps, the default is p+2.
  #   s0.k    : Scale matrix of the Inverse-Wishart prior for each Sigma_k.
  #   s0.eps  : Scale matrix of the Inverse-Wishart prior for Sigma.eps
  #   niter   : number of MCMC iteration
  #   burn    : desidered burn-in, set by default to 0.1 * niter
  #   ping    : logical, if TRUE a status message it's printed every iterations decile.
  #             Defaults to TRUE.
  #
  # Value:
  #   An object of class 'mbsts' which is a list with the following components
  #
  #   Smodel      : Resulting 'SSModel' (is it really useful? maybe I can save only the needed matrices, Z,R,T)
  #   eta.samples : 'niter' draws from the distribution of eta_k
  #   eps.samples : 'niter' draws from the distribution of eps
  #   states.samples : draws from p(\alpha_t | Y_{1:T})
  #   sigma.eta   : 'niter' draws from the posterior distribution of Sigma_eta
  #   sigma.eps   : 'niter' draws from the posterior distribution of Sigma_eps
  #   Z           : 'niter' x N selection matrix
  #   beta        : N x p x 'niter' array containing the draws from the posterior distribution
  #                 of the regression coefficient matrix, p(beta | Sigma_eps, z)
  #   ...         : maybe not needed

  ### Dimensionalities & other inputs
  p <- dim(y)[2]  # number of time series

  # set default H (Zellner's g-prior)
  if (!is.null(X) & is.null(H)) {
    H <- inv(crossprod(X))
  }

  # set default ping
  if (is.null(ping)) {
    ping <- 0.1 * niter
  }
  sequence <- seq(ping, niter, by = ping)

  # set default nu0.k and nu0.eps
  if (is.null(nu0.k)) {
    nu0.k <- p + 2
  }
  if (is.null(nu0.eps)) {
    nu0.eps <- p + 2
  }

  # model definition
  Smodel$H[, , 1] <- rInvWishart(1, df = nu0.eps, Sigma = s0.eps)[, , 1]
  Q <- list()
  for (i in 1:length(unique(attr(Smodel, "eta_types")))) {
    Q[[i]] <- rInvWishart(1, df = nu0.k, Sigma = s0.k)[, , 1]
  }
  Q <- block.m(Q)
  Smodel$Q[, , 1] <- Q

  # get dim
  y <- Smodel$y
  t <- dim(y)[1]  # number of time points
  K <- dim(Smodel$R)[2]  # tot number of state disurbances
  k <- K/p  # number of state disturbances for each time series
  M <- dim(Smodel$T)[1]  # tot number of states
  m <- M/p  # number of states for each time series

  ### Empty arrays to store MCMC iterations
  eta.samples <- array(NA, c(nrow(y), K, niter - burn))
  eta.names <- unique(c(t(outer(attr(Smodel, "eta_types"), attr(Smodel$y, "dimnames")[[2]], paste))))
  colnames(eta.samples) <- eta.names

  Sigma.states <- array(NA, c(K, K, niter - burn))
  colnames(Sigma.states) <- attr(Smodel, "eta_types")
  rownames(Sigma.states) <- attr(Smodel, "eta_types")

  eps.samples <- array(NA, c(nrow(y), p, niter - burn))
  Sigma.obs <- array(NA, c(p, p, niter - burn))

  states.samples <- array(NA, c(nrow(y), M, niter - burn))
  colnames(states.samples) <- dimnames(Smodel$T)[[1]]

  if( !is.null(X)){
    Z <- matrix(NA, niter - burn, dim(X)[2])
    BETA <- array(0, c(dim(X)[2], p, niter - burn))
  }

  ### MCMC

  for (i in 1:niter) {

    # Monitoring status
    if (i %in% sequence) {
      print(i)
    }

    ### STEP 1: Durbin & Koopman simulation smoother from KFAS package
    states <- simulateSSM(Smodel, type = "states")[, , 1]
    eta <- array(simulateSSM(Smodel, type = "eta"), c(t, p, k))
    eps <- simulateSSM(Smodel, type = "epsilon")[, , 1]

    ### STEP 2: Sampling each Sigma.k from its posterior, p(Sigma_k | eta) ~ IW (nu.k, s.k)

    # 2.1. posterior mean
    nu.k <- nu0.k + t

    # 2.2. posterior variance
    Sigma.k <- list()
    for (j in 1:k) {
      eta.k <- eta[, , j]
      s.k <- s0.k + crossprod(eta.k)
      Sigma.k[[j]] <- rInvWishart(1, nu.k, s.k)[, , 1]
    }

    # 2.3. the simulation smoother in Step 1 inherits the new posterior variance
    Sigma.k <- block.m(Sigma.k)
    Smodel$Q <- array(Sigma.k, c(K, K, 1))


    if(is.null(X)){

      ### STEP 3: Sampling Sigma.eps from its posterior, p(Sigma.eps | eps) ~ IW(nu.eps, s.eps )

      # 3.1. posterior mean
      nu.eps <- nu0.eps + t

      # 3.2. posterior variance
      s.eps <- s0.eps + crossprod(eps)
      Sigma.eps <- rInvWishart(1, nu.eps, s.eps)[,,1]

      # 3.3. the simulation smoother in Step 1 inherits the new posterior variance
      Smodel$H <- array(Sigma.eps, c(p, p, 1))

    } else {

      ### STEP 3: Sampling p(z)

      # 3.1. response variable with time series component subtracted out
      y.star <- y - tcrossprod(states, Smodel$Z[, , 1])

      # 3.2. data augmentation step, i.e. initializing a selection vector z
      #      with all elements equal 1 (meaning that at first all regressors are selected).
      #      Given this selection, the log-likelihood is computed.
      z <- rep.int(1, dim(X)[2])
      lpy.c <- lpy.X(y.star, X[, z == 1, drop = F], H = H, nu0.eps = nu0.eps, s0.eps = s0.eps)

      # 3.3. exploring the space of all possible models: at each step we change one element of z
      #      from 1 to 0 while the others z_-j are held fixed and compute the log-likelihood of the
      #      resulting model. Under the assumption that the prior probabilities are the same (see
      #      pdf) the difference between the two log-likelihood is the log odd, log(oj) and we can
      #      simulate p(z) from a bernulli distribution having \pi = 1/(1+oj^(-1)). Finally, if the
      #      resulting zj element has changed (z = zp) then lpy.p becomes the new likelihood and the
      #      cycle re-starts from the new z selection vectors.
      for (j in sample(1:dim(X)[2])) {
        zp <- z
        zp[j] <- 1 - zp[j]
        lpy.p <- lpy.X(y.star, X[, zp == 1, drop = F], H = H[zp == 1, zp == 1, drop = F], nu0.eps = nu0.eps,
                       s0.eps = s0.eps)
        r <- (lpy.p - lpy.c) * (-1)^(zp[j] == 0)
        z[j] <- rbinom(1, 1, 1/(1 + exp(-r)))
        if (z[j] == zp[j]) {
          lpy.c <- lpy.p
        }
      }

      ### STEP 4: Sampling Sigma.eps from its posterior, p(Sigma.eps | Yt, z)

      # 4.1. posterior mean and variance
      if (sum(z) == 0) {
        nu.eps <- nu0.eps + t
        s.eps <- s0.eps + crossprod(y.star)
      } else {
        W <- inv(crossprod(X[, z == 1, drop = F]) + inv(H[z == 1, z == 1, drop = F]))
        M <- tcrossprod(W, X[, z == 1, drop = F]) %*% y.star
        nu.eps <- nu0.eps + t
        s.eps <- s0.eps + crossprod(y.star) - crossprod(M, inv(W)) %*% M
      }


      # 4.2. the simulation smoother in Step 1 inherits the new posterior variance
      Sigma.eps <- rInvWishart(1, nu.eps, s.eps)[, , 1]
      Smodel$H <- array(Sigma.eps, c(p, p, 1))

      ### STEP 5: Sampling Beta from its posterior, p(Beta | Yt, Sigma, zeta)

      if (sum(z) == 0) {
        beta <- matrix(0, nrow = ncol(X), ncol = p)
      } else {
        beta <- rmatrixnorm(n = 1, mean = M, L = W, R = Sigma.eps, force = T)
      }
    }

    # Save results
    if (i > burn) {
      eta.samples[, , i - burn] <- matrix(eta, ncol = K)
      Sigma.states[, , i - burn] <- Sigma.k
      eps.samples[, , i - burn] <- eps
      Sigma.obs[, , i - burn] <- Sigma.eps
      states.samples[, , i - burn] <- states
      if(!is.null(X)){
        Z[i - burn, ] <- z
        BETA[z == 1, , i - burn] <- beta
      }
    }
  }

  if(!is.null(X)){
    list <- list(eta.samples = eta.samples, eps.samples = eps.samples, states.samples = states.samples,
                 Sigma.eta = Sigma.states, Sigma.eps = Sigma.obs, Z.beta = Z, beta = BETA, X = X, y = y, Z = Smodel$Z,
                 T = Smodel$T, R = Smodel$R, niter = niter, burn = burn)
  } else {
    list <- list(eta.samples = eta.samples, eps.samples = eps.samples, states.samples = states.samples,
                 Sigma.eta = Sigma.states, Sigma.eps = Sigma.obs, y = y, Z = Smodel$Z,
                 T = Smodel$T, R = Smodel$R, niter = niter, burn = burn)
  }

  class(list) <- "mbsts"
  return(list)
}


######################################################################################
## Marginal log-likelihood for a given model: p(Y | rho)
######################################################################################

lpy.X <- function(y, X, H, s0.eps, nu0.eps) {

  # get dimensions
  t <- dim(y)[1]
  p <- dim(y)[2]

  # log-likelihood
  if (dim(X)[2] == 0) {
    nu.eps <- nu0.eps + t
    s.eps <- s0.eps + crossprod(y)
    -(p * t/2) * log(pi) + (nu0.eps/2) * log(det(s0.eps)) + lmvgamma(nu.eps/2, p) - lmvgamma(nu0.eps/2,
                                                                                             p) - (nu.eps/2) * log(det(s.eps))
  } else {
    W <- inv(crossprod(X) + inv(H))
    M <- tcrossprod(W, X) %*% y
    nu.eps <- nu0.eps + t
    s.eps <- s0.eps + crossprod(y) - crossprod(M, solve(W, M))
    (p/2) * log(det(inv(H) %*% W)) - (p * t/2) * log(pi) + (nu0.eps/2) * log(det(s0.eps)) + lmvgamma(nu.eps/2,
                                                                                                     p) - lmvgamma(nu0.eps/2, p) - (nu.eps/2) * log(det(s.eps))
  }
}


######################################################################################
## Function for creating block matrices given a list of matrices
######################################################################################

block.m <- function(m) {
  matrix <- m[[1]]
  n <- length(m)
  ind <- if (n > 1) {
    2:n
  } else {
    NULL
  }
  for (i in ind) {
    matrix <- as.matrix(bdiag(matrix, m[[i]]))
  }
  matrix
}

######################################################################################
## Faster solve
######################################################################################

inv <- function(x) {
  chol2inv(chol(x))
}
