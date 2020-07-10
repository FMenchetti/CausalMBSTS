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
#' @param nu0.r Degrees of freedom of the Inverse-Wishart prior for each Sigma_r. Set by default to n0.r = d + 2 where d is the number of time series in the multivariate model.
#' @param s0.r Scale matrix of the Inverse-Wishart prior for each Sigma_r.
#' @param nu0.eps Degrees of freedom of the Inverse-Wishart prior for Sigma_eps. Set by default to d + 2.
#' @param s0.eps Scale matrix of the Inverse-Wishart prior for Sigma.eps.
#' @param niter Number of MCMC iteration.
#' @param burn Desired burn-in, set by default to 0.1 * niter.
#' @param ping A status message it's printed every 'ping' iteration, defaults to 0.1 * 'niter'.
#'
#' @return An object of class 'mbsts' which is a list with the following components
#' \describe{
#'   \item{eta.samples}{'niter' draws from the distribution of eta_r.}
#'   \item{eps.samples}{'niter' draws from the distribution of eps.}
#'   \item{states.samples}{draws from p(alpha_t | Y_{1:T}).}
#'   \item{Sigma.r}{'niter' draws from the posterior distribution of Sigma_r.}
#'   \item{sigma.eps}{'niter' draws from the posterior distribution of Sigma_eps.}
#'   \item{Z.beta}{('niter'- 'burn') x P matrix of the models selected at each iteration.}
#'   \item{beta}{ P x d x ('niter' - 'burn') ) array of the draws from the posterior distribution of the regression coefficient matrix.}
#'   \item{X}{Predictor matrix.}
#'   \item{y}{Matrix of observations.}
#'   \item{Z}{(1 x m) selection matrix of the observation equation.}
#'   \item{T}{(m x m) matrix of the state equation.}
#'   \item{R}{(1 x r) matrix selecting the state disturbances.}
#'   \item{niter}{Number of mcmc iterations.}
#'   \item{burn}{Burn-in.}
#' }
#' @export
#'
#' @examples
#' ## Example 1 : local level + seasonal (d = 3)
#' y <- cbind(seq(0.5,200,by=0.5)*0.1 + rnorm(400),
#'            seq(100.25,200,by=0.25)*0.05 + rnorm(400),
#'            rnorm(400, 5,1))
#' model.1 <- SSModel(y ~ SSMtrend(degree = 1, Q = matrix(NA)) + SSMseasonal(period=7, Q = matrix(NA)))
#' mcmc.1 <- mbsts.mcmc(model.1, s0.r = diag(3), s0.eps = diag(3), niter = 100, burn = 10)
#'
#' ## Example 2 : local level + seasonal + covariates (d = 2)
#' y <- cbind(rnorm(100), rnorm(100, 2, 3))
#' X <- cbind(rnorm(100, 0.5, 1) + 5, rnorm(100, 0.2, 2) - 2)
#' model.2 <- SSModel(y ~ SSMtrend(degree = 1, Q = matrix(NA,2,2)) + SSMseasonal(period=7, Q = matrix(NA,2,2)))
#' mcmc.2 <- mbsts.mcmc(model.2, X = X, s0.r = diag(2), s0.eps = diag(2), niter = 100, burn = 10)

mbsts.mcmc <- function(Smodel, X = NULL, H = NULL, nu0.r = NULL, s0.r, nu0.eps = NULL, s0.eps, niter,
    burn, ping = NULL) {

    # MCMC to sample from the joint posterior of model parameters in an mbsts model.
    # For now, the case without contemporaneous predictors variables is not allowed
    #
    # Args:
    #   Smodel  : a multivariate state space model of class 'SSModel'
    #   X       : an optional T x N matrix of predictors
    #   H       : N x N variance-covariance matrix between regression coefficients.
    #             The default is Zellner's g-prior, H = (X'X)^(-1)
    #   nu0.r   : degrees of freedom of the Inverse-Wishart prior for each Sigma_r.
    #             The default is the smallest integer such that the expectation of eta_r exists,
    #             that is, nu0.r = p + 2 where p is the number of time series in the
    #             multivariate model
    #   nu0.eps : degrees of freedom of the Inverse-Wishart prior for Sigma_eps, the default is p+2.
    #   s0.r    : Scale matrix of the Inverse-Wishart prior for each Sigma_r.
    #   s0.eps  : Scale matrix of the Inverse-Wishart prior for Sigma.eps
    #   niter   : number of MCMC iteration
    #   burn    : desired burn-in, set by default to 0.1 * niter
    #   ping    : logical, if TRUE a status message it's printed every iterations decile.
    #             Defaults to TRUE.
    #
    # Value:
    #   An object of class 'mbsts' which is a list with the following components
    #
    #   Smodel      : Resulting 'SSModel' (is it really useful? maybe I can save only the needed matrices, Z,R,T)
    #   eta.samples : 'niter' draws from the distribution of eta_r
    #   eps.samples : 'niter' draws from the distribution of eps
    #   states.samples : draws from p(\alpha_t | Y_{1:T})
    #   Sigma.r   : 'niter' draws from the posterior distribution of Sigma_r
    #   sigma.eps   : 'niter' draws from the posterior distribution of Sigma_eps
    #   Z           : 'niter' x N selection matrix
    #   beta        : P x d x 'niter' array containing the draws from the posterior distribution
    #                 of the regression coefficient matrix, p(beta | Sigma_eps, z)
    #   ...         : maybe not needed

    ### Dimensionalities & other inputs
    d <- dim(y)[2]  # number of time series

    # set default H (Zellner's g-prior)
    if (!is.null(X) & is.null(H)) {
        H <- inv(crossprod(X))
    }

    # set default ping
    if (is.null(ping)) {
        ping <- 0.1 * niter
    }
    sequence <- seq(ping, niter, by = ping)

    # set default nu0.r and nu0.eps
    if (is.null(nu0.r)) {
        nu0.r <- d + 2
    }
    if (is.null(nu0.eps)) {
        nu0.eps <- d + 2
    }

    # model definition
    Smodel$H[, , 1] <- rInvWishart(1, df = nu0.eps, Sigma = s0.eps)[, , 1]
    Smodel$Q[, , 1] <- stateVarianceDef(Smodel, nu = nu0.r , s = s0.r)

    # get dim
    y <- Smodel$y
    t <- dim(y)[1]  # number of time points
    rr <- dim(Smodel$R)[2]  # tot number of state disturbances
    r <- rr/d  # number of state disturbances for each time series
    M <- dim(Smodel$T)[1]  # tot number of states
    m <- M/d  # number of states for each time series

    ### Empty arrays to store MCMC iterations
    eta.samples <- array(NA, c(nrow(y), rr, niter - burn))
    eta.names <- etaNamesDef(Smodel) # that's really ugly...is it possible to find an alternative to degree.1 and degree.2?
    #eta.names <- unique(c(t(outer(attr(Smodel, "eta_types"), attr(Smodel$y, "dimnames")[[2]], paste))))
    colnames(eta.samples) <- eta.names

    Sigma.states <- array(NA, c(rr, rr, niter - burn))
    colnames(Sigma.states) <- attr(Smodel, "eta_types") # eta.names it's maybe better
    rownames(Sigma.states) <- attr(Smodel, "eta_types")

    eps.samples <- array(NA, c(nrow(y), d, niter - burn))
    Sigma.obs <- array(NA, c(d, d, niter - burn))

    states.samples <- array(NA, c(nrow(y), M, niter - burn))
    colnames(states.samples) <- dimnames(Smodel$T)[[1]]

    if (!is.null(X)) {
        Z <- matrix(NA, niter - burn, dim(X)[2])
        BETA <- array(0, c(dim(X)[2], d, niter - burn))
    }

    ### MCMC

    for (i in 1:niter) {

        # Monitoring status
        if (i %in% sequence) {
            print(i)
        }

        ### STEP 1: Durbin & Koopman simulation smoother from KFAS package
        set.seed(i)
        states <- simulateSSM(Smodel, type = "states")[, , 1]
        set.seed(i)
        eta <- simulateSSM(Smodel, type = "eta")
        set.seed(i)
        eps <- simulateSSM(Smodel, type = "epsilon")[, , 1]

        ### STEP 2: Sampling each Sigma.r from its posterior, p(Sigma_r | eta) ~ IW (nu.r, s.r)

        # 2.1. posterior mean
        nu.r <- nu0.r + t

        # 2.2. posterior variance
        s.r <- etaPosteriorScaleMatrix(eta, eta.names, Smodel, s0.r)
        Sigma.r <- stateVarianceDef(Smodel, nu.r, s.r)

        # 2.3. the simulation smoother in Step 1 inherits the new posterior variance
        Smodel$Q[,,1] <- Sigma.r


        if (is.null(X)) {

            ### STEP 3: Sampling Sigma.eps from its posterior, p(Sigma.eps | eps) ~ IW(nu.eps, s.eps )

            # 3.1. posterior mean
            nu.eps <- nu0.eps + t

            # 3.2. posterior variance
            s.eps <- s0.eps + crossprod(eps)
            Sigma.eps <- rInvWishart(1, nu.eps, s.eps)[, , 1]

            # 3.3. the simulation smoother in Step 1 inherits the new posterior variance
            Smodel$H <- array(Sigma.eps, c(d, d, 1))

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
            #      simulate p(z) from a Bernoulli distribution having \pi = 1/(1+oj^(-1)). Finally, if the
            #      resulting zj element has changed (z = zp) then lpy.p becomes the new likelihood and the
            #      cycle re-starts from the new z selection vectors.
            for (j in sample(1:dim(X)[2])) {
                zp <- z
                zp[j] <- 1 - zp[j]
                lpy.p <- lpy.X(y.star, X[, zp == 1, drop = F], H = H[zp == 1, zp == 1, drop = F],
                  nu0.eps = nu0.eps, s0.eps = s0.eps)
                re <- (lpy.p - lpy.c) * (-1)^(zp[j] == 0)
                z[j] <- rbinom(1, 1, 1/(1 + exp(-re)))
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
            Smodel$H <- array(Sigma.eps, c(d, d, 1))

            ### STEP 5: Sampling Beta from its posterior, p(Beta | Yt, Sigma, zeta)

            if (sum(z) == 0) {
                beta <- matrix(0, nrow = ncol(X), ncol = d)
            } else {
                beta <- rmatrixnorm(n = 1, mean = M, L = W, R = Sigma.eps, force = T)
            }
        }

        # Save results
        if (i > burn) {
            eta.samples[, , i - burn] <- eta
            Sigma.states[, , i - burn] <- Sigma.r
            eps.samples[, , i - burn] <- eps
            Sigma.obs[, , i - burn] <- Sigma.eps
            states.samples[, , i - burn] <- states
            if (!is.null(X)) {
                Z[i - burn, ] <- z
                BETA[z == 1, , i - burn] <- beta
            }
        }
    }

    if (!is.null(X)) {
        list <- list(eta.samples = eta.samples, eps.samples = eps.samples, states.samples = states.samples,
            Sigma.r = Sigma.states, Sigma.eps = Sigma.obs, Z.beta = Z, beta = BETA, X = X, y = y,
            Z = Smodel$Z, T = Smodel$T, R = Smodel$R, niter = niter, burn = burn)
    } else {
        list <- list(eta.samples = eta.samples, eps.samples = eps.samples, states.samples = states.samples,
            Sigma.r = Sigma.states, Sigma.eps = Sigma.obs, y = y, Z = Smodel$Z, T = Smodel$T, R = Smodel$R,
            niter = niter, burn = burn)
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
    d <- dim(y)[2]

    # log-likelihood
    if (dim(X)[2] == 0) {
        nu.eps <- nu0.eps + t
        s.eps <- s0.eps + crossprod(y)
        -(d * t/2) * log(pi) + (nu0.eps/2) * log(det(s0.eps)) + lmvgamma(nu.eps/2, d) - lmvgamma(nu0.eps/2,
            d) - (nu.eps/2) * log(det(s.eps))
    } else {
        W <- inv(crossprod(X) + inv(H))
        M <- tcrossprod(W, X) %*% y
        nu.eps <- nu0.eps + t
        s.eps <- s0.eps + crossprod(y) - crossprod(M, solve(W, M))
        (d/2) * log(det(inv(H) %*% W)) - (d * t/2) * log(pi) + (nu0.eps/2) * log(det(s0.eps)) + lmvgamma(nu.eps/2,
            d) - lmvgamma(nu0.eps/2, d) - (nu.eps/2) * log(det(s.eps))
    }
}

######################################################################################
## Function for the definition of the state variances
######################################################################################

stateVarianceDef<-function(Smodel, nu, s){
  d <- attr(Smodel,"p")
  comp <- attr(Smodel, "eta_types")

  if("level" %in% comp ){
    if( "slope" %in% comp ){ degree <- 2 } else { degree <- 1}

    Q.trend <- matrix(0, d*degree, d*degree)
    for(i in 1:degree){
      if(is.list(s)){
        Q <- rInvWishart(1 , nu, s[[paste("degree.",i,sep="")]])[,,1]
      } else {Q <- rInvWishart(1 , nu, s)[,,1]}
      Q.trend[seq(from = i, by = degree, length = d), seq(from = i, by = degree, length = d)] <- Q
    }
  } else { Q.trend <- NA}

  if("seasonal" %in% comp){
    if(is.list(s)){
      Q.seas <- rInvWishart(1, nu, s[["seasonal"]])[,,1]
    } else {Q.seas <- rInvWishart(1, nu, s)[,,1]}
  } else {Q.seas <- NA}

  if("cycle" %in% comp){
    Q.cycle<-matrix(0, 2*d, 2*d)
    for(i in 1:2){
      if(is.list(s)){
        Q<-rInvWishart(1, nu, s[[paste("cycle.",i,sep="")]])[,,1]
      } else {Q<-rInvWishart(1, nu, s)[,,1]}
      Q.cycle[seq(from=i, by=2, length=d), seq(from=i, by=2, length=d)]<-Q
    }
  } else {Q.cycle <- NA}

  Q<-block.m(list(Q.trend, Q.seas, Q.cycle))
  return(Q)
}

######################################################################################
## Function for generating the names of the state disturbances
######################################################################################

etaNamesDef<-function(Smodel){
  d <- attr(Smodel,"p")
  comp <- attr(Smodel, "eta_types")

  if("level" %in% comp ){
    if( "slope" %in% comp ){ degree <- 2 } else { degree <- 1}
    names_trend <-c()
    for(i in 1:degree){
      names_trend[i] <- paste("degree.",i, sep="")
    }
    names_trend <- paste0(names_trend, rep(attr(Smodel$y, "dimnames")[[2]], each = degree))
  } else {names_trend <- NULL}

  if("seasonal" %in% comp){
    names_seas <- "seasonal"
    names_seas <- paste0(names_seas, attr(Smodel$y, "dimnames")[[2]])
  } else {names_seas <- NULL}

  if("cycle" %in% comp){
    names_cycle <- c()
    for(i in 1:2){
      names_cycle[i] <- paste("cycle.",i,sep="")
    }
    names_cycle <- paste0(names_cycle, rep(attr(Smodel$y, "dimnames")[[2]], each = 2))
  } else {names_cycle <- NULL}

  eta_names <- c(names_trend, names_seas, names_cycle)
  return(eta_names)
}

######################################################################################
## Function for computing the posterior scale matrix of the state disturbances
######################################################################################

etaPosteriorScaleMatrix <- function(eta_sim, eta_names, Smodel, s0.r){

  # Step 1: eta_sim is returned from 'simulateSSM' as t x rr x 1;
  #         the following re-shapes eta_sim in an array t x d x r
  #         so to have a t x d matrix of disturbances for each state that has a disturbance term
  comp <- unique(attr(Smodel, "eta_type"))
  d<-attr(Smodel, "p")
  r <- length(unique(comp))
  comp<-gsub(comp, pattern = "slope", replacement = "degree.2")
  comp<-gsub(comp, pattern = "level", replacement = "degree.1")

  if("cycle" %in% comp){
    comp <- gsub(comp, pattern = "cycle", replacement = "cycle.1")
    comp <- c(comp, "cycle.2")
    r <- r+1
  }

  eta <- array(NA, c(nrow(eta_sim), d, r))

  for(i in 1:length(comp)){
    eta[,,i] <- eta_sim[, grep(eta_names, pattern = comp[i]),]
  }

  # Step 2 : computing the posterior scale matrix for each disturbance
  S.r <- list()
  for(i in 1:r){
    S.r[[i]] <- s0.r + crossprod(eta[,,i])
  }
  names(S.r) <- comp

  return(S.r)
}

######################################################################################
## Function for creating block matrices given a list of matrices
######################################################################################

block.m <- function(m) {
    m <- m[!is.na(m)]
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
