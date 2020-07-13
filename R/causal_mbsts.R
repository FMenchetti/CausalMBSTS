######################################################################################
######################################################################################
####  Author:           Fiammetta Menchetti                                       ####
####                                                                              ####
####  Date last update: 2020-06-23                                                ####
####                                                                              ####
####  Content:          Joint causal effect estimation for MBSTS models           ####
####                                                                              ####
####                                                                              ####
####  Main function :   causal.mbsts                                              ####
####  Dependencies:     predict.mbsts                                             ####
####                    mbsts.mcmc (lpy.X,block.m)                                ####
####                                                                              ####
######################################################################################
######################################################################################


#' Causal effect estimation in a multivariate setting
#'
#' It estimates the general effect of an intervention in a multivariate time series
#' setting. It uses MCMC to sample from the joint posterior distribution of the parameters
#' of an MBSTS model before the intervention/treatment occurred. Then, it uses the
#' covariates post intervention to predict the counterfactual potential outcomes.
#' The prediction is done by sampling from the posterior predictive
#' distribution (ppd). Then the causal effect is computed by taking the difference between
#' the observed outcome of each time series and the mean of the ppd (credible intervals are
#' computed accordingly).
#'
#' @param Smodel A multivariate state space model of class 'SSModel'.
#' @param X Optional t x N matrix of predictors.
#' @param y t x d matrix of observations.
#' @param dates a vector of dates.
#' @param int.date Date of the intervention.
#' @param holi Optional vector of the same length as the time series, specifying the dates (if any) that should be excluded from the computation of the causal effect in the post period. The elements of the vector must be either 0 (the corresponding date is retained) or 1 (the corresponding date is excluded).
#' @param horizon Optional, vector of dates. If provided, a causal effect is computed for any time horizon. It defaults to the date of the last observation.
#' @param H P x P variance-covariance matrix of the regression coefficients. Set by default to H = (X'X)^(-1).
#' @param nu0.r Degrees of freedom of the Inverse-Wishart prior for each Sigma_r. Set by default to n0.r = d + 2 where d is the number of time series in the multivariate model.
#' @param s0.r Scale matrix of the Inverse-Wishart prior for each Sigma_r.
#' @param nu0.eps Degrees of freedom of the Inverse-Wishart prior for Sigma_eps. Set by default to d + 2.
#' @param s0.eps Scale matrix of the Inverse-Wishart prior for Sigma.eps.
#' @param niter Number of MCMC iteration.
#' @param burn Desidered burn-in, set by default to 0.1 * niter.
#' @param ping A status message it's printed every 'ping' iteration, defaults to 0.1 * 'niter'.
#'
#' @return A list with the following components
#'   \item{mcmc}{An object of class 'mbsts'.}
#'   \item{predict}{A list with the same components as those produced by the function 'predict.mbsts'}
#'   \item{y}{Observations in the analysis period excluding 'holi' (if provided).}
#'   \item{dates}{Dates in the analysis period excluding 'holi' (if provided).}
#'   \item{general}{General causal effect for all iterations.}
#'   \item{general.effect}{the estimated average causal effect and a 95\% credible interval.}
#'   \item{mean.general}{pointwise effect (may be not necessary).}
#'   \item{lower.general}{lower bound of the pointwise effect (may be not necessary).}
#'   \item{upper.general}{Upper bound of the pointwise effect (may be not necessary).}
#' @export
#'
#' @examples
#' ## Example 1 (daily data, d = 3, local level + seasonal + covariates)
#' # Generating a panel of observations and a vector of dates
#' set.seed(1)
#' y <- cbind(seq(0.5,200,by=0.5)*0.1 + rnorm(400),
#'            seq(100.25,200,by=0.25)*0.05 + rnorm(400),
#'            seq(1,400,by=1)*(-0.01) + rnorm(400, 0, 0.5))
#' dates <- seq.Date(from = as.Date('2019-01-10'),by = "days", length.out = 400)
#'
#' # Adding a fictional intervention and four covariates (they should be related to the outcome but unaffected by the intervention). To illustrate the functioning of Bayesian model selection, one covariate is assumed to be unrelated to y.
#' int.date <- as.Date('2019-11-05')
#' y.new <- y; y.new[dates >= int.date, ] <- y.new[dates >= int.date, ]*1.3
#' x1 <- y[,1]*0.5 + y[,2]*0.3 + y[,3]*0.1
#' x2 <- y[,2]*0.1 + rnorm(dim(y)[1],0,0.5)
#' x3 <- y[,3]*1.2 + rnorm(dim(y)[1],0,0.5)
#' x4 <- rnorm(dim(y)[1], 5, 10)
#' X <- cbind(x1, x2, x3, x4)
#'
#' # Some plots
#' par(mfrow=c(1,3))
#' for(i in 1:dim(y.new)[2]){
#'   plot(y.new[,i], x = dates, type='l', col='cadetblue', xlab='', ylab='', main= bquote(Y[.(i)]))
#'   lines(y[,i], x = dates, col='orange')
#'   }
#' par(mfrow=c(1,4))
#' for(i in 1:dim(X)[2]){
#'   plot(X[,i], type='l', col = 'darkgreen', x = dates, xlab='', ylab='', main = bquote(x[.(i)]))
#'   }
#'
#' # Model definition
#' model.1 <- SSModel(y ~ SSMtrend(degree = 1, Q = matrix(NA)) + SSMseasonal(period=7, Q = matrix(NA)))
#' causal.1 <- causal.mbsts(model.1, X = X, y = y.new, dates = dates, int.date = int.date, s0.r = 0.1*diag(3), s0.eps = 0.1*diag(3), niter = 100, burn = 10, horizon = c('2019-12-05','2020-02-13'))
#' causal.1$general.effect
#'
#' ## Example 2 (weekly data, local level + cycle, d = 2)
#' set.seed(1)
#' t <- seq(from = 0,to = 4*pi, length.out=300)
#' y <- cbind(3*sin(2*t)+rnorm(300), 2*cos(2*t) + rnorm(300))
#' dates <- seq.Date(from = as.Date("2015-01-01"), by = "week", length.out=300)
#' int.date <- as.Date("2020-02-27")
#' y[dates >= int.date,] <- y[dates >= int.date,]+2
#'
#' # Some plots
#' plot(y = y[,1], x=dates, type="l", col="cadetblue")
#' lines(y = y[,2], x = dates, col = "orange")
#' abline(v=int.date, col="red")
#'
#' # Model definition
#' model.2 <- SSModel(y ~ SSMtrend(degree = 1, Q = matrix(NA)) + SSMcycle(period=75, Q = matrix(NA)))
#' causal.2 <- causal.mbsts(model.2, y = y, dates = dates, int.date = int.date, s0.r = 0.01*diag(2), s0.eps = 0.1*diag(2), niter = 100, burn = 10)
#' causal.2$general.effect



causal.mbsts <- function(Smodel, X = NULL, y, dates, int.date, holi = NULL, horizon = NULL, H = NULL,
    nu0.r = NULL, s0.r, nu0.eps = NULL, s0.eps, niter, burn = NULL, ping = NULL) {

    # It estimates the general effect of an intervention in a multivariate time series
    # setting. It uses MCMC to sample from the joint posterior distribution of the parameters
    # of an MBSTS model before the intervention/treatment occurred. Then, it uses the
    # covariates post intervention to predict the counterfactual potential outcomes.
    # The prediction is done by sampling from the posterior predictive
    # distribution (ppd). Then the causal effect is computed by taking the difference between
    # the observed outcome of each time series and the mean of the ppd (credible intervals are
    # computed accordingly).

    # Args:
    #   Smodel   : a multivariate state space model of class 'SSModel'
    #   X        : a t x P matrix of predictors
    #   y        : time series of observations
    #   dates    : a vector of dates
    #   int.date : date of the intervention
    #   holi     : Optional, 0-1 vector of the same length as the time series specifying whether 'dates[i]' should be included ('holi[i]' = 0) or excluded ('holi[i]' = 1)
    #   horizon  : Optional, vector of dates. If provided, a causal effect is computed for any time horizon. It defaults to the date of the last observation.
    #   H        : desidered P x P variance-covariance matrix between regression coefficients.
    #             The default is Zellner's g-prior, H = (X'X)^(-1)
    #   nu0.r    : degrees of freedom of the Inverse-Wishart prior for each Sigma_r.
    #             The default is the smallest integer such that the expectation of eta_r exists,
    #             that is, nu0.r = d + 2 where d is the number of time series in the
    #             multivariate model
    #   nu0.eps  : degrees of freedom of the Inverse-Wishart prior for Sigma_eps, the default is d+2.
    #   s0.r     : Scale matrix of the Inverse-Wishart prior for each Sigma_r.
    #   s0.eps   : Scale matrix of the Inverse-Wishart prior for Sigma.eps
    #   niter    : number of MCMC iteration
    #   burn     : desidered burn-in, set by default to 0.1 * niter
    #   ping     : logical, if TRUE a status message it's printed every iterations decile.
    #             Defaults to TRUE.
    #
    # Value: UPDATE
    #   joint.effect : the estimated average causal effect and a 95% credible interval
    #   joint        : joint causal effect for all iterations
    #   mean.joint   : pointwise effect
    #   lower.joint  : lower bound of the pointwise effect
    #   upper.joint  : upper bound of the pointwise effect
    #   mcmc         : an object of class 'mbsts'
    #   predict      : a list with the same components as those produced by the function 'predict.mbsts'
    #   adj.series   : observations in the analysis period excluding holidays
    #   adj.dates    : dates in the analysis period without holidays
    #   original.series : original time series (maybe not needed)
    #   original.dates : original dates (maybe not needed)

    ### STEP 1. Dividing pre and post periods

    ind <- dates < int.date
    X.pre <- X[ind, ]
    X.post <- X[!ind, ]
    if(is.null(dimnames(y)[[2]])) dimnames(y)[[2]] <- dimnames(Smodel$y)[[2]]
    y.pre <- y[ind, ]
    y.post <- y[!ind, ]

    # Estimating the model only in the pre-period
    Smodel$y <- y.pre
    attr(Smodel, "n") <- as.integer(nrow(y.pre))

    ### STEP 2. MCMC

    mbsts <- mbsts.mcmc(Smodel = Smodel, X = X.pre, H = H, nu0.r = nu0.r, s0.r = s0.r, nu0.eps = nu0.eps,
        s0.eps = s0.eps, niter = niter, burn = burn, ping = ping)

    ### STEP 3. In- and out-of-sample forecasts from the ppd
    predict <- predict(mbsts, steps.ahead = dim(y[!ind,])[1], X.post)

    ### STEP 4. Causal effect estimation
    # d <- dim(mbsts$y)[2]
    burn <- mbsts$burn
    y.post.rep <- array(y.post, c(nrow(y.post), ncol(y.post), niter - burn))
    y.diff <- y.post.rep - predict$post.pred.1

    # removing given dates
    holi <- holi[dates >= int.date]
    if (!is.null(holi)) {
        y.diff <- y.diff[holi == 0, , ]
        dates <- dates[holi == 0]
        y <- y[holi == 0, ]
    }

    # General causal effect (temporal average and cumulative sum)

    if (length(horizon) > 0) {
        mean.effect <- list()
        lower.bound <- list()
        upper.bound <- list()
        joint.effect <- list()
        for (i in 1:length(horizon)) {
            ind <- dates[dates >= int.date] <= horizon[i]
            mean.effect[[i]] <- apply(y.diff[ind, , ], c(1, 2), mean)
            lower.bound[[i]] <- apply(y.diff[ind, , ], c(1, 2), quantile, probs = 0.025)
            upper.bound[[i]] <- apply(y.diff[ind, , ], c(1, 2), quantile, probs = 0.975)
            joint.effect[[i]] <- cbind(mean = apply(colMeans(y.diff[ind, , ]), 1, mean),
                                       lower = apply(colMeans(y.diff[ind, , ]), 1, quantile, probs = 0.025),
                                       upper = apply(colMeans(y.diff[ind, , ]), 1, quantile, probs = 0.975),
                                       cum.sum = apply(colSums(y.diff[ind, , ]), 1, mean),
                                       cum.lower = apply(colSums(y.diff[ind, , ]), 1, quantile, probs = 0.025),
                                       cum.upper = apply(colSums(y.diff[ind, , ]), 1, quantile, probs = 0.975))
        }
    } else {
        mean.effect <- apply(y.diff, c(1, 2), mean)
        lower.bound <- apply(y.diff, c(1, 2), quantile, probs = 0.025)
        upper.bound <- apply(y.diff, c(1, 2), quantile, probs = 0.975)
        joint.effect <- cbind(mean = apply(colMeans(y.diff), 1, mean),
                              lower = apply(colMeans(y.diff), 1, quantile, probs = 0.025),
                              upper = apply(colMeans(y.diff), 1, quantile, probs = 0.975),
                              cum.sum = apply(colSums(y.diff), 1, mean),
                              cum.lower = apply(colSums(y.diff), 1, quantile, probs = 0.025),
                              cum.upper = apply(colSums(y.diff), 1, quantile, probs = 0.975))
    }

    list <- list(mcmc = mbsts, predict = predict, y = y, dates = dates, general = y.diff, general.effect = joint.effect,
        mean.general = mean.effect, lower.general = lower.bound, upper.general = upper.bound)
    class(list) <- "CausalMBSTS"
    return(list)
}

#--------------------------------------------------------------------------------------
