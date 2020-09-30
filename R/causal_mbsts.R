######################################################################################
######################################################################################
####  Author:           Fiammetta Menchetti                                       ####
####                                                                              ####
####  Date last update: 2020-06-23                                                ####
####                                                                              ####
####  Content:          Joint causal effect estimation for MBSTS models           ####
####                                                                              ####
####                                                                              ####
####  Main function :   CausalMBSTS                                               ####
####  Dependencies:     predict.mbsts                                             ####
####                    as.mcmc                                                   ####
####                                                                              ####
######################################################################################
######################################################################################


#' Causal effect estimation in a multivariate setting
#'
#' The function estimates the general effect of an intervention in a multivariate time series
#' setting. It uses MCMC to sample from the joint posterior distribution of the parameters
#' of an MBSTS model before the intervention/treatment occurred. Then, it uses the
#' post-intervention covariate values to predict the counterfactual potential outcomes.
#' The prediction is done by sampling from the posterior predictive
#' distribution (PPD). Then the causal effect is computed by taking the difference between
#' the observed outcome of each time series and the mean of the PPD (credible intervals are
#' computed accordingly).
#'
#' @param y t x d data.frame (or matrix) of observations, where d is the number
#'   of time series in the multivariate model.
#' @param components Character vector specifying the components of the
#'   multivariate structural time series model. Possible values are c("trend",
#'   "slope", "seasonal", "cycle").
#' @param seas.period Length of the seasonal pattern, if present.
#' @param cycle.period Length of the cycle pattern, if present.
#' @param X Optional t x N data frame (or matrix) of N predictors.
#' @param dates a vector of dates of length t (with elements of class
#'   \code{Date}) that correspond to observations in y.
#' @param int.date Date of the intervention (must be of class \code{Date}).
#' @param alpha Level of credible interval to report for the estimated causal
#'   effect. Default set to 0.05 (i.e., reporting a two-sided 95\% credible
#'   interval).
#' @param excl.dates Optional vector of length t, specifying the dates (if any)
#'   in the post period that should be excluded from the computation of the
#'   causal effect. The elements of the vector must be either 0 (the
#'   corresponding date is retained) or 1 (the corresponding date is excluded).
#'   The first part that corresponds to \code{dates < int.date} is ignored.
#' @param horizon Optional, vector of dates (with elements of class
#'   \code{Date}). If provided, a causal effect is computed for the time
#'   horizon(s) between \code{int.date} and each specified date. Defaults to the
#'   date of the last observation.
#' @param H P x P variance-covariance matrix of the regression coefficients. Set
#'   by default to H = c(X'X)^(-1) which is akin to the Zellner's g-prior. The
#'   value of the scaling factor is set to \code{c = 1}. Alternative priors
#'   could be H = c*diag((X'X)^(-1)) or H = c*I. See also Smith & Kohn, 1995
#'   that suggest setting \code{c} in the range [10,1000].
#' @param nu0.r Degrees of freedom of the Inverse-Wishart prior for each element
#'   of Sigma.r, a vector of errors for state r. Set by default to d + 2 (must
#'   be greater than d - 1).
#' @param s0.r Scale matrix of the Inverse-Wishart prior for each Sigma.r, a
#'   vector of errors for state r. Must be a (d x d) positive definite. Default
#'   set to the variance-covariance matrix of y multiplied by a scaling factor
#'   of 0.01.
#' @param nu0.eps Degrees of freedom of the Inverse-Wishart prior for Sigma.eps,
#'   a vector of observation errors for each time series. Set by default to d +
#'   2 (must be greater than d - 1).
#' @param s0.eps Scale matrix of the Inverse-Wishart prior for Sigma.eps, a
#'   vector of observation errors for each time series. Must be a (d x d)
#'   positive definite. Default set to the variance-covariance matrix of y
#'   multiplied by a scaling factor of 0.01.
#' @param niter Number of MCMC iterations.
#' @param burn Desired burn-in, set by default to 0.1 * \code{niter}.
#' @param ping A status message is printed every \code{ping} iteration. Default
#'   set to 0.1 * \code{niter}. Set to 0 to not track the status.
#'
#' @details {The assumed model is based on Normally distributed disturbance terms. The argument
#' \code{components} provides flexibility for model formulation, allowing to add
#' simultaneously up to four components that encapsulate the characteristics of a time series.
#'
#' The unknown parameters are the variance-covariance matrices of the error terms and,
#' if covariates are provided, the matrix of regression coefficients. Because of conjugacy, the
#' priors placed on the variance-covariance matrices of the error terms are Inverse-Wishart
#' distributions and the arguments (nu0.eps, s0.eps) and (nu0.r, s0.r) regulate their hyperparameters.
#'
#' The regression coeffiecients are assumed to follow a matrix-Normal prior and, to incorporate
#' a sparsity assumption, the prior mean is set to zero and a vector selecting the relevant
#' covariates is introduced with a data augmentation step.
#'
#' Sampling from the joint posterior distribution of the states and model parameters is performed
#' with a Gibbs sampler. The estimated model is then used to perform predictions of the counterfactual
#' potential outcomes in the period following the intervention date. In a final step, the predictions
#' are compared to the observed outcomes, thereby defining a causal effect at each time point (pointwise effect).
#'
#' The output component \code{general.effect} summarizes the estimates of two causal effects: the average
#' causal effect (temporal average of the pointwise effect) and the cumulative causal effect (cumulative
#' sum of the pointwise effect).
#'
#' Run \code{vignette("CausalMBSTS")} for a detailed example.
#'
#' For further details see Menchetti & Bojinov (2020).}
#'
#' @return A list with the following components:
#'   \item{mcmc}{An object of class \code{mbsts}.}
#'   \item{predict}{A list with the same components as those produced by the function \code{\link{predict.mbsts}}}
#'   \item{y}{Observations in the analysis period (excluding \code{excl.dates} if provided).}
#'   \item{dates}{Dates in the analysis period (excluding \code{excl.dates} if provided).}
#'   \item{general}{General causal effect for all iterations.}
#'   \item{general.effect}{Estimated average causal effect, cumulative causal effect and (1-\code{alpha})\%
#'   credible intervals. Returns a list if \code{horizon} is specified.}
#'   \item{mean.general}{Pointwise effect. Returns a list if \code{horizon} is specified.}
#'   \item{lower.general}{Lower bound of a (1-\code{alpha})\% credible interval of the pointwise effect.
#'   Returns a list if \code{horizon} is specified.}
#'   \item{upper.general}{Upper bound of a (1-\code{alpha})\% credible interval of the pointwise effect.
#'   Returns a list if \code{horizon} is specified.}
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
#' # Adding a fictional intervention and four covariates (they should be related
#' # to the outcome but unaffected by the intervention). To illustrate the
#' # functioning of Bayesian model selection, one covariate is assumed to be
#' # unrelated to y.
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
#' # Causal effect estimation
#' causal.1 <- CausalMBSTS(y.new, components = c("trend", "seasonal"), seas.period = 7, X = X,
#'                         dates = dates, int.date = int.date, s0.r = 0.1*diag(3), s0.eps = 0.1*diag(3),
#'                         niter = 100, burn = 10, horizon = as.Date(c('2019-12-05','2020-02-13')))
#' summary(causal.1)
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
#' # Causal effect estimation
#' causal.2 <- CausalMBSTS(y, components = c("trend", "cycle"), cycle.period = 75,
#'                         dates = dates, int.date = int.date, s0.r = 0.01*diag(2),
#'                         s0.eps = 0.1*diag(2), niter = 100, burn = 10)
#' summary(causal.2)



CausalMBSTS <- function(y, components, seas.period = NULL, cycle.period = NULL,
                        X = NULL, dates, int.date, alpha = 0.05, excl.dates = NULL,
                        horizon = NULL, H = NULL, nu0.r = NULL, s0.r, nu0.eps = NULL,
                        s0.eps, niter, burn = NULL, ping = NULL) {

    ### Parameter checks
    if(!is.matrix(y) && !is.data.frame(y)) stop("`y` must be a matrix or data.frame")
    if(!all(components %in% c("trend", "slope", "seasonal", "cycle")))
        stop("`components` must be one of 'trend', 'slop', 'seasonal', or 'cycle'")
    if(!missing(seas.period) && (!is.numeric(seas.period) ||length(seas.period) != 1))
        stop("`seas.period` must be a numeric vector of length one")
    if(!missing(cycle.period) && (!is.numeric(cycle.period) ||length(cycle.period) != 1))
        stop("`cycle.period` must be a numeric vector of length one")
    if(!missing(X)) {
        if(!is.matrix(X) && !is.data.frame(X)) stop("`X` must be a matrix or data.frame")
        if(nrow(X) != nrow (y)) stop("nrow(X) != nrow(X)")
    }
    if(!any(class(dates) %in% c("Date", "POSIXct", "POSIXlt", "POSIXt")))
        stop("`dates` must be a vector of class Date")
    if(length(dates) != nrow(y)) stop("length(dates) != nrow(y)")
    if(length(int.date) != 1 || !any(class(dates) %in% c("Date", "POSIXct", "POSIXlt", "POSIXt")))
        stop("`int.date` must be a Date of length 1")
    if(!missing(excl.dates)) {
        if(any(!as.integer(excl.dates) %in% c(0L, 1L))) stop("`excl.dates` must be 0/1 or TRUE/FALSE")
        if(length(excl.dates) != nrow(y)) stop("length(`excl.dates`) != nrow(`y`)")
    }
    if(!missing(horizon) && !any(class(horizon) %in% c("Date", "POSIXct", "POSIXlt", "POSIXt")))
        stop("`horizon` must be a Date object")
    if(missing(s0.r)){s0.r <- NULL}
    if(missing(s0.eps)){s0.eps <- NULL}

    ### STEP 1. Dividing pre and post periods
    ind <- dates < int.date
    if(all(!ind)) stop("All 'dates' are prior to 'int.date'.")
    X.pre <- X[ind, ]
    X.post <- X[!ind, ]
    if(is.null(dimnames(y)[[2]])) dimnames(y)[[2]] <- paste("y",seq(1,dim(y)[2],by=1), sep="")
    y.pre <- y[ind, ]
    y.post <- y[!ind, ]

### STEP 2. MCMC
    mbsts_args <- list(y = y.pre, components = components, seas.period = seas.period,
                       cycle.period = cycle.period, X = X.pre, H = H, nu0.r = nu0.r,
                       s0.r = s0.r, nu0.eps = nu0.eps, s0.eps = s0.eps, niter = niter,
                       burn = burn, ping = ping)

    mbsts_args <- mbsts_args[sapply(mbsts_args, function(x) !is.null(x))] # !is.na(x) &&
    mbsts <- do.call(as.mbsts, mbsts_args)

    ### STEP 3. In- and out-of-sample forecasts from the PPD
    predict <- predict(mbsts, steps.ahead = dim(y[!ind,])[1], X.post)

    ### STEP 4. Causal effect estimation
    burn <- mbsts$burn
    y.post.rep <- array(y.post, c(nrow(y.post), ncol(y.post), niter - burn))
    y.diff <- y.post.rep - predict$post.pred.1

    # removing given dates
    if (!is.null(excl.dates)) {
        excl.dates.post <- excl.dates[dates >= int.date]
        y.diff <- y.diff[excl.dates.post == 0, , ]
        dates <- dates[excl.dates == 0]
        y <- y[excl.dates == 0, ]
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
            lower.bound[[i]] <- apply(y.diff[ind, , ], c(1, 2), quantile, probs = alpha/2)
            upper.bound[[i]] <- apply(y.diff[ind, , ], c(1, 2), quantile, probs = 1-alpha/2)
            joint.effect[[i]] <- cbind(mean = apply(colMeans(y.diff[ind, , ]), 1, mean),
                                       lower = apply(colMeans(y.diff[ind, , ]), 1, quantile, probs = alpha/2),
                                       upper = apply(colMeans(y.diff[ind, , ]), 1, quantile, probs = 1-alpha/2),
                                       cum.sum = apply(colSums(y.diff[ind, , ]), 1, mean),
                                       cum.lower = apply(colSums(y.diff[ind, , ]), 1, quantile, probs = alpha/2),
                                       cum.upper = apply(colSums(y.diff[ind, , ]), 1, quantile, probs = 1-alpha/2))
        }
    } else {
        mean.effect <- apply(y.diff, c(1, 2), mean)
        lower.bound <- apply(y.diff, c(1, 2), quantile, probs = alpha/2)
        upper.bound <- apply(y.diff, c(1, 2), quantile, probs = 1-alpha/2)
        joint.effect <- cbind(mean = apply(colMeans(y.diff), 1, mean),
                              lower = apply(colMeans(y.diff), 1, quantile, probs = alpha/2),
                              upper = apply(colMeans(y.diff), 1, quantile, probs = 1-alpha/2),
                              cum.sum = apply(colSums(y.diff), 1, mean),
                              cum.lower = apply(colSums(y.diff), 1, quantile, probs = alpha/2),
                              cum.upper = apply(colSums(y.diff), 1, quantile, probs = 1-alpha/2))
    }

    list_res <- list(mcmc = mbsts, predict = predict, y = y, dates = dates, general = y.diff, general.effect = joint.effect,
        mean.general = mean.effect, lower.general = lower.bound, upper.general = upper.bound)
    class(list_res) <- "CausalMBSTS"
    return(list_res)
}

#--------------------------------------------------------------------------------------
