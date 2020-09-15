#' Plotting function for object of class CausalMBSTS
#'
#' Given an object of class 'CausalMBSTS', the function draws: i) the plot of
#' the estimated (pointwise) causal effect; ii) the original time series plotted
#' against the predicted counterfactual; iii) posterior predictive checks; iv)
#' regressor inclusion probabilities (only for models with a regression
#' component).
#'
#' @importFrom forecast Acf
#' @param x Object of class 'CausalMBSTS'
#' @param int.date Date of the intervention.
#' @param type A character string indicating the type of plot to be produced.
#'   Possible values in 'c('impact', 'forecast', 'ppchecks', 'inclusion.probs')'.
#'   See Details for further explanation.
#' @param prob Regressors inclusion probabilities above 'prob' are plotted.
#'   Optional, only required for type = 'inclusion.prob'.
#' @param ... Arguments passed to other methods (currently unused).
#'
#' @details
#' Option 'impact' for parameter \code{type} plots the general causal effect at every time points in the post
#'   period. Multiple plots will be generated, corresponding to each combination of time series and horizon (if specified).
#'   Option 'forecast' plots the observed time series against the predicted counterfactual, one plot per each
#'   combination of time series and horizon (if specified). 'ppchecks' draws posterior predictive checks for the model
#'   estimated in the pre-period. They include four plots generated for each time series (and horizon). The plots are
#'   (1) density of posterior mean vs. density of the data before intervention, (2) Histogram of maximum in-sample forecasts
#'   and Bayes p-value, (3) QQ-plot of residuals, and (4) ACF of residuals. Option 'inclusion.probs' plots the regressors'
#'   inclusion probabilities above 'prob'.
#'
#' @return NULL, invisibly.
#' @export
#'
#' @examples
#'
#' ## Example 1 (daily data, d = 3, local level + seasonal + covariates)
#' y <- cbind(seq(0.5,200,by=0.5)*0.1 + rnorm(400),
#'            seq(100.25,200,by=0.25)*0.05 + rnorm(400),
#'            seq(1,400,by=1)*(-0.01) + rnorm(400, 0, 0.5))
#' dates <- seq.Date(from = as.Date('2019-01-10'),by = "days", length.out = 400)
#'
#' # Adding a fictional intervention and four covariates. To illustrate the
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
#' # Model definition
#' causal.1 <- CausalMBSTS(y.new, components = c("trend", "seasonal"), seas.period = 7,
#'                         X = X, dates = dates, int.date = int.date,
#'                         s0.r = 0.1*diag(3), s0.eps = 0.1*diag(3), niter = 20,
#'                         burn = 5, horizon = as.Date(c('2019-12-05','2020-02-13')))
#'
#' ## Plotting
#' plot(causal.1, int.date = int.date, type = 'inclusion.probs', prob = 0.1)
#' # as expected, x4 is rarely included in the model
#' par(mar = c(2,2,2,2))
#' par(mfrow=c(2,3))
#' plot(causal.1, int.date = int.date, type = c('impact', 'forecast'))
#' par(mfrow=c(3,4))
#' plot(causal.1, type = 'ppchecks', int.date = int.date)
#'
#' ## Example 2
#' set.seed(1)
#' t <- seq(from = 0,to = 4*pi, length.out=300)
#' y <- cbind(3*sin(2*t)+rnorm(300), 2*cos(2*t) + rnorm(300))
#' dates <- seq.Date(from = as.Date("2015-01-01"), by = "week", length.out=300)
#' int.date <- as.Date("2020-02-27")
#' y[dates >= int.date,] <- y[dates >= int.date,]+2
#'
#' # Model definition
#' causal.2 <- CausalMBSTS(y, components = c("trend", "cycle"), cycle.period = 75,
#'                         dates = dates, int.date = int.date,
#'                         s0.r = 0.01*diag(2), s0.eps = 0.1*diag(2),
#'                         niter = 100, burn = 10)
#'
#' # Plotting
#' par(mfrow=c(2,4))
#' plot(causal.2, type = 'ppchecks', int.date = int.date)
#' par(mfrow=c(2,2))
#' plot(causal.2, type = c('impact','forecast'), int.date = int.date)


plot.CausalMBSTS <- function(x, int.date, type = c("impact", "forecast", "ppchecks"), prob = NULL, ...) {

    ## Parameter checks
    if(!all(type %in% c("impact", "forecast", "ppchecks", "inclusion.probs")))
        stop("allowed 'type' values are 'impact', 'forecast', 'ppchecks' and 'inclusion.probs'")

    ### Causal effect plot
    if ("impact" %in% type) {
        plotImpact(x, int.date = int.date)
    }

    ### Plot Observed vs Forecast
    if ("forecast" %in% type) {
        plotForecast(x, int.date = int.date)
    }

    ### Posterior predictive checks
    if ("ppchecks" %in% type) {
        plotChecks(x, int.date = int.date)
    }

    # Regressor index plot
    if ("inclusion.probs" %in% type) {
        plotInclusionProb(x, prob = prob)
    }

    return(invisible())
}

# ----------------------------------------------------------------------------------------

#' @import graphics
#'
plotImpact <- function(CausalMBSTS, int.date) {
    dates <- CausalMBSTS$dates
    dim <- dim(CausalMBSTS$mean.general)

    if (!is.null(dim)) {
        d <- dim[2]
        for (i in 1:d) {
            ylim <- c(min(CausalMBSTS$lower.general[, i]), max(CausalMBSTS$upper.general[, i]))
            x <- dates[dates >= int.date]
            main <- paste("Pointwise effect ", "Y", i, sep = "")
            plot(y = CausalMBSTS$mean.general[, i], x = x, type = "l", col = "blue", ylim = ylim,
                main = main, ylab = "", xlab = "")
            lines(y = CausalMBSTS$upper.general[, i], x = x, lty = 2)
            lines(y = CausalMBSTS$lower.general[, i], x = x, lty = 2)
        }
    } else {
        for (j in 1:length(CausalMBSTS$mean.general)) {
            dim <- dim(CausalMBSTS$mean.general[[j]])
            start <- which(dates == int.date)
            end <- start + dim[1]
            x <- dates[start:(end - 1)]
            d <- dim[2]
            for (i in 1:d) {
                ylim <- c(min(CausalMBSTS$lower.general[[j]][, i]), max(CausalMBSTS$upper.general[[j]][,i]))
                main <- paste("Pointwise effect ", "Y", i, ", horizon ", j, sep = "")
                plot(y = CausalMBSTS$mean.general[[j]][, i], x = x, type = "l", col = "blue", ylim = ylim,
                  main = main, ylab = "", xlab = "")
                lines(y = CausalMBSTS$upper.general[[j]][, i], x = x, lty = 2)
                lines(y = CausalMBSTS$lower.general[[j]][, i], x = x, lty = 2)
            }
        }
    }
}


#------------------------------------------------------------------------------------------------

#' @import graphics
#'
plotForecast <- function(CausalMBSTS, int.date) {
    dates <- CausalMBSTS$dates
    y <- CausalMBSTS$y
    start <- which(dates == int.date) - round(0.4 * sum(dates < int.date))
    post.mean <- apply(CausalMBSTS$predict$post.pred, c(1, 2), mean)

    if (!is.list(CausalMBSTS$mean.general)) {
        for (i in 1:dim(y)[2]) {
            end <- dim(y)[1]
            x <- dates[start:end]
            yi <- y[start:end, i]
            ylim <- c(min(yi, post.mean[start:end, i]), max(yi, post.mean[start:end, i]))
            main <- paste("Forecasted series ", "Y", i, sep = "")
            plot(y = yi, x = x, type = "l", ylim = ylim, ylab = "", xlab = "", main = main)
            lines(post.mean[start:end, i], col = "blue", x = x)
            abline(v = int.date, col = "red")
        }
    } else {
        for (j in 1:length(CausalMBSTS$mean.general)) {
            dim <- dim(CausalMBSTS$mean.general[[j]])
            end <- which(dates == int.date) + dim[1]
            x <- dates[start:(end - 1)]
            d <- dim[2]
            for (i in 1:d) {
                yi <- y[start:(end - 1), i]
                ylim <- c(min(yi, post.mean[start:(end - 1), i]), max(yi, post.mean[start:(end -
                  1), i]))
                main <- paste("Forecasted series ", "Y", i, ", horizon ", j, sep = "")
                plot(y = yi, x = x, type = "l", ylim = ylim, ylab = "", xlab = "", main = main)
                lines(post.mean[start:(end - 1), i], col = "blue", x = x)
                abline(v = int.date, col = "red")
            }
        }
    }
}

#--------------------------------------------------------------------------------------------

#' @import graphics stats
#'
plotChecks <- function(CausalMBSTS, int.date) {
    dates <- CausalMBSTS$dates
    ind <- dates < int.date
    y <- CausalMBSTS$y
    post.pred <- CausalMBSTS$predict$post.pred.0
    post.pred.mean <- apply(post.pred, c(1, 2), mean)

    for (i in 1:dim(y)[2]) {
        # Density of posterior mean vs density of the data before intervention
        plot(density(y[ind, i]), xlab = "", ylab = "",
             main = paste0("Density comparison, Y", i))
        lines(density(post.pred.mean[, i]), col = "blue")

        # Histograms & Bayesian p-value
        max.distrib <- apply(post.pred, c(2, 3), max)
        pvalue <- sum(max.distrib[i, ] >= max(y[ind, i]))/ncol(max.distrib)
        hist(max.distrib[i, ], 30, col = "lightblue", border = "grey", main = paste0("Bayesian p-value = ",
            round(pvalue, 2), ", Y", i), xlab = "Max. in-sample forecasts")
        abline(v = max(y[ind, i]), col = "darkblue", lwd = 3)

        # Residual plots
        y.rep <- matrix(y[ind, i], nrow(y[ind, ]), (CausalMBSTS$mcmc$niter - CausalMBSTS$mcmc$burn),
            byrow = F)
        res <- (y.rep - (post.pred[, i, ] - CausalMBSTS$mcmc$eps.samples[, i, ]))
        std.res <- t(apply(res, 1, FUN = "/", sqrt(CausalMBSTS$mcmc$Sigma.eps[i, i, ])))
        qqnorm(rowMeans(std.res), main = paste0("Residual QQ-plot, Y", i))
        qqline(rowMeans(std.res))
        Acf(rowMeans(std.res), main = paste0("Residual ACF, Y", i))
    }
}

#-----------------------------------------------------------------------------------------------

#' @import graphics
#'
plotInclusionProb <- function(CausalMBSTS, prob = .5) {
    if (is.null(prob))
        {
            prob <- 0.5
        }  # check here: if no reg above that prob threshold is present stop and write a message (e.g. no reg above that threshold, try to lower the threshold)
    X <- CausalMBSTS$mcmc$X
    select <- apply(CausalMBSTS$mcmc$Z.beta, 2, mean) > prob  # selecting regressors' with prob of inclusion bigger than some amount

    plot(apply(CausalMBSTS$mcmc$Z.beta, 2, mean)[select], type = "h", lwd = 2,
         ylab = expression(paste("Pr(", italic(z[j] == 1), "|", italic(y), ",X)", sep = "")),
         las = 2, xaxt = "none", xlab = "", main = "Regressors' inclusion probabilities")
    axis(1, seq(1, ncol(X[, select]), 1), las = 2, labels = colnames(X)[select])
}
