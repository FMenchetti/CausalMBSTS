
#' Summarize results of causal effect estimation with CausalMBSTS
#'
#' The method extracts and computes various results of the causal analysis with CausalMBSTS.
#'
#' @param x An object of class 'CausalMBSTS', a result of a call to \code{\link{CausalMBSTS}}.
#'
#' @return Returns an object of class \code{summary.CausalMBSTS}, which is a list of data frames corresponding to each
#' date provided in \code{horizon} (or its default value) with the following columns:
#'   \item{mean}{Estimated average causal effect.}
#'   \item{lower}{Lower bound of the two-sided (1-\code{alpha})\% credible interval.}
#'   \item{upper}{Upper bound of the two-sided (1-\code{alpha})\% credible interval.}
#'   \item{cum.sum}{Pointwise effect}
#'   \item{lower.general}{Lower bound of a (1-\code{alpha})\% credible interval of the pointwise effect.}
#'   \item{upper.general}{Upper bound of a (1-\code{alpha})\% credible interval of the pointwise effect.}
#'   \item{bayes.pval}{Bayesian p-value for the average causal effect.}
#'   \item{pct.causal.eff}{Probability of a causal effect (\%)}
#'
#' @examples
#' set.seed(1)
#' t <- seq(from = 0,to = 4*pi, length.out=300)
#' y <- cbind(3*sin(2*t)+rnorm(300), 2*cos(2*t) + rnorm(300))
#' dates <- seq.Date(from = as.Date("2015-01-01"), by = "week", length.out=300)
#' int.date <- as.Date("2020-02-27")
#' y[dates >= int.date,] <- y[dates >= int.date,]+2
#'
#' # Causal effect estimation
#' causal.2 <- CausalMBSTS(y, components = c("trend", "cycle"), cycle.period = 75,
#'                         dates = dates, int.date = int.date, s0.r = 0.01*diag(2),
#'                         s0.eps = 0.1*diag(2), niter = 100, burn = 10)
#'
#' sum.causal.2 <- summary(causal.2)
#' print(sum.causal.2, digits = 2)
#' sum.causal.2$horizon_default
#'
#' @export
summary.CausalMBSTS <- function(x){

  effect <- x$general.effect
  pred <- x$predict$post.pred
  niter <- dim(pred)[3]
  start <- dim(x$predict$post.pred.0)[1] +1 # intervention date
  summary_impact <- list()

  if (is.list(effect)){
    for(i in 1:length(effect)){

      # Bayesian p-value
      end <- start + dim(x$mean.general[[i]])[1] -1 # interval between the intervention and the time horizon
      true.sum <- colSums(x$y[start:end,])
      pred.sum <- apply(pred[start:end,,], c(2,3), sum) # cumulative causal effect at each iteration
      replicas <- cbind(pred.sum,true.sum) # cumulative sum of the observations in each replicated data set
      p <- apply(rbind(rowSums(replicas >= true.sum)/(niter + 1), rowSums(replicas <= true.sum)/(niter+1)),2,min)
      prob_causal_eff <- data.frame(bayes.pval = p, pct.causal.eff = (1-p)*100)

      # Outputting estimated effect(s) and p-values(s)
      output <- cbind(effect[[i]], prob_causal_eff)
      rownames(output) <-dimnames(x$y)[[2]]

      # saving output
      summary_impact[[i]] <- output
    }

    names(summary_impact) <- paste0("horizon",1:length(effect))

  } else {

    # Bayesian p-value
    end <- dim(x$y)[1]
    true.sum <- colSums(x$y[start:end,])
    pred.sum <- apply(pred[start:end,,], c(2,3), sum) # cumulative causal effect at each iteration
    replicas <- cbind(pred.sum,true.sum) # cumulative sum of the observations in each replicated data set
    p <- apply(rbind(rowSums(replicas >= true.sum)/(niter + 1), rowSums(replicas <= true.sum)/(niter+1)),2,min)
    prob_causal_eff <- data.frame(bayes.pval = p, pct.causal.eff = (1-p)*100)

    # Outputting estimated effect(s) and p-values(s)
    output <- cbind(effect, prob_causal_eff)
    rownames(output) <- dimnames(x$y)[[2]]

    # saving output
    summary_impact[["horizon_default"]] <- output
  }

  class(summary_impact) <- "summary.CausalMBSTS"
  summary_impact
}



#' Format and print the estimated causal effect(s), credible interval(s), and Bayesian p-value(s) into a clear output.
#'
#'
#' @param object An object of class 'summary.CausalMBSTS', a result of a call to \code{\link{summary.CausalMBSTS}}.
#' @param digits Number of digits to display.
#' @param ... Additional arguments passed to other methods.
#' @return Invisibly, \code{object}.
#' @export

print.summary.CausalMBSTS <- function(object, digits = max(3, getOption("digits") - 3), ...){

  for(i in 1:length(object)){
    object[[i]] <- round(object[[i]], digits = digits)
    cat("\nResults for ",names(object)[i],":\n\n", sep="")

    # effect
    average.ci <- paste("(", object[[i]][,"lower"], ",", object[[i]][,"upper"], ")" , sep = "")
    cumulative.ci <- paste("(", object[[i]][,"cum.lower"], ",", object[[i]][,"cum.upper"], ")" , sep = "")
    impact <- data.frame(object[[i]][,"mean"], average.ci, object[[i]][,"cum.sum"], cumulative.ci, object[[i]][,"bayes.pval"], object[[i]][,"pct.causal.eff"])
    rownames(impact) <- rownames(object[[i]])
    colnames(impact) <- c("Avg. Effect", "95% CI", "Cumulative Effect", "95% CI", "Bayesian p-value", "Prob. of a causal effect (%)")
    print(impact)
  }
  cat("\n")
  return(invisible(object))
}
