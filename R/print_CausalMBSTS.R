
#' Printing method for an object of class CausalMBSTS
#'
#' The method prints the results of a causal analysis in a nice and concise output format.
#'
#' @param x An object of class 'CausalMBSTS'.
#'
#' @return
#' @export
#'
#' @examples
#' ## Example 2
#' set.seed(1)
#' t <- seq(from = 0,to = 4*pi, length.out=300)
#' y <- cbind(3*sin(2*t)+rnorm(300), 2*cos(2*t) + rnorm(300))
#' dates <- seq.Date(from = as.Date("2015-01-01"), by = "week", length.out=300)
#' int.date <- as.Date("2020-02-27")
#' y[dates >= int.date,] <- y[dates >= int.date,]+2
#'
#' # Model definition
#' model.2 <- model(y, components = c("trend", "cycle"), cycle.period = 75)
#' causal.2 <- causal.mbsts(model.2, dates = dates, int.date = int.date, s0.r = 0.01*diag(2), s0.eps = 0.1*diag(2), niter = 100, burn = 10)
#'
#' # Print
#' print(causal.2)

print.CausalMBSTS <- function(x){

  effect <- x$general.effect
  pred <- x$predict$post.pred
  niter <- dim(pred)[3]
  start <- dim(x$predict$post.pred.0)[1] +1 # intervention date
  summary_impact <- list()

  if (is.list(effect)){
    effect <- lapply(effect, round, digits = 2)

    for(i in 1:length(effect)){

      # effect
      average <- effect[[i]][,"mean"]
      cumulative <- effect[[i]][,"cum.sum"]
      average.ci <- paste("[", effect[[i]][,"lower"], ", ", effect[[i]][,"upper"], "]" , sep = "")
      cumulative.ci <- paste("[", effect[[i]][,"cum.lower"], ", ", effect[[i]][,"cum.upper"], "]" , sep = "")
      impact <- data.frame( average, average.ci, cumulative, cumulative.ci)
      colnames(impact) <- c("Avg. Effect", "95% CI", "Cumulative Effect", "95% CI")

      # Bayesian p-value
      end <- start + dim(x$mean.general[[i]])[1] -1 # interval between the intervention and the time horizon
      true.sum <- colSums(x$y[start:end,])
      pred.sum <- apply(pred[start:end,,], c(2,3), sum) # cumulative causal effect at each iteration
      replicas <- cbind(pred.sum,true.sum) # cumulative sum of the observations in each replicated data set
      p <- apply(rbind(rowSums(replicas >= true.sum)/(niter + 1),rowSums(replicas <= true.sum)/(niter+1)),2,min)
      prob_causal_eff <- data.frame(round(p,4), (1-round(p,4))*100)
      rownames(prob_causal_eff) <- rownames(impact) <- dimnames(x$y)[[2]]
      colnames(prob_causal_eff) <- c("Bayesian p-value", "Prob. of a causal effect (%)")

      # saving output
      summary_impact[[i]] <- list(impact, prob_causal_eff)
      names(summary_impact[[i]]) <- c("Effect", "P-value")
      }

    names(summary_impact) <- paste("horizon",1:length(effect))

  } else {
    effect <- round(effect, digits = 2)

    # effect
    average <- effect[,"mean"]
    cumulative <- effect[,"cum.sum"]
    average.ci <- paste("[", effect[,"lower"], ", ", effect[,"upper"], "]" , sep = "")
    cumulative.ci <- paste("[", effect[,"cum.lower"], ", ", effect[,"cum.upper"], "]" , sep = "")
    impact <- data.frame( average, average.ci, cumulative, cumulative.ci)
    colnames(impact) <- c("Avg. Effect", "95% CI", "Cumulative Effect", "95% CI")

    # Bayesian p-value

    end <- dim(x$y)[1]
    true.sum <- colSums(x$y[start:end,])
    pred.sum <- apply(pred[start:end,,], c(2,3), sum) # cumulative causal effect at each iteration
    replicas <- cbind(pred.sum,true.sum) # cumulative sum of the observations in each replicated data set
    p <- apply(rbind(rowSums(replicas >= true.sum)/(niter + 1),rowSums(replicas <= true.sum)/(niter+1)),2,min)
    prob_causal_eff <- data.frame(round(p,4), (1-round(p,4))*100)
    rownames(prob_causal_eff) <- rownames(impact) <- dimnames(x$y)[[2]]
    colnames(prob_causal_eff) <- c("Bayesian p-value", "Prob. of a causal effect (%)")

    # saving output
    summary_impact<- list(impact, prob_causal_eff)
    names(summary_impact) <- c("Effect", "P-value")
  }

  print(summary_impact)
}
