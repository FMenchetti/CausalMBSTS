
#' Plotting function for object of class CausalMBSTS
#'
#' Function to plot the results of the causal analysis using 'causal.mbsts'
#'
#' "impact" plots the general causal effect at every time points in the post period.
#' "forecast" plots the observed time series against the predicted counterfactual.
#' "ppchecks" draws posterior predictive checks for the model estimated in the pre period.
#' "inclusion.probs" plots the regressors inclusion probabilities.
#'
#' @importFrom forecast Acf
#' @param CausalMBSTS Object of class 'CausalMBSTS'
#' @param int.date Date of the intervention.
#' @param type A character string indicating the type of plot to be produced. Possible values in 'type = c("impact","forecast","ppchecks", "inclusion.probs")'.
#'
#' @return
#' @export
#'
#' @examples
#' plot(causal)
plot.CausalMBSTS<-function(CausalMBSTS, int.date, type = c("impact", "forecast", "ppchecks")){

  # given the estimated models saved in RDS format, the function draws:
  # i) causal effect plots; ii) observed vs forecast plots;
  # iii) posterior predictive checks; iv) inclusion probabilities.

  # Args:
  #  ...


      ### Causal effect plot
      if("impact" %in% type){
        plotImpact(CausalMBSTS, int.date = int.date)
      }

      ### Plot Observed vs Forecast
      if("forecast %in% type"){
        plotForecast(CausalMBSTS, int.date = int.date)
      }

      ### Posterior predictive checks
      if("ppchecks %in% type"){
        plotChecks(CausalMBSTS, int.date = int.date)
      }

      # Regressor index plot
      if("inclusion.probs" %in% type){
        plotInclusionProb(CausalMBSTS)
      }

    }

# ----------------------------------------------------------------------------------------

plotImpact <- function(CausalMBSTS, int.date){
  dates<-CausalMBSTS$dates
  dim <- dim(CausalMBSTS$mean.general)

  if (!is.null(dim)){
    d <- dim[2]
    for(i in 1:d){
      ylim<-c(min(CausalMBSTS$lower.general[,i]), max(CausalMBSTS$upper.general[,i]))
      x <- dates[dates >= int.date]
      main<- paste("Pointwise impact ","Y",i, sep="")
      plot(y=CausalMBSTS$mean.general[,i], x=x,type="l",col="blue",ylim=ylim,
           main= main, ylab="", xlab="")
      lines(y=CausalMBSTS$upper.general[,i],x=x, lty = 2)
      lines(y=CausalMBSTS$lower.general[,i],x=x, lty = 2)
    }
  } else {
    for(j in 1:length(CausalMBSTS$mean.general)){
      dim <- dim(CausalMBSTS$mean.general[[j]])
      start<-which(dates==int.date)
      end<-start + dim[1]
      x<-dates[start:(end-1)]
      d <- dim[2]
      for(i in 1:d){
        ylim<-c(min(CausalMBSTS$lower.general[[j]][,i]), max(CausalMBSTS$upper.general[[j]][,i]))
        main<- paste("Pointwise impact ","Y",i, sep="")
        plot(y=CausalMBSTS$mean.general[[j]][,i], x=x,type="l",col="blue",ylim=ylim,
             main= main, ylab="", xlab="")
        lines(y=CausalMBSTS$upper.general[[j]][,i],x=x, lty = 2)
        lines(y=CausalMBSTS$lower.general[[j]][,i],x=x, lty = 2)
      }
    }
  }
}

#------------------------------------------------------------------------------------------------

plotForecast <- function(CausalMBSTS, int.date){
  dates<-CausalMBSTS$dates
  y<-CausalMBSTS$y
  start<-which(dates==int.date) - round(0.4*sum(dates < int.date))
  post.mean <- apply(CausalMBSTS$predict$post.pred, c(1,2), mean)

  if (!is.list(CausalMBSTS$mean.general)){
    for(i in 1:dim(y)[2]){
      end<-dim(y)[1]
      x<-dates[start:end]
      yi<-y[start:end,i]
      ylim<-c(min(yi,post.mean[start:end,i]), max(yi,post.mean[start:end,i]))
      plot(y=yi,x=x,type="l",ylim=ylim,ylab="")
      lines(post.mean[start:end,i],col="blue",x=x)
      abline(v=int.date,col="red")
    }
  } else {
    for(j in 1:length(CausalMBSTS$mean.general)){
      dim <- dim(CausalMBSTS$mean.general[[j]])
      end<-which(dates == int.date) + dim[1]
      x<-dates[start:(end-1)]
      d <- dim[2]
      for(i in 1:d){
        yi<-y[start:(end-1),i]
        ylim<-c(min(yi,post.mean[start:(end-1),i]), max(yi,post.mean[start:(end-1),i]))
        plot(y=yi,x=x,type="l",ylim=ylim,ylab="")
        lines(post.mean[start:(end-1),i],col="blue",x=x)
        abline(v=int.date,col="red")
      }
    }
  }
}

#--------------------------------------------------------------------------------------------

plotChecks<-function(CausalMBSTS, int.date){
  dates<-CausalMBSTS$dates
  ind<-dates<int.date
  y<-CausalMBSTS$y
  post.pred<-CausalMBSTS$predict$post.pred.0
  post.pred.mean<-apply(post.pred,c(1,2),mean)

  for(i in 1:dim(y)[2]){
    # Density of posterior mean vs density of the data before intervention
    plot(density(y[ind,i]), xlab="",ylab="",main="")
    lines(density(post.pred.mean[,i]),col="blue")

    # Histograms & Bayesian p-value
    max.distrib<-apply(post.pred,c(2,3),max)
    pvalue<-sum(max.distrib[i,] >= max(y[ind,i]))/ncol(max.distrib)
    hist(max.distrib[i,], 30, col="lightblue",border="grey", main = paste("p.value =",round(pvalue,2)))
    abline(v=max(y[ind,i]), col="darkblue",lwd=3)

    # Residual plots
    y.rep<-matrix(y[ind,i], nrow(y[ind,]), (CausalMBSTS$mcmc$niter-CausalMBSTS$mcmc$burn), byrow = F)
    res<-(y.rep - (post.pred[,i,]-CausalMBSTS$mcmc$eps.samples[,i,]))
    std.res<-t(apply(res,i,FUN="/",sqrt(CausalMBSTS$mcmc$Sigma.eps[i,i,])))
    qqnorm(rowMeans(std.res),main="")
    qqline(rowMeans(std.res))
    Acf(rowMeans(std.res),main="")
  }
}

#-----------------------------------------------------------------------------------------------

plotInclusionProb<-function(CausalMBSTS, prob = prob){
  if(is.null(prob)){prob <- 0.5} # check here: if no reg above that prob threshold is present stop and write a message (e.g. no reg above that threshold, try to lower the threshold)
  X<-CausalMBSTS$mcmc$X
  select<- apply(CausalMBSTS$mcmc$Z.beta, 2, mean)>prob # selecting regressors with prob of inclusion bigger than some amount
  plot(apply(CausalMBSTS$mcmc$Z.beta,2,mean)[select], type="h", lwd=2,
       ylab=expression(paste("Pr(",italic(z[j] == 1),"|",italic(y),",X)",sep="")),
       las=2, xaxt="none", xlab="", cex.lab=2)
  axis(1, seq(1,ncol(X[,select]),1),las=2, labels = gsub(colnames(X)[select],pattern = "pezzi.", replacement = ""), cex.axis=2)
}
