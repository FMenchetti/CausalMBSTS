######################################################################################
######################################################################################
####  Author:           Fiammetta Menchetti                                       ####
####                                                                              ####
####  Date last update: 2020-03-06                                                ####
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


causal.mbsts<-function(Smodel, X, y, dates, int.date, horizon = NULL, H = NULL, nu0.k = NULL, s0.k, 
                       nu0.eps = NULL, s0.eps, niter, burn = NULL, ping = NULL){
  
  # It estimates a specific joint effect in a multivariate time series
  # setting. It uses MCMC to sample from the joint posterior distribution of the parameters
  # of an MBSTS model before the intervention/treatment occurred. Then, it uses the
  # covariates post intervention to predict the counterfactual potential outcome 
  # The prediction is done by sampling from the posterior predictive 
  # distribution (ppd). Then the causal effect is computed by taking the difference between 
  # the observed outcome of each time series and the mean of the ppd (credible intervals are
  # computed accordingly).
  
  # Args:
  #   Smodel   : a multivariate state space model of class 'SSModel'
  #   X        : a T x N matrix of predictors
  #   dates    : a vector of dates
  #   int.date : date of the intervention
  #   H        : desidered N x N variance-covariance matrix between regression coefficients.
  #             The default is Zellner's g-prior, H = (X'X)^(-1) 
  #   nu0.k    : degrees of freedom of the Inverse-Wishart prior for each Sigma_k.
  #             The default is the smallest integer such that the expectation of eta_k exists,
  #             that is, nu0.k = p + 2 where p is the number of time series in the 
  #             multivariate model
  #   nu0.eps  : degrees of freedom of the Inverse-Wishart prior for Sigma_eps, the default is p+2.
  #   s0.k     : Scale matrix of the Inverse-Wishart prior for each Sigma_k.  
  #   s0.eps   : Scale matrix of the Inverse-Wishart prior for Sigma.eps      
  #   niter    : number of MCMC iteration
  #   burn     : desidered burn-in, set by default to 0.1 * niter
  #   ping     : logical, if TRUE a status message it's printed every iterations decile. 
  #             Defaults to TRUE. 
  #
  # Value:
  #   joint.effect : the estimated average causal effect and a 95% credible interval
  #   joint        : joint causal effect for all iterations
  #   mean.joint   : pointwise effect
  #   lower.joint  : lower bound of the pointwise effect
  #   upper.joint  : upper bound of the pointwise effect
  #   mbsts        : an object of class 'mbsts'
  #   predict      : a list with the same components as those produced by the function 'predict.mbsts'
  #   adj.series   : observations in the analysis period excluding holidays
  #   adj.dates    : dates in the analysis period without holidays
  
  ### STEP 1. Dividing pre and post periods
  
  ind<-dates < int.date
  X.pre  <-X[ind,]
  X.post <-X[!ind,]
  y.pre  <-y[ind,]
  y.post <-y[!ind,]
  
  # Estimating the model only in the pre-period
  Smodel$y<-y.pre  
  attr(Smodel,"n")<-as.integer(nrow(y.pre))
  
  ### STEP 2. MCMC 
  
  mbsts<-mbsts.mcmc(Smodel = Smodel, X = X.pre, H = NULL, nu0.k = nu0.k, s0.k = s0.k, nu0.eps = nu0.eps, 
                    s0.eps = s0.eps, niter=niter, burn=burn, ping=ping)
  
  ### STEP 3. In- and out-of-sample forecasts from the ppd
  predict<-predict.mbsts(mbsts, X.post)
 
 
  ### STEP 4. Causal effect estimation 
  p<-dim(mbsts$y)[2]
  burn<-mbsts$burn
  y.post.rep<-array(y.post, c(nrow(y.post),ncol(y.post),niter-burn)) 
  y.diff<-y.post.rep-predict$post.pred.1 
  
  # removing holidays 
  holi<-hol.dummy(dates[dates>=int.date]) 
  y.diff<-y.diff[holi==0,,] 
  adj.dates<-dates[hol.dummy(dates)==0]
  
  # Joint causal effect
  
  if(length(horizon)>0){ mean.effect<-list() ; lower.bound<-list() ; upper.bound<-list() ; joint.effect<-list()
    for(i in 1:length(horizon)){
      ind<-adj.dates[adj.dates>=int.date] <= horizon[i]
      mean.effect[[i]]<-apply(y.diff[ind,,],c(1,2),mean) 
      lower.bound[[i]]<-apply(y.diff[ind,,],c(1,2), quantile, probs = 0.025)
      upper.bound[[i]]<-apply(y.diff[ind,,],c(1,2), quantile, probs = 0.975)
      joint.effect[[i]]<-cbind(mean  = apply(colMeans(y.diff[ind,,]),1,mean), 
                          lower = apply(colMeans(y.diff[ind,,]),1,quantile, probs=0.025),
                          upper = apply(colMeans(y.diff[ind,,]),1,quantile, probs=0.975))
    }
  } else {
    mean.effect<-apply(y.diff,c(1,2),mean) 
    lower.bound<-apply(y.diff,c(1,2), quantile, probs = 0.025)
    upper.bound<-apply(y.diff,c(1,2), quantile, probs = 0.975)
    joint.effect<-cbind(mean  = apply(colMeans(y.diff),1,mean), 
                        lower = apply(colMeans(y.diff),1,quantile, probs=0.025),
                        upper = apply(colMeans(y.diff),1,quantile, probs=0.975))
  }
  
  return(list(mcmc = mbsts, predict = predict, adj.series = y[hol.dummy(dates)==0,], 
              adj.dates = adj.dates, joint = y.diff, joint.effect = joint.effect,  
              mean.joint = mean.effect, lower.joint = lower.bound, upper.joint = upper.bound,
              original.series = y, original.dates = dates))
}

#--------------------------------------------------------------------------------------
