######################################################################################
######################################################################################
####  Author:           Fiammetta Menchetti                                       ####
####                                                                              ####
####  Date last update: 2020-03-06                                                ####
####                                                                              ####
####  Content:          Prediction of a given MBSTS model                         ####
####                                                                              ####
####                                                                              ####
####  Main function :   predict.mbsts                                             ####
####  Dependencies:     mbsts.mcmc (lpy.X,block.m)                                ####
####                                                                              ####
######################################################################################
######################################################################################



predict.mbsts<-function(mbsts, newdata){
  
  # Given an object of class 'mbsts' and 'newdata' provided in matrix form, this function
  # provides in-sample forecasts and out-of-sample forecasts, both based on drawing
  # from the posterior predictive distribution. Note that NA values are not allowed in 
  # the new regressor matrix. 
  #
  # Args:
  #   mbsts   : an object of class 'mbsts'
  #   newdata : matrix of new data
  #
  # Value:
  #   post.pred.0 : T x p x niter array of in-sample forecasts
  #   post.pred.1 : S x p x niter array out-of-sample forecasts, where S is the number of
  #                 of forecasted periods (set to the length of provided new data)
  #   post.pred   : (T + S) x p x niter array combining in- and out-of-sample forecasts
  
  ### Dimensionalities & other objects
  
  # get dim
  t    <-dim(mbsts$y)[1]                  # number of time points
  K    <-dim(mbsts$R)[2]                  # tot number of state disurbances
  p    <-dim(mbsts$y)[2]                  # number of time series
  k    <-K/p                              # number of state disturbances for each time series
  M    <-dim(mbsts$T)[1]                  # tot number of states 
  step <-dim(newdata)[1]                  # number of obs to forecast
  niter<-mbsts$niter - mbsts$burn         # number of simulations
  last <- dim(mbsts$eta.samples)[1]       # last draw from the posterior p(alpha_t|Y_n)
  
  
  ### Empty arrays to store iterations
  eta.new<-array(NA, c(step,K,niter))
  
  states.new<-array(NA, c(step,M,niter))
  colnames(states.new)<-colnames(mbsts$states.samples)
  
  y.star<-array(NA, c(step,p,niter))
  
  post.pred.0<-array(NA, c(t,p,niter))           # storing in-sample draws  
  post.pred.1<-array(NA, c(step,p,niter))        # storing out-of-samples draws 
  post.pred  <-array(NA, c(t+step,p,niter))      # all together
  
  
  for(i in 1:niter){
    
    ### STEP 1: Get in-samples draws from the ppd.
    
    # Details:  let theta = (alpha, beta, Sigma.eps, z, Sigma.eta),
    #           during the MCMC we got samples from the joint posterior distribution p(theta|y). So
    #           now we can draw new values y.tilde from the posterior predictive distribution by simply
    #           taking the "niter" draws from the joint p(theta|y) and substitute them into model equations
    #           (we assume that y.new comes from the same distribution of y so that y.new is 
    #           independent of y given theta).
    
    post.pred.0[,,i]<-tcrossprod(mbsts$states.samples[,,i], mbsts$Z[,,1]) + mbsts$X %*% mbsts$beta[,,i]  + mbsts$eps.samples[,,i]
    
    ### STEP 2: Get new out-of-samples draws from the ppd.
    
    # Details: this time we should point out that an out-of sample draw y.new is no more independent
    #          of past y given theta. To see that, let's say we want to sample y_t+k | Y_t.
    #          This time, theta = (alpha_t+k,...,alpha_t+1, alpha_t, beta, Sigma.eps, z, Sigma.eta)
    #          and from our MCMC we have just the posterior of (alpha_t, beta, Sigma.eps, z, Sigma.eta)
    #          because we sampled from the full conditional p(alpha_t | Y_t, beta, Sigma.eps, z, Sigma.eta).
    #          Thus, to sample from y_t+k | y_t we should "recover" the missing dependence structure
    #          p(alpha_t+k,...,alpha_t+1 | theta') whre theta' is our "old" theta, that is 
    #          theta' = alpha_t, beta, Sigma.eps, z, Sigma.eta. More details in the pdf. 
    #          Conversely, y*_t+k =  y_t+k - Z alpha_t+k is independent of y*_t given theta',  
    #          because there's no more time component and we can just write y*. 
    #          During the MCMC we sampled from the joint posterior p(beta,Sigma.eps,z|y*),
    #          and as in step 1 we can use those draws to sample from the 
    #          posterior predictive distribution y.new* | y* by simply taking the draws from 
    #          p(beta,Sigma.eps,z|y*) and then draw from p(y.new*|beta,Sigma.eps,z).  
    #           
    
    # 2.1. Sampling new states
    #if(seed){set.seed(seed)}
    eta.new<-matrix(mvrnorm(1,rep(0,K),mbsts$Sigma.eta[,,i]),nrow=K)
    states.new[1,,i]<-mbsts$T[,,1] %*% matrix(mbsts$states.samples[last,,i],M) + mbsts$R[,,1] %*% eta.new
    
    ind<-if(step>1){2:step} else {NULL}
    for( j in ind){
      eta.new<-matrix(mvrnorm(1,rep(0,K),mbsts$Sigma.eta[,,i]),nrow=K)
      states.new[j,,i]<-mbsts$T[,,1] %*% matrix(states.new[j-1,,i],nrow=M) + mbsts$R[,,1]%*%eta.new 
    }
    
    # 2.2. Sampling y*.new.
    
    y.star[,,i]<-newdata %*% mbsts$beta[,,i] + mvrnorm(1, rep(0,p), mbsts$Sigma.eps[,,i])
    
    # 2.3. Out-of-sample draws
    
    post.pred.1[,,i]<-states.new[,,i] %*% t(mbsts$Z[,,1]) + y.star[,,i]
    
    ### STEP 3: Combining in-sample and out-of-samples draws
    
    post.pred[,,i]<-rbind(post.pred.0[,,i], post.pred.1[,,i])
  } 
  
  return(list(post.pred.0 = post.pred.0, post.pred.1 = post.pred.1, post.pred = post.pred ))  
}
