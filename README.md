# CausalMBSTS

## What is CausalMBSTS?
CausalMBSTS is an R package that can be used to infer a causal effect of an intervention 
on a multivariate time series. For example, if two items are frequently bought together (e.g., printer and ink) but only one of them gets a discount, what is the effect of the discount on the sales of both products? To answer this question, the dependence structure between the items must be considered. The package does so by using Multivariate Bayesian Structural Time Series models (MBSTS). 

## Why MBSTS models?
MBSTS models are flexible tools that also allow for a transparent way to deal with uncertainty. Flexibility comes from an option for a user to add sub-components (e.g., trend, seasonality, and cycle) that encapsulate the characteristics of a time series. Transparency comes from the Bayesian framework, in which our uncertainty surrounding model parameters is explicitly accounted for by setting prior distributions that represent our prior knowledge of the phenomenon. Then, by combining the priors with the data, we obtain posterior distributions that can be used for several purposes (e.g., compute a point estimate, credible intervals, perform prediction). The package implements the Gaussian linear model that, possibly after transformation of the dependent variable, provides an adequate representation for many time series that are analyzed in practice.     

## How does the package work?
The main function included in the package is <tt>CausalMBSTS()</tt>. The user is asked to provide the data, a vector of dates (including the intervention date) and to select the desired model components. The user can set custom values of the hyperparameters of the priors or use the default values. Then, the function uses a Gibbs sampler to draw from the joint posterior distribution of the model parameters before the intervention. The next step is the prediction of the counterfactual outcomes in the absence of intervention, which is done by sampling from the posterior predictive distribution (PPD). The causal effect at each point in time after the intervention is computed as the difference between the observed outcome and the PPD. Plotting and printing methods for the resulting object are provided. 

Aside from causal inference, the package can also be used for forecasting. This task can be completed with the function <tt>as.mbsts()</tt> coupled with the method <tt>predict.mbsts()</tt>: the former estimates the MBSTS model and outputs an object of class <tt>mbsts</tt>; the latter performs a prediction step for a given number of future periods. 

## Worked examples
The example below shows how to use CausalMBSTS to infer the effect of an intervention on a multivariate time series. First, we need to generate some random data and a fictional intervention date. To keep the interpretation simple, assume we have weekly sales counts from two related products (winter wool hats and gloves) and that, say, on February 27, 2020 the store manager introduced a big price discount on all gloves to boost sales. 

```{r}
# Get CausalMBSTS from GitHub
library(devtools)
install_github("FMenchetti/CausalMBSTS")

# Set seed & random data generation
set.seed(1)
t <- seq(from = 0,to = 4*pi, length.out=300)
y <- cbind(3*sin(2*t)+rnorm(300), 2*cos(2*t) + rnorm(300))
dates <- seq.Date(from = as.Date("2015-01-01"), by = "week", length.out=300)
int.date <- as.Date("2020-02-27")
y[dates >= int.date,] <- y[dates >= int.date,]+2

# Some plots
plot(y = y[,1], x=dates, type="l", col="cadetblue")
lines(y = y[,2], x = dates, col = "orange")
abline(v=int.date, col="red")
```

From the plots of the two time series it is possible to notice a clear cyclical pattern with a period of approximately 1.5 years, equivalent to 75 weeks. Thus, we can try to estimate a model with two subcomponents: trend and cycle. To speed up computations, we set <tt>niter = 100</tt> but we recommend to do at least 1000 iterations and check the convergence of the Markov chain. 

```{r}
# Causal effect estimation
causal.2 <- CausalMBSTS(y, components = c("trend", "cycle"), cycle.period = 75,
                        dates = dates, int.date = int.date, s0.r = 0.01*diag(2),
                        s0.eps = 0.1*diag(2), niter = 100, burn = 10)
summary(causal.2)

# Causal effect plot
par(mfrow=c(1,2))
plot(causal.2, int.date = int.date, type = "impact")

# Observed sales vs counterfactual sales plot
plot(causal.2, int.date = int.date, type = "forecast")
```

The results show a significant causal effect of the intervention on the sales of both products: on average, sales of gloves and hats increased of 2 units every week; at the end of the analysis period, the sales attributable to the discount amounted to approximately 74 units for each product. 

The Bayesian p-value reported in the summary measures to what extent this result is due to chance. It is computed as the proportion of times that the cumulative sum of sales exceeds the observed sum; thus, a huge number indicate that the counterfactual sales are often higher than the observed sales after the discount, which clearly contradicts the hypothesis of a positive causal effect. In this case, the p-values are very small (0.011 and 0.022), yielding a probability of a causal effect around 98 % for both the time series. 

Before reporting the results of any analysis based on a parametric model, it is recommended to check model adequacy. In a Bayesian framework, this can be done with posterior predictive checks. 

```{r}
# Posterior predictive checks
par(mar = c(2,2,2,2)) ; par(mfrow = c(2,4))
plot(causal.2, int.date = int.date, type = "ppchecks")
```
Starting from the left, the above plots show: the density of observed data vs. the posterior predictive mean (good overlap means that the model fits the data); ii) observed maximum compared to the distribution of the maximum from the posterior draws (extreme p-values indicate that the model produces extreme observations); iii) Normal QQ-Plot of standardized residuals; iv) autocorrelation function of standardized residuals. In this case, the posterior predictive checks are in favor of the model, but it is also important to check the convergence of the Markov chain with some trace plots (for illustration purposes, only the trace plots of the variance-covariance matrices of the observation and trend disturbances are displayed).

```{r}
# Trace plots
par(mar = c(2,2,2,2)) ; par(mfrow = c(4,3))
mcmc <- causal.2$mcmc
plot(mcmc$Sigma.eps[1,1,], type = "l", main = "Variance of obs residuals Y1")
plot(mcmc$Sigma.eps[2,2,], type = "l", main = "Variance of obs residuals Y2")
plot(mcmc$Sigma.eps[1,2,], type = "l", main = "Covariance of obs residuals")
plot(mcmc$Sigma.r[1,1,], type = "l", main = "Variance of trend residuals Y1")
plot(mcmc$Sigma.r[2,2,], type = "l", main = "Variance of trend residuals Y2")
plot(mcmc$Sigma.r[1,2,], type = "l", main = "Covariance of trend residuals")
```
The trace plots show that convergence is soon reached for the variance-covariance matrix of the observation disturbances, whereas the variance-covariance matrix of the trend disturbances doesn't seem to have converged. Thus, it would be better to repeat the above analysis with a higher number of iterations.

The second example shows how to use CausalMBSTS to forecast the evolution of a multivariate time series. In line with the above example, let's assume that at the end of August 2019 the store manager wants to predict future sales of gloves and hats for the entire winter season (approximately 26 weeks) before placing an order to the suppliers. First, we create an MBSTS model (the same as the one above) and then we use this model to forecast sales.

```{r}
# Set seed & random data generation
set.seed(1)
t <- seq(from = 0,to = 4*pi, length.out=222)
y <- cbind(3*sin(2*t)+rnorm(222), 2*cos(2*t) + rnorm(222))
dates <- seq.Date(from = as.Date("2015-06-01"), by = "week", length.out=222)

# MBSTS model definition
sales_model <- as.mbsts(y, components = c("trend", "cycle"), cycle.period = 75,
                        s0.r = 0.01*diag(2), s0.eps = 0.1*diag(2), niter = 100, burn = 10) 

# Prediction step
pred <- predict(sales_model, steps.ahead = 26)

# Some plots
par(mfrow = c(1,1))
y.mean <- apply(pred$post.pred,c(1,2),mean)
new.dates <- seq.Date(from = as.Date("2015-06-01"), by = "week", length.out=248)
plot(y = y.mean[,1], x = new.dates, type="l", col="cadetblue")
lines(y = y.mean[,2], x = new.dates, col = "orange")
abline( v = dates[222], col = "red")
```

## Further readings
Menchetti \& Bojinov (2020), Estimating causal effects in the presence of partial interference using multivariate Bayesian structural time series models <https://arxiv.org/abs/2006.12269>
