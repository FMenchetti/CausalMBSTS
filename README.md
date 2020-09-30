# CausalMBSTS

## What is CausalMBSTS?
CausalMBSTS is an R package that can be used to infer a causal effect of an intervention 
on a multivariate time series. For example, if two items are frequently bought together (e.g., printer and ink) but only one of them gets a discount, what is the effect of the discount on the sales of both products? To answer this question, the dependence structure between the items must be considered. The package does so by using Multivariate Bayesian Structural Time Series models (MBSTS). 

## Why MBSTS models?
MBSTS models are flexible tools that also allow for a transparent way to deal with uncertainty. Flexibility comes from an option for a user to add sub-components (e.g., trend, seasonality, and cycle) that encapsulate the characteristics of a time series. Transparency comes from the Bayesian framework, in which our uncertainty surrounding model parameters is explicitly accounted for by setting prior distributions that represent our prior knowledge of the phenomenon. Then, by combining the priors with the data, we obtain posterior distributions that can be used for several purposes (e.g., compute a point estimate, credible intervals, perform prediction). The package implements the Gaussian linear model that, possibly after transformation of the dependent variable, provides an adequate representation for many time series that are analyzed in practice.     

## How does the package work?
The main function included in the package is <tt>CausalMBSTS()</tt>. The user is asked to provide the data, a vector of dates (including the intervention date) and to select the desired model components. The user can set custom values of the hyperparameters of the priors or use the default values. Then, the function uses a Gibbs sampler to draw from the joint posterior distribution of the model parameters before the intervention. The next step is the prediction of the counterfactual outcomes in the absence of intervention, which is done by sampling from the posterior predictive distribution (PPD). The causal effect at each point in time after the intervention is computed as the difference between the observed outcome and the PPD. Plotting and printing methods for the resulting object are provided. 

Aside from causal inference, the package can also be used for forecasting. This task can be completed with the function <tt>as.mbsts()</tt> coupled with the method <tt>predict.mbsts()</tt>: the former estimates the MBSTS model and outputs an object of class <tt>mbsts</tt>; the latter performs a prediction step for a given number of future periods. 

## Installing the package from GitHub

```{r}
# Get CausalMBSTS from GitHub
devtools::install_github("FMenchetti/CausalMBSTS", build_vignettes = TRUE)
```

## Further readings
Menchetti \& Bojinov (2020), Estimating causal effects in the presence of partial interference using multivariate Bayesian structural time series models <https://arxiv.org/abs/2006.12269>
