# CausalMBSTS

## What is CausalMBSTS?
CausalMBSTS is an R package that can be used to infer the causal effect of an intervention 
on a multivariate time series. For example, if two items are frequentely bought together (e.g., printer and ink) but only one of them gets a discount, what is effect of the discount on the sales of both products? To answer this question, the dependence structure between the items must be taken into account. The package does so by using Multivariate Bayesian Structural Time Series models (MBSTS). 

## Why MBSTS models?
MBSTS models are flexible tools that also allow for a transparent way to deal with uncertainty. Flexibility comes from our ability to add sub-components (e.g., trend,
seasonality, and cycle) that encapsulate the characteristics of a time series. Transparency comes from the Bayesian framework, in which our uncertainty surrounding model parameters is explicitely accounted for by setting prior distributions that represent our prior knowledge of the phenomenon. Then, by combining the priors with the data, we obtain posterior distributions that can be used for a number of purposes (e.g., compute a point estimate, credible intervals, perform prediction). The package implements the Gaussian linear model that, sometimes after transformation of the dependent variable, provides an adequate representation for many time series that are analyzed in practice.     

## How does the package work?
The main function included in the package is \code{CausalMBSTS()}. The user is asked to provide the data, a vector of dates (including the intervention date) and to select the desired model components. The user can set custom values on the hyperparameters of the priors as well as use the default values. Then, the function uses a Gibbs sampler to draw from the joint posterior distribution of the model parameters before the intervention. The following step is the prediction of the counterfactual outcomes in the absence of intervention, which is done by sampling from the posterior predictive distribution (PPD). The causal effect at each point in time after the intervention is computed as the difference between the observed outcome and the PPD. Plotting and printing methods for the resulting object are provided. Aside from causal inference, the package can also be used for forecasting. This task can be completed with the function \code{as.mbsts} coupled with the method \code{predict.mbsts}: the former estimates the MBSTS model and outputs an object of class \code{mbsts}; the latter performs a prediction step for a given number of future periods. 

## Examples



## Further readings
