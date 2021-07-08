# DynBalancing
Dynamic Balancing for Estimating the Effects of Dynamic Treatments



The package can be used for causal inference with dynamic (time-varying) treatments.  The package works for balanced and unbalanced panels with high-dimensional (and time-varying) covariates.  Effects are estimated through penalized regression with dynamic balancing weights. 

## Installation

```
install.packages("remotes")
library(remotes) 
remotes::install_github("dviviano/DynBalancing")
library(DynBalancing)
```

## Implementation

The package contains the functions: 

- DynBalancing_ATE
- DynBalancing_History
- DynBalncing_Het_ATE


See the 
<a href="https://dviviano.github.io/blog/posts/package_illustration/index.html">  vignettes </a> for its implementation. Details on the arguments are available typing help(.) and the name of the function. 

## Reference 

Viviano, Davide, and Jelena Bradic. "Dynamic covariate balancing: estimating treatment effects over time." arXiv preprint arXiv:2103.01280 (2021)

## Support 

Davide Viviano: dviviano@ucsd.edu 

