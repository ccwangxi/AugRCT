
<!-- README.md is generated from README.Rmd. Please edit that file -->

# AugRCT

<!-- badges: start -->

<!-- badges: end -->

The `AugRCT` package is used to compare multiple methods in augmenting the current control arms with historical data in Randomized Clinical Trials. The methods can be divided into __four subgroups__, 

* naive frequentist methods (no borrowing and pooling, i.e. full borrowing), 
* PS matching only/plus Bayesian information borrowing, 
* PS stratification only/plus Bayesian information borrowing, 
* IPTW only/plus Bayesian information borrowing. 

We require binary outcome(e.g. overall response rate) and the hypothesis test is superiority test favoring a higher value. This package includes codes for both simulation and real data application.

We have two R functions for __simulation__. For each simulation setting, 

1. `Simu.Hybrid` implements all the hybrid approaches (PS plus Bayesian)

2. `Simu.Freq` implements all the frequentist approaches (naive and PS only). 

From both functions, a data frame with point estimate and 95% CI of treatment effect in risk difference of all the methods will be returned. 

Similarly, there are two R functions in __real data application__. With available real data, 

3. `RealData.Hybrid` implements all the hybrid approaches (PS plus Bayesian)

4. `RealData.Freq` implements all the frequentist approaches (naive and PS only). 

From both functions, a data frame with point estimate and 95% CI of treatment effect in risk difference of all the methods will be returned.  

## Installation

You can install AugRCT by opening the `AugRCT.Rproj` file. Then load the `AugRCT` package.

``` r
library(AugRCT)
```

## Example

Please run `AugRCT-vignette.Rmd` under `AugRCT\vignettes` to check the example for simulation on R server at Merck and real data application. 
