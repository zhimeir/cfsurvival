# cfsurvival
An R package that implements the conformalized survival analysis methodology.

## Overview
The goal of ***cfsurvival*** is to provide a lower predictive lower bound for the survival 
time of an individual, with the guarantee that with probability (1-$\alpha$) it is no larger
than the true survival time. 

## Installation
To install this package, run the following command in your R console:
```{r}
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("zhimeir/cfsurvival")
```
## Dependencies
The following [R](https://www.r-project.org/) packages (version 4.0.2) packages are required:
 - [quantreg](https://cran.r-project.org/web/packages/quantreg/index.html)
 - [grf](https://github.com/grf-labs/grf)
 - [quantregForest](https://cran.r-project.org/web/packages/quantregForest/index.html)
 - [randomForestSRC](https://cran.r-project.org/web/packages/randomForestSRC/index.html)
 - [survival](https://cran.r-project.org/web/packages/survival/index.html)
 - [tidyverse](https://www.tidyverse.org/)
 - [fishmethods](https://cran.r-project.org/web/packages/fishmethods/index.html)
 - [foreach](https://cran.r-project.org/web/packages/foreach/index.html)
 - [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)
 - [GauPro](https://cran.r-project.org/web/packages/GauPro/index.html)
 - [gbm](https://cran.r-project.org/web/packages/gbm/index.html)
 - [np](https://cran.r-project.org/web/packages/np/index.html)
 - [conTree](http://statweb.stanford.edu/~jhf/conTree/)

## Usage example 
```{r}

# Generate data
n <- 500
X <- runif(n,0,2)
T <- exp(X+rnorm(n,0,1))
C <- rexp(n,rate = 0.01)
event <- (T <= C)
time <- pmin(T,C)
data <- data.frame(X = X, C = C, event = event, censored_T = censored_T)

# Prediction point
n_test <- 10
X <- runif(n_test, 0, 2)

# Run cfsurv with c_0 = 
res <- cfsurv(x, c, X, C, event, time, alpha=0.1, model="cox")
```
## License 
This packakge is is distributed under the MIT license.

