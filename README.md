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
 - [conTree](https://github.com/bnaras/conTree)

## Usage example 
```{r}
# Library the package
library(cfsurvival)

# Generate data
set.seed(24601)
n <- 2000
X <- runif(n, 0, 2)
T <- exp(X + rnorm(n,0,1))
C <- rexp(n, rate = 0.05)
event <- (T <= C)
censored_T <- pmin(T, C)
data <- data.frame(X = X, C = C, event = event, censored_T = censored_T)

# Prediction point
n_test <- 1000
X_test <- runif(n_test, 0, 2)
T_test <- exp(X_test + rnorm(n,0,1))

# Running cfsurvival under completely independent censoring with c0 = 30 
c0 <- 30
pr_list <- rep(0.5, n)
pr_new_list <- rep(0.5, n_test)

# Use the Cox model
res <- cfsurv(x = X_test, c_list = c0, pr_list = pr_list, pr_new_list = pr_new_list,
             Xtrain = X, C = C, event = event, time = censored_T, 
             alpha = 0.1, model = "aft")

# Examine the result
cat(sprintf("The coverage is %.3f.\n", mean(res <= T_test)))
```
## License 
This packakge is is distributed under the MIT license.

