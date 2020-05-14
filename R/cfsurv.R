#' cfsurv
#'
#' The main function to generate a predictive conformal confidence interval for a unit's survival time.
#'
#' @param x the covariate for prediction point. (Currently only one-dimensional.)
#' @param r the censoring time for the prediction point.
#' @param data a data frame containing the training data. It should contain four columns: (X,R,event,censoredT). X is the covariate; R is the censoring tim;, event is the indicator if T<=R; censored_T is the censored survial time, i.e., the minimum of T and R.
#' @param alpha The miscoverage rate.
#' @param seed The random seed. Default is 24601.
#' @param model The model used to fit the quantile. Choices include "cox" and "randomforest". Default is "cox".
#' @param dist The distribution of T used in the cox model. Choices include "weibull", "exponential" and "gaussian". Default is "weibull".
#' @param h The bandwidth for the local confidence interval. Default is 1.
#' @export
#'
#' @examples
#' # Generate data
#' n <- 500
#' X <- runif(n,0,2)
#' T <- exp(X+rnorm(n,0,1)) 
#' R <- rexp(n,rate = 0.01)
#' event <- T<=R
#' censored_T <- pmin(T,R)
#' data <- data.frame(X=X,R=R,evemt=event,censored_T=censored_T)
#'
#' # Prediction point
#' x <- seq(0,2,by=.4)
#' r <- 2
#'
#' # Run cfsurv
#' res <- cfsurv(x,r,data,alpha=0.1,model="randomforest")
#'


# function to construct conformal confidence interval
cfsurv <- function(x,r,data,
                   alpha=0.05,
                   type="marginal",
                   seed = 24601,
                   model = "cox",
                   dist= "weibull",
                   h=1){
  ## Check if the required packages are installed
  ## Solution found from https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them  
  list.of.packages <- c("ggplot2",
                        "grf",
                        "quantregForest",
                        "randomForestSRC",
                        "survival",
                        "tidyverse",
                        "fishmethods",
                        "foreach",
                        "doParallel")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
  suppressPackageStartupMessages(res <- lapply(X=list.of.packages,FUN=require,character.only=TRUE))

  ## Process the input
  ## Check the length of x and r: only two cases are supported. length(r)=1, or length(r)=length(x)
  len_r <- length(r)
  len_x <- length(x)
  if(len_r>1 & len_r!=len_x){
    stop("The length of R is not compatible with that of X!")
  }

  ## Check the type of the model. Only "cox" and "randomforest" are supported
  if(model %in% c("cox","randomforest")==0) stop("The regression model is not supported.")

  ## Check the type of the confidence inteval
  if(type %in% c("marginal","local")==0) stop("The type of confidence interval is not supported.")

  ## Check the value of alpha
  if (alpha>=1 | alpha<=0) stop("The value of alpha is out of bound.")
  
  ## Check the columns of the data frame
  if(!is.data.frame(data))stop("data should be a data frame.")
  check_columnname <- "X"%in%names(data)& 
     "R"%in%names(data)&
     "event"%in%names(data)&
     "censored_T"%in%names(data)
  if(!check_columnname)stop("The columns in data are not correct. Please check the document.")

  ## set random seed
  set.seed(seed)

  ## Split the data into the training set and the calibration set
  n = dim(data)[1]
  n_train = n/2
  n_calib = n-n_train
  I_fit <- sample(1:n,n_train,replace = FALSE)
  data_fit <- data[I_fit,]
  data_calib <- data[-I_fit,] 
  
  ## Run the main function and gather resutls
  if(model == "cox"){
    res = cox_based(x,r,alpha,
                    data_fit,
                    data_calib,
                    type,
                    dist,
                    h)
   }
  if(model == "randomforest"){
    res = rf_based(x,r,alpha,
                   data_fit,
                   data_calib,
                   type,
                   h)
  }
  return(res)

  
}

