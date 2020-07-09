#' Predictive confidence interval for survival data
#'
#' The main function to generate a predictive conformal confidence interval for a unit's survival time.
#'
#' @param x a vector of the covariate for test point. 
#' @param r the censoring time for the test point.
#' @param Xtrain a n-by-p matrix of the covariate of the training data.
#' @param R a length n vector of the censoring time of the training data.
#' @param event a length n vector of indicators if the time observed is censored. TRUE corresponds to NOT censored, and FALSE censored.
#' @param time  a vevtor of length n, containing the observed survival time.
#' @param alpha a number between 0 and 1, speciifying the miscoverage rate.
#' @param seed an integer random seed (default: 24601).
#' @param model Options include "cox", "randomforest", "Powell", "Portnoy" and "PengHuang". This determines the model used to fit the condditional quantile (default: "cox").
#' @param dist either "weibull", "exponential" or "gaussian" (default: "weibull"). The distribution of T used in the cox model. 
#' @param h the bandwidth for the local confidence interval. Default is 1.
#'
#' @return low_ci a value of the lower bound for the survival time of the test point.
#' @return includeR 0 or 1, indicating if [r,inf) is included in the confidence interval.
#'
#' @examples
#' # Generate data
#' n <- 500
#' X <- runif(n,0,2)
#' T <- exp(X+rnorm(n,0,1))
#' R <- rexp(n,rate = 0.01)
#' event <- T<=R
#' time <- pmin(T,R)
#' data <- data.frame(X=X,R=R,event=event,censored_T=censored_T)
#' # Prediction point
#' x <- seq(0,2,by=.4)
#' r <- 2
#' # Run cfsurv
#' res <- cfsurv(x,r,X,R,event,time,alpha=0.1,model="cox")
#'
#' @export

# function to construct conformal confidence interval
cfsurv <- function(x,r,Xtrain,R,event,time,
                   alpha=0.05,
                   type="marginal",
                   seed = 24601,
                   model = "cox",
                   dist= "weibull",
                   I_fit = NULL,
                   h=1){
  ## Check if the required packages are installed
  ## Solution found from https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
  list.of.packages <- c("ggplot2",
                        "quantreg",
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
  X <- Xtrain
  len_r <- length(r)
  if(is.null(dim(x)[1])){
    len_x <- length(x)
    p <- 1
  }else{
    len_x <- dim(x)[1]
    p <- dim(x)[2]
  }
  
  
  if(is.null(dim(X)[1])){
    n <- length(X)
    pX <- 1
  }else{
    n <- dim(X)[1]
    pX <- dim(X)[2]
  }
  
  
  if(len_r>1 & len_r!=len_x){
    stop("The length of R is not compatible with that of X!")
  }

  ## Check the type of the model. Only "cox" and "randomforest" are supported
  if(model %in% c("cox","cox_scaled","randomforest","pow","portnoy","PengHuang")==0) stop("The regression model is not supported.")

  ## Check the type of the confidence inteval
  if(type %in% c("marginal","local")==0) stop("The type of confidence interval is not supported.")

  ## Check the value of alpha
  if (alpha>=1 | alpha<=0) stop("The value of alpha is out of bound.")

  ## Check the dimensions of the data 
  xnames <- paste0('X', 1:p)
  if(n != length(R))stop("The number of rows in X does not match the length of R.")
  if(length(R) != length(event))stop("The length of R does not match the length of event.")
  if(length(event) != length(time))stop("The length of event does not match the length of time.")
  if(p != pX) stop("The dimension of the test point does not match the dimension of the training point.")

  data <- as.data.frame(cbind(R,event,time,X))
  colnames(data) <- c("R","event","censored_T",xnames)
  ## set random seed
  set.seed(seed)

  ## Split the data into the training set and the calibration set
  n = dim(data)[1]
  n_train = n/2
  n_calib = n-n_train
  if(is.null(I_fit)){
    I_fit <- sample(1:n,n_train,replace = FALSE)
  }
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
  if(model == "cox_scaled"){
    res = cox_scaled(x,r,alpha,
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
  if(model == "pow"){
    res = pow_based(x,r,alpha,
                   data_fit,
                   data_calib,
                   type,
                   h)
  }
  if(model == "portnoy"){
    res = portnoy_based(x,r,alpha,
                   data_fit,
                   data_calib,
                   type,
                   h)
  }
  if(model == "PengHuang"){
    res = ph_based(x,r,alpha,
                   data_fit,
                   data_calib,
                   type,
                   h)
  }
  return(res)


}

