#' Conformalized survival analysis
#'
#' The main function to construct a predictive confidence interval 
#' for the survival time
#' @param x the covariate for prediction point
#' @param r the censoring time for the prediction point
#'

# function to construct conformal confidence interval
cfsurv <- function(x,r,data,alpha=0.05,seed = 24601,dist= "weibull"){
  fit <- quantreg_model(data,dist,seed)
  
  ## obtain final confidence interval
  res <- sapply(x,lower_ci,mdl=fit$mdl,r=r,alpha=alpha,data_calib = fit$data_calib)
  return(res)
}

## Internal functions
quantreg_model <- function(data,dist,seed){
  ## set random seed
  set.seed(seed)

  ## Split the data
  n = dim(data)[1]
  n_train = n/2
  n_calib = n-n_train
  I_fit <- sample(1:n,n_train,replace = FALSE)
  data_fit <- data[I_fit,]
  data_calib <- data[-I_fit,] 
  
  ## fitting the model
  mdl <- survreg(Surv(censored_T,event)~X,data=data_fit,dist=dist)

  ##   Return results
  return(list(mdl = mdl,data_calib=data_calib))
}

# construct nonconformity scores for one sided confidence interval with a lower bound
lower_score <- function(mdl,data_calib,alpha){
  res <- predict(mdl,
                 newdata = data.frame(X=data_calib$X,censored_T=data_calib$censored_T,event=data_calib$event),
                 type="quantile",
                 p=alpha)
  quant_low <-  res
  
  ## calibration
  score  <-  pmin(quant_low,data_calib$R) - data_calib$censored_T
  return(list(score = score))
}

# construct the one-sided confidence interval for T
lower_ci <- function(mdl,x,r,alpha,data_calib){
  res <- predict(mdl,
                 newdata = data.frame(X=data_calib$X),
                 type="quantile",
                 p=alpha)
  ##   browser()
  quant_low <-  res
  
  ## calibration
  score  <-  pmin(data_calib$R,quant_low)-data_calib$censored_T
 
  res <- predict(mdl,
                 newdata = data.frame(X=x),
                 type="quantile",
                 p=alpha)
  new_quant_low <-  res
  corr_term <- quantile(c(score,Inf),1-alpha)
  ci_low <- min(new_quant_low,r)-corr_term
  ci_low <- pmin(ci_low,r-1)
  ci_low <- pmax(ci_low,0)
  includeR <- ifelse(min(new_quant_low,r)-r<=corr_term,1,0)
  return(list(ci_low=ci_low,includeR = includeR,corr_term = corr_term))
}
