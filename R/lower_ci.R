#' lower_ci
#'
#' construct the one-sided confidence interval for a unit's survival time T
#'
#' @param mdl The fitted model
#' @param x The covariate of the prediction point
#' @param r The observe time of the prediction point
#' @param alpha The miscaverage rate.
#' @param data A data frame used for calibration. Should contain four columns: (X,R,event,censored_T). 
#' @param mdl The fitted model to estimate the conditional quantile. The default is NULL.
#' @param quant_lo The fitted conditional quantile for the calibration data.
#' @param new_quant_lo The fitted conditional quantile for the tprediction data.
#'
#' @export


lower_ci <- function(x,r,alpha,
                     data=data,
                     mdl=NULL,
                     quant_lo=NULL,
                     new_quant_lo=NULL){
  ## Add a bit of noise to the score to break the ties
  sigma_noise <- 1e-6
  len_x <- length(x)
  len_r <- length(r)
  n <- dim(data)[1]

  ## The nonconformty scores
  if(is.null(quant_lo)){
  res <- predict(mdl,
                 newdata = data.frame(X=data$X),
                 type="quantile",
                 p=alpha)
  quant_lo <-  res  
  }
  ## calibration
  score  <-  pmin(data$R,quant_lo)-data$censored_T+rnorm(n,0,sigma_noise)
  if(is.null(new_quant_lo)){
  res <- predict(mdl,
                 newdata = data.frame(X=x),
                 type="quantile",
                 p=alpha)
  new_quant_lo <-  res
  }
  corr_term <- quantile(c(score,Inf),1-alpha,type=1)
  extra_noise <- rnorm(len_x,0,sigma_noise)
  ci_low <- pmin(new_quant_lo,r)-corr_term+extra_noise
  ci_low <- pmin(ci_low,r)
  ci_low <- pmax(ci_low,0)
  includeR <- ifelse(pmin(new_quant_lo,r)-r+extra_noise<=corr_term,1,0)
  ##   return(list(ci_low=ci_low,includeR = includeR,corr_term = corr_term))
  return(list(ci_low=ci_low,includeR = includeR))
}
