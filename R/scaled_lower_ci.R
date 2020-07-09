#' One-sided predictive confidence interval (scaled)
#'
#' construct the one-sided confidence interval for a unit's survival time T
#'
#' @param x a vector of the covariate of the test data.
#' @param r the censoring time of the test data.
#' @param alpha a number betweeo 0 and 1, specifying the miscaverage rate.
#' @param data a data frame used for calibration, containing four columns: (X,R,event,censored_T). 
#' @param mdl The fitted model to estimate the conditional quantile (default is NULL).
#' @param quant_lo the fitted conditional quantile for the calibration data (default is NULL).
#' @param new_quant_lo the fitted conditional quantile for the test data (default is NULL).
#'
#' @return low_ci a value of the lower bound for the survival time of the test point.
#' @return includeR 0 or 1, indicating if [r,inf) is included in the confidence interval.
#'
#' @family confint
#'
#' @export


scaled_lower_ci <- function(x,r,alpha,
                     data,
                     quant_lo,
                     new_quant_lo,
                     weight,
                     new_weight
                     ){
  ## Add a bit of noise to the score to break the ties
  sigma_noise <- 1e-6
  len_r <- length(r)
  if(is.null(dim(x)[1])){
    len_x <- length(x)
    p <- 1
  }else{
    len_x <- dim(x)[1]
    p <- dim(x)[2]
  }
  if(len_r>1 & len_r!=len_x){
    stop("The length of R is not compatible with that of X!")
  }
  n <- dim(data)[1]
  xnames <- paste0("X",1:p)

  ## The nonconformity scores
  ## calibration
  score  <-  pmin(data$R,quant_lo)-data$censored_T+rnorm(n,0,sigma_noise)
  ## Normalize the score by the interquantile range
  score <- score/weight
  corr_term <- quantile(c(score,Inf),1-alpha,type=1)
  extra_noise <- rnorm(len_x,0,sigma_noise)
  ci_low <- pmin(new_quant_lo,r)-new_weight*corr_term+extra_noise
  ci_low <- pmin(ci_low,r)
  ci_low <- pmax(ci_low,0)
  includeR <- ifelse(pmin(new_quant_lo,r)-r+extra_noise<=new_weight*corr_term,1,0)
  ##   return(list(ci_low=ci_low,includeR = includeR,corr_term = corr_term))
  return(list(ci_low=ci_low,includeR = includeR))
}
