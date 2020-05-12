#' lower_ci
#' construct the one-sided confidence interval for T
#'
#' @param mdl The fitted model
#' @param x The covariate of the prediction point
#' @param r The observe time of the prediction point
#' @param alpha level
#' @param data data frame
#' @export


lower_ci <- function(mdl,x,r,alpha,data){
  ## Add a bit of noise to the score to break the ties
  sigma_noise <- 1e-6
  n <- dim(data)[1]

  ## The nonconformty scores
  res <- predict(mdl,
                 newdata = data.frame(X=data$X),
                 type="quantile",
                 p=alpha)
  quant_lo <-  res  
  ## calibration
  score  <-  pmin(data$R,quant_lo)-data$censored_T+rnorm(n,0,sigma_noise)

 
  res <- predict(mdl,
                 newdata = data.frame(X=x),
                 type="quantile",
                 p=alpha)
  new_quant_lo <-  res
  corr_term <- quantile(c(score,Inf),1-alpha,type=1)
  extra_noise <- rnorm(1,0,sigma_noise)
  ci_low <- min(new_quant_lo,r)-corr_term+extra_noise
  ci_low <- pmin(ci_low,r)
  ci_low <- pmax(ci_low,0)
  includeR <- ifelse(min(new_quant_lo,r)-r+extra_noise<=corr_term,1,0)
  ##   return(list(ci_low=ci_low,includeR = includeR,corr_term = corr_term))
  return(list(ci_low=ci_low,includeR = includeR))
}
