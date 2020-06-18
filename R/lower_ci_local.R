#' One-sided local predictive confidence interval 
#'
#' Construct a local one-sided confidence interval for a unit's survival time T
#'
#' @param x a vector of the ovariate of the test data. 
#' @param r the censoring time of the test data.
#' @param data a data frame used for calibration, containing four columns: (X,R,evebt,censored_T)
#' @param alpha a number between 0 and 1, specifying the miscoverage rate.
#' @param mdl The fitted model to estimate the conditional quantile (default is NULL).
#' @param quant_lo the fitted conditional quantile for the calibration data (default is NULL).
#' @param new_quant_lo the fitted conditional quantile for the test data (default is NULL).
#' @param h a number specifying thebandwidth (default: 1).
#'
#' @return low_ci a value of the lower bound for the survival time of the test point.
#' @return includeR 0 or 1, indicating if [r,inf) is included in the confidence interval.
#'
#' @family confint
#'
#' @export

lower_ci_local <- function(x,r,alpha,
                           data,
                           h,
                           mdl=NULL,
                           quant_lo=NULL,
                           new_quant_lo=NULL){
  n <- dim(data)[1]
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
  xnames <- paste0("X",1:p)
  data_test <- data.frame(x)
  colnames(data_test) <- xnames

  if(is.null(quant_lo)){
  res <- predict(mdl,
                 newdata = data,
                 type="quantile",
                 p=alpha)
  quant_lo <-  res
  }
  
  ## calibration
  sigma_noise <- 1e-6
  score <-  pmin(quant_lo,data$R)-data$censored_T+rnorm(n,0,sigma_noise)

  if(is.null(new_quant_lo)){
  res <- predict(mdl,
                 newdata = data_test,
                 type="quantile",
                 p=alpha)
  new_quant_lo <-  res
  }
  
  score_vec <- c(score,Inf)
  sort_score <- sort(score_vec)
  order_score <- order(score_vec)
  corr_term <- c()
  for(ind.r in 1:len_r){
    w <- get_weight_r(data,h,r[ind.r],r[ind.r])
    sort_w <- w[order_score]
    idxw <- min(which(cumsum(sort_w)>=1-alpha))
    corr_term <- c(corr_term,sort_score[idxw])
  }
  extra_noise <- rnorm(len_x,0,sigma_noise)
  ci_low <- pmin(new_quant_lo,r)-corr_term+extra_noise
  ci_low <- pmin(ci_low,r)
  ci_low <- pmax(ci_low,0)
  includeR <- ifelse(pmin(r,new_quant_lo)-r+extra_noise<=corr_term,1,0)
  ##   return(list(ci_low=ci_low,includeR = includeR,corr_term = corr_term))
  return(list(ci_low=ci_low,includeR = includeR))
}

get_weight_r <- function(data,h,r0,r){
  n <- dim(data)[1]
  R <- data$R
  w <- rep(0,n+1)
	w[1:n] <- exp(-(R-r0)^2/2/h^2)
	w[n+1] <- exp(-(r-r0)^2/2/h^2)
  w <- w/sum(w)
  return(w=w)
}
