#' lower_ci_local
#'
#' A function to compute the local lower bound of the survival time
#' @param mdl the fitted moodel
#' @param x covariate of the prediction point. Can be a point or a vector
#' @param r observe time of the prediction point. Can be point or a vector. When it is a vector, its length should match the length of x.
#' @param data data frame containing the calibration data
#' @param alpha construct a level 1-alpha confidence interval
#' @param h bandwidth (default: 1)
#' @export

lower_ci_local <- function(mdl,x,r,data,alpha,h=1){
  n <- dim(data)[1]
  res <- predict(mdl,
                 newdata = data.frame(X=data$X),
                 type="quantile",
                 p=alpha)
  quant_lo <-  res
  
  ## calibration
  sigma_noise <- 1e-6
  score <-  pmin(quant_lo,data$R)-data$censored_T+rnorm(n,0,sigma_noise)
 
  res <- predict(mdl,
                 newdata = data.frame(X=x),
                 type="quantile",
                 p=alpha)
  new_quant_lo <-  res

  w <- get_weight_r(data,h,r,r)
  score_vec <- c(score,Inf)
  sort_score <- sort(score_vec)
  order_score <- order(score_vec)
  sort_w <- w[order_score]
  idxw <- min(which(cumsum(sort_w)>=1-alpha))
  corr_term <- sort_score[idxw]
  extra_noise <- rnorm(1,0,sigma_noise)
  ci_low <- min(new_quant_lo,r)-corr_term+extra_noise
  ci_low <- pmin(ci_low,r)
  ci_low <- pmax(ci_low,0)
  includeR <- ifelse(min(r,new_quant_lo)-r+extra_noise<=corr_term,1,0)
  return(list(ci_low=ci_low,includeR = includeR,corr_term = corr_term))
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
