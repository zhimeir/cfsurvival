#' lower_score
#'
#' @param model The fitted model
#' @param data calibration data frame
#' @param alpha 1-alpha is the level of the confidence interval 
#' @export


lower_score <- function(mdl,data,alpha){
  res <- predict(mdl,
                 newdata = data.frame(X=data$X),
                 type="quantile",
                 p=alpha)
  quant_lo <-  res
  
  ## calibration
  score  <-  pmin(quant_lo,data$R) - data$censored_T
  return(list(score = score))
}

