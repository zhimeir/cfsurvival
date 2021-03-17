#' Confidence interval based on AFT model
#'
#' Construct conformal predictive interval based on the accelareted failure-time (AFT) model
#'
#' @param x a vector of the covariate of the test data.
#' @param c the censoring time cutoff.
#' @param alpha a number betweeo 0 and 1, specifying the miscaverage rate.
#' @param data_fit a data frame, containing the training data.
#' @param data_calib a data frame, containing the calibration data.
#' @param dist The distribution of T used in the aft model (default: weibull).
#' @param weight_calib The weight corresponding to the calibration data.
#' @param weight_new The weight corresponding to the test data.
#'
#' @return low_ci a value of the lower bound for the survival time of the test point.
#' @return includeR 0 or 1, indicating if [r,inf) is included in the confidence interval.
#'
#' @family model
#'
#' @export

aft_based <- function(x,c,alpha,
                      data_fit,
                      data_calib,
                      dist = "weibull",
                      weight_calib,
                      weight_new
                      ){
  ## Check the dimensionality of the input
  if(is.null(dim(x)[1])){
    len_x <- length(x)
    p <- 1
  }else{
    len_x <- dim(x)[1]
    p <- dim(x)[2]
  }

  ## Keep only the data points with C>=c
  ## Transform min(T,C) to min(T,c) 
  weight_calib <- weight_calib[data_calib$C>=c]
  data_calib <- data_calib[data_calib$C>=c,]
  data_calib$censored_T <- pmin(data_calib$censored_T,c)
  
  ## Fit the survival model
  xnames <- paste0("X",1:p)
  fmla <- as.formula(paste("Surv(censored_T, event) ~ ", paste(xnames, collapse= "+")))
  mdl <- survreg(fmla,data=data_fit,dist=dist)

  ## The fitted quantile for the calibration data
  res <- predict(mdl,
                newdata = data_calib,
                type="quantile",
                p=alpha)
  quant <-  res  
  score <- pmin(c,quant)-data_calib$censored_T
  
  ## The fitted quantile for the new data
  newdata <- data.frame(x)
  colnames(newdata) <- xnames
  res <- predict(mdl,
                  newdata = newdata,
                  type="quantile",
                  p=alpha)
  new_quant <-  res

  ## Compute the calibration term
  calib_term <- sapply(X=weight_new,get_calibration,score=score,
                        weight_calib=weight_calib,alpha=alpha)
  ## obtain final confidence interval
  lower_bnd <- pmin(new_quant,c)-calib_term
 
  lower_bnd <- pmax(lower_bnd,0)
  lower_bnd <- pmin(lower_bnd,c)
  return(lower_bnd)
}


