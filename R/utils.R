#' Fitting the censoring probability P(C>=c|X)
#'
#' @export

censoring_prob <- function(fit,calib,test=NULL,
                           xnames,c){

  ## Fitting P(-C<=-c_0|X) (since P(C>=c_0|X)=P(-C<=-c_0|X))
  fit$C <- -fit$C
  fmla <- with(fit,as.formula(paste("C ~ ", paste(xnames, collapse= "+"))))
  capture.output(bw <- npcdistbw(fmla),file ='NULL')

  ## Computing censoring scores for the calibration data
  newdata_calib <- calib
  newdata_calib$C <- -c
  pr_calib<- npcdist(bws=bw,newdata = newdata_calib)$condist

  ## Computing the censoring scores for the test data
  if(!is.null(test)){
    newdata <- cbind(test,C=-c)
    newdata <- data.frame(newdata)
    colnames(newdata) <- c(xnames,"C")
    pr_new <- npcdist(bws=bw,newdata=newdata)$condist
  }else{pr_new=NULL}
  return(list(pr_calib=pr_calib,pr_new=pr_new))
}


#' Computing the calibration term with covaraite shift
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


get_calibration <- function(score,weight_calib,weight_new,alpha){
  ## Check input format
  if(length(score)!=length(weight_calib)) stop("The length of score is not compatible with the length of weight!")

  if(!is.numeric(alpha)) stop("alpha should be a real number between 0 and 1!")
  if(alpha>1 | alpha<0) stop("alpha should be a real number between 0 and 1!")

  ## Computing the calibration term
  weight <- c(weight_calib,weight_new)
  weight <- weight/sum(weight)
  score_vec <- c(score,Inf)
  sort_score <- sort(score_vec)
  order_score <- order(score_vec)
  sort_w <- weight[order_score]
  idxw <- min(which(cumsum(sort_w)>=1-alpha))
  calib_term <- sort_score[idxw]

  return(calib_term)
}

 
