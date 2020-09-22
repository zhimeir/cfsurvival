#' Confidence interval based on Cox model and an Adaptive Score
#'
#' Construct conformal predictive interval based on cox model
#'
#' @param x a vector of the covariate of the test data.
#' @param r the censoring time of the test data.
#' @param alpha a number betweeo 0 and 1, specifying the miscaverage rate.
#' @param data_fit a data frame, containing the training data.
#' @param data_calib a data frame, containing the calibration data.
#' @param type either "marginal" or "local". Determines the type of confidence interval.
#' @param dist The distribution of T used in the cox model.
#'
#' @return low_ci a value of the lower bound for the survival time of the test point.
#' @return includeR 0 or 1, indicating if [r,inf) is included in the confidence interval.
#'
#' @family model
#'
#' @export

adaptive_cox_based <- function(x,r,alpha,
                      data_fit,
                      data_calib){
  ## Check the format of the input data 
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

  ## Keep only the data with R>=r and transform min(T,R) to min(T,r)
  data_fit <- data_fit[data_fit$R>=r,]
  data_fit$censored_T <- pmin(data_fit$censored_T,r)
  data_calib <- data_calib[data_calib$R>=r,]
  data_calib$censored_T <- pmin(data_calib$censored_T,r)

  ## Fit the model
  xnames <- paste0("X",1:p)
  fmla <- as.formula(paste("Surv(censored_T, event) ~ ", paste(xnames, collapse= "+")))
  mdl <- coxph(fmla,data=data_fit)

  #Extract the fitted survival function
  res <- survfit(mdl,
                 newdata=data_calib,
                 stype=1)
  fit_time <- res$time
  fit_surv <- res$surv
  score <- rep(0,dim(data_calib)[1])
  for(i in 1:dim(data_calib)[1]){
    score[i] <- extract_surv_prob(data_calib$censored_T[i],fit_time,fit_surv[,i])
  }
 
  #The fitted survival function for the new data
  newdata <- data.frame(x)
  colnames(newdata) <- xnames
  res <- survfit(mdl,
                 newdata=newdata,
                 stype=1)
  test_time <- res$time
  test_surv <- res$surv

  ## Compute the weight 
  ## Estimating the CDF of the censoring mechanism P(C>=c0)
  ##   browser()
  fmla <- with(data_fit,as.formula(paste("R ~ ", paste(xnames, collapse= "+"))))
  bw <- npcdistbw(fmla)
  cdf_calib<- npcdist(bws=bw,newdata = data_calib)
  newdata <- cbind(newdata,R = r)
  cdf_new <- npcdist(bws=bw,newdata=newdata)
  weight_calib <- cdf_calib$condist
  weight_new <- cdf_new$condist
  weight_calib <-  1/(1-weight_calib)
  weight_new <-  1/(1-weight_new)

  ## Obtain the calibration term
  calib_term <- sapply(X=weight_new,get_calibration,score=score,
                       weight_calib=weight_calib,alpha=alpha)

  ## Obtain the final confidence interval
  lower_bnd <- rep(0,len_x)
  for(i in 1:len_x){
    ind <- min(which(test_surv[,i]<=calib_term[i]))
    ## What if  the  index does not exist?
    lower_bnd[i] <- test_time[ind]
  }
  return(lower_bnd)

}

## Function to extract P(T>t) [to be changed to P(T>=t)] 
extract_surv_prob <- function(t,time,surv){
  if(t<min(time)){
    surv_prob <- 1
  }else{
    tau <- max(which(t>=time))
    surv_prob <- surv[tau]
  }
  return(surv_prob)
}

