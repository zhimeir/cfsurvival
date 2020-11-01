#' Confidence interval based on Cox model
#'
#' Construct conformal predictive interval based on cox model
#'
#' @param x a vector of the covariate of the test data.
#' @param c the censoring time of the test data.
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

cox_based <- function(x,c,alpha,
                      data_fit,
                      data_calib,
                      type,
                      dist,
                      weight_calib,
                      weight_new){
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
  
  if(type == "quantile"){
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
    ##     score <- pmin(data_calib$C,quant)-data_calib$censored_T
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
    }

  if(type == "percentile"){
    ## Fit the model
    xnames <- paste0("X",1:p)
    fmla <- as.formula(paste("Surv(censored_T, event) ~ ", paste(xnames, collapse= "+")))
    ##     mdl <- coxph(fmla,data=data_fit)
    mdl <- survreg(fmla,data=data_fit,dist="weibull")
    shape <- 1/mdl$scale
    scale <- exp(coef(mdl)[1])
    coef <- coef(mdl)[-1]
    score <- unlist(map2(data_calib$X1,data_calib$censored_T,
                  get_survival_fun,shape,scale,coef))

    ## Extract the fitted survival function
    ##     res <- survfit(mdl,
    ##                    newdata=data_calib,
    ##                    stype=1)
    ##     fit_time <- res$time
    ##     fit_surv <- res$surv
    ##     score <- rep(0,dim(data_calib)[1])
    ##     for(i in 1:dim(data_calib)[1]){
    ##       score[i] <- extract_surv_prob(data_calib$censored_T[i],fit_time,fit_surv[,i])
    ##     }

    #The fitted survival function for the new data
    ##     newdata <- data.frame(x)
    ##     colnames(newdata) <- xnames
    ##     res <- survfit(mdl,
    ##                   newdata=newdata,
    ##                   stype=1)
    ##     test_time <- res$time
    ##     test_surv <- res$surv

    ## Obtain the calibration term
    calib_term <- sapply(X=weight_new,get_calibration,score=score,
                         weight_calib=weight_calib,alpha=alpha)

    ## Obtain the final confidence interval
    lower_bnd <- rep(0,len_x)
    for(i in 1:len_x){
      ## ind <- min(which(test_surv[,i]<=calib_term[i]))
      ## What if the index does not exist?
      lower_bnd[i] <- scale*(-log(calib_term[i])*exp(-coef*data_test$X1[i]))^(1/shape)
      }
    }
 
  lower_bnd <- pmax(lower_bnd,0)
  lower_bnd <- pmin(lower_bnd,c)
  return(lower_bnd)
}


get_survival_fun <- function(x,t,sh,sc,coef){
  s <- -(t/sc)^sh*exp(coef*x)
  s <- exp(s)
  return(s)
}
