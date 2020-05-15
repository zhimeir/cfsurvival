#' cox_based
#'
#' Construct conformal prediction model based on cox model
#'
#' @param x the covariate of the prediction point.
#' @param r the censoring time of the prediction point.
#' @param alpha Miscoverage rate.
#' @param data_fit traing data.
#' @param data_calib calibration data.
#' @param type type of confidence interval: marginal or local.
#' @param dist The distribution of T used in the cox model.
#'
#' @return ci_low The lower bound of survival time.
#' @return includeR An indicator suggesting if [r,inf) is included in the confidence interval.
#'
#' @export

cox_based <- function(x,r,alpha,
                      data_fit,
                      data_calib,
                      type,
                      dist,
                      h){
  
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
  fmla <- as.formula(paste("Surv(censored_T, event) ~ ", paste(xnames, collapse= "+")))
  mdl <- survreg(fmla,data=data_fit,dist=dist)

  ## obtain final confidence interval
  if(type == "marginal"){
    res <- lower_ci(x,
                mdl=mdl,
                r=r,
                alpha=alpha,
                data = data_calib)
    return(res)
  }
    if(type == "local"){
        res <- lower_ci_local(x,
                              mdl=mdl,
                              r=r,
                              alpha=alpha,
                              data=data_calib,
                              h=h)
        return(res)
      }
}
