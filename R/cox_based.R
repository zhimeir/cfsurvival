#' cox_based
#'
#' Construct conformal prediction model based on cox model
#'
#' @param x the covariate of the prediction point
#' @param r the observing time of the prediction point
#' @param alpha level of confidence. The target is to contruct a level 1-alpha confidence interval
#' @param data_fit traing data
#' @param data_calib calibration data
#' @param type type of confidence interval: marginal or local
#' @param dist distribution of T
#' @export

cox_based <- function(x,r,alpha,
                      data_fit,
                      data_calib,
                      type,
                      dist,
                      h){
  
  len_x <- length(x)
  len_r <- length(r)
  mdl <- survreg(Surv(censored_T,event)~X,data=data_fit,dist=dist)
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
