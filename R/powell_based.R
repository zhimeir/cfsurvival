#' Confidence interval based on Powell (1986)
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

pow_based <- function(x,r,alpha,
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
  mdl <- crq(fmla,data=data_fit,method = "Pow")

  ## obtain final confidence interval
  if(type == "marginal"){
    res <- lower_ci(x,
                r=r,
                alpha=alpha,
                data = data_calib,
                quant_lo = quant_lo,
                new_quant_lo = new_quant_lo
    )
    return(res)
  }
    if(type == "local"){
        res <- lower_ci_local(x,
                              r=r,
                              alpha=alpha,
                              data=data_calib,
                              quant_lo = quant_lo,
                              new_quant_lo = new_quant_lo,
                              h=h)
        return(res)
      }
}
