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
  fmla <- as.formula(paste("Curv(censored_T, R,ctype=\"right\") ~ ", paste(xnames, collapse= "+")))
  mdl <- crq(fmla,data=data_fit,method = "Pow",taus = alpha)
  res <- predict(mdl,
        newdata = data_calib,
        type = "quantile")
  quant_lo <- res

  newdata <- data.frame(x)
  colnames(newdata) <- xnames
  res <- predict(mdl,
        newdata = newdata,
        type = "quantile")
  new_quant_lo <- res
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
      if(len_r == 1){
        res <- lower_ci_local(x,
                              r=r,
                              alpha=alpha,
                              data=data_calib,
                              quant_lo = quant_lo,
                              new_quant_lo = new_quant_lo,
                              h=h)
        return(res)
      }
      ## When there are multiple pairs of (x,r)
      if(len_x == len_r & len_r>1){
        res <- pmap(list(x=x,r=r,new_quant_lo = new_quant_lo),
                    lower_ci_local,
                    alpha=alpha,
                    quant_lo = quant_lo,
                    data=data_calib,
                    h=h)
        return(res)
      }

      }
}
