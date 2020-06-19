#' Confidence interval based on Portnoy (2003)
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

portnoy_based <- function(x,r,alpha,
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
  mdl <- crq(fmla,data=data_fit,method = "Portnoy")
  mdl_coef <- coef(mdl,taus = alpha)

  ## A function to get result
  extract_res <- function(x,mdl_coef){
    res <- mdl_coef[1] + x%*%mdl_coef[-1]
    return(res)
  }
  ## Extract scores
  if(p==1){
    quant_lo <- sapply(data_calib[,colnames(data_calib)%in%xnames],extract_res,mdl_coef = mdl_coef)%>%t()
    new_quant_lo <- sapply(x,extract_res,mdl_coef=mdl_coef)%>%t()
  }else{
    quant_lo <- apply(data_calib[,colnames(data_calib)%in%xnames],1,extract_res,mdl_coef = mdl_coef)%>%t()
    new_quant_lo <- apply(x,1,extract_res,mdl_coef=mdl_coef)%>%t()
  }
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
      if(len_x == len_r){
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


## A function to get result
extract_res <- function(x){
  res <- mdl_coef[1,] + x%*%mdl_coef[-1,]
  return(res)
}
