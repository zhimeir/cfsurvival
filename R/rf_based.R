#' rf_based
#'
#' Construct a conformal prediction confidence interval based on random forest
#'
#' @param x covariate of the prediction point
#' @export

rf_based <- function(x,r,
                     alpha,
                     data_fit,
                     data_calib,
                     type,
                     h){
  ## Parameters
  n_calib <- dim(data_calib)[1]
  if(is.null(dim(x)[1])){
    len_x <- length(x)
  }else{
    len_x <- dim(x)[1]
  }
  len_r <- length(r)
  p <- dim(x)[2]
  ntree <- 1000
  nodesize <- 80
  xnames <- paste0("X",1:p)
  data_test <- rbind(data_calib$X,x)
  colnames(data_test) <- xnames
  
  fmla <- as.formula(paste("censored_T ~ ",paste(xnames,collapse="+")))
  ## Fit the model
  mdl <- crf.km(fmla, ntree = ntree, 
                 nodesize = nodesize,
                 data_train = data_fit[,names(data_fit)%in%
                                       c(xnames,"censored_T","event")], 
                 data_test = data_test, 
                 yname = 'censored_T', 
                 iname = 'event',
                 tau = alpha,
                 method = "grf")
  quant_lo <- mdl$predicted[1:n_calib]
  new_quant_lo <- tail(mdl$predicted,-n_calib)

  ## obtain final confidence interval
  if(type == "marginal"){
    res <- lower_ci(x,
                  r=r,
                  alpha=alpha,
                  data=data_calib,
                  quant_lo =quant_lo,
                  new_quant_lo=new_quant_lo)
    return(res)
    } 
    if(type == "local"){
      ## When there is only one value for r
      if(len_r==1){
        res <- lower_ci_local(x=x,
                              r=r,
                              alpha=alpha,
                              data=data_calib,
                              quant_lo =quant_lo,
                              new_quant_lo=new_quant_lo,
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
