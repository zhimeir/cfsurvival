#' Fitting the censoring probability P(C>=c|X)
#'
#' @keywords internal

censoring_prob <- function(data_fit,
                           data_calib,
                           data_test,
                           xnames,c){
  fmla <- with(data_fit,as.formula(paste("C ~ ", paste(xnames, collapse= "+"))))
  bw <- npcdistbw(fmla)
  newdata_calib <- data_calib
  newdata_calib$C <- c
  pr_calib<- 1-npcdist(bws=bw,newdata = newdata_calib)$condist
  newdata <- cbind(data_test,C = c)
  newdata <- data.frame(newdata)
  colnames(newdata) <- c(xnames,"C")
  pr_new <- 1-npcdist(bws=bw,newdata=newdata)$condist

  return(list(pr_calib=pr_calib,pr_new=pr_new))
}

#' Function to extract P(T>t) [to be changed to P(T>=t)] 
#'
#' @keywords internal
extract_surv_prob <- function(t,time,surv){
  if(t<min(time)){
    surv_prob <- 1
  }else{
    tau <- max(which(t>=time))
    surv_prob <- surv[tau]
  }
  return(surv_prob)
}


 
