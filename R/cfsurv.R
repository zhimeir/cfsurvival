#' cfsurv
#'
#' The main function to generate a predictive conformal confidence interval for a unit's survival time.
#'
#' @param x the covariate for prediction point
#' @param r the censoring time for the prediction point
#' @param data the data frame (X,R,censoredT)
#' @export


# function to construct conformal confidence interval
cfsurv <- function(x,r,data,
                   alpha=0.05,
                   type="marginal",
                   seed = 24601,
                   model = "cox",
                   dist= "weibull",
                   h=NULL){
  ## Process the input
  len_r <- length(r)
  len_x <- length(x)
  if(len_r>1 & len_r!=len_x){
    stop("The length of R is not compatible with that of X!")
  }

  ## set random seed
  set.seed(seed)

  ## Split the data into the training set and the calibration set
  n = dim(data)[1]
  n_train = n/2
  n_calib = n-n_train
  I_fit <- sample(1:n,n_train,replace = FALSE)
  data_fit <- data[I_fit,]
  data_calib <- data[-I_fit,] 

  mdl <- fit_model(data_fit,
                   model=model, 
                   dist = dist)
  
  ## obtain final confidence interval
  if(type == "marginal"){
    if(len_r == 1){
    res <- sapply(x,lower_ci,
                mdl=mdl,
                r=r,
                alpha=alpha,
                data = data_calib)
    return(res)
    }else{
    if(len_r == len_x){
      res <- map2(x,
                  r,
                  lower_ci,
                  mdl=mdl,
                  alpha=alpha,
                  data=data_calib)
    return(res)
    } 
   }
  }
    if(type == "local"){
      if(len_r == 1){
        res <- sapply(x,lower_ci_local,
                      mdl=mdl,
                      r=r,
                      alpha=alpha,
                      data=data_calib,
                      h=h)
        return(res)
      }else{
      res <- map2(x,
                  r,
                  lower_ci_local,
                  mdl=mdl,
                  alpha=alpha,
                  data=data_calib,
                  h=h)
      return(res)
      }
    }
}

