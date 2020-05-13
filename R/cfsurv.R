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
  ## Check if the required packages are installed
  ##   list.of.packages <- c("ggplot2", "grf", "quantregForest", "randomForestSRC", "survival")
  ##   new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  ##   if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

  ## Process the input
  len_r <- length(r)
  len_x <- length(x)
  if(len_r>1 & len_r!=len_x){
    stop("The length of R is not compatible with that of X!")
  }

  if(model %in% c("cox","randomforest")==0) stop("The regression method is not supported.")

  ## set random seed
  set.seed(seed)

  ## Split the data into the training set and the calibration set
  n = dim(data)[1]
  n_train = n/2
  n_calib = n-n_train
  I_fit <- sample(1:n,n_train,replace = FALSE)
  data_fit <- data[I_fit,]
  data_calib <- data[-I_fit,] 
  
  ## Run the main function and gather resutls
  if(model == "cox"){
    res = cox_based(x,r,alpha,
                    data_fit,
                    data_calib,
                    type,
                    dist,
                    h)
   }
  if(model == "randomforest"){
    res = rf_based(x,r,alpha,
                   data_fit,
                   data_calib,
                   type,
                   h)
  }
  return(res)

  
}

