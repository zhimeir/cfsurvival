#' Predictive confidence interval for survival data
#'
#' The main function to generate a predictive conformal confidence interval for a unit's survival time.
#'
#' @param x a vector of the covariate for test point. 
#' @param c_list the cutoff for the censoring time (a positive number or a vector).
#' @param pr_list the censoring probability at the cutoff for the fitting and calibration datasets.
#' @param pr_new_list the censoring probability at the cutoff for the test points.
#' @param Xtrain a n-by-p matrix of the covariate of the training data.
#' @param C a length n vector of the censoring time of the training data.
#' @param event a length n vector of indicators if the time observed is censored. TRUE corresponds to NOT censored, and FALSE censored.
#' @param time  a vevtor of length n, containing the observed survival time.
#' @param alpha a number between 0 and 1, speciifying the miscoverage rate.
#' @param seed an integer random seed (default: 24601).
#' @param model Options include "aft", "np", "randomforest", "Powell", "Portnoy", "PengHuang", "distBoost", "quantBoost", "gpr". This determines the model used to fit the conformal score (default: "cox").
#' @param dist either "weibull", "exponential" or "gaussian" (default: "weibull"). The distribution of T used in the AFT model. 
#' @param I_fit The set of indices corresponding to the fitting data set. (Default: NULL)
#'
#' @examples
#' # Generate data
#' set.seed(24601)
#' n <- 2000
#' X <- runif(n, 0, 2)
#' T <- exp(X + rnorm(n, 0, 1))
#' C <- rexp(n, rate = 0.05)
#' event <- (T <= C)
#' censored_T <- pmin(T, C)
#' data <- data.frame(X = X, C = C, event = event, censored_T = censored_T)
#'
#' # Prediction point
#' n_test <- 1000
#' X_test <- runif(n_test, 0, 2)
#' T_test <- exp(X_test + rnorm(n,0,1))
#'
#' # Running cfsurvival with c0 = 30
#' c0 <- 30
#' pr_list <- rep(0.5, n)
#' pr_new_list <- rep(0.5, n_test)
#' res <- cfsurv(x = X_test, c_list = c0, pr_list = pr_list, pr_new_list = pr_new_list,
#'              Xtrain = X, C = C, event = event, time = censored_T, 
#'              alpha = 0.1, model = "aft")
#'
#' # Examine the result
#' cat(sprintf("The coverage is %.3f.\n", mean(res <= T_test)))
#'
#' @export

# function to construct conformal confidence interval
cfsurv <- function(x,c_list=NULL,
                   pr_list=NULL,
                   pr_new_list=NULL,
                   Xtrain,C,event,time,
                   alpha=0.05,
                   seed = 24601,
                   model = "aft",
                   dist= "weibull",
                   I_fit = NULL,
                   ftol=.1,tol=.1,
                   n.tree=500
                   ){
  ## Check if the required packages are installed
  ## Solution found from https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
  list.of.packages <- c("quantreg",
                        "grf",
                        "quantregForest",
                        "randomForestSRC",
                        "survival",
                        "tidyverse",
                        "fishmethods",
                        "foreach",
                        "doParallel",
                        "GauPro",
                        "gbm",
                        "np", 
                        "conTree")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
  suppressPackageStartupMessages(res <- lapply(X=list.of.packages,FUN=require,character.only=TRUE))
  ## Process the input
  ## Check the length of x and c: only two cases are supported. length(r)=1, or length(r)=length(x)
  X <- Xtrain
  if(is.null(dim(x)[1])){
    len_x <- length(x)
    p <- 1
  }else{
    len_x <- dim(x)[1]
    p <- dim(x)[2]
  }
  
  
  if(is.null(dim(X)[1])){
    n <- length(X)
    pX <- 1
  }else{
    n <- dim(X)[1]
    pX <- dim(X)[2]
  }
  

  ## Check the type of the model. Only "cox" and "randomforest" are supported
  if(model %in% c("aft","randomforest","pow","portnoy","PengHuang",
                  "distBoost","gpr", "quantBoost")==0) 
    stop("The regression model is not supported.")


  ## Check the value of alpha
  if (alpha>=1 | alpha<=0) stop("The value of alpha is out of bound.")

  ## Check the dimensions of the data 
  xnames <- paste0('X', 1:p)
  if(n != length(C))stop("The number of rows in X does not match the length of R.")
  if(length(C) != length(event))stop("The length of R does not match the length of event.")
  if(length(event) != length(time))stop("The length of event does not match the length of time.")
  if(p != pX) stop("The dimension of the test point does not match the dimension of the training point.")

  data <- as.data.frame(cbind(C,event,time,X))
  colnames(data) <- c("C","event","censored_T",xnames)

  ## set random seed
  set.seed(seed)

  ## Split the data into the training set and the calibration set
  n = dim(data)[1]
  n_train = n/2
  n_calib = n-n_train
  if(is.null(I_fit)){
    I_fit <- sample(1:n,n_train,replace = FALSE)
  }
  data_fit <- data[I_fit,]
  data_calib <- data[-I_fit,]
  newdata <- data.frame(x)
  colnames(newdata) <- xnames
  
  ## If c is not specified, select c automatically 
  if(is.null(c_list)){
    ref_length <- 10
    c_list <- seq(min(data_fit$C),max(data_fit$C),length=ref_length)
  }
  
  if(length(c_list)==1){
    c <- c_list
    if(is.null(pr_list) | is.null(pr_new_list)){
      res <- censoring_prob(fit = data_fit, calib = data_calib,
                            test = newdata, xnames = xnames,
                            c = c, method = "distBoost",n.tree = 500)
      pr_calib <- res$pr_calib
      pr_new <- res$pr_new
    }else{
      pr_calib <- pr_list[-I_fit]
      pr_new <- pr_new_list
    }
  }else{
    if(is.null(pr_list) | is.null(pr_new_list)){
      res <- selection_c(X=data_fit[,colnames(data_fit)%in%xnames],
                         C=data_fit$C,
                         event=data_fit$event,
                         time=data_fit$censored_T,
                         weight_ref=NULL,
                         alpha,c_ref=c_list,
                         type=type,dist=dist)
      c <- res$c_opt
      res <- censoring_prob(data_fit,data_calib,newdata,xnames,c,ftol,tol)
      pr_calib <- res$pr_calib
      pr_new <- res$pr_new
    }else{
      weight_ref <- 1/pr_list[I_fit,]
      res <- selection_c(X=data_fit[,colnames(data_fit)%in%xnames],
                         C=data_fit$C,
                         event=data_fit$event,
                         time=data_fit$censored_T,
                         alpha,c_ref=c_list,
                         weight_ref=weight_ref,
                         type=type,dist=dist)
      c <- res$c_opt
      pr_calib <- pr_list[-I_fit,c_list==c] 
      pr_new <- pr_new_list[,c_list==c]
    }
  }
  ## Computing the weight for the calibration data and the test data
  weight_calib <- 1/pr_calib
  weight_new <- 1/pr_new
 
  ## Run the main function and gather resutls
  if(model == "distBoost"){
    res = distBoost_based(x,c,alpha,
                    data_fit,
                    data_calib,
                    weight_calib,
                    weight_new,
                    n.tree)
   }

  if(model == "quantBoost"){
    res = quantBoost_based(x,c,alpha,
                    data_fit,
                    data_calib,
                    weight_calib,
                    weight_new,
                    n.tree) 
  }

  if(model == "gpr"){
    res = gpr_based(x,c,alpha,
                    data_fit,
                    data_calib,
                    weight_calib,
                    weight_new)
   }
  
  if(model == "aft"){
    res = aft_based(x,c,alpha,
                    data_fit,
                    data_calib,
                    dist,
                    weight_calib,
                    weight_new)
   }
  
  if(model == "np"){
    res = np_based(x,c,alpha,
                    data_fit,
                    data_calib,
                    weight_calib,
                    weight_new,
                    ftol,
                    tol)
   }
  if(model == "randomforest"){
    res = rf_based(x,c,alpha,
                   data_fit,
                   data_calib,
                   weight_calib,
                   weight_new)
  }
  
  if(model == "pow"){
    res = pow_based(x,c,alpha,
                   data_fit,
                   data_calib,
                   weight_calib,
                   weight_new)
  }
  if(model == "portnoy"){
    res = portnoy_based(x,c,alpha,
                   data_fit,
                   data_calib,
                   weight_calib,
                   weight_new)
  }
  if(model == "PengHuang"){
    res = ph_based(x,c,alpha,
                   data_fit,
                   data_calib,
                   weight_calib,
                   weight_new)
  }
  return(res)


}

