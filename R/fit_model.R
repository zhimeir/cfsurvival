#' Fit the model to estimate the quantile of survival time
#'
#' @param data data frame used to fit the quantile regression model. It should contain (X,R,censoredT)
#' @param model the model used to fit the quantile regression. Takes value in ("cox","randomforest"). The default is "cox".
#' @return The fitted model

fit_model <- function(data,
                      alpha=NULL,
                      model = "cox",
                      dist = "weibull"){
  if(model == "cox"){
  mdl <- survreg(Surv(censored_T,event)~X,data=data,dist=dist)
  }

  if(model == "randomforest"){
    ##     Set up parameters
    ntree <- 1000
    nodesize <- 80

    mdl <- crf.km(fmla, ntree = ntree, nodesize = nodesize, data_train = data, data_test = data, 
               yname = 'T', iname = 'event', tau = alpha, method = "grf", calibrate_taus = alpha)
    ##     Yc <- Yc$predicted
    }

  return(mdl)
}

