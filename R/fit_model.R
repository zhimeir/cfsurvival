#' Fit the model to estimate the quantile of survival time
#'
#' @param data data frame used to fit the quantile regression model. It should contain (X,R,censoredT)
#' @param model the model used to fit the quantile regression. Takes value in ("cox","randomforest"). The default is "cox".
#' @return The fitted model

fit_model <- function(data,
                      model = "cox",
                      dist = "weibull"){
  if(model == "cox"){
  mdl <- survreg(Surv(censored_T,event)~X,data=data,dist=dist)
  }

  return(mdl)
}

