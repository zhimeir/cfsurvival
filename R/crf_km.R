#' crf.km
#'
#' A function from https://github.com/AlexanderYogurt/censored_ExtremelyRandomForest
#' to perform censored quantile random forest. It is slightly modified.
#'
#' @export

## source("metrics.R")
## source("help_functions.R")
## source("Csurv.R")
## source("crf.R")
## source("cranger.R")
## source("cgrf.R")
## library(fishmethods)
## 
## list.of.packages <- c("survival", "randomForest", "ranger", "grf")
## new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
## if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
## 
## library(survival)
## library(randomForest)
## library(ranger)
## library(grf)

crf.km <- function(fmla, ntree, 
                   nodesize,
                   data_train,
                   data_test,
                   yname, 
                   iname,
                   tau,
                   xname = NULL, 
                   calibrate_taus = c(0.1, 0.5, 0.9),
                   honesty = TRUE, 
                   method = "grf", 
                   splitrule = "extratrees", 
                   nnb = FALSE, reg.split = FALSE) {
  # build Forest model
  if (method == "randomForest") {
    rf <- randomForest(fmla, data = data_train, ntree = ntree, nodesize=nodesize)
    # get proximity matrix
    proxMtx <- rf.getWeights(rf, data_train, data_test, yname, iname)
  } else if (method == "ranger") {
    if (is.null(xname)) {
      rf <- ranger(fmla, data = data_train, num.trees = ntree, min.node.size=nodesize, splitrule=splitrule)
    } else {
      rf <- ranger(fmla, data = data_train[,c(xname,yname)], num.trees = ntree, min.node.size=nodesize, splitrule=splitrule)
    }
    # get proximity matrix
    proxMtx <- ranger.getWeights(rf, data_train, data_test, yname, iname, xname)
  } else if (method == "grf") {
    if (is.null(xname)) {
      rf <- quantile_forest(X=data_train[ ,!(names(data_train) %in% c(yname, iname)), drop=F],
                            Y=data_train[,yname],
                            quantiles = calibrate_taus, 
                            num.trees = ntree,
                            min.node.size = nodesize,
                            honesty = honesty,
                            regression.splitting = reg.split)
    } else {
      rf <- quantile_forest(X=data_train[ ,xname], Y=data_train[,yname], quantiles = calibrate_taus, 
                            num.trees = ntree, min.node.size = nodesize, honesty = honesty,
                            regression.splitting = reg.split)
    }
    proxMtx <- as.matrix(grf.getWeights(rf, data_test, yname, iname, xname))
  }
  
  # censor forest
  n <- nrow(data_test)
  ntrain <- nrow(data_train)
  Yc <- rep(NA, n) # to store new predictions
  Ytrain <- data_train[[yname]]
  censorInd <- data_train[[iname]]
  right_censoring <- TRUE
  
  # find minimum
  for (r in 1:NROW(data_test)) {
    # C survival estimate
    if (nnb) {
      boot.idx <- proxMtx[r, ] > 1/ntrain
      C.boot <- Ytrain[boot.idx]
      i.boot <- 1 - censorInd[boot.idx]
      C.km <- survfit(Surv(C.boot, i.boot) ~ 1, type = 'kaplan-meier')
      C.surv <- stepfun(C.km$time, c(1, C.km$surv))
    } else {
      boot.idx <- proxMtx[r, ] > 0
    }

    max.uncensored <- max(Ytrain[boot.idx & censorInd==1])
    candidateY <- sort(Ytrain[boot.idx & Ytrain<=max.uncensored], decreasing = FALSE)

    Yc[r] <- candidateY[1]
    min_loss <- 1000000
    denom <- sapply(Ytrain, function(x) {1*(Ytrain >= x)%*%proxMtx[r, ]})
    denom[denom==0] <- 1
    base <- 1 - (proxMtx[r, ]/denom)
    base <- base^(1-censorInd)
    nn <- 1
    for (lambda in candidateY) {
      if (right_censoring) {
        if (nnb) {
          kappa <- C.surv(lambda)
        } else {
          kappa <- C.surv(lambda, Ytrain, base)
          #print(kappa)
        }
        loss1 <- proxMtx[r, ]%*%(1*(Ytrain > lambda))
        #print(kappa)
        loss <- abs((1-tau)*kappa - loss1)
      } else {
        #
      }
      if (loss < min_loss) {
        Yc[r] <- lambda
        min_loss <- loss
      } else if (loss == min_loss) {
        Yc[r] <- (Yc[r]*nn + lambda) / (nn+1)
        nn <- nn + 1
      }
    }
  }
  return(
    list(
      'predicted'=Yc,
      'proxMtx'=proxMtx
    )
  )
}
