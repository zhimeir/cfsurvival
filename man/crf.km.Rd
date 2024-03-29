% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crf_km.R
\name{crf.km}
\alias{crf.km}
\title{crf.km}
\usage{
crf.km(
  fmla,
  ntree,
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
  nnb = FALSE,
  reg.split = FALSE
)
}
\description{
The function is based on the scripts from https://github.com/AlexanderYogurt/censored_ExtremelyRandomForest to perform censored quantile random forest.
}
