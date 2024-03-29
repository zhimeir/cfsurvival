% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{censoring_prob}
\alias{censoring_prob}
\title{Fitting the censoring probability P(C>=c|X)}
\usage{
censoring_prob(
  fit,
  calib,
  test = NULL,
  method = "distBoost",
  xnames,
  c,
  ftol = 0.1,
  tol = 0.1,
  n.tree = 40
)
}
\description{
Fitting the censoring probability P(C>=c|X)
}
