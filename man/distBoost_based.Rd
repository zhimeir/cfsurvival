% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distBoost_based.R
\name{distBoost_based}
\alias{distBoost_based}
\title{Confidence interval based on distributional boosting}
\usage{
distBoost_based(
  x,
  c,
  alpha,
  data_fit,
  data_calib,
  weight_calib,
  weight_new,
  n.tree = 100
)
}
\arguments{
\item{x}{a vector of the covariate of the test data.}

\item{c}{the censoring time of the test data.}

\item{alpha}{a number betweeo 0 and 1, specifying the miscaverage rate.}

\item{data_fit}{a data frame, containing the training data.}

\item{data_calib}{a data frame, containing the calibration data.}

\item{type}{either "marginal" or "local". Determines the type of confidence interval.}

\item{dist}{The distribution of T used in the cox model.}
}
\value{
low_ci a value of the lower bound for the survival time of the test point.

includeR 0 or 1, indicating if [r,inf) is included in the confidence interval.
}
\description{
Construct conformal predictive interval based on distributional boosting
}
\seealso{
Other model: 
\code{\link{aft_based}()},
\code{\link{gpr_based}()},
\code{\link{np_based}()},
\code{\link{ph_based}()},
\code{\link{portnoy_based}()},
\code{\link{pow_based}()},
\code{\link{quantBoost_based}()},
\code{\link{rf_based}()}
}
\concept{model}
