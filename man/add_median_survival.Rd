% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_median_survival.R
\name{add_median_survival}
\alias{add_median_survival}
\title{Calculate median survival time}
\usage{
add_median_survival(
  data,
  time,
  status,
  var,
  round = 1,
  statistics = TRUE,
  weights = NULL,
  conf.type = "log-log"
)
}
\arguments{
\item{data}{data.frame or data.table containing survival data.}

\item{time}{The time interval from start of observation until date of event (e.g. disease progression or death)
or censoring.}

\item{status}{Variable specifying if event occured or data has been censored. Default behaviour inherited from the \code{surv_fit} function in the \code{survival} package, with 0
indicating censored data and 1 indicating event.}

\item{var}{Variable tested for Influence on outcome.}

\item{round}{rounds the results to the specified number of decimal places (default 1)}

\item{statistics}{Logical value. If TRUE pvalue is printed. Default is TRUE. Statistical test is log-rank test.}

\item{weights}{character variable specifying the name of the weights column. Weights have to be added to the original dataframe in order to be applied correctly.}

\item{conf.type}{Method to calculate confidence intervals. Log-log method is the default in SAS.}
}
\description{
This function fits a survival curve for each subgroup defined by var and extracts the median survival times.
Additionally the median survival time for the whole sample is calculated.
}
