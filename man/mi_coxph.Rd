% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mi_coxph.R
\name{mi_coxph}
\alias{mi_coxph}
\title{Coxph with multiple imputed values}
\usage{
mi_coxph(data, time, status, vars, prop.var = NULL, m = 5, ...)
}
\arguments{
\item{data}{data.frame or data.table containing survival data.}

\item{time}{the time interval from start of observation until date of event (e.g. disease progression or death)
or censoring.}

\item{status}{variable specifying if event occured or data has been censored.}

\item{vars}{variables tested for Influence on outcome.}

\item{prop.var}{variable for which propensity scores should be calculated. If no value is provided (prop.var = NULL), no weights are used in coxph. Default is NULL.}

\item{m}{Number of multiple imputations. The default is \code{m=5}.}

\item{...}{additional arguments to be passed on to coxph function}
}
\description{
This code generates a coxph model with multiple imputed values for missing data. Imputation is conducted with the \code{\link[mice:mice]{mice::mice()}} function.
Pooled results obtained with \code{\link[mice:pool]{mice::pool()}} and respective pvalues and 95\% CIs are presented as results.
}
