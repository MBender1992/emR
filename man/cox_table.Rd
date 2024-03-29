% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cox_table.R
\name{cox_table}
\alias{cox_table}
\title{Print results from univariate or multivariate cox regression}
\usage{
cox_table(
  data,
  time,
  status,
  vars,
  rgroup = NULL,
  footnote = NULL,
  printHTML = TRUE,
  univariate = TRUE,
  weights = NULL,
  ...
)
}
\arguments{
\item{data}{data.frame or data.table containing survival data.}

\item{time}{the time interval from start of observation until date of event (e.g. disease progression or death)
or censoring.}

\item{status}{variable specifying if event occured or data has been censored.}

\item{vars}{one or more variables defined as character strings to be included in the table}

\item{rgroup}{a vector of character strings containing headings for row groups.}

\item{footnote}{character string passed on to \code{tfoot} argument of \code{htmlTable}}

\item{printHTML}{Logical value. If TRUE output is printed as htmlTable. Default is TRUE.}

\item{univariate}{Logical value. If TRUE output of univariate cox regression is printed. Else output of multivariate
cox regression is printed. Default is TRUE.}

\item{weights}{character variable specifying the name of the weights column. Weights have to be added to the original dataframe in order to be applied correctly.}

\item{...}{additional arguments to be passed on to \code{cox_output}}
}
\description{
This function combines the results from multiple calls of \code{cox_output} and prints the output
as an htmlTable generated with the \code{htmlTable} package or the output of multivariate cox regression.
}
