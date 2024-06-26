% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ate_weights.R
\name{ate_weights}
\alias{ate_weights}
\title{Calculate average treatment effect using inverse propensity score weighting}
\usage{
ate_weights(data, vars, prop.var)
}
\arguments{
\item{data}{Character vector specifying rownames of the table (empty columns should be named with "").}

\item{vars}{variables tested for Influence on outcome. NAs within vars should be replaced with a pseudocategory, e.g. "missing". A more sophisticated
approach with multiple imputation of missing values followed by propensity score calculation is provided with the \code{mi_PS} function of this package.}

\item{prop.var}{variable for which propensity scores should be calculated}
}
\description{
Propensity scores (PS) are calculated using fitted values obtained in a logistic regression. The inverse propensity score weighting is subsequently conducted
by 1/PS for the factor level that corresponds to 1 and by 1/(1-PS) for the factor level that corresponds to 0.
}
\examples{
# example weights calculation with the lung dataset
df <- survival::lung
df$sex <- factor(df$sex)
df$ph.ecog <- ifelse(is.na(df$ph.ecog), "Missing", df$ph.ecog)
df$ph.karno <- ifelse(is.na(df$ph.karno), "Missing", df$ph.karno)
df$weights.ate <- ate_weights(data = df, vars = c("age", "ph.ecog", "ph.karno"), prop.var = "sex")
}
