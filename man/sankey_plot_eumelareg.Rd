% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sankey_plot_eumelareg.R
\name{sankey_plot_eumelareg}
\alias{sankey_plot_eumelareg}
\title{Draw Sankey Plot in EuMelaReg style}
\usage{
sankey_plot_eumelareg(data, ...)
}
\arguments{
\item{data}{data.frame containing the different strata with one row for each patient}

\item{...}{Additional arguments for the \code{sankeyNetwork} function}
}
\description{
This function draws a sankey plot to display the directed flow between strata (from left to right).
A dataframe with one row per patient is used as input variable. This dataframe is transformed to a frequency table
and links between the strata are calculated. The frequencies of each flow from one stratum to another are added and
duplicate links are removed. The nodes are named after the variables within the strata combined with
a number indicating the stratum. The sankey plot is then drawn using the \code{networkD3} package. Additional
arguments can be passed on to the \code{sankeyNetwork} function from said package.
}
