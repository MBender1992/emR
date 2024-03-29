% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ROCplot_eumelareg.R
\name{ROCplot_eumelareg}
\alias{ROCplot_eumelareg}
\title{Plot ROC curve with youden index and auc.}
\usage{
ROCplot_eumelareg(
  predictions,
  labels,
  print.auc = TRUE,
  print.cutoff = TRUE,
  cutoff.method = c("youden", "closest.topleft"),
  auc.label.x = 0.7,
  auc.label.y = 0.1,
  cutoff.label.x = 0.15,
  cutoff.label.y = 0.9,
  textsize = 12
)
}
\arguments{
\item{predictions}{vector containing predicted probabilities obtained from a prediction model.}

\item{labels}{observed labels of two outcome possibilities. Preferred are binomial values with 0 = control and 1 = case.}

\item{print.auc}{logical value to decide whether auc should be printed within plot. Default is TRUE.}

\item{print.cutoff}{logical value to decide whether cutoff should be printed within plot. Default is TRUE.}

\item{cutoff.method}{method to be used to calculate cutoff. Allowed values are "youden" and "closest.topleft".}

\item{auc.label.x}{x-coordinate of auc label.}

\item{auc.label.y}{y-coordinate of auc label.}

\item{cutoff.label.x}{x-coordinate of cutoff label.}

\item{cutoff.label.y}{y-coordinate of cutoff label}

\item{textsize}{argument to define size of axes title, axes text and plot text.}
}
\description{
This function plots a default EuMelaReg ROC curve using the package \code{\link[=plotROC]{plotROC()}}. By default Youden-Index and AUC are printed and displayed in the graph.
Additional ggplot layers can be added to the graph.
}
