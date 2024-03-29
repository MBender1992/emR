% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_stats.R
\name{get_stats}
\alias{get_stats}
\title{Get description stats for survival analysis}
\usage{
get_stats(
  data,
  strat,
  outcome,
  digits = 1,
  add_total_col = "last",
  show_all_values = TRUE,
  hrzl_prop = FALSE,
  statistics = FALSE,
  html = TRUE,
  header_count = TRUE,
  continuous_fn = Gmisc::describeMedian,
  ...
)
}
\arguments{
\item{data}{data.frame or data.table}

\item{strat}{variable used for stratification.}

\item{outcome}{outcome variable to be summarized.}

\item{digits}{the number of decimals used.}

\item{add_total_col}{This adds a total column to the resulting table.
You can also specify if you want the total column "first" or "last"
in the column order.}

\item{show_all_values}{Show all values in proportions. For factors with only two values
it is most sane to only show one option as the other one will just be a complement
to the first, i.e. we want to convey a proportion. For instance sex - if you know
gender then automatically you know the distribution of the other sex as it's 100 \% - other \%.
To choose which one you want to show then set the \code{default_ref} parameter.}

\item{hrzl_prop}{This is default FALSE and indicates
that the proportions are to be interpreted in a vertical manner.
If we want the data to be horizontal, i.e. the total should be shown
and then how these differ in the different groups then set this to TRUE.}

\item{statistics}{Add statistics, fisher test for proportions and Wilcoxon
for continuous variables. See details below for more customization.}

\item{html}{If HTML compatible output should be used. If \code{FALSE}
it outputs LaTeX formatting}

\item{header_count}{Set to \code{TRUE} if you want to add a header count,
e.g. Smoking; No. 25 observations, where there is a new line after the
factor name. If you want a different text for the second line you can
specifically use the \code{\link[base]{sprintf}} formatting, e.g. "No. \%s patients".}

\item{continuous_fn}{The method to describe continuous variables. The
default is \code{\link[Gmisc]{describeMean}}.}

\item{...}{additional arguments passed on to \code{getDescriptionStatsBy} function from \code{Gmisc} package}
}
\description{
Calculate descriptive statistics grouped by strata defined by the variable \code{strata} for
different outcome variables.
}
