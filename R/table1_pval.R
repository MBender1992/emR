#' Print table1 with pvalues
#'
#' This function calculates statistics and prints results in a table generated with [table1()]. For categorical data
#' chi² test is used, for numerical data Welch ANOVA is conducted
#' @param data data.frame or data.table containing patient characteristics and demographics
#' @param strat strata to group by
#' @param vars character vector specifying columns that shall be included in the table
#' @inheritParams table1::table1
#' @inheritParams rstatix::anova_test
#' @param html logical indicating whether output should include html characters or not
#' @param ... additional arguments passed on to [table1()] function.
#' @export

table1_pval <- function(data, strat, vars, footnote = NULL, html = TRUE, white.adjust = TRUE, ...){

  if (any(is.na(data[[strat]]))) {
    warning("There are missing values in the grouping variable. These are excluded from statistical analysis.")
  }
  # calculate pvalue for each variable grouped by strata
  rndr <- function(x, name, ...) {
    if (length(x) == 0) {
      y <- data[[name]]
      ind <- !is.na(y)
      y <- y[ind]
      kw_formula <- as.formula(paste(name, strat, sep = "~"))
      s <- rep("", length(render.default(x=y, name=name, ...)))
      if (is.numeric(y)) {
        p <- suppressWarnings(rstatix::anova_test(data, formula = kw_formula)$p)
      } else {
        p <- chisq.test(table(y, droplevels(data[[strat]][ind])))$p.value
      }
      if(html == TRUE){
        s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
      } else {
        s[2] <- format.pval(p, digits=3, eps=0.001)
      }
      s
    } else {
      render.default(x=x, name=name, ...)
    }
  }

  # add names to new strata
  rndr.strat <- function(label, n, ...) {
    ifelse(n==0, label, render.strat.default(label, n, ...))
  }

  # add column with pvalues
  data[[strat]] <- factor(data[[strat]], levels = c(levels(data[[strat]]),2) , labels = c(levels(data[[strat]]),"P-value"))

  formula <- as.formula(paste("~", paste(vars, collapse = "+")  , "|", strat))
  table1(formula, data, droplevels = F, render=rndr,
         #render.missing = NULL,
         render.strat=rndr.strat, footnote=footnote, ...)

}



