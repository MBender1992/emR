#' Generate outcome table
#'
#' Wrapper around \code{add_median_survival} to transform the data to htmltable format in order to be used in display
#' with \code{htmlTable::htmlTable}. If time and status of
#' @param data data.frame or data.table
#' @param time Time for survival
#' @param status censor variable
#' @param surv_names names of the survival to be displayed in the table (e.g "Median PFS (95% CI)")
#' @param var Variable tested for Influence on outcome
#' @param bestres column containing data for best response
#' @param ORR column containing data for over all response rate
#' @param DCR column containing data for disease control rate
#' @param weights character variable specifying the name of the weights column. Weights have to be added to the original dataframe in order to be applied correctly.
#' @param footnote add footnote
#' @param font font style for the table
#' @param statistics Logical value. If TRUE pvalue is printed. Default is TRUE. Default test statistics are wilcoxon (or anova if n > 2)
#' @param ... add additional css styling arguments to addHtmlTableStyle from the \code{htmlTable} package
#' for numerical data and fisher exact test for categorical data.
#' @export

outcome_table_survival <- function(data, time, status, surv_names, var, bestres = NULL, weights = NULL,
                                   ORR = NULL, DCR = NULL,statistics = TRUE, footnote = NULL, font = "calibri", ...){

  input <- data.frame(time = time,
                      status = status,
                      rownames = surv_names)

  tmp <- mapply(add_median_survival, time = input[,1], status = input[,2], MoreArgs = list(data = data, var = var, statistics = statistics, weights = weights))
  med_surv <- t(rbind.data.frame(tmp))
  if(statistics == TRUE){
    colnames(med_surv) <- c(sort(as.character(unique(data[[var]]))), "Total", "pvalue")
  } else {
    colnames(med_surv) <- c(sort(as.character(unique(data[[var]]))), "Total")
  }
  rownames(med_surv) <- input$rownames

  table_data <- list()
  if (!is.null(bestres)) table_data[["Best response"]] <- get_stats(data = data, strat = var, outcome = bestres, statistics = statistics)
  if (!is.null(ORR)) table_data[["ORR"]] <- get_stats(data = data, strat = var, outcome = ORR, statistics = statistics)
  if (!is.null(DCR)) table_data[["DCR"]] <- get_stats(data = data, strat = var, outcome = DCR, statistics = statistics)
  if (!is.null(med_surv)) table_data[["Survival"]] <- med_surv

  rgroup <- c()
  n.rgroup <- c()
  output_data <- NULL
  for (varlabel in names(table_data)){
    output_data <- rbind(output_data,
                         table_data[[varlabel]])
    rgroup <- c(rgroup,
                varlabel)
    n.rgroup <- c(n.rgroup,
                  nrow(table_data[[varlabel]]))
  }

  tmp <-sapply(colnames(output_data), gsub, pattern = "No. ", replacement = "(N = ")
  tmp <- gsub(",",".",tmp)
  if (statistics == TRUE) {
    colnames(output_data) <- c(sapply(tmp[-length(tmp)], paste, ")", sep = ""), "P-value")
  } else {
    colnames(output_data) <- sapply(tmp, paste, ")", sep = "")
  }

  output_data_style <- addHtmlTableStyle(output_data, ...)
  htmlTable(output_data_style, align="cccc",
            rgroup=rgroup, n.rgroup=n.rgroup,
            rgroupCSSseparator="",
            rowlabel="",
            tfoot=footnote,
            ctable= TRUE
            )
}



