#' Print table with univariate, full mutlivariable and reduced (backwards selection) multivariable cox regression
#'
#' This function combines the results of a univariable, full multivariable model and a backwards selection model.
#' @param data data.frame or data.table containing survival data.
#' @param time the time interval from start of observation until date of event (e.g. disease progression or death)
#' or censoring.
#' @param time2 ending time of the interval for interval censored or counting process data only.
#' Intervals are assumed to be open on the left and closed on the right, (start, end].
#' For counting process data, event indicates whether an event o
#' @param status variable specifying if event occured or data has been censored.
#' @param vars one or more variables defined as character strings to be included in the table
#' @param rgroup a vector of character strings containing headings for row groups.
#' @param p.thres pvalue threshold for backwards selection model.
#' @param footnote character string passed on to \code{tfoot} argument of \code{htmlTable} or \code{add_footer_lines} of \code{flextable}
#' @param fixed.var specifies fixed variables to be included in the cox model.
#' @param output string specifying the type of output. Allowed values are "html" generating a htmlTable, "flex" generating a flextable or "df" returning a data.frame
#' @param ... additional arguments to be passed on to \code{cox_table_combined}
#' @export

cox_table_combined <- function(data, time, time2 = NULL, status, vars, fixed.var = NULL, rgroup = NULL, p.thres = 0.1,
                               footnote = NULL, output = "html", ...){
  if(!is.null(time2)) {
    time1 <- time
    univ <- cox_table_time_dependent_cov(data = data, time1 = time1, time2 = time2, status = status,
                      vars = vars, footnote = NULL, printHTML = FALSE, univariate = TRUE)
    fullmodel <- cox_table_time_dependent_cov(data = data, time1 = time1, time2 = time2, status = status,
                           vars = vars, footnote = NULL, printHTML = FALSE, univariate = FALSE)
    backwardsmodel <- cox_table_time_dependent_cov(data = data, time1 = time1, time2 = time2, status = status,  vars = vars, footnote = NULL, fixed.var = fixed.var,
                                printHTML = FALSE, univariate = FALSE, modeltype = "backwards", p.thres = p.thres)
  } else{
    univ <- cox_table(data = data, time = time, status = status,
                      vars = vars, footnote = NULL, printHTML = FALSE, univariate = TRUE)
    fullmodel <- cox_table(data = data, time = time, status = status,
                           vars = vars, footnote = NULL, printHTML = FALSE, univariate = FALSE)
    backwardsmodel <- cox_table(data = data, time = time, status = status,  vars = vars, footnote = NULL, fixed.var = fixed.var,
                                printHTML = FALSE, univariate = FALSE, modeltype = "backwards", p.thres = p.thres)
  }

  tmp <- dplyr::left_join(tibble::rownames_to_column(univ$res), tibble::rownames_to_column(fullmodel$res), by = "rowname")
  tmp <- textshape::column_to_rownames(tmp, "rowname")
  out <- dplyr::left_join(tibble::rownames_to_column(tmp), tibble::rownames_to_column(backwardsmodel$res), by = "rowname")
  out <- textshape::column_to_rownames(out, "rowname")
  if(is.null(rgroup)){
    rgroup <- fullmodel$rgroup
  }

  if(output == "html"){
    htmlTable(out, rgroup = rgroup, n.rgroup = fullmodel$n.group, cgroup = c("Univariate cox regression", "Multivariate cox regression <br />(full model)",
                                                                             "Multivariate cox regression <br />(backwards selection&#42;)"), n.cgroup = c(2,2,2),
              tfoot = paste(footnote, " &#42;: Threshold for backwards selection was p<",p.thres, ".", sep = ""), ...)
  } else if(output == "flex") {
    tbl <- tibble::rownames_to_column(out, "level")
    ind <- max.col(sapply(rgroup, grepl, tbl$level))
    tbl$groups <- rgroup[ind]
    tbl$level <- ifelse(tbl$level != tbl$groups,  substr(tbl$level, nchar(tbl$groups)+1, nchar(tbl$level)), tbl$level)
    colnames(tbl) <- c("level","HR1", "pval1", "HR2", "pval2", "HR3", "pval3", "groups")
    tbl <- as_grouped_data(tbl, groups = "groups")

    # tbl <- tibble::rownames_to_column(tbl,"var")
    ft <- as_flextable(tbl, hide_grouplabel = TRUE)
    ft <- set_header_labels(ft,level = "", HR1 = "HR (95% CI)", HR2 = "HR (95% CI)", HR3 = "HR (95% CI)",
                            pval1 = "pvalue", pval2 = "pvalue", pval3 = "pvalue")
    ft <- fontsize(ft, size = 10)
    ft <-add_header_row(ft, colwidths = c(1,2, 2,2),
                        values = c("","Univariate cox regression", "Multivariate cox regression (full model)",
                                   "Multivariate cox regression (backwards selection*)"))
    ft <- bold(ft, part = "header")
    ft <- bold(ft, i = ~ !is.na(groups))
    ft <- padding(ft,i = ~ is.na(groups), padding.left=15)
    ft <- add_footer_lines(ft, paste(footnote, "*: Threshold for backwards selection was p<",p.thres, ".", sep = ""))

    # print flextable
    FitFlextableToPage(ft, pgwidth = 7)
  } else if(output == "df"){
    out
  } else {
    stop("Please specify output. Allowed values are html for a htmltable, flex for a flextable or df for a data.frame")
  }
}
