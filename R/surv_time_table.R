#' Calculate x-months survival and show results in a table
#'
#' This function calculates x-months survival grouped by an optional stratification factor and displays the results in a Table.
#' @param data data.frame or data.table containing survival data.
#' @param time The time interval from start of observation until date of event (e.g. disease progression or death)
#' or censoring.
#' @param status Variable specifying if event occured or data has been censored. Default behaviour inherited from the \code{surv_fit} function in the \code{survival} package, with 0
#' indicating censored data and 1 indicating event.
#' @param var Variable tested for Influence on outcome.
#' @param times Time in months for which survival should be calculated.
#' @export

surv_time_table <- function(data, time, status, times, var, print.html = TRUE){
  ls_total <- lapply(times, function(x){
    tmp <- mapply(survival_time, time = time, status = status,
                  MoreArgs = list(times = x, data = data))
    as.data.frame(t(dplyr::bind_rows(tmp)))
  })
  res_total <- do.call(rbind, ls_total)

  ls_var <- lapply(times, function(x){
    tmp <- mapply(survival_time, time = time, status = status,
                  MoreArgs = list(var = var, times = x, data = data))
    as.data.frame(t(dplyr::bind_rows(tmp)))
  })
  res_var <- do.call(rbind, ls_var)
  colnames(res_var) <- levels(data[[var]])
  res_var$Total <- res_total$V1
  htmlTable::htmlTable(res_var)
}
