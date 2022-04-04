#' Calculate survival fraction at given timepoints
#'
#' This function calculates x-months survival grouped by an optional stratification argument.
#' @param data data.frame or data.table containing survival data.
#' @param time The time interval from start of observation until date of event (e.g. disease progression or death)
#' or censoring.
#' @param status Variable specifying if event occured or data has been censored. Default behaviour inherited from the \code{surv_fit} function in the \code{survival} package, with 0
#' indicating censored data and 1 indicating event.
#' @param var Variable tested for Influence on outcome.
#' @param times Time in months for which survival should be calculated.
#' @export

survival_time <- function(data, time, status, var = NULL, times){
  if(!is.null(var)){
    tmp <- summary(survival::survfit(Surv(eval(parse(text = time)), eval(parse(text = status))) ~ eval(parse(text = var)), data = data), times = times, extend = TRUE)
    res <- sapply(1:length(levels(data[[var]])), function(x){
      paste(round(tmp$surv[x], 3) * 100,"% (", round(tmp$lower[x], 3) * 100,"%-", round(tmp$upper[x], 3) * 100,"%)", sep = "")
    })
    res <- data.frame(res)
    rownames(res) <- levels(dat[[var]])
    colnames(res) <- paste(times, "-months survival", sep = "")
    return(res)
  } else {
    tmp <- summary(survival::survfit(Surv(eval(parse(text = time)), eval(parse(text = status))) ~ 1, data = data), times = times, extend = TRUE)
    res <- paste(round(tmp$surv, 3) * 100,"% (", round(tmp$lower, 3) * 100,"%-", round(tmp$upper, 3) * 100,"%)", sep = "")
    res <- data.frame(res)
    colnames(res) <- paste(times, "-months survival", sep = "")
    rownames(res) <- "Total"
    return(res)
  }
}
