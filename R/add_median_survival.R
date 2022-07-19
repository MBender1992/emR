#' Calculate median survival time
#'
#' This function fits a survival curve for each subgroup defined by var and extracts the median survival times.
#' Additionally the median survival time for the whole sample is calculated.
#' @param data data.frame or data.table containing survival data.
#' @param time The time interval from start of observation until date of event (e.g. disease progression or death)
#' or censoring.
#' @param statistics Logical value. If TRUE pvalue is printed. Default is TRUE. Statistical test is log-rank test.
#' @param status Variable specifying if event occured or data has been censored. Default behaviour inherited from the \code{surv_fit} function in the \code{survival} package, with 0
#' indicating censored data and 1 indicating event.
#' @param var Variable tested for Influence on outcome.
#' @param round rounds the results to the specified number of decimal places (default 1)
#' @param weights character variable specifying the name of the weights column. Weights have to be added to the original dataframe in order to be applied correctly.
#' @export


add_median_survival <- function(data, time, status, var, round = 1, statistics = TRUE, weights = NULL){

  if(!is.null(weights)){
    if(length(levels(data[[var]])) >2){
      stop("IPS weighted pvalues can only be calculated for factors with exactly 2 levels.")
    }
    weights <-  data[[weights]]
  }

  fit <- surv_fit(Surv(eval(parse(text = time)), eval(parse(text = status))) ~ eval(parse(text = var)), data = data, weights = weights)
  surv_med <- surv_median(fit)

  if(!is.null(weights)){
    dat_logrank <- data[!is.na(data[[time]])]
    pval <- logrank_IPSW_RISCA(dat_logrank[[time]], dat_logrank[[status]], ifelse(dat_logrank[[var]] == levels(dat_logrank[[var]])[1], 0, 1), weights)$p.value
  } else {
    pval <- surv_pvalue(fit)$pval
  }

  pval <- ifelse(pval < 0.0001, "< 0.0001", round(pval,3))
  tbl <- data.frame(sapply(1:length(surv_med$median),function(x){
    paste(round(surv_med$median[x],round), " (", round(surv_med$lower[x],round),"-", round(surv_med$upper[x],round),")", sep = "")
  }))

  fit <- surv_fit(Surv(eval(parse(text = time)), eval(parse(text = status))) ~ 1, data = data, weights = weights)
  surv_med <- surv_median(fit)
  tmp <- data.frame(sapply(1:length(surv_med$median),function(x){
    paste(round(surv_med$median[x],round), " (", round(surv_med$lower[x],round),"-", round(surv_med$upper[x],round),")", sep = "")
  }))

  if (statistics == TRUE){
    res <-rbind(tbl, tmp, pval)
    rownames(res) <- c(sort(as.character(unique(data[[var]]))), "Total", "pvalue")
  } else {
    res <-rbind(tbl, tmp)
    rownames(res) <- c(sort(as.character(unique(data[[var]]))), "Total")
  }

  colnames(res) <- "Median (95% CI)"
  return(res)
}
