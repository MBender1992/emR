#' Calculate median survival time
#'

#' Fitting a survival fit for each subgroup defined by var and extracting median survival times.
#' Additionally the total survival time for the whole sample is calculated.

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
#' @param weights character variable specifying the name of the weights column. Weights (numerical vector) have to be added to the original dataframe in order to be applied correctly.
#' @param conf.type Method to calculate confidence intervals. Log-log method is the default in SAS and was therefore set as the default here too.
#' @export
#' @examples
#' # example survival calculation with the lung dataset
#' df <- survival::lung
#' add_median_survival(data = df, time = "time", status = "status", var = "ph.ecog")
#' # calculate median survival with inverse propensity score weighting for sex
#' df$sex <- factor(df$sex)
#' df$ph.ecog <- ifelse(is.na(df$ph.ecog), "Missing", df$ph.ecog)
#' df$ph.karno <- ifelse(is.na(df$ph.karno), "Missing", df$ph.karno)
#' df$weights.ate <- ate_weights(data = df, vars = c("age", "ph.ecog", "ph.karno"), prop.var = "sex")
#' # the internal logrank test (logrank_IPSW_RISCA) needs a status variable of 0 and 1
#' df$status <- ifelse(df$status == 2, 1, 0)
#' add_median_survival(data = df, time = "time", status = "status",
#' var = "sex", weights = "weights.ate")

add_median_survival <- function(data, time, status, var, round = 1, statistics = TRUE, weights = NULL, conf.type = "log-log"){

  if(!is.null(weights)){
     weights <-  data[[weights]]
  }

  surv_formula <- as.formula(paste("Surv(", time, ",", status,")~",var, sep = ""))
  fit <- surv_fit(surv_formula, data = data, weights = weights, conf.type = conf.type)

  surv_med <- surv_median(fit)

  if(!is.null(weights)){
    ind <- which(!is.na(data[[time]]))
    dat_logrank <- data[ind,]
    pval <- logrank_IPSW_RISCA(dat_logrank[[time]], dat_logrank[[status]], ifelse(dat_logrank[[var]] == levels(dat_logrank[[var]])[1], 0, 1), weights)$p.value
  } else {
    pval <- surv_pvalue(fit)$pval
  }

  pval <- ifelse(pval < 0.0001, "< 0.0001", round(pval,4))
  tbl <- data.frame(sapply(1:length(surv_med$median),function(x){
    paste(round(surv_med$median[x],round), " (", round(surv_med$lower[x],round),"-", round(surv_med$upper[x],round),")", sep = "")
  }))

  surv_formula <- as.formula(paste("Surv(", time, ",", status,")~ 1", sep = ""))
  fit <- surv_fit(surv_formula, data = data, weights = weights, conf.type = conf.type)

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

