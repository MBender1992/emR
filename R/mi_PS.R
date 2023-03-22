#' Calculate Propensity Scores with multiple imputation
#'
#' This code generates propensity scores with multiple imputed values for missing data. Imputation is conducted with the [mice::mice()] function.
#' Results are averaged for each row of entries representing single iterations.
#' @param data data.frame or data.table containing survival data.
#' @param vars variables tested for Influence on outcome.
#' @inheritParams mice::mice
#' @param prop.var variable for which propensity scores should be calculated. If no value is provided (prop.var = NULL), no weights are used in coxph. Default is NULL.
#' @param ... additional arguments to be passed on to coxph function
#' @export


mi_PS <- function(data, vars, prop.var = NULL,  m = 5, ...){
  weights_ate <- NULL

  # define data and variables
  dat <- as.data.frame(data)
  dat <- dat[,c(vars)]
  vars_input <- paste(vars, collapse = " + ")

  # impute missing data
  set.seed(9)
  imp <- mice(dat, m = m)
  imp_comp <- complete(imp, "long")

  # calculate coxph and frequency of factor levels for each iteration of the multiple imputation
  ls_PS <- lapply(1:m, function(x){
    tmp <- imp_comp[imp_comp$.imp ==1,]
    if(!is.null(prop.var)) tmp$weights.ate <- ate_weights(tmp, vars, prop.var = prop.var)
    return(tmp)
  })
  ls_PS_vec <- lapply(1:m, function(x){
    ls_PS[[x]]$weights.ate
  })
  res <-as.data.frame(do.call(cbind, ls_PS_vec))
  rowMeans(res)
}
