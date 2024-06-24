#' Meta analysis forest plot with bootstrapping CIs
#'
#' This code generates a forestplot from a meta analysis coxph model.
#' @param data data
#' @param time the time interval from start of observation until date of event (e.g. disease progression or death)
#' or censoring.
#' @param status variable specifying if event occured or data has been censored.
#' @param vars variables tested for Influence on outcome.
#' @param var grouping variable
#' @param n.boot number of bootstrapping resamples
#' @param n.sample number of samples drawn. Default is the length of the input data.frame
#' @param meta.group variable for which meta analysis should be conducted. Usually the outcome of interest (e.g. treatment).
#' @param univariate Logical value. If TRUE output of univariate cox regression is printed. Else output of multivariate
#' cox regression is printed. Default is FALSE.
#' @param weights character variable specifying the name of the weights column. Weights have to be added to the original dataframe in order to be applied correctly.
#' @param ... additional arguments passed on to coxph
#' @export

bt_coxph_meta_analysis <- function(data, time, status, vars, var, meta.group, weights = NULL,
                                   univariate = FALSE, n.boot = 10000, n.sample = nrow(data), ...){
  vars_input <- NULL
  warning("This function is outdated. The correct function is included in the Remission script.")
  if (!is.factor(data[[meta.group]]))
    stop("Grouping variable has to be a factor.")
  if (length(levels(data[[meta.group]])) != 2)
    stop("Grouping factor must have exactly two levels")
  if (is.numeric(data[[var]]))
    stop("Numeric variables are not allowed as input.")
  if (univariate == FALSE) {
    vars_coxph <- c(vars[vars != var], meta.group)
    vars_input <- paste(vars_coxph, collapse = " + ")
  } else {
    vars_input <- meta.group
  }
  res <- pbapply::pblapply(1:length(levels(data[[var]])), function(x) {
    dat <- data[data[[var]] == levels(data[[var]])[x], ]
    scaling_factor <- nrow(dat)/nrow(data)
    if(!is.null(weights)) weights <- dat[[weights]]
    fit <- coxph(as.formula(paste("Surv(", time, ", ", status,
                                  ") ~ ", vars_input, sep = "")), data = dat, weights = weights)
    df <- as.data.frame(broom::tidy(fit, conf.int = TRUE))
    bt_fit <- pbapply::pbreplicate(n = n.boot, {
      dat_boot <- dat[sample(nrow(dat), size = round(n.sample * scaling_factor), replace = T), ]
      if(!is.null(weights)) weights <- dat_boot[[weights]]
      surv <- tryCatch(coxph(as.formula(paste("Surv(", time, ", ", status,
                                              ") ~ ", vars_input, sep = "")), data = dat_boot, weights = weights), error=function(e) NULL)
      if(!is.null(surv)){
        surv$coefficients
      } else rep(NA, length(levels(data[[var]])))
    })
    if(is.list(bt_fit)) bt_fit <- t(do.call(rbind, bt_fit))
    bt_df <- data.frame(matrixStats::rowQuantiles(bt_fit, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
    colnames(bt_df) <- c("conf.low", "estimate", "conf.high")
    df$conf.low <- bt_df$conf.low
    df$conf.high <- bt_df$conf.high
    df$estimate <- bt_df$estimate
    df$p.value <- sapply(1:nrow(df), function(i) ci_pval(est = df[i,
    ]$estimate, l = df[i, ]$conf.low, u = df[i, ]$conf.high,
    log.trans = FALSE))
    df$pvalue <- df$var <- var
    df$level <- levels(data[[var]])[x]
    df$N <- dim(dat)[1]
    df
  })
  return(res)
}
