#' Meta analysis forest plot
#'
#' This code generates a forestplot from a meta analysis coxph model.
#' @inheritParams survminer::ggforest
#' @param data data
#' @param time the time interval from start of observation until date of event (e.g. disease progression or death)
#' or censoring.
#' @param status variable specifying if event occured or data has been censored.
#' @param vars variables tested for Influence on outcome.
#' @param meta.group variable for which meta analysis should be conducted. Usually the outcome of interest (e.g. treatment).
#' @param univariate Logical value. If TRUE output of univariate cox regression is printed. Else output of multivariate
#' @param varnames Character vector specifying rownames of the table (empty columns should be named with "").
#' @param point_size Size of mean points.
#' @param line_size Size of errorbar line.
#' @param weights Vector containing sample weights.
#' @param n.boot number of bootstrapping resamples
#' @param n.sample number of samples drawn. Default is the length of the input data.frame
#' @param vjust_text vertical adjustment of text containing information about events, global pvalue, AIC and concordance index
#' @param y_breaks argument to supply manual y_breaks as a numerical vector. Default is NULL and breaks are set automatically within the function.
#' @param ylim argument to supply manual y limits as numerical vector of length 2. Default is NULL and limits are set automatically within the function.
#' @export

bt_forestplot_meta_eumelareg <- function (data, time, status, vars, meta.group, univariate = TRUE, weights = NULL,
                                       main = "Hazard ratio for disease progression or death (95% CI)",
                                       y_breaks = NULL, cpositions = c(0, 0.1, 0.3), refLabel =  "reference",
                                       point_size = 4, fontsize = 1,line_size = 0.9, vjust_text = 1.2, noDigits = 2,
                                       n.boot = 10000, n.sample = nrow(data),  varnames = NULL, ylim = NULL){

    if(!is.null(weights))  data$weights_cox <- weights
    ls <- lapply(vars, bt_coxph_meta_analysis, data = data, time = time, weights = "weights_cox",
                 status = status, vars = vars, meta.group = meta.group,
                 univariate = univariate, n.boot = n.boot, n.sample = n.sample)

    cox_formula <- as.formula(paste("Surv(", time, ", ", status,") ~ ",meta.group, "+", paste(vars, collapse = "+"), sep = ""))
    fit_total <- coxph(cox_formula, weights = data[["weights_cox"]], data = data)
    df_total <- as.data.frame(broom::tidy(fit_total, conf.int = TRUE))
    perform bootstrapping
    bt_fit <- pbapply::pbreplicate(n = n.boot, {
      dat_boot <- data[sample(nrow(data), size = n.sample, replace=T), ]
      surv <- coxph(cox_formula, data = dat_boot)
      surv$coefficients
    })

    # convert to df and format dataframe
    bt_df <- data.frame(matrixStats::rowQuantiles(bt_fit, probs = c(0.025, 0.5, 0.975)))
    colnames(bt_df) <- c("conf.low", "estimate", "conf.high")
    df_total$conf.low <- bt_df$conf.low
    df_total$conf.high <- bt_df$conf.high
    df_total$estimate <- bt_df$estimate
    df_total$p.value <- sapply(1: nrow(df_total), function(i) ci_pval(est = df_total[i,]$estimate, l = df_total[i,]$conf.low, u = df_total[i,]$conf.high, log.trans = FALSE))
    df_total <- df_total[stringr::str_detect(df_total$term, meta.group),]
    df_total$var <- "Overall"
    df_total$level <- ""
    df_total$N <- dim(data)[1]
    rownames(df_total) <- "Overall HR"

  toShow <- lapply(1:length(ls), function(x) {
    toShow <- if(is.data.frame(ls[[x]][[1]]))  do.call(rbind, ls[[x]])  else ls[[x]]
    toShow <- toShow[toShow$term == paste(meta.group, levels(data[[meta.group]])[2],sep = ""), -c(1, 4)]
    rownames(toShow) <- paste(toShow$var, toShow$level, sep = "")
    toShow <- toShow[, c("var", "level", "N", "p.value", "estimate", "conf.low", "conf.high")]
    return(toShow)
  })
  toShow <- do.call(rbind, toShow)
  df_total <- df_total[, c("var", "level", "N", "p.value", "estimate", "conf.low", "conf.high")]
  toShow <- rbind(toShow, df_total)

  # plot the Forestplot
  forest_plotFUN(toShow = toShow, main = main,  y_breaks = y_breaks, cpositions = cpositions, point_size = point_size, varnames = varnames,
                 fontsize = fontsize, line_size = line_size, vjust_text = vjust_text, refLabel = refLabel, noDigits = noDigits, ylim =ylim)

}


