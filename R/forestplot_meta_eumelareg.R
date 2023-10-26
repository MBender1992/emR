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
#' @param weights Logical variable specifying whether inverse propensity score weighting should be applied.
#' @param imputation Logical variable specifying whether missing should be replaced by multiple imputation
#' @param impIter number of iterations for multiple imputation.
#' @param vjust_text vertical adjustment of text containing information about events, global pvalue, AIC and concordance index
#' @param y_breaks argument to supply manual y_breaks as a numerical vector. Default is NULL and breaks are set automatically within the function.
#' @param ylim argument to supply manual y limits as numerical vector of length 2. Default is NULL and limits are set automatically within the function.
#' @export

forestplot_meta_eumelareg <- function (data, time, status, vars, meta.group, univariate = TRUE, weights = NULL, imputation = FALSE,
                                       main = "Hazard ratio for disease progression or death (95% CI)",
                                       y_breaks = NULL, cpositions = c(0, 0.1, 0.3), impIter = 25, refLabel =  "reference",
                                       point_size = 4, fontsize = 1,line_size = 0.9, vjust_text = 1.2, noDigits = 2,
                                       varnames = NULL, ylim = NULL){

  if(imputation == TRUE){
    ls <- lapply(vars, mi_coxph_meta_analysis, data = data, time = time, weights = weights, m = impIter,
                 status = status, vars = vars, meta.group = meta.group)
  } else{
    if(!is.null(weights)){
      data$weights_cox <- weights
    } #else {
      #data$weights_cox <- NULL
    #} remove the comment if function throws error
    ls <- lapply(vars, coxph_meta_analysis, data = data, time = time, weights = "weights_cox",
                 status = status, vars = vars, meta.group = meta.group,
                 univariate = univariate)
  }
  toShow <- lapply(1:length(ls), function(x) {
    toShow <- if(is.data.frame(ls[[x]][[1]]))  do.call(rbind, ls[[x]])  else ls[[x]]
    toShow <- toShow[toShow$term == paste(meta.group, levels(data[[meta.group]])[2],sep = ""), -c(1, 4)]
    rownames(toShow) <- paste(toShow$var, toShow$level, sep = "")
    toShow <- toShow[, c("var", "level", "N", "p.value", "estimate", "conf.low", "conf.high")]
    return(toShow)
  })
  toShow <- do.call(rbind, toShow)

  # plot the Forestplot
  forest_plotFUN(toShow = toShow, main = main,  y_breaks = y_breaks, cpositions = cpositions, point_size = point_size, varnames = varnames,
                 fontsize = fontsize, line_size = line_size, vjust_text = vjust_text, refLabel = refLabel, noDigits = noDigits, ylim =ylim)

}
