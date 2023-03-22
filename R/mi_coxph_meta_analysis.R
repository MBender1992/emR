#' COmparison of two treatments with cox regression
#'
#' This code generates the input for a forestplot comparing two treatment options.
#' @param data data
#' @param time the time interval from start of observation until date of event (e.g. disease progression or death)
#' or censoring.
#' @param status variable specifying if event occured or data has been censored.
#' @param vars variables tested for Influence on outcome.
#' @param var grouping variable
#' @param meta.group variable for which meta analysis should be conducted. Usually the outcome of interest (e.g. treatment).
#' @param weights Logical value indicating wether inverse propensity score weighting should be applied
#' @param m number of iterations for multiple imputation
#' @export

mi_coxph_meta_analysis <- function(data, time, status, vars, var, meta.group, weights = FALSE, m = 25){

  vars_input <- NULL

  if(!is.factor(data[[meta.group]])) stop("Grouping variable has to be a factor.")
  if(length(levels(data[[meta.group]])) != 2) stop("Grouping factor must have exactly two levels")
  if(length(vars) <= 1) stop("Multiple imputation is only possible with more than one input variable")
  if(is.numeric(data[[var]])) stop(paste("Please categorize variable", var, "in order to perform multiple imputation"))

  lapply(1:length(levels(data[[var]])), function(x) {
      vars_coxph <- c(vars[vars != var], meta.group)
      dat <- data[data[[var]] == levels(data[[var]])[x],]
      for(i in vars_coxph){
        dat[[i]] <- droplevels(dat[[i]])
        if(length(levels(dat[[i]])) < 2){
          vars_coxph <- vars_coxph[vars_coxph != i]
        }
      }
      fit_mi <- mi_coxph(data = dat, time = time, status = status, vars = vars_coxph, prop.var = if(weights == TRUE) meta.group,  m = m)
      df <- as.data.frame(fit_mi$fit)
      df$robust.se <- df$std.error
      df <- df[,c(1,2,3,8,4,5,6,7)]
      df$var <- var
      df$level <- levels(data[[var]])[x]
      df$N <- dim(dat)[1]
      df})
}





