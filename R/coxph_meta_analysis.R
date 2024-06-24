#' Meta analysis forest plot
#'
#' This code generates a forestplot from a meta analysis coxph model.
#' @param data data
#' @param time the time interval from start of observation until date of event (e.g. disease progression or death)
#' or censoring.
#' @param status variable specifying if event occured or data has been censored.
#' @param vars variables tested for Influence on outcome.
#' @param var grouping variable
#' @param meta.group variable for which meta analysis should be conducted. Usually the outcome of interest (e.g. treatment).
#' @param univariate Logical value. If TRUE output of univariate cox regression is printed. Else output of multivariate
#' cox regression is printed. Default is FALSE.
#' @param weights character variable specifying the name of the weights column. Weights have to be added to the original dataframe in order to be applied correctly.
#' @param na.rm logical indicating whether missing values (Specified as "missing" or "Missing") should be removed
#' @param ... additional arguments passed on to coxph
#' @export

coxph_meta_analysis <- function(data, time, status, vars, var, meta.group, weights = NULL, univariate = FALSE, na.rm = FALSE, ...){

  vars_input <- NULL

  if(!is.factor(data[[meta.group]])) stop("Grouping variable has to be a factor.")
  if(length(levels(data[[meta.group]])) != 2) stop("Grouping factor must have exactly two levels")

  if(!is.numeric(data[[var]])){
    fct_lvls <- if(na.rm == TRUE){
      levels(data[[var]])[!stringr::str_detect(levels(data[[var]]), "[Mm]issing")]
    } else {
      fct_lvls <- levels(data[[var]])
    }
    res <- lapply(fct_lvls, function(x) {
      if (univariate == FALSE) {
        vars_coxph <- c(vars[vars != var], meta.group)
        vars_input <- paste(vars_coxph, collapse = " + ")
      }
      else {
        vars_input <- meta.group
      }
      dat <- data[data[[var]] == x,] # levels(data[[var]])[x]
      fit <- coxph(as.formula(paste("Surv(", time, ", ",
                                    status, ") ~ ", vars_input, sep = "")),
                   data = dat, weights = if(!is.null(weights)) dat[[weights]])
      df <- as.data.frame(broom::tidy(fit, conf.int = TRUE))
      df$var <- var
      df$level <- x
      df$N <- dim(dat)[1]
      df
    })
   return(res)
  } else {
    dat <- data
    fit <- coxph(as.formula(paste("Surv(", time, ", ",
                                  status, ") ~ ", vars_input, sep = "")),
                 data = dat, weights = if(!is.null(weights)) dat[[weights]])
    df <- as.data.frame(broom::tidy(fit, conf.int = TRUE))
    df$var <- var
    df$level <- ""
    df$N <- dim(dat)[1]
    df
  }
}




