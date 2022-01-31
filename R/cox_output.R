#' Calculate cox regression and extract Hazard Ratio and pvalues
#'
#' This function calculates a univariate cox regression with the \code{coxph} function from the \code{survival}
#' package if the argument vars is length 1; otherwise a multivariate cox regression is calculated. Hazard ratios and pvalues are extracted and shown in a table format.
#' If a variable in univariate cox regression is divided in more than 2 groups the wald test statistic
#' for overall significance is also shown.
#' @param data data.frame or data.table containing survival data.
#' @param time the time interval from start of observation until date of event (e.g. disease progression or death)
#' or censoring.
#' @param status variable specifying if event occured or data has been censored.
#' @param vars variables tested for Influence on outcome. Defines type of cox regression: if length is 1 the output
#' is a univariate cox regression, for length > 1 the output is a multivariate cox regression
#' @param fixed.var specifies fixed variables to be included in the cox model.
#' @param modeltype character value. Allowed values include: "full" and "backwards".
#' Specify if full model or backwards selection model should be used for multivariate cox regression.
#' @param p.thres pvalue threshold for backwards selection model.
#' @param niter number of iterations for backwards selection model.
#' @param weights character variable specifying the name of the weights column. Weights have to be added to the original dataframe in order to be applied correctly.
#' @param output Defines the output object of the function.  Default is \"table\". Another possible option is \"fit\" to return a coxph object.
#' @export

cox_output <- function(data, time, status, vars, fixed.var = NULL, output = "table",  modeltype = "full", p.thres = 0.1, niter = 10, weights = NULL){


  weights <- if(!is.null(weights))data[[weights]]

  if (length(vars) > 1) {
    vars_input <- paste(vars, collapse = " + ")
  } else {
    vars_input <- vars
  }

  formula <- paste("Surv(",time,", ", status,") ~ ", vars_input, sep = "")
  fit_cox <- coxph(as.formula(formula), data = data, weights = weights)
  res_cox <- summary(fit_cox)
  res_ci <- res_cox$conf.int
  names <- rownames(res_ci)

  if (length(vars) > 1) {
    if (modeltype == "backwards") {
      for (i in 1:niter){

        if (!is.null(fixed.var)) {
          ind <- rownames(res_cox$coefficients)[c(which(stringr::str_detect(rownames(res_cox$coefficients), fixed.var)), which(res_cox$coefficients[,5] < p.thres))]
        } else {
          ind <- rownames(res_cox$coefficients)[which(res_cox$coefficients[,5] < p.thres)]
            }

        vars_backwards <- vars[!is.na(charmatch(vars,ind))]
        vars_input <- paste(vars_backwards, collapse = " + ")


        formula <- paste("Surv(",time,", ", status,") ~ ", vars_input, sep = "")
        fit_cox <- coxph(as.formula(formula), data = data, weights = weights)
        res_cox <- summary(fit_cox)
        res_ci <- res_cox$conf.int
        names <- rownames(res_ci)
      }
    } else if (modeltype == "full") {
      res_ci <- res_ci
    } else {
      stop("Please specify model type. Allowed values are \"full\" for a complete model or \"backwards\" for a model
           with backwards selection")
    }
  }

  HR <- paste(round(res_ci[,1],2)," (", round(res_ci[,3],2), "-", round(res_ci[,4],2), ")", sep = "")
  pval.coef <- ifelse(res_cox$coefficients[,5] < 0.001,"<0.001",round(res_cox$coefficients[,5],3))
  pval.wald <- ifelse(res_cox$waldtest[3] < 0.001, "<0.001",round(res_cox$waldtest[3],3))

  out <- data.frame(HR = HR, p = as.character(pval.coef))
  if (length(vars) == 1){
    if (!is.numeric(data[[vars]])){
      out <- rbind(c(NA,as.character(pval.wald)), out)
      rownames(out) <- c(vars, names)
    }  else {
      rownames(out) <- names
    }
  } else {
    rownames(out) <- names
  }
  colnames(out) <- c("HR (95% CI)", "pvalue")

  if (output == "table") {
    return(out)
  } else if (output == "fit") {
    return(fit_cox)
  } else {
    stop("Please specify output. Allowed values include \"table\" to output a dataframe and \"fit\" to return a coxph object")
  }

}




