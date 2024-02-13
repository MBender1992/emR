#' Coxph with multiple imputed values
#'
#' This code generates a coxph model with multiple imputed values for missing data. Imputation is conducted with the [mice::mice()] function.
#' Pooled results obtained with [mice::pool()] and respective pvalues and 95% CIs are presented as results.
#' @param data data.frame or data.table containing survival data.
#' @param time the time interval from start of observation until date of event (e.g. disease progression or death)
#' or censoring.
#' @param status variable specifying if event occured or data has been censored.
#' @param vars variables tested for Influence on outcome.
#' @inheritParams mice::mice
#' @param prop.var variable for which propensity scores should be calculated. If no value is provided (prop.var = NULL), no weights are used in coxph. Default is NULL.
#' @param ... additional arguments to be passed on to coxph function
#' @export

mi_coxph <-
  function(data,
           time,
           status,
           vars,
           prop.var = NULL,
           m = 5,
           ...) {
    weights.ate <- NULL

    # define data and variables
    dat <- as.data.frame(data)
    dat <- dat[, c(vars, time, status)]
    vars_input <- paste(vars, collapse = " + ")

    # impute missing data
    set.seed(9)
    imp <- mice(dat, m = m)
    imp_comp <- complete(imp, "long")

    # calculate coxph and frequency of factor levels for each iteration of the multiple imputation
    fitm <- lapply(1:m, function(x) {
      tmp <- imp_comp[imp_comp$.imp == x, ]
      terms <-
        as.character(utils::tail(data.frame(rbind(
          lapply(data, class)
        )[, vars]), n = 1))
      names(terms) <- vars
      allTerms <- lapply(seq_along(terms), function(i) {
        var <- names(terms)[i]
        if (terms[i] %in% c("factor", "character")) {
          adf <- as.data.frame(table(tmp[, var]))
          cbind(var = var, adf, pos = 1:nrow(adf))
        }
        else if (terms[i] %in% c("numeric", "integer")) {
          data.frame(
            var = var,
            Var1 = "",
            Freq = nrow(tmp),
            pos = 1
          )
        }
      })
      if (!is.null(prop.var))
        tmp$weights.ate <- ate_weights(tmp, vars, prop.var = prop.var)
      fit <-
        coxph(as.formula(paste(
          "Surv(", time, ", ", status, ") ~ ", vars_input, sep = ""
        )), weights = weights.ate, data = tmp)
      list(allTerms = allTerms, fit = fit)
    })

    # extract frequency of factor levels and calculate mean
    allTerms_ls <- lapply(1:m, function(x) {
      fitm[[x]]$allTerms
    })
    n <- nested_list_mean(allTerms_ls, vars)

    # pool results from coxph function
    estimates <- lapply(1:m, function(x)
      fitm[[x]]$fit)
    res <- summary(pool(estimates), conf.int = TRUE)
    res$df <- NULL
    colnames(res) <-
      c("term",
        "estimate",
        "std.error",
        "statistic",
        "p.value",
        "conf.low",
        "conf.high")
    if (!is.null(prop.var))
      message("Inverse propensity score weighting was used within cox regression.")

    # return results
    list(fit = res, nfit = n)
  }
