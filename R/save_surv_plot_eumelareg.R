#' Save survival curve generated with \code{survplot_eumelareg}
#'
#' This function saves a default EuMelaReg survival plot (Kaplan-Meier plot) produced with \code{survminer}.
#' @param data dataframe containing survival data.
#' @param time The time interval from start of observation until date of event (e.g. disease progression or death)
#' or censoring.
#' @param status Variable specifying if event occured or data has been censored.
#' @param var Variable tested for Influence on outcome.
#' @param legend.labs Character vector specifying legend labels. Used to replace the names of the strata from the fit.
#' @param regID Name of column containing registry information. Default is "REGISTR" as defined by SAS output.
#' @param registry Name of EuMelaReg registry. If \code{registry = "EUMELAREG"} data from all registries is plotted.
#' @param path File path where the plot shall be stored. By default a folder "Results/survival_curve" is generated
#' within the working directory and objects are stored in this subfolder.
#' @param ... Additional arguments passed to \code{\link[survminer]{ggsurvplot}}.
#' @export

save_surv_plot_eumelareg <- function(data,time, status, var,legend.labs=NULL, regID = "REGISTR", registry = "EUMELAREG",path = NULL,...){
  warning("This function is outdated. Please use save_surv_plot_EMR.")
  LOT <- stringr::str_extract(time, "\\d")
  tmp <- stringr::str_extract(status, stringr::fixed(c("PFS","OS"), ignore_case=TRUE))
  survival <- tmp[-which(is.na(tmp))]

  if (is.null(path)) {
    if (!dir.exists(paste("Results/survival_curve/",registry, sep = ""))) {
      dir.create(paste("Results/survival_curve/",registry, sep = ""), recursive = TRUE)
    }
    path <- paste("Results/survival_curve/",registry, sep = "")
  }

  if (ifelse(is.na(stringr::str_extract(var,"\\d$") != LOT),TRUE, stringr::str_extract(var,"\\d$") != LOT)){
    warning("Please check if grouping variable corresponds to the correct line of treatment")
  }

  if (registry != "EUMELAREG"){
    data <-  data[eval(parse(text=regID)) == registry]
    data[[var]] <- droplevels(data[[var]])
  }

  filename <- paste("/survival_curve_",var,"_",registry,"_",survival,LOT,".png", sep = "")
  png(paste(path, filename, sep = ""), units="in", width=12, height=6.2, res=600)
  if (!is.null(legend.labs)){

    t <- try(parse(text = legend.labs), silent = TRUE)
    if("try-error" %in% class(t)){
      labs <- legend.labs
    } else {
      labs <- eval(parse(text = legend.labs))
    }

    if(survival == "PFS"){
      p <- survplot_eumelareg(data = data, time = time, status = status,
                              legend.labs = labs,
                              var = var,ylab = "Probability of Progression Free Survival",...)
    } else if(survival == "OS"){
      p <- survplot_eumelareg(data = data, time = time, status = status,
                              legend.labs = labs,
                              var = var,  ylab = "Probability of Overall Survival",...)
    } else {
      stop("Please specify type of survival")
    }
  } else {
    if(survival == "PFS"){
      p <- survplot_eumelareg(data = data, time = time, status = status,
                              var = var,ylab = "Probability of Progression Free Survival",...)
    } else if(survival == "OS"){
      p <- survplot_eumelareg(data = data, time = time, status = status,
                              var = var,  ylab = "Probability of Overall Survival",...)
    } else {
      stop("Please specify type of survival")
    }
  }
  print(p)
  dev.off()
}


