#' Save survival plot generated with survplot_eumelareg as png
#'
#' This function saves the output of  survplot_eumelareg as a png file with a default resolution of 600dpi in a 12x8 in format.
#' @param data data.frame or data.table containing survival data.
#' @param time the time interval from start of observation until date of event (e.g. disease progression or death)
#' or censoring.
#' @param status variable specifying if event occured or data has been censored.
#' @param var variable tested for Influence on outcome.
#' @param weights character variable specifying the name of the weights column. Weights have to be added to the original dataframe in order to be applied correctly.
#' @param path path where the files should be saved
#' @param width width of the png file in inches
#' @param height height of the png file in inches
#' @param ... additional arguments to be passed on to survplot_eumelareg() function
#' @export

save_surv_plot_EMR <- function(data, time, status, var,  weights = NULL, path, width = 12, height = 8, ...){
  if(!is.null(weights)){
    filename <- paste(time, var, "IPSW", sep = "_")
  } else {
    filename <- paste(time, var, sep = "_")
  }
  p <- survplot_eumelareg(data, time = time, status = status, var = var, weights = weights,
                          legend.labs = sort(unique(data[[var]])), ...)
  if (is.list(p)) {
    png(paste(path, filename,"_plot.png", sep = ""), units="in", width=width, height=height, res=600)
    print(p$plot)
    dev.off()
  } else {
    png(paste(path, filename,".png", sep = ""), units="in", width=width, height=height, res=600)
    print(p)
    dev.off()
  }
}
