#' Save Pie Chart generated with \code{pie_chart_eumelareg}
#'
#' This code saves a pie chart to show distribution of variables grouped by EuMelaReg registries.
#' @param data dataframe.
#' @param y Name of column to be shown.
#' @param regID Name of column containing information about registries.
#' @param registry Name of EuMelaReg registry.
#' @param path File path where the plot shall be stored. By default a folder "Results/pie_charts" is generated
#' within the working directory and objects are stored in this subfolder.
#' @param title Plot title.
#' @param ... Additional parameters passed on to \code{pie_chart_eumelareg}
#' @export

save_pie_chart_eumelareg <- function(data,y,regID,registry = "EUMELAREG",title=NULL, path = NULL,...){

  if (registry != "EUMELAREG"){
    data <- data[eval(parse(text=regID)) == registry]
  }

  if (is.null(path)) {
    if (!dir.exists("Results/pie_charts")) {
      dir.create("Results/pie_charts", recursive = TRUE)
    }
    path <-  "Results/pie_charts"
  }

  filename <- paste("/pie_chart_",y,"_",registry,".png", sep = "")
  png(paste(path,filename,sep=""), units="in", width=5, height=4, res=600)
  p <- pie_chart_eumelareg(data[[y]], text.size = 4,
                           title = paste(title, registry),...)
  print(p)
  dev.off()
}



