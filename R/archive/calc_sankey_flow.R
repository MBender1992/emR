#' Calculate data flow
#'
#' This function calculates the flow from 1 stratum to the next.
#' @param x number of column containing source data
#' @param y number of column containing target data
#' @param data data.frame containing the different strata with one row for each patient
#' @export

calc_sankey_flow <- function(x,y,data){
  dat <- as.data.frame(table(data.frame(source = data[,x], target = data[,y])))
  dat$source <-paste0(dat$source, "_",x)
  dat$target <- paste0(dat$target,"_",y)
  return(dat)
}
