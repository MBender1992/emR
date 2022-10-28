#' Print results from stratified univariate cox regression
#'
#' This function combines the results from multiple calls of \code{cox_output} and prints the output
#' as an htmlTable generated with the \code{htmlTable} package stratified by a factor.
#' @param data data.frame or data.table containing survival data.
#' @param time the time interval from start of observation until date of event (e.g. disease progression or death)
#' or censoring.
#' @param status variable specifying if event occured or data has been censored.
#' @param vars one or more variables defined as character strings to be included in the table
#' @param strata factor to stratidfy by
#' @param footnote character string passed on to \code{tfoot} argument of \code{htmlTable}
#' @param printHTML Logical value. If TRUE output is printed as htmlTable. Default is TRUE.
#' @param univariate Logical value. If TRUE output of univariate cox regression is printed. Else output of multivariate
#' cox regression is printed. Not defined yet. TBD.
#' @param ... additional arguments to be passed on to \code{cox_output}
#' @export

strat_cox_table <- function(data, time, status, vars, footnote = NULL, strata,
                            printHTML = TRUE, univariate = TRUE, ...){

  if(univariate == FALSE) stop("Function not defined yet for multivariable cox regression")

  out_tmp <- lapply(1:length(levels(data[[strata]])), function(x){
    ind <- data[[strata]] == levels(data[[strata]])[x]
    input <- data[ind]

    n <- length(vars)
    tmp <- lapply(vars, cox_output, data = input, time = time, status = status,...)
    res <- lapply(1:n, function(x){
      if(dim(tmp[[x]])[1] == 1) {
        tmp[[x]]
      }
      else {
        tmp[[x]][-1,]
      }
    })
    dplyr::bind_rows(res)[1]
  })

  out <- dplyr::bind_cols(out_tmp)
  rownames(out) <- gsub(">=", "&#8805", rownames(out))
  colnames(out) <- levels(data[[strata]])

  if (printHTML == TRUE){
    htmlTable::htmlTable(out, tfoot = footnote)
  } else{
    return(out)
  }
}
