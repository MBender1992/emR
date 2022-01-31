#' Print results from univariate or multivariate cox regression
#'
#' This function combines the results from multiple calls of \code{cox_output} and prints the output
#' as an htmlTable generated with the \code{htmlTable} package or the output of multivariate cox regression.
#' @param data data.frame or data.table containing survival data.
#' @param time the time interval from start of observation until date of event (e.g. disease progression or death)
#' or censoring.
#' @param status variable specifying if event occured or data has been censored.
#' @param vars one or more variables defined as character strings to be included in the table
#' @param rgroup a vector of character strings containing headings for row groups.
#' @param footnote character string passed on to \code{tfoot} argument of \code{htmlTable}
#' @param printHTML Logical value. If TRUE output is printed as htmlTable. Default is TRUE.
#' @param weights character variable specifying the name of the weights column. Weights have to be added to the original dataframe in order to be applied correctly.
#' @param univariate Logical value. If TRUE output of univariate cox regression is printed. Else output of multivariate
#' cox regression is printed. Default is TRUE.
#' @param ... additional arguments to be passed on to \code{cox_output}
#' @export


cox_table <- function(data, time, status, vars, rgroup = NULL, footnote = NULL,
                      printHTML = TRUE, univariate = TRUE, weights = NULL, ...){

  if(is.null(rgroup)){
    rgroup <- vars
  }

  n <- length(vars)
  tmp <- lapply(vars, cox_output, data = data, time = time, status = status, weights = weights, ...)
  res <- lapply(1:n, function(x){
    if(dim(tmp[[x]])[1] == 1) {
      tmp[[x]]
      }
    else {
      tmp[[x]][-1,]
      }
    })
  out <- dplyr::bind_rows(res)
  n.rgroup <- unlist(lapply(1:n, function(x) dim(res[[x]])[1]))
  pvals <- unlist(lapply(1:n, function(x){tmp[[x]][1,]$pvalue}))

  # set attribute argument of rgroup to add global pvalues to the end of the column (wald-test)
  attr_pval <- sapply(1:n, function(x){list(pvals[x])})
  names(attr_pval) <- 1:n
  attr(rgroup, "add") <- attr_pval

  if (univariate == FALSE){
    out <- cox_output(data = data, time = time, status = status, vars = vars, ...)
    attr(rgroup, "add") <- NULL
  }

  rownames(out) <- gsub(">=", "&#8805", rownames(out))
  if (printHTML == TRUE){
    htmlTable::htmlTable(out, tfoot = footnote, rgroup = rgroup, n.rgroup = n.rgroup)
  } else{
    list(res = out, rgroup = rgroup, n.group = n.rgroup)
  }
}



