#' Calculate mean for each list element in a nested list
#'
#' This code takes each element i of nested list j and calculates the mean of the list elements with the same index.
#' @param list nested list
#' @param vars variables tested for Influence on outcome.
#' @export

nested_list_mean <- function(list, vars){
  lapply(1:length(vars), function(x){
    tmp <- lapply(list, `[[`, x)
    varlist <-sapply(tmp, "[[", 3)
    tmp[[1]]$Freq <- sapply(1:ifelse(is.null(dim(varlist)), 1, dim(varlist)[1]), function(y){
      if (!is.null(dim(varlist))) round(mean(varlist[y,]))
      else round(mean(varlist))
    })

    n <- tmp[[1]]
    return(n)
  })
}
