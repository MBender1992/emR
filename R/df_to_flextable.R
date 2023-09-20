#' Transfrom dataframe into flextable object
#'
#' @param tbl dataframe containing table data
#' @param data original data which was used to generate the table
#' @param vars_tbl variables used to generate the table
#' @param indent logical indicating whether factor levels should be indented. Default is true
#' @export

df_to_flextable <- function(tbl, data, vars_tbl, indent = TRUE){
  levels <- unlist(sapply(vars_tbl, function(var){levels(data[[var]])}))
  if(indent == TRUE) levels <- paste("  ", levels, sep = "")
  levels <- c(levels, "  Mean (SD)", "  Median [Min, Max]")
  ind <- !tbl[,1] %in% levels
  flx_tbl <- flextable(tbl)
  flx_tbl <- bold(flx_tbl, i = ind, j = 1)
  flx_tbl <- theme_flextable(flx_tbl)
  # flx_tbl <- bold(flx_tbl, bold = TRUE, part = "header")
  flx_tbl <- padding(flx_tbl, i=!ind, j=1, padding.left=10)
  FitFlextableToPage(flx_tbl, pgwidth = 7)
}


