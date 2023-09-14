#' Expand flextable to fit on word document
#'
#' @param ft flextable object
#' @param pgwidth width of the desired word document. Default is 6 inches
#' @export

FitFlextableToPage <- function(ft, pgwidth = 6){

  ft_out <- autofit(ft)

  ft_out <- width(ft_out, width = dim(ft_out)$widths*pgwidth /(flextable_dim(ft_out)$widths))
  return(ft_out)
}
