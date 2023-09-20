#' Transform table1 to flextable object
#'
#' @param tbl1 table1 object
#' @param pgwidth width of the table in word document
#' @param textsize fontsize of text body
#' @param headsize fontsize of header
#' @export

tbl1_to_flextable <- function(tbl1, pgwidth = 7, textsize = 9, headsize = 10){
  tbl1 <- table1::t1flex(tbl1)
  tbl1 <- theme_flextable(tbl1, textsize = textsize, headsize = headsize)
  FitFlextableToPage(tbl1, pgwidth = pgwidth)
}