#' Define common formatting for flextables
#'
#' @param ft flextable object
#' @param textsize fontsize of text body
#' @param headsize fontsize of header
#' @export

theme_flextable <- function(ft, textsize = 9, headsize = 10){
  tmp <- flextable::bold(ft, bold = TRUE, part = "header")
  tmp <- fontsize(tmp,  size = textsize, part = "body")
  tmp <- padding(tmp, padding = 0, part = "all")
  tmp <- line_spacing(tmp, space = 1, part = "all")
  out <- fontsize(tmp,  size = headsize, part = "header")
  return(out)
}
