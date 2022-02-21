#' Standardize date format of clinical studies
#'
#' Missing months are replaced by June, missing days by 15th of the month. The date format is then transferred to a
#' yyyymmdd
#' format.
#' @param x vector containing dates.
#' @export

convert_date <- function(x){
  x <- ifelse(stringr::str_detect(x, "^\\d{4}$"), paste(x, "-06", sep = ""), x)
  x <- ifelse(stringr::str_detect(x, "^\\d{4}-\\d{2}$"), paste(x, "-15", sep = ""), x)
  ymd <- lubridate::ymd(x)
  dmy <- lubridate::dmy(x)
  dmy[is.na(dmy)] <- ymd[is.na(dmy)]
  return(dmy)
}



