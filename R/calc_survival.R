#' Calculate survival interval
#'
#' This function calculates the survival interval between start of therapy and time of event (e.g. progression, death, etc.)
#' in months.
#' @param startDate start of therapy.
#' @param eventDate date of event.
#' @export

# calculate PFS
calc_survival <- function(startDate, eventDate){
  time <- eventDate - startDate + 1
  round(lubridate::time_length(time,unit="months"),2)
}
