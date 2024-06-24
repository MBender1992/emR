#' Part of helper functions used in the emR package
#'
#' This function rounds values to the closest number, which can be divided by five and adjust to different orders of
#' magnitude by using different exponentials of five depending on the difference between the maximum and minimum value in x.
#'
#' @param x Vector of values to be rounded.
#' @export
roundUp <- function(x) 5^ceiling(log(abs(x),5))

#' Part of helper functions used in this package
#'
#' Creates breaks for \code{scale_y_continuous} between 0 and the maximum value of y for graphs with positive y, rounded with the \code{roundUp}
#' function. This ensures that the space between y-axis breaks is always a multiple of five.
#'
#' @param y Vector used to construct the sequence between 0 and the maximum value. If used as \code{scale_y_continuous(breaks = breakFun)},
#' y is inherited from the ggplot object.
#' @seealso [breakFun2()] for negative y.
#' @export

breakFun <- function(y) seq(0, ceiling(max(y)), by = roundUp(max(y)/10))

#' Part of helper functions used in this package
#'
#' Creates breaks for \code{scale_y_continuous} between 0 and the minimum value of y for graphs with negative y, rounded with the \code{roundUp}
#' function. This ensures that the space between y-axis breaks is always a multiple of five.
#'
#' @param y Vector used to construct the sequence between 0 and the minimum value. If used as \code{scale_y_continuous(breaks = breakFun)},
#' y is inherited from the ggplot object.
#' @seealso [breakFun()] for positive y.
#' @export

breakFun2 <- function(y) seq(0, floor(min(y)), by = -roundUp(max(y)/10))



