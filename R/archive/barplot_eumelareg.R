#' Barplot with points showing the data distribution in EuMelaReg style
#'
#'
#' This function plots a default EuMelaReg Barplot produced with \code{ggplot2}. Similar to a ggplot2 object,
#' additional custom layers can be added with the \code{+} operator.
#' @param data data.frame or data.table.
#' @param x Independent variable
#' @param y Dependent variable
#' @param n.vjust Offset of n over the maximum value for each group
#' @param bar.border.size Size of the border surrounding the bars.
#' @param bar.alpha Alpha value of the barplot fill argument. If < 1 a semi-transparent plot is generated and
#'                  the graph cannot be saved as metafile anymore.
#' @param bar.width Width of the bars.
#' @param errorbar.fun Function used for calculation of the errrorbars. Default is standard deviation.
#' @param errorbar.width Width of the errorbar whiskers.
#' @param errorbar.size Thickness of errorbar line.
#' @param textsize Size of the axes text defined by the function \code{theme_eumelareg_barplot}
#' @param dotplot.binwidth Indirectly regulates size of the dots. Relative to the plot size.
#' @param brewer.palette Color palette from package \code{RColorBrewer}. Can only be specified if custom.palette is NULL.
#' @param custom.palette Custom color palette used in \code{scale_fill_manual()}. Can only be specified if brewer.palette is NULL.
#' @param ylab Label for y-axis.
#' @param y.breaks Argument to define breaks on y-axis. Default is \code{breakFun}, which displays the breaks as
#'                 a multiple of 5.
#' @examples
#' # generate data with 6 random groups with 500 observations each
#' library(data.table)
#' dt <- data.table(group = factor(rep(1:6,500)), value = rnorm(3000, 50,15))
#'
#' # Barplot with standard deviation and points showing the distribution of the data
#' barplot_eumelareg(data = dt, x="group", y = "value", n.vjust =.07,
#'                   dotplot.binwidth = 0.5, brewer.palette = "Blues")
#' @export

barplot_eumelareg <- function(data,x,y,n.vjust = 0.05, bar.border.size = 1, bar.alpha = 0.5, bar.width = 0.5, y.breaks = breakFun,
                              errorbar.fun = "mean_sdl", errorbar.width = 0.3, errorbar.size = 1.3, textsize = 14,
                              brewer.palette = NULL,custom.palette = NULL, ylab = "Value", dotplot.binwidth = 1){
  n <- y.pos <- NULL
  data <- data.table::data.table(data)

  n <- aggregate(as.formula(paste(y,"~",x)), data = data, FUN = length)
  n$y.pos <- data[data[, .I[eval(parse(text=y)) == max(eval(parse(text=y)))], by=eval(parse(text=x))]$V1][[y]] * (1+n.vjust)
  lim <- max(data[[y]])* (1+2*n.vjust)


  p <- ggplot(data, aes_string(x, y, fill = x)) +
    geom_bar(stat = "summary", fun = "mean", color = "black", size = bar.border.size,
             alpha = bar.alpha, width = bar.width)+
    geom_errorbar(stat = "summary", fun.data = errorbar.fun,
                  fun.args = list(mult = 1), width  = errorbar.width, size = errorbar.size) +
    geom_dotplot(binaxis = "y", stackdir = "center", binwidth = dotplot.binwidth)+ # oder geom_jitter als Alternative?
    geom_text(data = n, aes(label = eval(parse(text=y)), y = y.pos),  vjust = -1)+
    theme_eumelareg_barplot(textsize = textsize) +
    scale_y_continuous(expand = c(0, 0),breaks = y.breaks,  limits  = c(0, lim))+
    ylab(ylab)

  if(!is.null(brewer.palette) & is.null(custom.palette)) {
    p + scale_fill_brewer(palette = brewer.palette)
  } else if (!is.null(custom.palette) & is.null(brewer.palette)) {
    p + scale_fill_manual(values = custom.palette)
  } else  if (!is.null(brewer.palette) & !is.null(custom.palette)) {
    stop("Please only specify one of the two arguments \"brewer palette\" and \"custom palette\" ")
  } else {
    p
  }

}

