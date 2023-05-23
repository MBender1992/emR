#' Base function to display a forestplot
#'
#' This code generates a helper function used in the wrappers "forestplot_eumelareg" and "forestplot_meta_eumelareg".
#' @inheritParams survminer::ggforest
#' @param toShow Results of a cox regression containing variable names, variable levels, hazards, confidence intervals and p-values, which will be converted to a forestplot.
#' @param point_size Size of mean points.
#' @param line_size Size of errorbar line.
#' @param vjust_text vertical adjustment of text containing information about events, global pvalue, AIC and concordance index
#' @param y_breaks argument to supply manual y_breaks as a numerical vector. Default is NULL and breaks are set automatically within the function.
#' @param ylim argument to supply manual y limits as numerical vector of length 2. Default is NULL and limits are set automatically within the function.
#' @export



forest_plotFUN <- function(toShow, main, y_breaks, cpositions, point_size, fontsize, line_size, vjust_text, refLabel, noDigits, ylim, varnames){
  if(!is.null(varnames)) toShow$var <- stringi::stri_replace_all_fixed(toShow$var, pattern = vars, replacement = varnames, vectorize_all = FALSE)
  toShowExp <- toShow[, 5:7]
  toShowExp[is.na(toShowExp)] <- 0
  toShowExp <- format(exp(toShowExp), digits = noDigits)
  toShowExpClean <- data.frame(toShow, pvalue = signif(toShow[, 4], noDigits + 1), toShowExp)
  toShowExpClean$stars <- paste0(round(toShowExpClean$p.value,noDigits + 1), " ",
                                 ifelse(toShowExpClean$p.value < 0.05, "*", ""),
                                 ifelse(toShowExpClean$p.value < 0.01, "*", ""),
                                 ifelse(toShowExpClean$p.value < 0.001, "*", ""))
  toShowExpClean$ci <- paste0("(", toShowExpClean[, "conf.low.1"], " - ", toShowExpClean[, "conf.high.1"], ")")
  toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] = "<0.001 ***"
  toShowExpClean$stars[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$ci[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$estimate[is.na(toShowExpClean$estimate)] = 0
  toShowExpClean$var = as.character(toShowExpClean$var)
  toShowExpClean$var[duplicated(toShowExpClean$var)] = ""
  toShowExpClean$N <- paste0("(N=", toShowExpClean$N, ")")
  toShowExpClean$levelN <- paste(toShowExpClean$level, toShowExpClean$N)
  toShowExpClean$estimateCI <- paste(toShowExpClean$estimate.1,  toShowExpClean$ci)
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1,]
  toShowExpClean$estimate <- ifelse(toShowExpClean$estimate == 0, NA, toShowExpClean$estimate)
  rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high,   na.rm = TRUE)
  if(is.null(y_breaks)) breaks <- grDevices::axisTicks(rangeb/2, log = TRUE, nint = 7) else breaks <- y_breaks
  rangeplot <- rangeb
  rangeplot[1] <- rangeplot[1] - diff(rangeb)
  rangeplot[2] <- rangeplot[2] + 0.15 * diff(rangeb)
  if (!is.null(ylim)) {
    rangeplot <- log(ylim)
    if (any(1.3*log(utils::tail(breaks, n = 1)) < toShowExpClean[!is.na(toShowExpClean$conf.high),]$conf.high)) message("Some upper confidence intervals have been cut in favor of better display.")
    if (any(1.3*log(breaks[1]) >  toShowExpClean[!is.na(toShowExpClean$conf.low),]$conf.low)) message("Some lower confidence intervals have been cut in favor of better display.")
    toShowExpClean$estimate <- ifelse(toShowExpClean$estimate < log(ylim[1]), NA, toShowExpClean$estimate)
    toShowExpClean$conf.high <- ifelse(toShowExpClean$estimate < log(ylim[1]), NA, toShowExpClean$conf.high)
    toShowExpClean$conf.low <- ifelse(toShowExpClean$estimate > log(ylim[2]), NA, toShowExpClean$conf.low)
    toShowExpClean$conf.high <- ifelse(1.3*log(utils::tail(breaks, n = 1)) < toShowExpClean$conf.high, 1.3*log(utils::tail(breaks, n = 1)), toShowExpClean$conf.high)
    toShowExpClean$conf.low <- ifelse(1.3*log(breaks[1]) > toShowExpClean$conf.low,  1.3*log(breaks[1]), toShowExpClean$conf.low)
  }
  width <- diff(rangeplot)
  y_variable <- rangeplot[1] + cpositions[1] * width
  y_nlevel <- rangeplot[1] + cpositions[2] * width
  y_cistring <- rangeplot[1] + cpositions[3] * width
  y_stars <- rangeb[2]
  x_annotate <- seq_len(nrow(toShowExpClean))
  annot_size_mm <- fontsize * as.numeric(grid::convertX(unit(theme_get()$text$size,  "pt"), "mm"))

  p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) +
    geom_rect(aes(xmin = seq_along(var) - 0.5, xmax = seq_along(var) +  0.5, ymin = exp(rangeplot[1]),
                  ymax = exp(rangeplot[2]),  fill = ordered(seq_along(var)%%2 + 1))) +
    # color of the rectangles
    scale_fill_manual(values = c("#FFFFFF33", "grey95"), guide = "none") +
    geom_errorbar(aes(ymin = exp(conf.low), ymax = exp(conf.high)), size = line_size, width = 0) +
    geom_point(pch = 16, size = point_size, color = "#009AA6") +
    geom_hline(yintercept = 1, linetype = 2) + coord_flip(ylim = exp(rangeplot)) +
    ggtitle(main) +
    theme_light() +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = "none",
          panel.border = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = fontsize*13),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = fontsize *13, hjust = 0.5)) +
    xlab("") +
    annotate(geom = "text", x = x_annotate, y = exp(y_variable), label = toShowExpClean$var,
             fontface = "bold",  hjust = 0, size = annot_size_mm) +
    annotate(geom = "text", x = x_annotate, y = exp(y_nlevel), hjust = 0,
             label = toShowExpClean$levelN, size = annot_size_mm) +
    annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
             label = toShowExpClean$estimateCI, size = annot_size_mm) +
    annotate(geom = "text", x = x_annotate, y = if (!is.null(ylim))   ylim[2] - 0.4 * ylim[2]  else exp(y_stars),
             label = toShowExpClean$stars, size = annot_size_mm,   hjust = -0.2, fontface = "italic")

  if (!is.null(y_breaks)) {
    p <- p + scale_y_log10(name = "", expand = c(0.02,  0.02), breaks = breaks)
  } else {
    p <- p + scale_y_log10(name = "", labels = sprintf("%g",breaks), expand = c(0.02, 0.02), breaks = breaks)
  }

  gt <- suppressWarnings(ggplot_gtable(ggplot_build(p)))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  ggpubr::as_ggplot(gt)
}




