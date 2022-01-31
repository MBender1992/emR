#' Survival curve in EuMelaReg style
#'
#' This function plots a default EuMelaReg survival plot (Kaplan-Meier plot) produced with \code{survminer}.
#' @param data data.frame or data.table containing survival data.
#' @param time the time interval from start of observation until date of event (e.g. disease progression or death)
#' or censoring.
#' @param status variable specifying if event occured or data has been censored.
#' @param var variable tested for Influence on outcome.
#' @param xlim x axis limit.
#' @param xlab X-axis label.
#' @param ylab y-axis label.
#' @param break.y.by interval of breaks on the y-axis.
#' @param break.time.by Interval of breaks on the x-axis (time axis).
#' @param ggtheme function, ggplot2 theme name. Default value is theme_eumelareg_surv_plot. Allowed values include ggplot2 official themes: see \code{theme}.
#' @param tables.theme function, ggplot2 theme name. Default value is theme_eumelareg_surv_plot. Allowed values include ggplot2 official themes: see \code{theme}.
#' @param axes.offset logical value. If TRUE the space between the plot origin and the axes is removed.
#' @inheritParams survminer::ggsurvplot
#' @param risk.table.title the title to be used for the risk table. Default is no title.
#' @param risk.table.title.size size of risk table title.
#' @param risk.table.text.size size of risk table text.
#' @param risk.table.width relative width of the risk table.
#' @param merge logical value. If TRUE survival curve and median survival table are plotted in the same graph. Else
#' two separate figures are generated. Default is FALSE.
#' @param plot.width relative width of the survival plot
#' @param plot.height relative height of the survival plot. The risk table is adjusted accordingly.
#' @param plot.margin.left numerical. Used to adjust the plot horizontally.
#' @param plot.axes.text.size text size of x-axis and y-axis text.
#' @param plot.axes.title.size text size of x-axis and y-axis titles.
#' @param legend.labs character vector specifying legend labels. Used to replace the names of the strata from the fit.
#' Should be given in the same order as those strata.
#' @param legend.position he position of legends ("none", "left", "right", "bottom", "top", or two-element numeric vector)
#' @param legend.size size of the legend title and legend text.
#' @param legend.title name of legend title.
#' @param pval.coord coords of pvalue within plot.
#' @param pval.size size of pvalue within plot
#' @param weights character variable specifying the name of the weights column. Weights have to be added to the original dataframe in order to be applied correctly.
#' @details Further arguments can be obtained from the [ggsurvplot()] function.
#' @export


survplot_eumelareg <- function (data, time = "time", status = "status", var = NULL,
                                xlab = "Time in months", ylab = "Probability of Survival", axes.offset = FALSE,
                                break.y.by = 0.1, break.time.by = 3,  xlim = c(0, 48),
                                ggtheme = theme_eumelareg_surv_plot(), tables.theme = theme_eumelareg_surv_table(),
                                plot.width = 0.838, plot.height = 0.7, plot.margin.left = 20,
                                plot.axes.text.size =12, plot.axes.title.size = 12, weights = NULL,
                                risk.table.width = 0.92, risk.table.title = NULL,
                                risk.table.title.size = 12, risk.table.text.size = 12,
                                legend.position = "top",legend.title = NULL, legend.labs = NULL, legend.size = 12,
                                pval = TRUE, pval.coord = c(1,0.1), pval.size = 12,
                                merge = FALSE, palette = "jco",  ...)
{

  ## data preprocessing, filter out data with missing values of the target variable and assign legend labels based on factor levels
  if(!is.null(var)){
    data <- data[!which(is.na(var))]
    if (is.null(legend.labs)) {
      legend.labs <- sort(unique(data[[var]]))
      legend.labs.risk.table <- gsub(">", "&gt;", legend.labs)
    } else {
      legend.labs.risk.table <- gsub(">", "&gt;", legend.labs)
      data[[var]] <- factor(data[[var]])
    }
  } else {
    if (is.null(legend.labs)) {
      legend.labs <- "All"
      legend.labs.risk.table <- ""
    }
  }
  # calculate height of table based on presence of risk table title
  table.height <- 1 - plot.height
  if (is.null(risk.table.title)) {
    risk.table.title <- ""
    table.height <- 1- plot.height + 0.05
  }

  ## fit survival
  if(!is.null(var)){fit <- surv_fit(Surv(eval(parse(text = time)), eval(parse(text = status))) ~ eval(parse(text = var)), data = data, weights = weights)
  } else {
    fit <- surv_fit(Surv(eval(parse(text = time)), eval(parse(text = status))) ~ 1, data = data, weights = weights)
    pval <- FALSE
  }


  ## plot survival curve
  ggsurv <- ggsurvplot(fit, data = data, xlab = xlab, ylab = ylab,
                       pval = pval, pval.size = pval.size/2.835,  xlim = xlim,
                       break.y.by = break.y.by, break.time.by = break.time.by,
                       ggtheme = ggtheme, tables.theme = tables.theme, axes.offset = axes.offset,
                       font.tickslab = plot.axes.text.size, font.x = plot.axes.title.size,
                       font.y = plot.axes.title.size, font.legend = legend.size,
                       legend.labs = legend.labs, legend.title = legend.title, palette = palette, pval.coord = pval.coord,...)

  ## adjust survival curve
  ggsurv$plot <- ggsurv$plot + theme(legend.position = legend.position,
                                     plot.margin = unit(c(5.5, 5.5, 5.5, plot.margin.left), "points"))

  ## draw risk table and remove unnecessary lines and text
  risk_table <- ggrisktable(fit, data = data, risk.table.title = risk.table.title, xlim = xlim,
                            fontsize = risk.table.text.size/2.835, break.time.by = break.time.by, # size/2.835 from points to mm
                            legend.labs = legend.labs.risk.table, ...) +
    theme(axis.line.y = element_blank(), axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), axis.line.x = element_blank(),
          axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(), plot.title = element_text(face = "bold", size = risk.table.title.size))

  ## calculate median survival and draw as table
  if(!is.null(var)){
    surv_med <- surv_median(fit)
    tbl <- as.data.frame(table(data[!is.na(eval(parse(text = time))),
                                    droplevels(eval(parse(text = var)))]))
    tbl$median <- sapply(1:length(surv_med$median), function(x) {
      paste(surv_med$median[x], " (", surv_med$lower[x],
            "-", surv_med$upper[x], ")", sep = "")
    })
    rownames(tbl) <- tbl$Var1
    tbl$Var1 <- NULL
    colnames(tbl) <- c("No. of patients", "Median  (95% CI)")
    tblGrob <- gridExtra::tableGrob(tbl, theme = gridExtra::ttheme_minimal())

  }

  ## define blank plot to adjust the ggarrange panel
  blankPlot <- ggplot() + geom_blank(aes(1, 1)) + theme(plot.background = element_blank(),
                                                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                        panel.border = element_blank(), panel.background = element_blank(),
                                                        axis.title.x = element_blank(), axis.title.y = element_blank(),
                                                        axis.text.x = element_blank(), axis.text.y = element_blank(),
                                                        axis.ticks = element_blank(), axis.line = element_blank())


  # arrange survplot and survtable
  p1 <- cowplot::ggdraw() +
    cowplot::draw_plot(risk_table, x = 0, y = 0, width = risk.table.width, height = table.height) +
    cowplot::draw_plot(ggsurv$plot, x = 0.04, y = 1 - plot.height, width = plot.width, height = plot.height)

  # arrange survplot+table with median survival table
  if(!is.null(var)){
    p2 <- ggpubr::ggarrange(tblGrob, blankPlot, nrow = 2)
    if (merge == TRUE) {
      ggpubr::ggarrange(p1, p2, ncol = 2, widths = c(2, 1))
    }
    else {
      list(plot = p1, table = ggpubr::ggarrange(tblGrob))
    }
  } else {
    return(p1)
  }
}




