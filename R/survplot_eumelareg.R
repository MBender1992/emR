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
#' @param conf.type Method to calculate confidence intervals. Log-log method is the default in SAS.
#' @param text.size size of plot text.
#' @param risk.table.width relative width of the risk table.
#' @param merge logical value. If TRUE survival curve and median survival table are plotted in the same graph. Else
#' two separate figures are generated. Default is FALSE.
#' @param plot.width relative width of the survival plot
#' @param plot.height relative height of the survival plot. The risk table is adjusted accordingly.
#' @param plot.margin.left numerical. Used to adjust the plot horizontally.
#' @param legend.labs character vector specifying legend labels. Used to replace the names of the strata from the fit.
#' Should be given in the same order as those strata.
#' @param landmarks times for which landmark survival data shall be provided.
#' Should be given in the same order as those strata.
#' @param legend.position he position of legends ("none", "left", "right", "bottom", "top", or two-element numeric vector)
#' @param legend.title name of legend title.
#' @param pval.coord coords of pvalue within plot.
#' @param weights character variable specifying the name of the weights column. Weights have to be added to the original dataframe in order to be applied correctly.
#' @param extract.legend logical indicating whether only the legend should be extracted from the plot to build manual arrangements of plots and legend with [ggarrange()]
#' @details Further arguments can be obtained from the [ggsurvplot()] function.
#' @export

survplot_eumelareg <- function (data, time = "time", status = "status", var = NULL,
                                xlab = "Time in months", ylab = "Probability of Survival", axes.offset = FALSE,
                                break.y.by = 0.1, break.time.by = 3,  xlim = c(0, 48),
                                ggtheme = theme_eumelareg_surv_plot(), tables.theme = theme_eumelareg_surv_table(),
                                plot.width = 0.838, plot.height = 0.88, plot.margin.left = NULL,
                                text.size = 12, weights = NULL, landmarks = c(12,24), conf.type = "log-log",
                                risk.table.width = 0.9, risk.table.title = NULL,
                                legend.position = "top",legend.title = "", legend.labs = NULL,
                                pval = TRUE, pval.coord = c(1,0.1), merge = FALSE, palette = "jco", extract.legend = FALSE,  ...)
{

  ## data preprocessing, filter out data with missing values of the target variable and assign legend labels based on factor levels
  if(!is.null(var)){
    if (is.data.table(data)) {
      data <- data[which(!is.na(data[[var]]))]
    }
    else {
      data <- data[which(!is.na(data[[var]])), ]
    }
    levels <- levels(data[[var]])
   if(is.null(plot.margin.left)) plot.margin.left <- 12 * sqrt(max(nchar(levels)))
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
  if(!is.null(var)){
    fit <- surv_fit(Surv(eval(parse(text = time)), eval(parse(text = status))) ~ eval(parse(text = var)), data = data, weights = weights, conf.type = conf.type)
    if(!is.null(weights)){
      fit_table <- surv_fit(Surv(eval(parse(text = time)), eval(parse(text = status))) ~ eval(parse(text = var)), data = data, weights = NULL, conf.type = conf.type)
      pval <- FALSE
    } else {
      fit_table <- fit
    }

  } else {
    fit <- surv_fit(Surv(eval(parse(text = time)), eval(parse(text = status))) ~ 1, data = data, weights = weights, conf.type = conf.type)
    pval <- FALSE
  }

  ## plot survival curve
  ggsurv <- ggsurvplot(fit, data = data, xlab = xlab, ylab = ylab,
                       pval = pval, pval.size = text.size/2.835,  xlim = xlim,
                       break.y.by = break.y.by, break.time.by = break.time.by,
                       ggtheme = ggtheme, tables.theme = tables.theme, axes.offset = axes.offset,
                       font.tickslab = text.size, font.x = text.size,
                       font.y = text.size, font.legend = text.size,
                       legend.labs = legend.labs, legend.title = legend.title, palette = palette, pval.coord = pval.coord, ...)

  ## adjust survival curve position
  ggsurv$plot <- ggsurv$plot + theme(legend.position = legend.position,
                                     plot.margin = unit(c(5.5, 5.5, 5.5, plot.margin.left), "points"))

  # calculate and display propensity score weighted pvalue
  if(!is.null(weights)){
    dat_logrank <- data[!is.na(data[[time]])]
    pval.ipsw <- logrank_IPSW_RISCA(dat_logrank[[time]], dat_logrank[[status]], ifelse(dat_logrank[[var]] == levels(dat_logrank[[var]])[1], 0, 1), weights[!is.na(data[[time]])])$p.value

    if(pval.ipsw < 0.0001){
      pval.ipsw <- "< 0.0001"
      xjust <- 2.8*(text.size+1)/12
      plabel <- paste("p ", pval.ipsw, sep = "")
    } else if (pval.ipsw < 0.001 & pval.ipsw >= 0.0001){
      pval.ipsw <- format(pval.ipsw, scientific = F, digits = 1)
      xjust <- 2.8*(text.size+1)/12
      plabel <- paste("p = ", pval, sep = "")
    }else {
      pval.ipsw <- round(pval.ipsw, 3)
      xjust <- 2.4*(text.size+1)/12
      plabel <- paste("p = ", pval.ipsw, sep = "")
    }
    ggsurv$plot <- ggsurv$plot + annotate(geom = "text", x = pval.coord[1]+xjust, y = pval.coord[2], label = plabel, size = text.size/2.835)

  }

  ## draw risk table and remove unnecessary lines and text
  risk_table <- ggrisktable(fit_table, data = data, risk.table.title = risk.table.title, xlim = c(0, xlim[2]-1),
                            fontsize = text.size/2.835, break.time.by = break.time.by, # size/2.835 from points to mm
                            legend.labs = legend.labs.risk.table, ...) +
    theme(axis.line.y = element_blank(), axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), axis.line.x = element_blank(),
          axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(), plot.title = element_text(face = "bold", size = text.size))

  ## calculate median survival and draw as table
  if(!is.null(var)){
    surv_med <- surv_median(fit)
    if(is.data.table(data)){
      tbl <- as.data.frame(table(data[!is.na(eval(parse(text = time))),
                                      if(is.factor(data[[var]])) droplevels(eval(parse(text = var)))]))
    } else {
      tbl <- as.data.frame(table(data[!is.na(data[[time]]), var]))
    }

    tbl$median <- sapply(1:length(surv_med$median), function(x) {
      paste(round(surv_med$median[x],2), " (", round(surv_med$lower[x],2),
            "-", round(surv_med$upper[x],2), ")", sep = "")
    })
    rownames(tbl) <- tbl$Var1
    tbl$Var1 <- NULL
    colnames(tbl) <- c("No. of patients", "Median  (95% CI)")
    tblGrob <- gridExtra::tableGrob(tbl, theme = gridExtra::ttheme_minimal())

  }

  if(!is.null(var)){
    lndmrk <- lapply(landmarks, function(t){
      survival_time(data = data, time = time, status = status, var = var, times = t)
    })
    lndmrk <- rlist::list.cbind(lndmrk)
    lndmrkGrob <- gridExtra::tableGrob(lndmrk, theme = gridExtra::ttheme_minimal())
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

  # functionality to only extract the legend of the plot
  if(extract.legend == TRUE){
    legend <- ggpubr::get_legend(ggsurv$plot)
    return(legend)
  }

  # arrange survplot+table with median survival table
  if(!is.null(var)){
    p2 <- ggpubr::ggarrange(tblGrob,lndmrkGrob, blankPlot, nrow =3)
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


