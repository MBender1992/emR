#' Survival curve in EuMelaReg style
#'
#' This function plots a default EuMelaReg survival plot (Kaplan-Meier plot) produced with \code{survminer}.
#' @param fit A survfit object generated with the \code{survfit()} function from the \code{survival} package.
#' @param data data.frame or data.table containing the same data as used in the \code{survfit()} function.
#' @param time The time interval from start of observation until date of event (e.g. disease progression or death)
#' or censoring.
#' @param status Variable specifying if event occured or data has been censored.
#' @param var Variable tested for Influence on outcome.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param break.y.by Interval of breaks on the y-axis.
#' @param break.time.by Interval of breaks on the x-axis (time axis).
#' @param ggtheme function, ggplot2 theme name. Default value is theme_eumelareg_surv_plot. Allowed values include ggplot2 official themes: see \code{theme}.
#' @param tables.theme function, ggplot2 theme name. Default value is theme_eumelareg_surv_plot. Allowed values include ggplot2 official themes: see \code{theme}.
#' @param axes.offset Logical value. If TRUE the space between the plot origin and the axes is removed.
#' @inheritParams survminer::ggsurvplot
#' @param risk.table.y.text Logical value. Default value is TRUE. If FALSE, the y axis tick labels of tables will be hidden.
#' @param risk.table.title The title to be used for the risk table
#' @param table.margin.left Numerical. Used to adjust the risk table horizontally.
#' @param legend.labs Character vector specifying legend labels. Used to replace the names of the strata from the fit.
#' Should be given in the same order as those strata.
#' @examples
#' library(survival)
#' data("lung")
#' lung$sex <- ifelse(lung$sex == 1, "Male", "Female")
#' fit <- survfit(Surv(time, status) ~ sex, data = lung)
#' survplot_eumelareg(fit, lung, var = "sex",table.margin.left = 0,  legend.labs = c("Male", "Female"))
#' @seealso [ggsurvplot()]
#' @export

survplot_eumelareg <- function(fit,data,time = "time", status = "status", var, xlab = "Time in months",
                               ylab = "Probability of Overall Survival",pval = TRUE, break.y.by = 0.1,
                               break.time.by = 3, ggtheme = theme_eumelareg_surv_plot(),
                               tables.theme = theme_eumelareg_surv_table(), axes.offset = TRUE,
                               risk.table = "absolute", risk.table.y.text = TRUE,risk.table.title = "No. at Risk",
                               table.margin.left = 0, legend.labs = NULL, palette = "RdYlBu"){

  # plot survival curve
  ggsurv <- ggsurvplot(fit,data = data, xlab = xlab,  ylab = ylab, pval = pval,
                       break.y.by = break.y.by, break.time.by = break.time.by, ggtheme = ggtheme,
                       tables.theme = tables.theme, axes.offset = axes.offset, risk.table = risk.table,
                       risk.table.y.text = risk.table.y.text, risk.table.title = risk.table.title,
                       legend.labs =  legend.labs, palette =  palette)

  # adjust position of risk table
  ggsurv$table <- ggsurv$table +
    theme(plot.margin = unit(c(5.5, 5.5, 5.5, table.margin.left), "points"))

  # define table with Median survival (displayed on the right of the figure)
  surv_med <- surv_median(fit)

  tbl <- as.data.frame(table(data[[var]]))
  tbl$median <- sapply(1:length(surv_med$median),function(x){
    paste(surv_med$median[x], " (", surv_med$lower[x],"-", surv_med$upper[x],")", sep = "")
  })
  rownames(tbl) <- tbl$Var1
  tbl$Var1 <- NULL
  colnames(tbl) <- c("No. of patients", "Median  (95% CI)")
  tblGrob <- gridExtra::tableGrob(tbl, theme = gridExtra::ttheme_minimal())

  # calculate cox regression model
  # cox.fit <- coxph(as.formula(paste("Surv(",time,",",status,") ~ ", var)), data =  data)
  # res <- summary(cox.fit)
  #
  # # extract coefficient and 95% CI from cox regression
  # HR <- res$coefficients[2]
  # lower <- res$conf.int[3]
  # upper <- res$conf.int[4]
  #
  # # define text below table on the right of the figure
  # HR_string <- paste("Hazard Ratio for \n Disease Progression or Death ", format(round(HR,2), nsmall = 2),
  #                    "\n (95% CI, ", format(round(lower,2), nsmall = 2),"-",
  #                    format(round(upper,2), nsmall = 2),")", sep = "")
  # p_string <- ifelse(res$waldtest[3] < 0.001, paste("Wald test: P<0.001", sep = ""),
  #                    paste("Wald test: P=",round(res$waldtest[3],3), sep = ""))
  # txt <- paste(HR_string, p_string, sep = "\n")

  # add blank plot for arranging
  blankPlot <- ggplot()+geom_blank(aes(1,1))+
    theme(plot.background = element_blank(), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),
          panel.border = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(),
          axis.title.y = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(),
          axis.ticks = element_blank(), axis.line = element_blank()
    )

  # arrange plot, table and text
  p1 <- ggpubr::ggarrange(ggsurv$plot, ggsurv$table, ncol = 1, heights = c(3,1))
  # p2 <- ggpubr::ggarrange(tblGrob,blankPlot,ggpubr::text_grob(txt), ncol = 2, nrow = 4)
  p2 <- ggpubr::ggarrange(tblGrob,blankPlot,  nrow = 2)
  ggpubr::ggarrange(p1, p2, ncol = 2, widths = c(2,1))
}


