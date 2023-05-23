#' Forest Plot for Cox Proportional Hazards Model in EuMelaReg style
#'
#' This code generates a forest plot from a coxph model.
#' @inheritParams survminer::ggforest
#' @param fit an object of class coxph or a list containing a mipo object and a list with number patients for each factor level
#' @param varnames Character vector specifying rownames of the table (empty columns should be named with "").
#' @param vars variables that were used in the coxph model
#' @param point_size Size of mean points.
#' @param line_size Size of errorbar line.
#' @param vjust_text vertical adjustment of text containing information about events, global pvalue, AIC and concordance index
#' @param subset Display only specified variables as a subset of cox regression. If multivariable cox model was calculated the results shown are adjusted for all input variables but only the specified subset will be displayed.
#' @param y_breaks argument to supply manual y_breaks as a numerical vector. Default is NULL and breaks are set automatically within the function.
#' @param ylim argument to supply manual y limits as numerical vector of length 2. Default is NULL and limits are set automatically within the function.
#' @export

forestplot_eumelareg <- function (fit, data = NULL, vars = NULL, main = "Hazard ratio for disease progression or death (95% CI)", y_breaks = NULL,
                                  cpositions = c(0, 0.1, 0.3), point_size = 4, fontsize = 1,line_size = 0.9, vjust_text = 1.2, subset = NULL,
                                  refLabel = "reference", noDigits = 2, varnames = NULL, ylim = NULL){

  conf.high <- conf.low <- estimate <- var <-  NULL

  if(any(class(fit) %in% "list")){
    model <- fit$fit
  } else {
    model <- fit
  }
  insight::get_data(model)
  if(any(class(model) %in% "mipo.summary")){
    if(is.null(data) | is.null(vars)) stop("Please provide data and variables argument.")
    data <- as.data.frame(data)
    terms <- as.character(utils::tail(data.frame(rbind(lapply(data, class))[, vars]), n = 1))
    names(terms) <- vars
    coef <- model
    message("Using results from multiple imputation.")
  } else {
    stopifnot(inherits(model, "coxph"))
    data <- as.data.frame(survminer:::.get_data(model, data = data))
    terms <- attr(model$terms, "dataClasses")[-1]
    coef <- as.data.frame(broom::tidy(model, conf.int = TRUE))
    gmodel <- broom::glance(model)
  }
  if(!is.null(fit$nfit)){
    allTerms <- fit$nfit
  } else {
    allTerms <- lapply(seq_along(terms), function(i) {
      var <- names(terms)[i]
      if (terms[i] %in% c("factor", "character")) {
        adf <- as.data.frame(table(data[, var]))
        cbind(var = var, adf, pos = 1:nrow(adf))
      }
      else if (terms[i] == "numeric") {
        data.frame(var = var, Var1 = "", Freq = nrow(data),
                   pos = 1)
      }
      else {
        vars = grep(paste0("^", var, "*."), coef$term,
                    value = TRUE)
        data.frame(var = vars, Var1 = "", Freq = nrow(data),
                   pos = seq_along(vars))
      }
    })
  }

  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N","pos")
  inds <- apply(allTermsDF[, 1:2], 1, paste0, collapse = "")
  rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
  toShow <- cbind(allTermsDF, coef[inds, ])[, c("var",
                                                "level", "N", "p.value", "estimate",
                                                "conf.low", "conf.high", "pos")]
    # only display subset of variables
  if(!is.null(subset)) toShow <- toShow[toShow$var %in% subset,]

  # plot the Forestplot
  forest_plotFUN(toShow = toShow, main = main,  y_breaks = y_breaks, cpositions = cpositions, point_size = point_size, varnames = varnames,
                 fontsize = fontsize, line_size = line_size, vjust_text = vjust_text, refLabel = refLabel, noDigits = noDigits, ylim =ylim)
  # toShowExp <- toShow[, 5:7]
  # toShowExp[is.na(toShowExp)] <- 0
  # toShowExp <- format(exp(toShowExp), digits = noDigits)
  # toShowExpClean <- data.frame(toShow, pvalue = signif(toShow[,4], noDigits + 1), toShowExp)
  # toShowExpClean$stars <- paste0(round(toShowExpClean$p.value,
  #                                      noDigits + 1), " ", ifelse(toShowExpClean$p.value <  0.05, "*", ""),
  #                                ifelse(toShowExpClean$p.value <   0.01, "*", ""),
  #                                ifelse(toShowExpClean$p.value < 0.001, "*", ""))
  # toShowExpClean$ci <- paste0("(", toShowExpClean[, "conf.low.1"],
  #                             " - ", toShowExpClean[, "conf.high.1"], ")")
  # toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] = refLabel
  # toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] = "<0.001 ***"
  # toShowExpClean$stars[is.na(toShowExpClean$estimate)] = ""
  # toShowExpClean$ci[is.na(toShowExpClean$estimate)] = ""
  # toShowExpClean$estimate[is.na(toShowExpClean$estimate)] = 0
  # toShowExpClean$var = as.character(toShowExpClean$var)
  # toShowExpClean$var[duplicated(toShowExpClean$var)] = ""
  # toShowExpClean$N <- paste0("(N=", toShowExpClean$N, ")")
  # toShowExpClean$levelN <- paste(toShowExpClean$level, toShowExpClean$N) #neu
  # toShowExpClean$estimateCI <- paste(toShowExpClean$estimate.1, toShowExpClean$ci) # neu
  # toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1,]
  # toShowExpClean$estimate <- ifelse(toShowExpClean$estimate == 0, NA, toShowExpClean$estimate)
  # toShowExpClean <- toShowExpClean[toShowExpClean$var != "(weights)",] #
  # rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high, na.rm = TRUE)
  # if(is.null(y_breaks)) breaks <- grDevices::axisTicks(rangeb/2, log = TRUE, nint = 7) else breaks <- y_breaks
  # rangeplot <- rangeb
  # rangeplot[1] <- rangeplot[1] - diff(rangeb)
  # rangeplot[2] <- rangeplot[2] + 0.15 * diff(rangeb)
  # if (!is.null(ylim)) {
  #   rangeplot <- log(ylim)
  #   toShowExpClean$conf.high <- ifelse(log(ylim[2]) < toShowExpClean$conf.high, NA, toShowExpClean$conf.high)
  #   toShowExpClean$conf.low <- ifelse(log(ylim[1]) > toShowExpClean$conf.low, NA, toShowExpClean$conf.low)
  # }
  # width <- diff(rangeplot)
  # y_variable <- rangeplot[1] + cpositions[1] * width
  # y_nlevel <- rangeplot[1] + cpositions[2] * width
  # y_cistring <- rangeplot[1] + cpositions[3] * width
  # y_stars <- rangeb[2]
  # x_annotate <- seq_len(nrow(toShowExpClean))
  # annot_size_mm <- fontsize * as.numeric(grid::convertX(unit(theme_get()$text$size,"pt"), "mm"))
  #
  # p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) +
  #   geom_rect(aes(xmin = seq_along(var) - 0.5, xmax = seq_along(var) +
  #                   0.5, ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
  #                 fill = ordered(seq_along(var)%%2 + 1))) +
  #   # color of the rectangles
  #   scale_fill_manual(values = c("#FFFFFF33","grey95"), guide = "none") +
  #   # show confidence intervals
  #   geom_errorbar(aes(ymin = exp(conf.low), ymax = exp(conf.high)), size = line_size, width =0) +
  #   # plot mean points
  #   geom_point(pch = 16, size = point_size, color = "#009AA6") +
  #   # add no effect line at 1
  #   geom_hline(yintercept = 1, linetype = 2) +
  #   coord_flip(ylim = exp(rangeplot)) +
  #   ggtitle(main) +
  #   theme_light() +
  #   theme(panel.grid.minor.y = element_blank(),
  #         panel.grid.minor.x = element_blank(),
  #         panel.grid.major.y = element_blank(),
  #         legend.position = "none",
  #         panel.border = element_blank(),
  #         axis.text.x = element_text(size = fontsize *13),
  #         axis.title.y = element_blank(),
  #         axis.text.y = element_blank(),
  #         axis.ticks.y = element_blank(),
  #         plot.title = element_text(size = fontsize *13, hjust = 0.5)) +
  #   xlab("") +
  #   annotate(geom = "text", x = x_annotate,
  #            y = exp(y_variable), label = toShowExpClean$var, fontface = "bold",
  #            hjust = 0, size = annot_size_mm) +
  #   annotate(geom = "text", x = x_annotate, y = exp(y_nlevel), hjust = 0, label = toShowExpClean$levelN,
  #            # vjust = -0.1,
  #            size = annot_size_mm) +
  #   # Annotate mean HR
  #   annotate(geom = "text",  x = x_annotate, y = exp(y_cistring), label = toShowExpClean$estimateCI,
  #            size = annot_size_mm) +
  #   # Annotate stars
  #   annotate(geom = "text", x = x_annotate, y = if(!is.null(ylim)) ylim[2]-0.4*ylim[2] else exp(y_stars),
  #            label = toShowExpClean$stars, size = annot_size_mm,
  #            hjust = -0.2, fontface = "italic") #+
  #
  # if(inherits(model, "coxph")){
  #   p <- p +  annotate(geom = "text", x = 0.5, y = exp(y_variable),
  #                      label = paste0("# Events: ",   gmodel$nevent, "; Global p-value (Log-Rank): ",
  #                                     format.pval(gmodel$p.value.log, eps = ".001"),
  #                                     " \nAIC: ", round(gmodel$AIC, 2), "; Concordance Index: ",
  #                                     round(gmodel$concordance, 2)), size = annot_size_mm,
  #                      hjust = 0, vjust = vjust_text, fontface = "italic")
  # }
  #
  #
  # if(!is.null(y_breaks)){
  #   p <- p + scale_y_log10(name = "", expand = c(0.02, 0.02), breaks = y_breaks)
  # } else {
  #   p <- p + scale_y_log10(name = "", labels = sprintf("%g", breaks), expand = c(0.02, 0.02), breaks = breaks)
  # }
  #
  # gt <-suppressWarnings(ggplot_gtable(ggplot_build(p)))
  # gt$layout$clip[gt$layout$name == "panel"] <- "off"
  # ggpubr::as_ggplot(gt)
}



