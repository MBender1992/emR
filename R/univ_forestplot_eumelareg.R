#' Forest Plot for Cox Proportional Hazards Model in EuMelaReg style
#'
#' This code generates a forest plot from a coxph model.
#' @inheritParams survminer::ggforest
#' @param time the time interval from start of observation until date of event (e.g. disease progression or death)
#' or censoring.
#' @param status variable specifying if event occured or data has been censored.
#' @param varnames Character vector specifying rownames of the table (empty columns should be named with "").
#' @param vars variables that were used in the coxph model
#' @param point.size Size of mean points.
#' @param line.size Size of errorbar line.
#' @param vjust_text vertical adjustment of text containing information about events, global pvalue, AIC and concordance index
#' @param y_breaks argument to supply manual y_breaks as a numerical vector. Default is NULL and breaks are set automatically within the function.
#' @param ylim argument to supply manual y limits as numerical vector of length 2. Default is NULL and limits are set automatically within the function.
#' @export

univ_forestplot_eumelareg <- function (data, time, status, vars, main = "Hazard ratio for disease progression or death (95% CI)", y_breaks = NULL,
                                  cpositions = c(0, 0.1, 0.3), point.size = 3, fontsize = 0.8,line.size = 0.7, vjust_text = 1.2,
                                  refLabel = "reference", noDigits = 2, varnames = NULL, ylim = NULL){

    conf.high <- conf.low <- estimate <- var <-  NULL

    data <- as.data.frame(data)
    ls <- lapply(1:length(vars), function(x){
      tmp <- cox_output(vars[x], data = data, time = time, status = status, output = "fit")
      terms <- attr(tmp$terms, "dataClasses")[-1]
      coef <- as.data.frame(broom::tidy(tmp, conf.int = TRUE))
      list(terms = terms, coef = coef)
    })

    terms <- unlist(lapply(1:length(ls), function(x){ls[[x]]$terms}))
    coef <- do.call(rbind, lapply(1:length(ls), function(x){ls[[x]]$coef}))

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


    allTermsDF <- do.call(rbind, allTerms)
    colnames(allTermsDF) <- c("var", "level", "N","pos")
    inds <- apply(allTermsDF[, 1:2], 1, paste0, collapse = "")
    rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
    toShow <- cbind(allTermsDF, coef[inds, ])[, c("var",
                                                  "level", "N", "p.value", "estimate",
                                                  "conf.low", "conf.high", "pos")]
    if (!is.null(varnames)) toShow$var <- varnames


    toShowExp <- toShow[, 5:7]
    toShowExp[is.na(toShowExp)] <- 0
    toShowExp <- format(exp(toShowExp), digits = noDigits)
    toShowExpClean <- data.frame(toShow, pvalue = signif(toShow[,4], noDigits + 1), toShowExp)
    toShowExpClean$stars <- paste0(round(toShowExpClean$p.value,
                                         noDigits + 1), " ", ifelse(toShowExpClean$p.value <  0.05, "*", ""),
                                   ifelse(toShowExpClean$p.value <   0.01, "*", ""),
                                   ifelse(toShowExpClean$p.value < 0.001, "*", ""))
    toShowExpClean$ci <- paste0("(", toShowExpClean[, "conf.low.1"],
                                " - ", toShowExpClean[, "conf.high.1"], ")")
    toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] = refLabel
    toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] = "<0.001 ***"
    toShowExpClean$stars[is.na(toShowExpClean$estimate)] = ""
    toShowExpClean$ci[is.na(toShowExpClean$estimate)] = ""
    toShowExpClean$estimate[is.na(toShowExpClean$estimate)] = 0
    toShowExpClean$var = as.character(toShowExpClean$var)
    toShowExpClean$var[duplicated(toShowExpClean$var)] = ""
    toShowExpClean$N <- paste0("(N=", toShowExpClean$N, ")")
    toShowExpClean$levelN <- paste(toShowExpClean$level, toShowExpClean$N) #neu
    toShowExpClean$estimateCI <- paste(toShowExpClean$estimate.1, toShowExpClean$ci) # neu
    toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1,]
    toShowExpClean$estimate <- ifelse(toShowExpClean$estimate == 0, NA, toShowExpClean$estimate)
    toShowExpClean <- toShowExpClean[toShowExpClean$var != "(weights)",] #
    rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high, na.rm = TRUE)
    if(is.null(y_breaks)) breaks <- grDevices::axisTicks(rangeb/2, log = TRUE, nint = 7) else breaks <- y_breaks
    rangeplot <- rangeb
    rangeplot[1] <- rangeplot[1] - diff(rangeb)
    rangeplot[2] <- rangeplot[2] + 0.15 * diff(rangeb)
    if (!is.null(ylim)) {
      rangeplot <- log(ylim)
      toShowExpClean$conf.high <- ifelse(log(ylim[2]) < toShowExpClean$conf.high, NA, toShowExpClean$conf.high)
      toShowExpClean$conf.low <- ifelse(log(ylim[1]) > toShowExpClean$conf.low, NA, toShowExpClean$conf.low)
    }
    width <- diff(rangeplot)
    y_variable <- rangeplot[1] + cpositions[1] * width
    y_nlevel <- rangeplot[1] + cpositions[2] * width
    y_cistring <- rangeplot[1] + cpositions[3] * width
    y_stars <- rangeb[2]
    x_annotate <- seq_len(nrow(toShowExpClean))
    annot_size_mm <- fontsize * as.numeric(grid::convertX(unit(theme_get()$text$size,"pt"), "mm"))

    p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) +
      geom_rect(aes(xmin = seq_along(var) - 0.5, xmax = seq_along(var) +
                      0.5, ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
                    fill = ordered(seq_along(var)%%2 + 1))) +
      # color of the rectangles
      scale_fill_manual(values = c("#FFFFFF33","grey95"), guide = "none") +
      # show confidence intervals
      geom_errorbar(aes(ymin = exp(conf.low), ymax = exp(conf.high)),size = line.size, width =0) +
      # plot mean points
      geom_point(pch = 16, size = point.size, color = "#009AA6") +
      # add no effect line at 1
      geom_hline(yintercept = 1, linetype = 2) +
      coord_flip(ylim = exp(rangeplot)) +
      ggtitle(main) +
      theme_light() +
      theme(panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            legend.position = "none",
            panel.border = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      xlab("") +
      annotate(geom = "text", x = x_annotate,
               y = exp(y_variable), label = toShowExpClean$var, fontface = "bold",
               hjust = 0, size = annot_size_mm) +
      annotate(geom = "text", x = x_annotate, y = exp(y_nlevel), hjust = 0, label = toShowExpClean$levelN,
               # vjust = -0.1,
               size = annot_size_mm) +
      # Annotate mean HR
      annotate(geom = "text",  x = x_annotate, y = exp(y_cistring), label = toShowExpClean$estimateCI,
               size = annot_size_mm) +
      # Annotate stars
      annotate(geom = "text", x = x_annotate, y = if(!is.null(ylim)) ylim[2]-0.4*ylim[2] else exp(y_stars),
               label = toShowExpClean$stars, size = annot_size_mm,
               hjust = -0.2, fontface = "italic") #+

    if(!is.null(y_breaks)){
      p <- p + scale_y_log10(name = "", expand = c(0.02, 0.02), breaks = y_breaks)
    } else {
      p <- p + scale_y_log10(name = "", labels = sprintf("%g", breaks), expand = c(0.02, 0.02), breaks = breaks)
    }

    gt <-suppressWarnings(ggplot_gtable(ggplot_build(p)))
    gt$layout$clip[gt$layout$name == "panel"] <- "off"
    ggpubr::as_ggplot(gt)
}



