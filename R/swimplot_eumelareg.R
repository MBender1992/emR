#' Swimmer plot in EuMelaReg style
#'
#' This function plots a default EuMelaReg swimmerplot using \code{ggplot}.
#' @param data data.frame or data.table containing response data. Each event of response has to be defined in an own row.
#' Data for Response type, durable response, ongoing treatment and treatmend end have to entered in the respective columns.
#' @param ID column name that contains subject IDs.
#' @param end column name that contains time after which the subject died or was lost to follow-up.
#' @param response.start column that contains time after index date (e.g. diagnosis, treatment start) with documented response.
#' @param response.end column that contains time after index date with documented end of response.
#' @param response.type column that contains information whether patient had complete or partial response.
#' @param continued.response column that contains (logical) value whether response is ongoing or not.
#' @param durable.response column containing (logical) value whether response was durable (equal to or longer than six months).
#' @param strat grouping factor. Color of bars are colored accordingly.
#' @param stack.groups logical value indicating whether subjects should only be sorted by time or also stacked into groups defined by the strat argument. Default is TRUE.
#' @param symbol.size size of symbols indicating response types
#' @param xlab label of x-axis
#' @param ylab label of y-axis
#' @export

swimplot_eumelareg <- function(data, ID, end, response.start = NULL, response.end = NULL, response.type = NULL, stack.groups = TRUE,
                               continued.response = NULL, durable.response = NULL, strat = NULL, symbol.size = 5,
                               xlab = NULL, ylab = NULL){
  responseType <- NULL
  # define symbols for response types
  unicode <-  list(triangle=sprintf('\u25B2'), circle=sprintf('\u25CF'), square=sprintf('\u25A0'), arrow=sprintf('\u2794'))

  if(!is.null(response.end) & !is.null(continued.response)){
    data[[response.end]] <- ifelse(is.na(data[[response.end]] & data[[continued.response]] == FALSE), data[[end]], data[[response.end]])
  }
  # define name of column names
  column_names <- c(ID, "time", "responseType")

  # define df for durable responses
  if(!is.null(durable.response)){
    durable <- data[data[[durable.response]] == TRUE, c(ID, durable.response)]
    durable$responseType <- "Durable"
    colnames(durable) <- column_names
    durable$time <- -0.25
  } else durable <- NULL

  # define df for continued responses
  if(!is.null(continued.response)){
    continued <- data[data[[continued.response]] == TRUE, c(ID, end)]
    continued$responseType <- "Continued Treatment"
    continued$End <- continued$End + 0.25 # adjust arrow position
    colnames(continued) <- column_names
  } else continued <- NULL

  # define df for CR and PR Start
  if(!is.null(response.start) & !is.null(response.type)){
    response_start <- data[!is.na(data[[response.start]]), c(ID, response.type, response.start)]
    response_start$responseType <- paste(response_start$Response, "start")
    response_start <- response_start[,c(1,3,4)]
    colnames(response_start) <- column_names
  } else response_start <- NULL

  # define df for CR and PR End
  if(!is.null(response.end)){
    response_end <- data[!is.na(data[[response.end]]), c(ID, response.end)]
    response_end$responseType <- "Response end"
    colnames(response_end) <- column_names
  } else response_end <- NULL

  # bind all dfs into one
  dat_shapes <- rbind(response_start, response_end, durable, continued)
  responseLevels <- c("Complete response start", "Partial response start", "Response end", "Durable", "Continued Treatment")
  if(!is.null(response.start) & !is.null(response.type)){
    dat_shapes$responseType <-  factor(dat_shapes$responseType, levels=responseLevels)
    dat_shapes <- dplyr::arrange(dat_shapes, dplyr::desc(responseType))
  }

  # define data used for bars
  if(!is.null(strat)){
    dat_bar <- data[, c(ID, strat, end)]
  } else {
    dat_bar <- data[, c(ID,  end)]
  }

  dat_bar <- dplyr::distinct(dat_bar)
  dat_bar$SUBJID <- forcats::fct_reorder(.f=dat_bar[[ID]], .x=as.numeric(dat_bar[[end]]), .desc = FALSE)
  if(!is.null(strat) & stack.groups == TRUE){
    dat_bar$SUBJID <- forcats::fct_reorder2(.f=dat_bar[[ID]], .x= as.numeric(dat_bar[[end]]), .y = as.numeric(dat_bar[[strat]]), .desc = TRUE)
  }

  # define max value for ylimits
  max_val <- max(dat_bar[[end]])

  # plot
  if(!is.null(strat)){
    p <- ggplot(dat_bar,aes_string(x= ID, y = end)) +
      geom_bar(stat = "identity", aes_string(fill = strat), color = "black", width = 0.7) +
      ggsci::scale_fill_jco(alpha = 0.7)
  } else {
    p <- ggplot(dat_bar,aes_string(x= ID, y = end)) +
      geom_bar(stat = "identity", fill = "grey", color = "black", width = 0.7)
  }

  p <- p +
    coord_flip() +
    scale_colour_manual(values = c("Complete response start" = if(!is.null(response.start) & !is.null(response.type)) "#E41A1C",
                                   "Partial response start" = if(!is.null(response.start) & !is.null(response.type)) "#377EB8",
                                   "Response end" = if(!is.null(response.end)) "black",
                                   "Durable" = if(!is.null(durable.response)) "black",
                                   "Continued Treatment" =  if(!is.null(continued.response)) "black")) +

    scale_shape_manual(values = c("Complete response start" =if(!is.null(response.start) & !is.null(response.type)) unicode[["triangle"]],
                                  "Partial response start" = if(!is.null(response.start) & !is.null(response.type)) unicode[["triangle"]],
                                  "Response end" = if(!is.null(response.end)) unicode[["circle"]],
                                  "Durable" = if(!is.null(durable.response)) unicode[["square"]],
                                  "Continued Treatment" = if(!is.null(continued.response)) unicode[["arrow"]])) +
    scale_y_continuous(limits = c(-0.5, max_val + 3.5), breaks = seq(0, max_val + 3, by = 3)) +
    labs(fill = strat, colour = "Symbol Key", shape = "Symbol Key",  x = xlab, y = ylab, title = "Swimmer Plot",
         caption ="Durable defined as subject with six months or more of confirmed response") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(size=7, hjust=0),
          panel.grid = element_blank())

  if(!is.null(dat_shapes)){
    p <- p + geom_point(data=dat_shapes, aes_string(ID, "time", colour="responseType", shape="responseType"), size= symbol.size)
  }
  print(p)

}


