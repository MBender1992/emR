#' Convert table1 into data.frame
#'
#' @param x table1 object
#' @param ... additional arguments (see table1 documentation)
#' @export

as_data_frame_table1 <- function(x, ...) {
  obj <- attr(x, "obj")
  with(obj, {
    rlh <- if (is.null(rowlabelhead) || rowlabelhead=="") " " else rowlabelhead
    z <- lapply(contents, function(y) {
      y <- as.data.frame(y, stringsAsFactors=F)
      y2 <- data.frame(x=paste0(c("", rep("  ", nrow(y) - 1)), rownames(y)), stringsAsFactors=F)
      y <- cbind(setNames(y2, rlh), y)
      y
    })
    df <- do.call(rbind, z)
    df <- rbind(c("", ifelse(is.na(headings[2,]), "", sprintf("(N=%s)", headings[2,]))), df)
    colnames(df) <- c(rlh, headings[1,])
    rownames(df) <- NULL
    noquote(df)
  })
}
