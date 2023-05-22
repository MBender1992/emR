#' Calculate bivariate residuals in a latent class analysis
#'
#' This code is adapted from
#' Author: Daniel Oberski
#' Date: 2017-08-01
#' Bivariate residual statistic for latent class analysis
#' Calculate the BVR for poLCA objects
#' @param fit fit object of class "poLCA" generated with the poLCA function.
#' @export

bvr <- function(fit) {
  stopifnot(class(fit) == "poLCA")

  ov_names <- names(fit$predcell)[1:(ncol(fit$predcell) - 2)]
  ov_combn <- combn(ov_names, 2)

  get_bvr <- function(ov_pair) {
    form_obs <- as.formula(paste0("observed ~ ", ov_pair[1], " + ", ov_pair[2]))
    form_exp <- as.formula(paste0("expected ~ ", ov_pair[1], " + ", ov_pair[2]))

    counts_obs <- xtabs(form_obs, data = fit$predcell)
    counts_exp <- xtabs(form_exp, data = fit$predcell)

    bvr <- sum((counts_obs - counts_exp)^2 / counts_exp)

    bvr
  }

  bvr_pairs <- apply(ov_combn, 2, get_bvr)
  names(bvr_pairs) <- apply(ov_combn, 2, paste, collapse = "<->")
  attr(bvr_pairs, "class") <- "dist"
  attr(bvr_pairs, "Size") <- length(ov_names)
  attr(bvr_pairs, "Labels") <- ov_names
  attr(bvr_pairs, "Diag") <- FALSE
  attr(bvr_pairs, "Upper") <- FALSE

  bvr_pairs
}
