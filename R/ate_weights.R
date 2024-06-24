#' Calculate average treatment effect using inverse propensity score weighting
#'
#' Propensity scores (PS) are calculated using fitted values obtained in a logistic regression. The inverse propensity score weighting is subsequently conducted
#' by 1/PS for the factor level that corresponds to 1 and by 1/(1-PS) for the factor level that corresponds to 0.
#' @param data Character vector specifying rownames of the table (empty columns should be named with "").
#' @param vars variables tested for Influence on outcome. NAs within vars should be replaced with a pseudocategory, e.g. "missing". A more sophisticated
#' approach with multiple imputation of missing values followed by propensity score calculation is provided with the \code{mi_PS} function of this package.
#' @param prop.var variable for which propensity scores should be calculated
#' @export
#' @examples
#' # example weights calculation with the lung dataset
#' df <- survival::lung
#' df$sex <- factor(df$sex)
#' df$ph.ecog <- ifelse(is.na(df$ph.ecog), "Missing", df$ph.ecog)
#' df$ph.karno <- ifelse(is.na(df$ph.karno), "Missing", df$ph.karno)
#' df$weights.ate <- ate_weights(data = df, vars = c("age", "ph.ecog", "ph.karno"), prop.var = "sex")

ate_weights <- function(data, vars, prop.var){
  data <- as.data.frame(data)
  vars <- vars[!vars %in% prop.var]
  vars_input <- paste(vars, collapse = " + ")
  formula <- as.formula(paste(prop.var,"~", vars_input, sep = ""))
  ps_model <- glm(formula, family = binomial, data = data)
  data$propensityScore <- predict(ps_model, type = "response")
 ifelse((data[[prop.var]] == levels(data[[prop.var]])[2]),(1/data$propensityScore), (1/(1-data$propensityScore)))
}


