#' @keywords internal
#' @aliases GLLAMMR-package
"_PACKAGE"

#' @useDynLib GLLAMMR, .registration = TRUE
#' @importFrom stats lm glm coef fitted residuals predict simulate logLik vcov
#'   model.matrix model.frame model.response terms as.formula update.formula
#'   nlminb optim plogis qlogis pnorm qnorm dnorm rnorm rbinom rpois runif
#'   sd var quantile aggregate na.omit complete.cases setNames pchisq qchisq
#'   AIC BIC binomial gaussian poisson optimHess
#' @importFrom graphics plot lines points legend abline text par hist barplot
#'   axis box grid mtext matplot polygon segments
#' @importFrom grDevices rainbow adjustcolor
#' @importFrom utils modifyList head tail packageVersion
#' @importFrom methods is
NULL
