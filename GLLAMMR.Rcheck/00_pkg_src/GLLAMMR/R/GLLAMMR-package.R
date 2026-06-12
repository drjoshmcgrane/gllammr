#' @keywords internal
#' @aliases GLLAMMR-package
"_PACKAGE"

#' @useDynLib GLLAMMR, .registration = TRUE
#' @importFrom stats lm glm coef fitted residuals predict simulate logLik vcov
#'   cooks.distance influence cor dpois lowess median na.fail optimize
#'   printCoefmat qqline qqnorm rexp rgamma shapiro.test update
#'   model.matrix model.frame model.response terms as.formula update.formula
#'   nlminb optim plogis qlogis pnorm qnorm dnorm rnorm rbinom rpois runif
#'   sd var quantile aggregate na.omit complete.cases setNames pchisq qchisq
#'   AIC BIC binomial gaussian poisson optimHess
#' @importFrom graphics plot lines points legend abline text par hist barplot
#'   axis box grid mtext matplot polygon segments boxplot curve image rect rug
#' @importFrom grDevices rainbow adjustcolor colorRampPalette
#' @importFrom utils modifyList head tail packageVersion data
#' @importFrom methods is
NULL

# Dataset names loaded via utils::data() inside the validation suite
utils::globalVariables(c("LSAT", "Science", "carcinoma", "grouseticks",
                         "sleepstudy", "toenail", "wine"))
