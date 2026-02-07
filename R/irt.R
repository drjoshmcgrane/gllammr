#' Fit Item Response Theory Models
#'
#' Fit Rasch, 2PL, or 3PL IRT models using GLLAMM framework
#'
#' @param response_matrix Matrix of item responses (persons x items), coded 0/1
#' @param model Type of IRT model: "Rasch", "2PL", or "3PL"
#' @param start Optional starting values
#' @param control Control parameters for optimization
#'
#' @return An object of class \code{gllamm_irt}
#'
#' @examples
#' \dontrun{
#' # Simulate Rasch data
#' set.seed(123)
#' n_persons <- 500
#' n_items <- 20
#' theta <- rnorm(n_persons, 0, 1)
#' difficulty <- rnorm(n_items, 0, 1)
#'
#' # Generate responses
#' responses <- matrix(NA, n_persons, n_items)
#' for (i in 1:n_persons) {
#'   for (j in 1:n_items) {
#'     p <- plogis(theta[i] - difficulty[j])
#'     responses[i, j] <- rbinom(1, 1, p)
#'   }
#' }
#'
#' # Fit Rasch model
#' fit_rasch <- fit_irt(responses, model = "Rasch")
#' summary(fit_rasch)
#'
#' # Fit 2PL model
#' fit_2pl <- fit_irt(responses, model = "2PL")
#' summary(fit_2pl)
#' }
#'
#' @export
fit_irt <- function(response_matrix, model = c("Rasch", "2PL", "3PL"),
                    start = NULL, control = list()) {

  model <- match.arg(model)
  model_code <- switch(model, Rasch = 1L, "2PL" = 2L, "3PL" = 3L)

  # Convert to long format
  n_persons <- nrow(response_matrix)
  n_items <- ncol(response_matrix)

  y_long <- as.vector(t(response_matrix))
  person_id <- rep(1:n_persons, each = n_items) - 1L  # 0-indexed
  item_id <- rep(1:n_items, times = n_persons) - 1L   # 0-indexed

  # Remove missing values
  complete_cases <- !is.na(y_long)
  y_long <- y_long[complete_cases]
  person_id <- person_id[complete_cases]
  item_id <- item_id[complete_cases]

  n_obs <- length(y_long)

  # Prepare TMB data
  tmb_data <- list(
    y = as.numeric(y_long),
    person_id = as.integer(person_id),
    item_id = as.integer(item_id),
    n_persons = as.integer(n_persons),
    n_items = as.integer(n_items),
    n_obs = as.integer(n_obs),
    model_type = as.integer(model_code)
  )

  # Initialize parameters
  if (is.null(start)) {
    # Initialize abilities to zero
    theta_init <- rep(0, n_persons)

    # Initialize difficulties from marginal item proportions
    item_means <- colMeans(response_matrix, na.rm = TRUE)
    difficulty_init <- -qlogis(pmin(pmax(item_means, 0.01), 0.99))

    # Initialize discriminations
    discrimination_init <- rep(1, n_items)

    # Initialize guessing parameters
    guessing_init <- rep(0.1, n_items)

    tmb_params <- list(
      theta = theta_init,
      difficulty = difficulty_init,
      discrimination = discrimination_init,
      guessing = guessing_init,
      log_sigma_theta = log(1.0)
    )
  } else {
    tmb_params <- start
  }

  # Create TMB object
  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    random = "theta",  # Integrate out person abilities
    DLL = "gllamm_irt",
    silent = TRUE
  )

  # Optimize
  control_defaults <- list(eval.max = 2000, iter.max = 1000, trace = 0)
  control <- modifyList(control_defaults, control)

  opt <- nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = control
  )

  # Get standard errors
  sdr <- try(TMB::sdreport(obj), silent = TRUE)

  # Extract parameters
  par_full <- obj$env$last.par.best

  difficulty_hat <- par_full[names(par_full) == "difficulty"]
  names(difficulty_hat) <- paste0("Item", 1:n_items)

  if (model_code >= 2) {
    discrimination_hat <- par_full[names(par_full) == "discrimination"]
    names(discrimination_hat) <- paste0("Item", 1:n_items)
  } else {
    discrimination_hat <- rep(1, n_items)
    names(discrimination_hat) <- paste0("Item", 1:n_items)
  }

  if (model_code == 3) {
    guessing_hat <- par_full[names(par_full) == "guessing"]
    names(guessing_hat) <- paste0("Item", 1:n_items)
  } else {
    guessing_hat <- NULL
  }

  theta_hat <- par_full[names(par_full) == "theta"]
  names(theta_hat) <- paste0("Person", 1:n_persons)

  sigma_theta_hat <- exp(par_full[names(par_full) == "log_sigma_theta"])

  # Construct result object
  result <- list(
    model = model,
    item_parameters = data.frame(
      difficulty = difficulty_hat,
      discrimination = discrimination_hat,
      guessing = if (!is.null(guessing_hat)) guessing_hat else NA
    ),
    person_abilities = theta_hat,
    ability_sd = sigma_theta_hat,
    logLik = -opt$objective,
    AIC = 2 * opt$objective + 2 * length(obj$par),
    BIC = 2 * opt$objective + log(n_persons) * length(obj$par),
    convergence = list(
      converged = (opt$convergence == 0),
      message = opt$message
    ),
    n_persons = n_persons,
    n_items = n_items,
    tmb_obj = obj,
    tmb_opt = opt,
    tmb_sdr = sdr
  )

  class(result) <- c("gllamm_irt", "gllamm")

  return(result)
}


#' @export
print.gllamm_irt <- function(x, ...) {
  cat("IRT Model (", x$model, ")\n\n", sep = "")
  cat("Number of persons:", x$n_persons, "\n")
  cat("Number of items:", x$n_items, "\n\n")

  cat("Item parameters:\n")
  print(round(x$item_parameters, 3))

  cat("\nAbility distribution:\n")
  cat("  Mean:", mean(x$person_abilities), "\n")
  cat("  SD:", sd(x$person_abilities), "\n")
  cat("  Estimated SD:", x$ability_sd, "\n\n")

  cat("Log-likelihood:", round(x$logLik, 2), "\n")
  cat("AIC:", round(x$AIC, 2), "\n")
  cat("BIC:", round(x$BIC, 2), "\n")

  invisible(x)
}


#' @export
summary.gllamm_irt <- function(object, ...) {
  print(object)

  cat("\nAbility quartiles:\n")
  print(quantile(object$person_abilities, c(0, 0.25, 0.5, 0.75, 1)))

  invisible(object)
}
