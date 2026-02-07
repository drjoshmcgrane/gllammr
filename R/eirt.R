#' Fit Explanatory Item Response Theory Models
#'
#' Fit IRT models where item parameters are modeled as functions of item covariates
#'
#' @param response_matrix Matrix of item responses (persons x items)
#' @param item_data Data frame of item-level covariates (must have n_items rows)
#' @param difficulty_formula Formula for difficulty regression (e.g., ~ word_freq + length)
#' @param discrimination_formula Formula for discrimination regression (e.g., ~ item_type)
#' @param model IRT model type: "Rasch", "2PL", or "GRM" (for polytomous)
#' @param start Optional starting values
#' @param control Control parameters for optimization
#'
#' @return An object of class \code{gllamm_eirt}
#'
#' @examples
#' \dontrun{
#' # Create item-level covariates
#' item_chars <- data.frame(
#'   word_frequency = rnorm(20),
#'   item_length = rpois(20, 5)
#' )
#'
#' # Fit EIRT model
#' fit <- fit_eirt(
#'   responses,
#'   item_data = item_chars,
#'   difficulty_formula = ~ word_frequency + item_length,
#'   discrimination_formula = ~ 1,  # Constant discrimination
#'   model = "2PL"
#' )
#'
#' summary(fit)
#' }
#'
#' @export
fit_eirt <- function(response_matrix,
                     item_data,
                     difficulty_formula = ~ 1,
                     discrimination_formula = ~ 1,
                     model = c("Rasch", "2PL", "GRM", "PCM", "GPCM"),
                     start = NULL,
                     control = list()) {

  model <- match.arg(model)

  # Validate inputs
  n_persons <- nrow(response_matrix)
  n_items <- ncol(response_matrix)

  if (nrow(item_data) != n_items) {
    stop("item_data must have ", n_items, " rows (one per item). Found: ", nrow(item_data))
  }

  # Determine if polytomous
  unique_vals <- unique(as.vector(response_matrix[!is.na(response_matrix)]))
  is_polytomous <- (model %in% c("GRM", "PCM", "GPCM") || length(unique_vals) > 2)

  if (is_polytomous) {
    # Model types for polytomous: 1=GRM, 2=PCM, 3=GPCM
    # But for EIRT, we use 2PL-like structure (estimated discrimination)
    model_type <- switch(model,
                         "GRM" = 2L,
                         "PCM" = 1L,  # PCM has fixed discrimination=1
                         "GPCM" = 2L)
    n_categories_per_item <- apply(response_matrix, 2, function(x) {
      length(unique(x[!is.na(x)]))
    })
    max_categories <- max(n_categories_per_item)
  } else {
    model_type <- switch(model, Rasch = 1L, "2PL" = 2L)
    n_categories_per_item <- rep(2, n_items)
    max_categories <- 2
  }

  # Create design matrices for item covariates
  W_diff <- model.matrix(difficulty_formula, data = item_data)
  W_disc <- model.matrix(discrimination_formula, data = item_data)

  p_diff <- ncol(W_diff)
  p_disc <- ncol(W_disc)

  # Convert to long format
  y_long <- as.vector(t(response_matrix))
  person_id <- rep(1:n_persons, each = n_items) - 1L
  item_id <- rep(1:n_items, times = n_persons) - 1L

  # Remove missing
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
    W_difficulty = as.matrix(W_diff),
    W_discrimination = as.matrix(W_disc),
    n_persons = as.integer(n_persons),
    n_items = as.integer(n_items),
    n_obs = as.integer(n_obs),
    model_type = as.integer(model_type),
    is_polytomous = as.integer(is_polytomous),
    n_categories_per_item = as.integer(n_categories_per_item),
    max_categories = as.integer(max_categories)
  )

  # Initialize parameters
  if (is.null(start)) {
    theta_init <- rep(0, n_persons)
    gamma_init <- rep(0, p_diff)
    delta_init <- rep(0, p_disc)
    epsilon_b_init <- rep(0, n_items)
    epsilon_a_init <- rep(0, n_items)

    threshold_resid_init <- matrix(0, n_items, max_categories - 1)

    tmb_params <- list(
      theta = theta_init,
      gamma = gamma_init,
      delta = delta_init,
      epsilon_b = epsilon_b_init,
      epsilon_a = epsilon_a_init,
      log_sigma_epsilon_b = log(0.5),
      log_sigma_epsilon_a = log(0.5),
      log_sigma_theta = log(1.0),
      threshold_resid = threshold_resid_init
    )
  } else {
    tmb_params <- start
  }

  # Create TMB object
  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    random = c("theta", "epsilon_b", "epsilon_a"),
    DLL = "gllamm_eirt",
    silent = TRUE
  )

  # Optimize
  control_defaults <- list(eval.max = 3000, iter.max = 2000, trace = 0)
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

  gamma_hat <- par_full[names(par_full) == "gamma"]
  names(gamma_hat) <- colnames(W_diff)

  delta_hat <- par_full[names(par_full) == "delta"]
  names(delta_hat) <- colnames(W_disc)

  difficulty_hat <- par_full[names(par_full) == "difficulty"]
  names(difficulty_hat) <- paste0("Item", 1:n_items)

  discrimination_hat <- par_full[names(par_full) == "discrimination"]
  names(discrimination_hat) <- paste0("Item", 1:n_items)

  theta_hat <- par_full[names(par_full) == "theta"]
  names(theta_hat) <- paste0("Person", 1:n_persons)

  sigma_epsilon_b_hat <- exp(par_full[names(par_full) == "log_sigma_epsilon_b"])
  sigma_epsilon_a_hat <- exp(par_full[names(par_full) == "log_sigma_epsilon_a"])
  sigma_theta_hat <- exp(par_full[names(par_full) == "log_sigma_theta"])

  # Construct result
  result <- list(
    model = model,
    regression_coefficients = list(
      difficulty = gamma_hat,
      discrimination = delta_hat
    ),
    item_parameters = list(
      difficulty = difficulty_hat,
      discrimination = discrimination_hat
    ),
    person_abilities = theta_hat,
    ability_sd = sigma_theta_hat,
    residual_sd = list(
      difficulty = sigma_epsilon_b_hat,
      discrimination = sigma_epsilon_a_hat
    ),
    logLik = -opt$objective,
    AIC = 2 * opt$objective + 2 * length(obj$par),
    BIC = 2 * opt$objective + log(n_persons) * length(obj$par),
    convergence = list(
      converged = (opt$convergence == 0),
      message = opt$message
    ),
    n_persons = n_persons,
    n_items = n_items,
    formulas = list(
      difficulty = difficulty_formula,
      discrimination = discrimination_formula
    ),
    item_data = item_data,
    tmb_obj = obj,
    tmb_opt = opt,
    tmb_sdr = sdr
  )

  class(result) <- c("gllamm_eirt", "gllamm")

  return(result)
}


#' Print EIRT model results
#'
#' @param x Object of class gllamm_eirt
#' @param ... Additional arguments
#'
#' @export
print.gllamm_eirt <- function(x, ...) {
  cat("Explanatory IRT Model (", x$model, ")\n\n", sep = "")
  cat("Number of persons:", x$n_persons, "\n")
  cat("Number of items:", x$n_items, "\n\n")

  cat("Difficulty regression:\n")
  cat("  Formula:", deparse(x$formulas$difficulty), "\n")
  cat("  Coefficients:\n")
  print(round(x$regression_coefficients$difficulty, 3))
  cat("  Residual SD:", round(x$residual_sd$difficulty, 3), "\n\n")

  cat("Discrimination regression:\n")
  cat("  Formula:", deparse(x$formulas$discrimination), "\n")
  cat("  Coefficients:\n")
  print(round(x$regression_coefficients$discrimination, 3))
  cat("  Residual SD:", round(x$residual_sd$discrimination, 3), "\n\n")

  cat("Ability distribution:\n")
  cat("  Mean:", round(mean(x$person_abilities), 3), "\n")
  cat("  SD:", round(sd(x$person_abilities), 3), "\n")
  cat("  Estimated SD:", round(x$ability_sd, 3), "\n\n")

  cat("Log-likelihood:", round(x$logLik, 2), "\n")
  cat("AIC:", round(x$AIC, 2), "\n")
  cat("BIC:", round(x$BIC, 2), "\n")

  invisible(x)
}


#' Summary of EIRT model
#'
#' @param object Object of class gllamm_eirt
#' @param ... Additional arguments
#'
#' @export
summary.gllamm_eirt <- function(object, ...) {
  print(object)

  cat("\n\nFitted Item Parameters (first 10 items):\n")
  n_show <- min(10, object$n_items)
  item_params_df <- data.frame(
    difficulty = object$item_parameters$difficulty[1:n_show],
    discrimination = object$item_parameters$discrimination[1:n_show]
  )
  print(round(item_params_df, 3))

  if (object$n_items > 10) {
    cat("  ... (", object$n_items - 10, " more items)\n")
  }

  invisible(object)
}
