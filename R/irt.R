#' Fit Item Response Theory Models
#'
#' Fit dichotomous (Rasch, 2PL, 3PL) or polytomous (GRM, PCM, GPCM, NRM) IRT models
#'
#' @param response_matrix Matrix of item responses (persons x items).
#'   For dichotomous models: coded 0/1.
#'   For polytomous models: coded 1, 2, ..., K (ordered categories).
#' @param model Type of IRT model:
#'   Dichotomous: "Rasch", "2PL", "3PL"
#'   Polytomous: "GRM" (Graded Response), "PCM" (Partial Credit),
#'               "GPCM" (Generalized Partial Credit), "NRM" (Nominal Response)
#' @param start Optional starting values
#' @param control Control parameters for optimization
#'
#' @return An object of class \code{gllamm_irt}
#'
#' @examples
#' \dontrun{
#' # Dichotomous example (Rasch)
#' set.seed(123)
#' n_persons <- 500
#' n_items <- 20
#' theta <- rnorm(n_persons, 0, 1)
#' difficulty <- rnorm(n_items, 0, 1)
#'
#' # Generate binary responses
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
#' # Polytomous example (GRM)
#' # Generate 5-category responses
#' responses_poly <- matrix(NA, n_persons, n_items)
#' thresholds <- matrix(seq(-2, 2, length.out = 4), n_items, 4, byrow = TRUE)
#' for (i in 1:n_persons) {
#'   for (j in 1:n_items) {
#'     probs <- c(plogis(theta[i] - thresholds[j, 1]),
#'                diff(plogis(theta[i] - thresholds[j, ])),
#'                1 - plogis(theta[i] - thresholds[j, 4]))
#'     responses_poly[i, j] <- sample(1:5, 1, prob = probs)
#'   }
#' }
#'
#' # Fit GRM model
#' fit_grm <- fit_irt(responses_poly, model = "GRM")
#' summary(fit_grm)
#' }
#'
#' @export
fit_irt <- function(response_matrix,
                    model = c("Rasch", "2PL", "3PL", "GRM", "PCM", "GPCM", "NRM"),
                    start = NULL, control = list()) {

  model <- match.arg(model)

  # Check if polytomous model
  is_polytomous <- model %in% c("GRM", "PCM", "GPCM", "NRM")

  # Detect response type from data
  unique_vals <- unique(as.vector(response_matrix[!is.na(response_matrix)]))
  n_categories <- length(unique_vals)

  # Validate response coding
  if (!is_polytomous && n_categories > 2) {
    stop("Dichotomous models (Rasch/2PL/3PL) require binary responses (0/1). ",
         "Found ", n_categories, " categories. Use polytomous models (GRM/PCM/GPCM) instead.")
  }

  if (is_polytomous && n_categories <= 2) {
    warning("Polytomous model specified but data appears dichotomous. ",
            "Consider using Rasch/2PL/3PL models instead.")
  }

  # Dispatch to appropriate function
  if (is_polytomous) {
    return(fit_irt_polytomous(response_matrix, model, start, control))
  } else {
    return(fit_irt_dichotomous(response_matrix, model, start, control))
  }
}


#' Internal function for dichotomous IRT models
#' @keywords internal
fit_irt_dichotomous <- function(response_matrix, model, start, control) {

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
  cat("Number of items:", x$n_items, "\n")

  # Check if polytomous model (item_parameters is a list)
  is_polytomous <- is.list(x$item_parameters) && "thresholds" %in% names(x$item_parameters)

  if (is_polytomous) {
    cat("Model type: Polytomous (max", x$max_categories, "categories)\n\n")

    cat("Item discriminations:\n")
    print(round(x$item_parameters$discrimination, 3))

    cat("\nItem thresholds (first 5 items):\n")
    for (i in 1:min(5, x$n_items)) {
      cat("  Item", i, ":", paste(round(x$item_parameters$thresholds[[i]], 3), collapse = ", "), "\n")
    }
    if (x$n_items > 5) {
      cat("  ... (", x$n_items - 5, " more items)\n")
    }
  } else {
    cat("Model type: Dichotomous\n\n")

    cat("Item parameters:\n")
    print(round(x$item_parameters, 3))
  }

  cat("\nAbility distribution:\n")
  cat("  Mean:", round(mean(x$person_abilities), 3), "\n")
  cat("  SD:", round(sd(x$person_abilities), 3), "\n")
  cat("  Estimated SD:", round(x$ability_sd, 3), "\n\n")

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


#' Internal function for polytomous IRT models
#' @keywords internal
fit_irt_polytomous <- function(response_matrix, model, start, control) {

  # Model codes for polytomous IRT
  model_code <- switch(model,
                       GRM = 1L,   # Graded Response Model
                       PCM = 2L,   # Partial Credit Model
                       GPCM = 3L,  # Generalized Partial Credit Model
                       NRM = 4L)   # Nominal Response Model

  # Dimensions
  n_persons <- nrow(response_matrix)
  n_items <- ncol(response_matrix)

  # Determine number of categories per item
  n_categories_per_item <- apply(response_matrix, 2, function(x) {
    length(unique(x[!is.na(x)]))
  })
  max_categories <- max(n_categories_per_item)

  # Validate response coding (should be 1, 2, ..., K)
  for (j in 1:n_items) {
    item_vals <- unique(response_matrix[!is.na(response_matrix[, j]), j])
    if (min(item_vals) != 1 || max(item_vals) != n_categories_per_item[j]) {
      stop("Item ", j, " responses must be coded 1 to ", n_categories_per_item[j],
           ". Found range: ", min(item_vals), " to ", max(item_vals))
    }
  }

  # Convert to long format
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
    n_categories_per_item = as.integer(n_categories_per_item),
    max_categories = as.integer(max_categories),
    n_persons = as.integer(n_persons),
    n_items = as.integer(n_items),
    n_obs = as.integer(n_obs),
    model_type = as.integer(model_code)
  )

  # Initialize parameters
  if (is.null(start)) {
    # Initialize abilities to zero
    theta_init <- rep(0, n_persons)

    # Initialize thresholds from marginal category proportions
    threshold_raw_init <- matrix(0, n_items, max_categories - 1)

    for (j in 1:n_items) {
      K <- n_categories_per_item[j]
      if (K > 1) {
        # Get category proportions
        item_responses <- response_matrix[, j]
        cat_props <- table(item_responses) / sum(!is.na(item_responses))

        # Compute cumulative proportions
        cum_props <- cumsum(cat_props)

        # Inverse logit to get thresholds
        thresholds <- qlogis(cum_props[-K])

        # Transform to ordered parameterization (log-differences)
        threshold_raw_init[j, 1] <- thresholds[1]
        if (K > 2) {
          for (k in 2:(K-1)) {
            threshold_raw_init[j, k] <- log(max(thresholds[k] - thresholds[k-1], 0.1))
          }
        }
      }
    }

    # Initialize discriminations
    # For PCM, these will be constrained to 1 in the template
    discrimination_init <- rep(1, n_items)

    tmb_params <- list(
      theta = theta_init,
      threshold_raw = threshold_raw_init,
      discrimination = discrimination_init,
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
    DLL = "gllamm_irt_poly",
    silent = TRUE
  )

  # Optimize
  control_defaults <- list(eval.max = 3000, iter.max = 1500, trace = 0)
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

  # Extract discrimination parameters
  discrimination_hat <- par_full[names(par_full) == "discrimination"]
  names(discrimination_hat) <- paste0("Item", 1:n_items)

  # Extract threshold parameters and reconstruct ordered thresholds
  threshold_raw_hat <- matrix(par_full[names(par_full) == "threshold_raw"],
                              n_items, max_categories - 1)

  # Reconstruct ordered thresholds for each item
  ordered_thresholds <- vector("list", n_items)
  for (j in 1:n_items) {
    K <- n_categories_per_item[j]
    if (K > 1) {
      tau <- numeric(K - 1)
      tau[1] <- threshold_raw_hat[j, 1]
      if (K > 2) {
        for (k in 2:(K-1)) {
          tau[k] <- tau[k-1] + exp(threshold_raw_hat[j, k])
        }
      }
      ordered_thresholds[[j]] <- tau
    }
  }
  names(ordered_thresholds) <- paste0("Item", 1:n_items)

  # Extract abilities
  theta_hat <- par_full[names(par_full) == "theta"]
  names(theta_hat) <- paste0("Person", 1:n_persons)

  sigma_theta_hat <- exp(par_full[names(par_full) == "log_sigma_theta"])

  # Construct result object
  result <- list(
    model = model,
    item_parameters = list(
      discrimination = discrimination_hat,
      thresholds = ordered_thresholds,
      n_categories = n_categories_per_item
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
    n_categories = n_categories_per_item,
    max_categories = max_categories,
    tmb_obj = obj,
    tmb_opt = opt,
    tmb_sdr = sdr
  )

  class(result) <- c("gllamm_irt", "gllamm")

  return(result)
}
