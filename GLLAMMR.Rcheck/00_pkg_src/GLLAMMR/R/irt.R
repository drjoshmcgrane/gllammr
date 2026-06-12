#' Fit Item Response Theory Models
#'
#' Fit dichotomous (Rasch, 2PL, 3PL) or polytomous (GRM, PCM, GPCM, NRM) IRT models.
#' Optionally include multi-level random effects for hierarchical or clustered data.
#'
#' @param response_matrix Matrix of item responses (persons x items).
#'   For dichotomous models: coded 0/1.
#'   For polytomous models: coded 1, 2, ..., K (ordered categories).
#' @param model Type of IRT model:
#'   Dichotomous: "Rasch", "2PL", "3PL"
#'   Polytomous: "GRM" (Graded Response), "PCM" (Partial Credit),
#'               "GPCM" (Generalized Partial Credit), "NRM" (Nominal Response)
#' @param person_data Optional data frame with person-level variables for multi-level models.
#'   Must have one row per person (rows correspond to rows in response_matrix).
#'   Used to specify grouping variables in the random effects formula.
#' @param random Optional random effects formula using lme4-style syntax.
#'   Examples: \code{~ (1 | class)}, \code{~ (1 | school/class)},
#'   \code{~ (1 | student) + (1 | time)}.
#'   Requires person_data to contain the grouping variables.
#' @param weights Optional vector of case weights. Length must equal number of persons.
#' @param method Estimation method. "auto" (default) uses "em" for
#'   single-level models and "laplace" whenever multi-level structure
#'   (\code{random}) or standard errors (\code{se = TRUE}) require it.
#'   "em" is Bock-Aitkin marginal maximum likelihood with fixed
#'   Gauss-Hermite quadrature (the mirt/TAM algorithm; typically 20-50x
#'   faster and evaluates the marginal likelihood more exactly than the
#'   Laplace approximation). "laplace" is the TMB path, required for
#'   multi-level models. EM person abilities are EAP scores; Laplace
#'   abilities are posterior modes.
#' @param quad_points Number of quadrature nodes for method = "em"
#'   (default 61)
#' @param se Compute parameter standard errors via TMB::sdreport (default
#'   FALSE, matching the default behavior of mirt; SE computation roughly
#'   doubles the fitting time for large person samples)
#' @param mc_items For 3PL model only: which items have guessing parameters.
#'   Can be: NULL (default, all items have guessing), logical vector (length = n_items),
#'   or integer vector (indices of MC items). Non-MC items use 2PL likelihood (no guessing).
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
#'
#' # 3PL with selective guessing (mixed MC and non-MC items)
#' # Assessment: 20 items, first 15 are MC, last 5 are open-ended
#' fit_3pl <- fit_irt(responses, model = "3PL", mc_items = 1:15)
#' # Only items 1-15 get guessing parameters
#' # Items 16-20 use 2PL likelihood (no guessing)
#'
#' # Multi-level IRT: students nested in classes
#' person_data <- data.frame(
#'   person_id = 1:n_persons,
#'   class_id = rep(1:10, each = 50)
#' )
#' fit_multilevel <- fit_irt(responses, model = "2PL",
#'                            person_data = person_data,
#'                            random = ~ (1 | class_id))
#' # theta_i = theta_0i + u_class[class[i]]
#' }
#'
#' @export
fit_irt <- function(response_matrix,
                    model = c("Rasch", "2PL", "3PL", "GRM", "PCM", "GPCM", "NRM"),
                    person_data = NULL,
                    random = NULL,
                    weights = NULL,
                    mc_items = NULL,
                    method = c("auto", "em", "laplace"),
                    quad_points = 61,
                    se = FALSE,
                    start = NULL, control = list()) {

  model <- match.arg(model)
  method <- match.arg(method)

  if (method == "auto") {
    method <- if (!is.null(random) || isTRUE(se)) "laplace" else "em"
  }
  if (method == "em" && !is.null(random)) {
    stop("method = \"em\" supports single-level models; ",
         "use method = \"laplace\" for multi-level IRT")
  }
  if (method == "em" && isTRUE(se)) {
    warning("Standard errors are not yet available under method = \"em\"; ",
            "ignoring se = TRUE. Use method = \"laplace\" for SEs.")
  }

  # Validate multi-level parameters
  has_random <- !is.null(random)

  if (has_random && is.null(person_data)) {
    stop("person_data must be provided when random effects are specified")
  }

  if (!is.null(person_data)) {
    n_persons <- nrow(response_matrix)
    if (nrow(person_data) != n_persons) {
      stop("person_data must have same number of rows as response_matrix (",
           nrow(person_data), " vs ", n_persons, ")")
    }
  }

  # Parse random effects if specified
  re_info <- NULL
  if (has_random) {
    # Parse the random effects formula
    re_terms <- parse_random_formula(random, person_data)

    # Create grouping matrix
    re_info <- create_grouping_matrix(re_terms, person_data)
  }

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

  # Validate weights if provided
  n_persons <- nrow(response_matrix)
  if (!is.null(weights)) {
    if (length(weights) != n_persons) {
      stop("Length of weights (", length(weights), ") must match number of persons (", n_persons, ")")
    }
    if (any(weights < 0, na.rm = TRUE)) {
      stop("All weights must be non-negative")
    }
    if (any(is.na(weights))) {
      stop("weights cannot contain missing values")
    }
  }

  # Validate mc_items (only for 3PL)
  if (!is.null(mc_items) && model != "3PL") {
    warning("mc_items parameter is only used for 3PL model. Ignoring.")
    mc_items <- NULL
  }

  # Dispatch to appropriate function
  if (is_polytomous) {
    if (!is.null(mc_items)) {
      warning("mc_items parameter is not applicable to polytomous models. Ignoring.")
    }
    if (method == "em") {
      return(fit_irt_em(response_matrix, model = model, weights = weights,
                        quad_points = quad_points, control = control))
    }
    return(fit_irt_polytomous(response_matrix, model, weights, re_info, se, start, control))
  } else {
    if (method == "em") {
      return(fit_irt_em(response_matrix, model = model, weights = weights,
                        mc_items = mc_items, quad_points = quad_points,
                        control = control))
    }
    return(fit_irt_dichotomous(response_matrix, model, weights, mc_items, re_info, se, start, control))
  }
}


#' Internal function for dichotomous IRT models
#' @keywords internal
fit_irt_dichotomous <- function(response_matrix, model, weights, mc_items, re_info, se, start, control) {

  model_code <- switch(model, Rasch = 1L, "2PL" = 2L, "3PL" = 3L)
  n_items <- ncol(response_matrix)

  # Process mc_items into indicator vector (for 3PL only)
  if (model == "3PL") {
    if (is.null(mc_items)) {
      # Default: all items have guessing
      mc_indicator <- rep(1L, n_items)
    } else if (is.logical(mc_items)) {
      # Logical vector
      if (length(mc_items) != n_items) {
        stop("mc_items length (", length(mc_items), ") must equal number of items (", n_items, ")")
      }
      mc_indicator <- as.integer(mc_items)
    } else if (is.numeric(mc_items)) {
      # Integer vector of indices
      if (any(mc_items < 1 | mc_items > n_items)) {
        stop("mc_items indices must be between 1 and ", n_items)
      }
      mc_indicator <- rep(0L, n_items)
      mc_indicator[mc_items] <- 1L
    } else {
      stop("mc_items must be NULL, logical vector, or integer vector")
    }
  } else {
    # Not 3PL: dummy vector (not used)
    mc_indicator <- rep(0L, n_items)
  }

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

  # Expand weights from person-level to observation-level
  if (!is.null(weights)) {
    weights_long <- rep(weights, each = n_items)[complete_cases]
  } else {
    weights_long <- rep(1.0, sum(complete_cases))
  }

  n_obs <- length(y_long)

  # Determine if using multi-level model
  has_random <- !is.null(re_info)

  # Prepare TMB data
  tmb_data <- list(
    y = as.numeric(y_long),
    person_id = as.integer(person_id),
    item_id = as.integer(item_id),
    n_persons = as.integer(n_persons),
    n_items = as.integer(n_items),
    n_obs = as.integer(n_obs),
    weights = as.numeric(weights_long),
    model_type = as.integer(model_code),
    mc_items = as.integer(mc_indicator)
  )

  # Add multi-level data if applicable
  if (has_random) {
    tmb_data$has_random <- 1L
    tmb_data$n_random_effects <- as.integer(re_info$n_re)
    tmb_data$group_ids <- as.matrix(re_info$group_ids)
    tmb_data$n_groups <- as.integer(re_info$n_groups)
    tmb_data$max_n_groups <- as.integer(max(re_info$n_groups))
    tmb_data$re_design <- as.matrix(re_info$re_design)
  } else {
    tmb_data$has_random <- 0L
    tmb_data$n_random_effects <- 0L
    tmb_data$group_ids <- matrix(0L, 0, 0)
    tmb_data$n_groups <- integer(0)
    tmb_data$max_n_groups <- 0L
    tmb_data$re_design <- matrix(0, 0, 0)
  }

  # Initialize parameters
  if (is.null(start)) {
    # Initialize abilities to zero (theta_0 for multi-level)
    if (has_random) {
      theta_0_init <- rep(0, n_persons)
    } else {
      theta_init <- rep(0, n_persons)
    }

    # Initialize difficulties from marginal item proportions
    item_means <- colMeans(response_matrix, na.rm = TRUE)
    difficulty_init <- -qlogis(pmin(pmax(item_means, 0.01), 0.99))

    # Initialize discriminations
    discrimination_init <- rep(1, n_items)

    # Initialize guessing parameters
    guessing_init <- rep(0.1, n_items)

    if (has_random) {
      # Multi-level model parameters
      max_n_groups <- max(re_info$n_groups)
      u_random_init <- matrix(0, max_n_groups, re_info$n_re)
      log_sigma_random_init <- rep(log(0.5), re_info$n_re)

      tmb_params <- list(
        theta_0 = theta_0_init,
        difficulty = difficulty_init,
        discrimination = discrimination_init,
        guessing = guessing_init,
        log_sigma_theta = log(1.0),
        u_random = u_random_init,
        log_sigma_random = log_sigma_random_init
      )
    } else {
      # Standard IRT parameters
      tmb_params <- list(
        theta = theta_init,
        difficulty = difficulty_init,
        discrimination = discrimination_init,
        guessing = guessing_init,
        log_sigma_theta = log(1.0)
      )
    }
  } else {
    tmb_params <- start
  }

  # Identification and dead-parameter maps:
  # - Rasch: discrimination/guessing never read; sigma_theta free (Rasch metric)
  # - 2PL/3PL: sigma_theta fixed at 1 (standard IRT identification; a free
  #   discrimination and a free ability SD are jointly unidentified)
  # - 3PL: guessing only for designated multiple-choice items
  map_list <- list()
  if (model_code == 1) {
    map_list$discrimination <- factor(rep(NA, n_items))
    map_list$guessing <- factor(rep(NA, n_items))
  } else if (model_code == 2) {
    map_list$guessing <- factor(rep(NA, n_items))
    map_list$log_sigma_theta <- factor(NA)
  } else {
    map_list$guessing <- factor(ifelse(mc_indicator == 1, seq_len(n_items), NA))
    map_list$log_sigma_theta <- factor(NA)
  }
  if (!is.null(map_list$log_sigma_theta)) {
    tmb_params$log_sigma_theta <- 0   # exp(0) = 1 when fixed
  }

  # Create TMB object with appropriate model and random effects
  if (has_random) {
    # Multi-level model: integrate out theta_0 and u_random
    tmb_data$model_name <- "irt_multilevel"
    obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = tmb_params,
      random = c("theta_0", "u_random"),
      map = map_list,
      DLL = "GLLAMMR",
      silent = TRUE
    )
  } else {
    # Standard IRT: integrate out theta only
    tmb_data$model_name <- "irt"
    obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = tmb_params,
      random = "theta",
      map = map_list,
      DLL = "GLLAMMR",
      silent = TRUE
    )
  }

  # Optimize. Box constraints keep discrimination positive and bounded
  # (unbounded a can diverge under the Laplace approximation) and guessing
  # inside its probability range.
  control_defaults <- list(eval.max = 2000, iter.max = 1000, trace = 0)
  control <- modifyList(control_defaults, control)

  par_names <- names(obj$par)
  lower <- rep(-Inf, length(par_names))
  upper <- rep(Inf, length(par_names))
  lower[par_names == "discrimination"] <- 0.05
  upper[par_names == "discrimination"] <- 10
  lower[par_names == "guessing"] <- 1e-3
  upper[par_names == "guessing"] <- 0.5

  opt <- nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    lower = lower,
    upper = upper,
    control = control
  )

  # Standard errors on request only (sdreport roughly doubles fit time)
  sdr <- if (se) try(TMB::sdreport(obj), silent = TRUE) else NULL

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
    # Guessing is mapped off for non-MC items; estimated values fill the
    # MC positions, others have no guessing parameter
    guessing_free <- par_full[names(par_full) == "guessing"]
    guessing_hat <- rep(0, n_items)
    guessing_hat[mc_indicator == 1] <- guessing_free
    names(guessing_hat) <- paste0("Item", 1:n_items)
  } else {
    guessing_hat <- NULL
  }

  # Extract abilities (different parameter name for multi-level)
  if (has_random) {
    theta_0_hat <- par_full[names(par_full) == "theta_0"]
    names(theta_0_hat) <- paste0("Person", 1:n_persons)
    theta_hat <- theta_0_hat  # Store as theta_hat for compatibility
  } else {
    theta_hat <- par_full[names(par_full) == "theta"]
    names(theta_hat) <- paste0("Person", 1:n_persons)
  }

  # Fixed at 1 (mapped) for 2PL/3PL identification
  log_sigma_free <- par_full[names(par_full) == "log_sigma_theta"]
  sigma_theta_hat <- if (length(log_sigma_free) == 0) 1 else exp(unname(log_sigma_free))

  # Extract random effects if present
  if (has_random) {
    u_random_hat <- matrix(par_full[names(par_full) == "u_random"],
                           tmb_data$max_n_groups, re_info$n_re)
    sigma_random_hat <- exp(par_full[names(par_full) == "log_sigma_random"])
    names(sigma_random_hat) <- re_info$group_names
  }

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
    mc_items = if (model_code == 3) mc_indicator else NULL,  # Store for 3PL
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

  # Add random effects information if multi-level model
  if (has_random) {
    # Compute composite abilities (theta_0 + random effects)
    composite_theta <- theta_0_hat  # Start with person deviations
    for (p in 1:n_persons) {
      for (re in 1:re_info$n_re) {
        group <- re_info$group_ids[p, re]
        if (group >= 0) {  # -1 indicates NA (partial nesting)
          composite_theta[p] <- composite_theta[p] + u_random_hat[group + 1, re]  # +1 for R indexing
        }
      }
    }

    # Compute ICCs (Intraclass Correlations)
    var_random <- sigma_random_hat^2
    var_person <- sigma_theta_hat^2
    var_total <- sum(var_random) + var_person + pi^2/3  # logistic variance

    icc_values <- c(var_random, var_person) / var_total
    names(icc_values) <- c(re_info$group_names, "Person")

    result$random_effects <- list(
      u_random = u_random_hat,
      sigma_random = sigma_random_hat,
      group_names = re_info$group_names,
      n_groups = re_info$n_groups,
      icc = icc_values,
      composite_theta = composite_theta
    )

    # Update class to indicate multi-level
    class(result) <- c("gllamm_irt_multilevel", "gllamm_irt", "gllamm")
  } else {
    result$random_effects <- NULL
    class(result) <- c("gllamm_irt", "gllamm")
  }

  return(result)
}


#' @export
print.gllamm_irt <- function(x, ...) {
  # Check if multi-level
  is_multilevel <- inherits(x, "gllamm_irt_multilevel")

  if (is_multilevel) {
    cat("Multi-Level IRT Model (", x$model, ")\n\n", sep = "")
  } else {
    cat("IRT Model (", x$model, ")\n\n", sep = "")
  }

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
    cat("Model type: Dichotomous\n")

    # Show MC items info for 3PL
    if (x$model == "3PL" && !is.null(x$mc_items)) {
      n_mc <- sum(x$mc_items)
      if (n_mc < x$n_items) {
        cat("MC items with guessing:", n_mc, "/", x$n_items, "\n")
      }
    }
    cat("\n")

    cat("Item parameters:\n")
    print(round(x$item_parameters, 3))
  }

  # Print random effects information for multi-level models
  if (is_multilevel) {
    cat("\nRandom Effects:\n")
    cat("  Grouping variables:", paste(x$random_effects$group_names, collapse = ", "), "\n")
    cat("  Number of groups:", paste(x$random_effects$n_groups, collapse = ", "), "\n\n")

    cat("Variance Components:\n")
    variance_table <- data.frame(
      Groups = c(x$random_effects$group_names, "Person", "Residual"),
      Variance = round(c(x$random_effects$sigma_random^2, x$ability_sd^2, pi^2/3), 4),
      Std.Dev = round(c(x$random_effects$sigma_random, x$ability_sd, sqrt(pi^2/3)), 4)
    )
    print(variance_table, row.names = FALSE)

    cat("\nIntraclass Correlations:\n")
    icc_df <- data.frame(
      Level = names(x$random_effects$icc),
      ICC = round(x$random_effects$icc, 4)
    )
    print(icc_df, row.names = FALSE)
    cat("\n")
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


#' Validate and auto-recode polytomous response matrices (shared by the
#' Laplace and EM estimation paths)
#'
#' @return list(response_matrix, n_categories_per_item, max_categories)
#' @keywords internal
validate_poly_responses <- function(response_matrix) {
  n_items <- ncol(response_matrix)
  # Determine number of categories per item
  n_categories_per_item <- apply(response_matrix, 2, function(x) {
    length(unique(x[!is.na(x)]))
  })
  max_categories <- max(n_categories_per_item)

  # Validate and auto-recode response coding (should be 1, 2, ..., K)
  needs_recoding <- FALSE
  recode_items <- integer(0)
  unobserved_cat_items <- integer(0)
  constant_items <- integer(0)

  for (j in 1:n_items) {
    item_vals <- unique(response_matrix[!is.na(response_matrix[, j]), j])
    min_val <- min(item_vals)
    max_val <- max(item_vals)
    n_cats <- n_categories_per_item[j]

    # Check if coded 0-based (dichotomous style)
    if (min_val == 0 && max_val == n_cats - 1) {
      needs_recoding <- TRUE
      recode_items <- c(recode_items, j)
    }
    # Properly coded 1-based with all categories observed
    else if (min_val == 1 && max_val == n_cats) {
      # OK
    }
    # 1-based coding but some intermediate categories never observed:
    # treat the item as having max_val categories
    else if (min_val == 1 && max_val > n_cats) {
      unobserved_cat_items <- c(unobserved_cat_items, j)
    }
    # Constant item (every respondent chose the same category)
    else if (n_cats == 1 && min_val >= 1) {
      constant_items <- c(constant_items, j)
    }
    else {
      stop("Item ", j, " has invalid response coding. ",
           "Found range [", min_val, ", ", max_val, "] but expected [1, ", n_cats, "] ",
           "for ", n_cats, "-category item. ",
           "Polytomous models require responses coded as 1, 2, ..., K.")
    }
  }

  # Auto-recode 0-based items to 1-based
  if (needs_recoding) {
    message("Note: Auto-recoding ", length(recode_items), " binary item(s) from 0/1 to 1/2 coding.\n",
            "  Items ", paste(recode_items, collapse = ", "), " detected as 0-based.\n",
            "  Polytomous models require 1-based coding (1, 2, ..., K).")
    for (j in recode_items) {
      response_matrix[, j] <- response_matrix[, j] + 1
    }
    # Recalculate n_categories_per_item after recoding
    n_categories_per_item <- apply(response_matrix, 2, function(x) {
      length(unique(x[!is.na(x)]))
    })
  }

  # Items with unobserved categories or no variance: use the maximum
  # observed value as the number of categories
  if (length(unobserved_cat_items) > 0) {
    warning("Item(s) ", paste(unobserved_cat_items, collapse = ", "),
            " have response categories that were never observed. ",
            "Treating the maximum observed value as the number of categories; ",
            "thresholds for unobserved categories are weakly identified.")
  }
  if (length(constant_items) > 0) {
    warning("Item(s) ", paste(constant_items, collapse = ", "),
            " have no response variance (all respondents chose the same ",
            "category); their item parameters are not identified.")
  }
  for (j in c(unobserved_cat_items, constant_items)) {
    n_categories_per_item[j] <- as.integer(max(response_matrix[, j], na.rm = TRUE))
  }
  max_categories <- max(n_categories_per_item)

  list(response_matrix = response_matrix,
       n_categories_per_item = n_categories_per_item,
       max_categories = max_categories)
}


#' Internal function for polytomous IRT models
#' @keywords internal
fit_irt_polytomous <- function(response_matrix, model, weights, re_info, se, start, control) {

  # Model codes for polytomous IRT
  model_code <- switch(model,
                       GRM = 1L,   # Graded Response Model
                       PCM = 2L,   # Partial Credit Model
                       GPCM = 3L,  # Generalized Partial Credit Model
                       NRM = 4L)   # Nominal Response Model

  # Dimensions
  n_persons <- nrow(response_matrix)
  n_items <- ncol(response_matrix)

  vp <- validate_poly_responses(response_matrix)
  response_matrix <- vp$response_matrix
  n_categories_per_item <- vp$n_categories_per_item
  max_categories <- vp$max_categories

  # Convert to long format
  y_long <- as.vector(t(response_matrix))
  person_id <- rep(1:n_persons, each = n_items) - 1L  # 0-indexed
  item_id <- rep(1:n_items, times = n_persons) - 1L   # 0-indexed

  # Remove missing values
  complete_cases <- !is.na(y_long)
  y_long <- y_long[complete_cases]
  person_id <- person_id[complete_cases]
  item_id <- item_id[complete_cases]

  # Expand weights from person-level to observation-level
  if (!is.null(weights)) {
    weights_long <- rep(weights, each = n_items)[complete_cases]
  } else {
    weights_long <- rep(1.0, sum(complete_cases))
  }

  n_obs <- length(y_long)

  # Determine if using multi-level model
  has_random <- !is.null(re_info)

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
    weights = as.numeric(weights_long),
    model_type = as.integer(model_code)
  )

  # Add multi-level data if applicable
  if (has_random) {
    tmb_data$has_random <- 1L
    tmb_data$n_random_effects <- as.integer(re_info$n_re)
    tmb_data$group_ids <- as.matrix(re_info$group_ids)
    tmb_data$n_groups <- as.integer(re_info$n_groups)
    tmb_data$max_n_groups <- as.integer(max(re_info$n_groups))
    tmb_data$re_design <- as.matrix(re_info$re_design)
  } else {
    tmb_data$has_random <- 0L
    tmb_data$n_random_effects <- 0L
    tmb_data$group_ids <- matrix(0L, 0, 0)
    tmb_data$n_groups <- integer(0)
    tmb_data$max_n_groups <- 0L
    tmb_data$re_design <- matrix(0, 0, 0)
  }

  # Initialize parameters
  if (is.null(start)) {
    # Initialize abilities to zero (theta_0 for multi-level)
    if (has_random) {
      theta_0_init <- rep(0, n_persons)
    } else {
      theta_init <- rep(0, n_persons)
    }

    # Initialize thresholds from marginal category proportions
    threshold_raw_init <- matrix(0, n_items, max_categories - 1)

    for (j in 1:n_items) {
      K <- n_categories_per_item[j]
      if (K > 1) {
        # Get smoothed category proportions over the full 1..K range
        # (unobserved categories get a small positive mass so that the
        # initial thresholds stay finite)
        item_responses <- response_matrix[, j]
        cat_counts <- table(factor(item_responses, levels = 1:K))
        n_j <- sum(!is.na(item_responses))
        cat_props <- (as.numeric(cat_counts) + 0.5) / (n_j + 0.5 * K)

        # Compute cumulative proportions
        cum_props <- cumsum(cat_props)

        # Inverse logit to get thresholds (clamped away from 0/1)
        thresholds <- qlogis(pmin(pmax(cum_props[-K], 1e-3), 1 - 1e-3))

        if (model %in% c("PCM", "GPCM")) {
          # PCM/GPCM: threshold_raw stores free step difficulties directly
          # Initialize as evenly spaced step difficulties
          threshold_raw_init[j, 1:(K-1)] <- seq(-1, 1, length.out = K-1)
        } else {
          # GRM: ordered parameterization (log-differences for spacing)
          threshold_raw_init[j, 1] <- thresholds[1]
          if (K > 2) {
            for (k in 2:(K-1)) {
              threshold_raw_init[j, k] <- log(max(thresholds[k] - thresholds[k-1], 0.1))
            }
          }
        }
      }
    }

    # Initialize discriminations
    # For PCM, these will be constrained to 1 in the template
    discrimination_init <- rep(1, n_items)

    if (has_random) {
      # Multi-level model parameters
      max_n_groups <- max(re_info$n_groups)
      u_random_init <- matrix(0, max_n_groups, re_info$n_re)
      log_sigma_random_init <- rep(log(0.5), re_info$n_re)

      tmb_params <- list(
        theta_0 = theta_0_init,
        threshold_raw = threshold_raw_init,
        discrimination = discrimination_init,
        log_sigma_theta = log(1.0),
        u_random = u_random_init,
        log_sigma_random = log_sigma_random_init
      )
    } else {
      # Standard polytomous IRT parameters
      tmb_params <- list(
        theta = theta_init,
        threshold_raw = threshold_raw_init,
        discrimination = discrimination_init,
        log_sigma_theta = log(1.0)
      )
    }
  } else {
    tmb_params <- start
  }

  # Identification and dead-parameter maps:
  # - PCM: discrimination never read; sigma_theta free (Rasch-family metric)
  # - GRM/GPCM/NRM: sigma_theta fixed at 1 (free discrimination and free
  #   ability SD are jointly unidentified)
  map_list <- list()
  if (model == "PCM") {
    map_list$discrimination <- factor(rep(NA, n_items))
  } else {
    map_list$log_sigma_theta <- factor(NA)
    tmb_params$log_sigma_theta <- 0   # exp(0) = 1 when fixed
  }

  # Create TMB object with appropriate model and random effects
  if (has_random) {
    # Multi-level polytomous model
    tmb_data$model_name <- "irt_poly_multilevel"
    obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = tmb_params,
      random = c("theta_0", "u_random"),
      map = map_list,
      DLL = "GLLAMMR",
      silent = TRUE
    )
  } else {
    # Standard polytomous IRT
    tmb_data$model_name <- "irt_poly"
    obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = tmb_params,
      random = "theta",
      map = map_list,
      DLL = "GLLAMMR",
      silent = TRUE
    )
  }

  # Optimize, with discrimination kept positive and bounded
  control_defaults <- list(eval.max = 3000, iter.max = 1500, trace = 0)
  control <- modifyList(control_defaults, control)

  par_names <- names(obj$par)
  lower <- rep(-Inf, length(par_names))
  upper <- rep(Inf, length(par_names))
  lower[par_names == "discrimination"] <- 0.05
  upper[par_names == "discrimination"] <- 10

  opt <- nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    lower = lower,
    upper = upper,
    control = control
  )

  # Standard errors on request only (sdreport roughly doubles fit time)
  sdr <- if (se) try(TMB::sdreport(obj), silent = TRUE) else NULL

  # Extract parameters
  par_full <- obj$env$last.par.best

  # Extract discrimination parameters (mapped off, fixed at 1, for PCM)
  discrimination_hat <- par_full[names(par_full) == "discrimination"]
  if (length(discrimination_hat) == 0) {
    discrimination_hat <- rep(1, n_items)
  }
  names(discrimination_hat) <- paste0("Item", 1:n_items)

  # Extract threshold parameters and reconstruct ordered thresholds
  threshold_raw_hat <- matrix(par_full[names(par_full) == "threshold_raw"],
                              n_items, max_categories - 1)

  # Reconstruct ordered thresholds for each item
  ordered_thresholds <- vector("list", n_items)
  for (j in 1:n_items) {
    K <- n_categories_per_item[j]
    if (K > 1) {
      if (model %in% c("PCM", "GPCM")) {
        # PCM/GPCM: threshold_raw stores free step difficulties directly
        tau <- threshold_raw_hat[j, 1:(K-1)]
      } else {
        # GRM: reconstruct ordered thresholds from log-difference parameterization
        tau <- numeric(K - 1)
        tau[1] <- threshold_raw_hat[j, 1]
        if (K > 2) {
          for (k in 2:(K-1)) {
            tau[k] <- tau[k-1] + exp(threshold_raw_hat[j, k])
          }
        }
      }
      ordered_thresholds[[j]] <- tau
    }
  }
  names(ordered_thresholds) <- paste0("Item", 1:n_items)

  # Extract abilities (different parameter name for multi-level)
  if (has_random) {
    theta_0_hat <- par_full[names(par_full) == "theta_0"]
    names(theta_0_hat) <- paste0("Person", 1:n_persons)
    theta_hat <- theta_0_hat  # Store as theta_hat for compatibility
  } else {
    theta_hat <- par_full[names(par_full) == "theta"]
    names(theta_hat) <- paste0("Person", 1:n_persons)
  }

  # Fixed at 1 (mapped) for GRM/GPCM/NRM identification
  log_sigma_free <- par_full[names(par_full) == "log_sigma_theta"]
  sigma_theta_hat <- if (length(log_sigma_free) == 0) 1 else exp(unname(log_sigma_free))

  # Extract random effects if present
  if (has_random) {
    u_random_hat <- matrix(par_full[names(par_full) == "u_random"],
                           tmb_data$max_n_groups, re_info$n_re)
    sigma_random_hat <- exp(par_full[names(par_full) == "log_sigma_random"])
    names(sigma_random_hat) <- re_info$group_names
  }

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

  # Add random effects information if multi-level model
  if (has_random) {
    # Compute composite abilities (theta_0 + random effects)
    composite_theta <- theta_0_hat  # Start with person deviations
    for (p in 1:n_persons) {
      for (re in 1:re_info$n_re) {
        group <- re_info$group_ids[p, re]
        if (group >= 0) {  # -1 indicates NA (partial nesting)
          composite_theta[p] <- composite_theta[p] + u_random_hat[group + 1, re]  # +1 for R indexing
        }
      }
    }

    # Compute ICCs (Intraclass Correlations)
    var_random <- sigma_random_hat^2
    var_person <- sigma_theta_hat^2
    var_total <- sum(var_random) + var_person + pi^2/3  # logistic variance

    icc_values <- c(var_random, var_person) / var_total
    names(icc_values) <- c(re_info$group_names, "Person")

    result$random_effects <- list(
      u_random = u_random_hat,
      sigma_random = sigma_random_hat,
      group_names = re_info$group_names,
      n_groups = re_info$n_groups,
      icc = icc_values,
      composite_theta = composite_theta
    )

    # Update class to indicate multi-level
    class(result) <- c("gllamm_irt_multilevel", "gllamm_irt", "gllamm")
  } else {
    result$random_effects <- NULL
    class(result) <- c("gllamm_irt", "gllamm")
  }

  return(result)
}
