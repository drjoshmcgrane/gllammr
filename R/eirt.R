#' Fit Explanatory Item Response Theory Models
#'
#' Fit IRT models where item parameters are modeled as functions of item covariates.
#' For polytomous models, supports GRM (cumulative logit), PCM (adjacent-categories),
#' and GPCM (adjacent-categories with discrimination). PCM can include optional
#' threshold-level predictors via threshold_formula.
#'
#' @param response_matrix Matrix of item responses (persons x items)
#' @param item_data Data frame of item-level covariates (must have n_items rows)
#' @param difficulty_formula Formula for item location/difficulty regression (e.g., ~ word_freq).
#'   Model: b_i = W_diff \%*\% gamma + epsilon_b (if item_residuals = TRUE)
#' @param discrimination_formula Formula for discrimination regression (e.g., ~ item_type).
#'   Applies to 2PL, GRM, and GPCM models.
#'   Model: log(a_i) = W_disc \%*\% delta + epsilon_a (if item_residuals = TRUE)
#' @param threshold_formula Formula for threshold-specific regression (e.g., ~ abstractness).
#'   Only used for PCM/GPCM models. When specified, enables threshold-difficulty regression
#'   (Kim & Wilson 2019 LPCM framework): delta_im = b_i + sum_k xi_km * x_ik + e_im
#' @param weights Optional vector of case weights (fweights or pweights).
#'   Length must equal number of persons. Default: all observations weighted equally.
#'   Weights are expanded to observation-level (person-item) internally.
#' @param model IRT model type: "Rasch" or "2PL" (dichotomous), or
#'   "GRM", "PCM", "GPCM" (polytomous)
#' @param item_residuals Logical. If TRUE (default), includes item-specific residuals
#'   (LLTM + error). If FALSE, uses pure LLTM where item parameters are exactly
#'   predicted by covariates with no residuals.
#' @param start Optional starting values list
#' @param control Control parameters for nlminb optimization
#'
#' @return An object of class \code{gllamm_eirt}
#'
#' @details
#' **Dichotomous models:**
#' - **Rasch**: P(Y=1) = logit^{-1}(theta - b_i), where b_i = W_diff \%*\% gamma [+ epsilon_b]
#' - **2PL**: P(Y=1) = logit^{-1}(a_i * (theta - b_i)), where log(a_i) = W_disc \%*\% delta [+ epsilon_a]
#'
#' **Polytomous models (Kim & Wilson 2019; De Boeck & Wilson 2004 framework):**
#' - **GRM**: Cumulative logit with ordered thresholds. b_i provides item location;
#'   step_param provides threshold spacing. Supports discrimination_formula.
#' - **PCM**: Adjacent-categories logit, two-fold parameterization (MFRM approach):
#'   delta_im = b_i + s_im, where sum(s_im) = 0 across steps for each item.
#'   When threshold_formula is specified, uses threshold regression:
#'   delta_im = b_i + sum_k xi_km * x_ik + e_im (Kim & Wilson LPCM framework)
#' - **GPCM**: As PCM with item-specific discrimination a_i.
#'   Supports both discrimination_formula and threshold_formula.
#'
#' **Item residuals:** When item_residuals = TRUE (default), item parameters include
#' residual terms epsilon_b and epsilon_a (LLTM + error). When FALSE, uses pure LLTM
#' where parameters are exactly determined by covariates.
#'
#' @examples
#' \dontrun{
#' # Dichotomous EIRT with item and discrimination predictors
#' item_chars <- data.frame(
#'   word_frequency = rnorm(20),
#'   item_length = rpois(20, 5),
#'   item_type = factor(sample(c("concrete", "abstract"), 20, replace = TRUE))
#' )
#' fit_2pl <- fit_eirt(responses, item_data = item_chars,
#'                     difficulty_formula = ~ word_frequency + item_length,
#'                     discrimination_formula = ~ item_type,
#'                     model = "2PL")
#'
#' # Pure LLTM (no residuals)
#' fit_lltm <- fit_eirt(responses, item_data = item_chars,
#'                      difficulty_formula = ~ word_frequency,
#'                      model = "Rasch",
#'                      item_residuals = FALSE)
#'
#' # Polytomous PCM (adjacent-categories logit, Rasch family)
#' fit_pcm <- fit_eirt(poly_responses, item_data = item_chars,
#'                     difficulty_formula = ~ abstractness,
#'                     model = "PCM")
#'
#' # PCM with threshold predictors (LPCM framework)
#' fit_pcm_thresh <- fit_eirt(poly_responses, item_data = item_chars,
#'                            difficulty_formula = ~ abstractness,
#'                            threshold_formula = ~ cognitive_level,
#'                            model = "PCM")
#'
#' # GPCM with all predictors
#' fit_gpcm <- fit_eirt(poly_responses, item_data = item_chars,
#'                      difficulty_formula = ~ word_frequency,
#'                      discrimination_formula = ~ item_type,
#'                      threshold_formula = ~ cognitive_level,
#'                      model = "GPCM")
#' }
#'
#' @param person_data Optional data frame with person-level variables for multi-level models
#' @param random Optional random effects formula (e.g., ~ (1 | class))
#'
#' @export
fit_eirt <- function(response_matrix,
                     item_data,
                     difficulty_formula = ~ 1,
                     discrimination_formula = ~ 1,
                     threshold_formula = NULL,
                     person_data = NULL,
                     random = NULL,
                     weights = NULL,
                     model = c("Rasch", "2PL", "GRM", "PCM", "GPCM"),
                     item_residuals = TRUE,
                     start = NULL,
                     control = list()) {

  model <- match.arg(model)

  # Multi-level: parse person-level random effects (mirrors fit_irt)
  has_random <- !is.null(random)
  if (has_random && is.null(person_data)) {
    stop("person_data must be provided when random effects are specified")
  }
  re_info <- NULL
  if (has_random) {
    if (nrow(person_data) != nrow(response_matrix)) {
      stop("person_data must have same number of rows as response_matrix (",
           nrow(person_data), " vs ", nrow(response_matrix), ")")
    }
    re_terms <- parse_random_formula(random, person_data)
    re_info <- create_grouping_matrix(re_terms, person_data)
  }

  # Validate item_residuals
  if (!is.logical(item_residuals) || length(item_residuals) != 1) {
    stop("item_residuals must be a single logical value (TRUE or FALSE)")
  }

  # Validate inputs
  n_persons <- nrow(response_matrix)
  n_items <- ncol(response_matrix)

  # Validate weights
  if (!is.null(weights)) {
    if (length(weights) != n_persons) {
      stop("weights length (", length(weights), ") must equal number of persons (", n_persons, ")")
    }
    if (any(weights < 0, na.rm = TRUE)) {
      stop("weights must be non-negative")
    }
    if (any(is.na(weights))) {
      stop("weights cannot contain NA values")
    }
  }

  if (nrow(item_data) != n_items) {
    stop("item_data must have ", n_items, " rows (one per item). Found: ", nrow(item_data))
  }

  # Determine if polytomous
  is_polytomous <- model %in% c("GRM", "PCM", "GPCM")

  if (is_polytomous) {
    # Determine polytomous model type
    # For PCM/GPCM: use threshold regression (type 4) if threshold_formula provided
    if (model == "PCM" && !is.null(threshold_formula)) {
      poly_model_type <- 4L  # PCM with threshold regression (LPCM framework)
    } else if (model == "GPCM" && !is.null(threshold_formula)) {
      poly_model_type <- 4L  # GPCM with threshold regression
      # Note: GPCM with threshold regression uses same likelihood as LPCM but with discrimination
    } else {
      poly_model_type <- switch(model, "GRM" = 1L, "PCM" = 2L, "GPCM" = 3L)
    }

    model_type <- 2L  # Placeholder (only used in dichotomous branch)
    n_categories_per_item <- apply(response_matrix, 2, function(x) {
      length(unique(x[!is.na(x)]))
    })
    max_categories <- max(n_categories_per_item)
  } else {
    poly_model_type <- 0L
    model_type <- switch(model, "Rasch" = 1L, "2PL" = 2L)
    n_categories_per_item <- rep(2L, n_items)
    max_categories <- 2L
  }

  # Create design matrices
  W_diff <- model.matrix(difficulty_formula, data = item_data)
  W_disc <- model.matrix(discrimination_formula, data = item_data)
  p_diff <- ncol(W_diff)
  p_disc <- ncol(W_disc)

  # Threshold regression design matrix (for PCM/GPCM with threshold_formula)
  if (!is.null(threshold_formula) && model %in% c("PCM", "GPCM")) {
    W_threshold <- model.matrix(threshold_formula, data = item_data)
  } else {
    W_threshold <- matrix(0, n_items, 1)
  }
  p_thresh <- ncol(W_threshold)

  # Validate and auto-recode responses for polytomous models
  if (is_polytomous) {
    needs_recoding <- FALSE
    recode_items <- integer(0)

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
      # Check if properly coded 1-based
      else if (min_val != 1 || max_val != n_cats) {
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
  }

  # Convert to long format
  y_long <- as.vector(t(response_matrix))
  person_id <- rep(1:n_persons, each = n_items) - 1L
  item_id <- rep(1:n_items, times = n_persons) - 1L

  # Remove missing
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

  # TMB data
  tmb_data <- list(
    y = as.numeric(y_long),
    person_id = as.integer(person_id),
    item_id = as.integer(item_id),
    W_difficulty = as.matrix(W_diff),
    W_discrimination = as.matrix(W_disc),
    W_threshold = as.matrix(W_threshold),
    n_persons = as.integer(n_persons),
    n_items = as.integer(n_items),
    n_obs = as.integer(n_obs),
    weights = as.numeric(weights_long),
    model_type = as.integer(model_type),
    is_polytomous = as.integer(is_polytomous),
    poly_model_type = as.integer(poly_model_type),
    item_residuals = as.integer(item_residuals),
    uses_discrimination = as.integer((!is_polytomous && model == "2PL") ||
                                       (is_polytomous && poly_model_type %in% c(1L, 3L))),
    n_categories_per_item = as.integer(n_categories_per_item),
    max_categories = as.integer(max_categories)
  )

  # Multi-level structure
  if (has_random) {
    tmb_data$has_random <- 1L
    tmb_data$n_random_effects <- as.integer(re_info$n_re)
    tmb_data$group_ids <- as.matrix(re_info$group_ids)
    tmb_data$n_groups <- as.integer(re_info$n_groups)
  } else {
    tmb_data$has_random <- 0L
    tmb_data$n_random_effects <- 0L
    tmb_data$group_ids <- matrix(0L, 0, 0)
    tmb_data$n_groups <- integer(0)
  }

  # Parameter dimensions
  n_step_cols <- max(1L, as.integer(max_categories) - 1L)

  # Initialize parameters
  if (is.null(start)) {
    step_param_init <- matrix(0, n_items, n_step_cols)
    if (model == "GRM" && max_categories > 2) {
      # Evenly-spaced initial thresholds around item location
      step_param_init[, 1] <- 0
      for (k in seq(2, n_step_cols)) {
        step_param_init[, k] <- log(2)
      }
    }

    # For GPCM, start with tighter discrimination variance to avoid instability
    init_log_sigma_a <- if (model == "GPCM") log(0.2) else log(0.5)

    tmb_params <- list(
      theta = rep(0, n_persons),
      gamma = rep(0, p_diff),
      delta = rep(0, p_disc),
      epsilon_b = rep(0, n_items),
      epsilon_a = rep(0, n_items),
      log_sigma_epsilon_b = log(0.5),
      log_sigma_epsilon_a = init_log_sigma_a,
      log_sigma_theta = log(1.0),
      step_param = step_param_init,
      xi = matrix(0, p_thresh, n_step_cols),
      e_step = matrix(0, n_items, n_step_cols),
      log_sigma_e_step = log(0.5),
      u_random = if (has_random) {
        matrix(0, max(re_info$n_groups), re_info$n_re)
      } else {
        matrix(0, 1, 1)
      },
      log_sigma_random = if (has_random) rep(log(0.5), re_info$n_re) else 0
    )
  } else {
    tmb_params <- start
  }

  # Build map to fix parameters the chosen model never reads; leaving them
  # free creates flat likelihood directions and singular Hessians
  map_list <- list()

  if (poly_model_type == 4L) {
    # LPCM: fix step_param (unused), but keep difficulty_formula active
    map_list$step_param <- factor(rep(NA, n_items * n_step_cols))
  } else {
    # All non-LPCM models: threshold-regression machinery unused
    map_list$xi <- factor(rep(NA, p_thresh * n_step_cols))
    map_list$e_step <- factor(rep(NA, n_items * n_step_cols))
    map_list$log_sigma_e_step <- factor(NA)
    if (!is_polytomous) {
      map_list$step_param <- factor(rep(NA, n_items * n_step_cols))
    }
  }

  # Discrimination machinery is only read by 2PL, GRM, and GPCM
  uses_discrimination <- tmb_data$uses_discrimination == 1L
  if (!uses_discrimination) {
    map_list$delta <- factor(rep(NA, p_disc))
    map_list$epsilon_a <- factor(rep(NA, n_items))
    map_list$log_sigma_epsilon_a <- factor(NA)
  }

  # Item residuals off: epsilon parameters unused
  if (!item_residuals) {
    map_list$epsilon_b <- factor(rep(NA, n_items))
    map_list$log_sigma_epsilon_b <- factor(NA)
    if (uses_discrimination) {
      map_list$epsilon_a <- factor(rep(NA, n_items))
      map_list$log_sigma_epsilon_a <- factor(NA)
    }
  }

  # Multi-level parameters are dead unless random effects are present
  if (!has_random) {
    map_list$u_random <- factor(rep(NA, length(tmb_params$u_random)))
    map_list$log_sigma_random <- factor(rep(NA, length(tmb_params$log_sigma_random)))
  }

  # Latent variables to integrate out (only those actually in the model)
  random_effects <- "theta"
  if (item_residuals) {
    random_effects <- c(random_effects, "epsilon_b")
    if (uses_discrimination) {
      random_effects <- c(random_effects, "epsilon_a")
    }
  }
  if (poly_model_type == 4L) {
    random_effects <- c(random_effects, "e_step")
  }
  if (has_random) {
    random_effects <- c(random_effects, "u_random")
  }

  # Create TMB object
  tmb_data$model_name <- "eirt"
  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    random = random_effects,
    map = map_list,
    DLL = "GLLAMMR",
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

  # Standard errors
  sdr <- try(TMB::sdreport(obj), silent = TRUE)

  # Extract fixed-effect parameters
  par_full <- obj$env$last.par.best

  gamma_hat <- par_full[names(par_full) == "gamma"]
  names(gamma_hat) <- colnames(W_diff)

  # delta is mapped off for models that never read discrimination
  delta_hat <- par_full[names(par_full) == "delta"]
  if (length(delta_hat) == 0) {
    delta_hat <- rep(NA_real_, p_disc)
  }
  names(delta_hat) <- colnames(W_disc)

  # Extract ADREPORT quantities
  xi_hat <- NULL
  sigma_e_step_hat <- NA

  if (!inherits(sdr, "try-error")) {
    sdr_summary <- summary(sdr, "report")

    difficulty_hat <- sdr_summary[rownames(sdr_summary) == "difficulty", "Estimate"]
    names(difficulty_hat) <- paste0("Item", seq_len(n_items))

    discrimination_hat <- sdr_summary[rownames(sdr_summary) == "discrimination", "Estimate"]
    names(discrimination_hat) <- paste0("Item", seq_len(n_items))

    theta_hat <- sdr_summary[rownames(sdr_summary) == "theta", "Estimate"]
    names(theta_hat) <- paste0("Person", seq_len(n_persons))

    sigma_epsilon_b_hat <- sdr_summary[rownames(sdr_summary) == "sigma_epsilon_b", "Estimate"]
    sigma_epsilon_a_hat <- sdr_summary[rownames(sdr_summary) == "sigma_epsilon_a", "Estimate"]
    sigma_theta_hat <- sdr_summary[rownames(sdr_summary) == "sigma_theta", "Estimate"]

    if (poly_model_type == 4L) {
      xi_rows <- rownames(sdr_summary) == "xi"
      if (any(xi_rows)) {
        xi_hat <- matrix(sdr_summary[xi_rows, "Estimate"],
                         nrow = p_thresh, ncol = max_categories - 1L)
        rownames(xi_hat) <- colnames(W_threshold)
        colnames(xi_hat) <- paste0("Threshold", seq_len(max_categories - 1L))
      }
      se_rows <- rownames(sdr_summary) == "sigma_e_step"
      if (any(se_rows)) {
        sigma_e_step_hat <- sdr_summary[se_rows, "Estimate"]
      }
    }
  } else {
    # Fallback: extract from raw parameters
    difficulty_hat <- rep(NA_real_, n_items)
    discrimination_hat <- rep(NA_real_, n_items)
    theta_hat <- par_full[names(par_full) == "theta"]
    sigma_epsilon_b_hat <- exp(par_full[names(par_full) == "log_sigma_epsilon_b"])
    sigma_epsilon_a_hat <- exp(par_full[names(par_full) == "log_sigma_epsilon_a"])
    sigma_theta_hat <- exp(par_full[names(par_full) == "log_sigma_theta"])
  }

  # Construct result
  result <- list(
    model = model,
    regression_coefficients = list(
      difficulty = gamma_hat,
      discrimination = delta_hat,
      threshold = xi_hat
    ),
    item_parameters = list(
      difficulty = difficulty_hat,
      discrimination = discrimination_hat
    ),
    person_abilities = theta_hat,
    ability_sd = sigma_theta_hat,
    residual_sd = list(
      difficulty = sigma_epsilon_b_hat,
      discrimination = sigma_epsilon_a_hat,
      threshold = sigma_e_step_hat
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
    n_categories_per_item = n_categories_per_item,
    max_categories = max_categories,
    is_polytomous = is_polytomous,
    poly_model_type = poly_model_type,
    formulas = list(
      difficulty = difficulty_formula,
      discrimination = discrimination_formula,
      threshold = threshold_formula
    ),
    item_data = item_data,
    tmb_obj = obj,
    tmb_opt = opt,
    tmb_sdr = sdr
  )

  class(result) <- c("gllamm_eirt", "gllamm")

  # Multi-level: group-level random effects, composite abilities, ICCs
  if (has_random) {
    max_n_groups <- max(re_info$n_groups)
    u_random_hat <- matrix(par_full[names(par_full) == "u_random"],
                           max_n_groups, re_info$n_re)
    sigma_random_hat <- exp(par_full[names(par_full) == "log_sigma_random"])
    names(sigma_random_hat) <- re_info$group_names

    # Composite ability: person deviation + group effects
    composite_theta <- theta_hat
    for (p in seq_len(n_persons)) {
      for (re in seq_len(re_info$n_re)) {
        group <- re_info$group_ids[p, re]
        if (group >= 0) {  # -1 indicates NA (partial nesting)
          composite_theta[p] <- composite_theta[p] + u_random_hat[group + 1, re]
        }
      }
    }

    # ICCs on the latent logistic scale (same convention as fit_irt)
    var_random <- sigma_random_hat^2
    var_person <- sigma_theta_hat^2
    var_total <- sum(var_random) + var_person + pi^2 / 3
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
    class(result) <- c("gllamm_eirt_multilevel", "gllamm_eirt", "gllamm")
  }

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

  # Threshold regression (for PCM/GPCM with threshold_formula)
  if (!is.null(x$regression_coefficients$threshold)) {
    n_thresholds <- x$max_categories - 1L
    cat("Threshold regression (", n_thresholds, " thresholds):\n", sep = "")
    cat("  Formula:", deparse(x$formulas$threshold), "\n")
    cat("  Coefficients (rows = predictors, cols = thresholds):\n")
    print(round(x$regression_coefficients$threshold, 3))
    cat("  Residual SD:", round(x$residual_sd$threshold, 3), "\n\n")
  }

  cat("Ability distribution:\n")
  cat("  Mean:", round(mean(x$person_abilities), 3), "\n")
  cat("  SD:", round(sd(x$person_abilities), 3), "\n")
  cat("  Estimated SD:", round(x$ability_sd, 3), "\n\n")

  if (!is.null(x$random_effects)) {
    cat("Group-level variance components:\n")
    for (g in seq_along(x$random_effects$sigma_random)) {
      cat("  ", x$random_effects$group_names[g],
          ": SD = ", round(x$random_effects$sigma_random[g], 3),
          " (", x$random_effects$n_groups[g], " groups)\n", sep = "")
    }
    cat("  ICC:", paste(names(x$random_effects$icc),
                        round(x$random_effects$icc, 3),
                        sep = " = ", collapse = ", "), "\n\n")
  }

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
    difficulty = object$item_parameters$difficulty[seq_len(n_show)],
    discrimination = object$item_parameters$discrimination[seq_len(n_show)]
  )
  print(round(item_params_df, 3))

  if (object$n_items > 10) {
    cat("  ... (", object$n_items - 10, " more items)\n")
  }

  invisible(object)
}
