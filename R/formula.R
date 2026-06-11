#' Parse GLLAMM formula
#'
#' Parses a formula into fixed effects and random effects components.
#' Supports lme4-style syntax without requiring lme4.
#'
#' @param formula A formula object with syntax: y ~ x + (terms | group)
#' @param data A data frame containing the variables
#'
#' @return A list with components:
#'   \item{fixed_formula}{Formula for fixed effects}
#'   \item{random_terms}{List of random effects specifications}
#'   \item{response_name}{Name of response variable}
#'
#' @keywords internal
parse_formula <- function(formula, data) {

  # Extract terms
  formula_terms <- terms(formula, data = data)

  # Find random effects terms (those with |)
  all_terms <- attr(formula_terms, "term.labels")

  # Initialize
  random_terms <- list()
  fixed_terms <- character(0)

  # Split terms into fixed and random; nested grouping (a/b/c) expands
  # lme4-style into (a), (a:b), (a:b:c)
  for (term in all_terms) {
    if (grepl("\\|", term)) {
      parsed_term <- parse_random_term(term, data)
      random_terms <- c(random_terms, expand_nested_random_term(parsed_term))
    } else {
      # This is a fixed effect
      fixed_terms <- c(fixed_terms, term)
    }
  }

  # Get response variable
  response_name <- as.character(formula[[2]])

  # Create fixed effects formula
  if (length(fixed_terms) > 0) {
    fixed_formula <- as.formula(paste(response_name, "~", paste(fixed_terms, collapse = " + ")))
  } else {
    fixed_formula <- as.formula(paste(response_name, "~ 1"))
  }

  list(
    fixed_formula = fixed_formula,
    random_terms = random_terms,
    response_name = response_name,
    original_formula = formula
  )
}


#' Parse a single random effects term
#'
#' @param term Character string with format "(terms | group)"
#' @param data Data frame for checking variables
#'
#' @return List with parsed random effects structure
#' @keywords internal
parse_random_term <- function(term, data) {

  # Remove whitespace
  term <- gsub("\\s+", "", term)

  # Check for parentheses and remove if present
  # Note: terms() may have already stripped them
  if (grepl("^\\(.*\\)$", term)) {
    # Remove outer parentheses
    term_inner <- gsub("^\\((.*)\\)$", "\\1", term)
  } else {
    # Already stripped by terms()
    term_inner <- term
  }

  # Check for || (uncorrelated random effects) BEFORE splitting
  uncorrelated <- grepl("\\|\\|", term_inner)

  # Split on | (for || this gives 3 parts, for | gives 2)
  if (uncorrelated) {
    # For ||, split and remove empty middle part
    parts <- strsplit(term_inner, "\\|")[[1]]
    parts <- parts[parts != ""]  # Remove empty strings
  } else {
    # For single |
    parts <- strsplit(term_inner, "\\|")[[1]]
  }

  if (length(parts) != 2) {
    stop("Random effects syntax must be: (term | group) or (term || group)")
  }

  re_terms <- parts[1]
  grouping <- parts[2]

  # Parse nested grouping (e.g., "school/class")
  if (grepl("/", grouping)) {
    grouping_vars <- strsplit(grouping, "/")[[1]]
    nested <- TRUE
  } else {
    grouping_vars <- grouping
    nested <- FALSE
  }

  # Parse random effects terms
  if (re_terms == "1") {
    # Random intercept only
    re_formula <- "~ 1"
  } else {
    # Random slopes (possibly with intercept)
    re_formula <- paste("~", re_terms)
  }

  list(
    formula = as.formula(re_formula),
    grouping = grouping_vars,
    nested = nested,
    uncorrelated = uncorrelated,
    original = term
  )
}


#' Expand a nested random-effects term into one term per level
#'
#' (x | a/b/c) becomes the lme4-equivalent (x | a) + (x | a:b) + (x | a:b:c).
#' Each expanded term carries an interaction grouping over the levels above it.
#'
#' @param parsed_term Output of parse_random_term()
#' @return List of parsed terms (length 1 for non-nested input)
#' @keywords internal
expand_nested_random_term <- function(parsed_term) {
  if (!isTRUE(parsed_term$nested) || length(parsed_term$grouping) <= 1) {
    parsed_term$grouping_vars <- parsed_term$grouping
    return(list(parsed_term))
  }

  vars <- parsed_term$grouping
  lapply(seq_along(vars), function(d) {
    t_d <- parsed_term
    t_d$grouping_vars <- vars[seq_len(d)]   # interaction over levels 1..d
    t_d$grouping <- paste(vars[seq_len(d)], collapse = ":")
    t_d$nested <- FALSE
    t_d$original <- paste0("(", deparse(parsed_term$formula[[2]]), " | ",
                           t_d$grouping, ")")
    t_d
  })
}


#' Extract model matrices
#'
#' Create design matrices for fixed and random effects
#'
#' @param parsed_formula Parsed formula object from parse_formula()
#' @param data Data frame
#'
#' @return List with X (fixed effects), Z (random effects), and grouping info
#' @keywords internal
make_model_matrices <- function(parsed_formula, data) {

  # Determine whether the response is available in `data` (it may be
  # absent when building matrices for prediction on new data)
  ff <- parsed_formula$fixed_formula
  resp_vars <- all.vars(ff[[2L]])
  has_response <- length(resp_vars) == 0 || all(resp_vars %in% names(data))

  # ---- Listwise deletion over ALL formula variables, up front ----
  # model.frame()/model.matrix() drop NA rows for X and y on their own,
  # but the random-effects design and grouping factor are built from
  # `data` directly: filtering once here keeps every piece aligned.
  used_vars <- if (has_response) all.vars(ff) else all.vars(ff[[3L]])
  for (rt in parsed_formula$random_terms) {
    used_vars <- union(used_vars, all.vars(rt$formula))
    used_vars <- union(used_vars, rt$grouping_vars %||% rt$grouping)
  }
  used_vars <- intersect(used_vars, names(data))
  n_original <- nrow(data)
  complete_idx <- which(stats::complete.cases(data[, used_vars,
                                                   drop = FALSE]))
  if (length(complete_idx) < nrow(data)) {
    if (has_response) {
      warning("Removing ", nrow(data) - length(complete_idx),
              " rows with missing values in model variables (listwise)")
    }
    data <- data[complete_idx, , drop = FALSE]
  }

  if (has_response) {
    # Fixed effects design matrix
    X <- model.matrix(ff, data = data)

    # Response
    y <- model.response(model.frame(ff, data = data))
  } else {
    tt <- stats::delete.response(stats::terms(ff, data = data))
    mf <- model.frame(tt, data = data)
    X <- model.matrix(tt, mf)
    y <- NULL
  }

  # Random effects design matrices
  Z_list <- list()
  groups_list <- list()
  n_random_coefs <- integer(0)

  for (i in seq_along(parsed_formula$random_terms)) {
    rt <- parsed_formula$random_terms[[i]]

    # Grouping factor: interaction over the term's grouping variables
    # (expanded nested terms carry the full path in grouping_vars)
    group_vars <- rt$grouping_vars %||% rt$grouping
    if (length(group_vars) > 1) {
      group_factor <- interaction(data[, group_vars], drop = TRUE)
    } else {
      group_factor <- factor(data[[group_vars]])
    }

    # Random effects design matrix for this term
    Z_formula <- rt$formula
    Z_i <- model.matrix(Z_formula, data = data)

    # Store number of random coefficients
    n_random_coefs <- c(n_random_coefs, ncol(Z_i))

    # Store Z matrix and grouping
    Z_list[[i]] <- Z_i
    groups_list[[i]] <- as.integer(group_factor) - 1L  # 0-indexed for C++
  }

  list(
    X = X,
    y = y,
    Z = Z_list,
    groups = groups_list,
    n_obs = nrow(X),
    n_fixed = ncol(X),
    n_random_terms = length(Z_list),
    n_random_coefs = n_random_coefs,
    n_groups = sapply(groups_list, function(g) length(unique(g))),
    complete_idx = complete_idx,
    n_original = n_original
  )
}


#' Align observation weights with listwise-deleted model data
#'
#' make_model_matrices() drops rows with missing values in any model
#' variable; user-supplied weights validated against the original data
#' length must be subset to the retained rows. Level-specific weight
#' lists are passed through untouched (they reference data columns and
#' are resolved downstream).
#'
#' @keywords internal
align_weights <- function(weights, model_data) {
  if (is.null(weights) || is.list(weights)) return(weights)
  if (length(weights) > model_data$n_obs) {
    return(weights[model_data$complete_idx])
  }
  weights
}


#' Validate formula
#'
#' Check that formula is properly specified
#'
#' @param formula Formula object
#' @param data Data frame
#'
#' @return TRUE if valid, stops with error otherwise
#' @keywords internal
validate_formula <- function(formula, data) {

  # Check formula is a formula
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }

  # Check for response variable
  if (length(formula) < 3) {
    stop("Formula must have a response variable: y ~ x + (1|group)")
  }

  # Check all variables exist in data
  vars <- all.vars(formula)
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  # Check for at least one random effect
  formula_str <- deparse(formula)
  if (!grepl("\\|", formula_str)) {
    warning("No random effects specified. Consider using glm() for fixed effects only models.")
  }

  TRUE
}


#' Drop the intercept column from a fixed-effects design matrix
#'
#' Cumulative-link (ordinal) models absorb the location into the thresholds;
#' keeping a free intercept alongside free thresholds leaves the model
#' unidentified (only their difference enters the likelihood).
#'
#' @param X Design matrix from model.matrix()
#' @return X without its "(Intercept)" column
#' @keywords internal
drop_intercept_column <- function(X) {
  ic <- which(colnames(X) == "(Intercept)")
  if (length(ic) > 0) {
    X <- X[, -ic, drop = FALSE]
  }
  X
}
