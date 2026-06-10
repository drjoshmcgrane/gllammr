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

  # Split terms into fixed and random
  for (term in all_terms) {
    if (grepl("\\|", term)) {
      # This is a random effects term
      random_terms <- c(random_terms, list(parse_random_term(term, data)))
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

  # Fixed effects design matrix
  X <- model.matrix(parsed_formula$fixed_formula, data = data)

  # Response
  y <- model.response(model.frame(parsed_formula$fixed_formula, data = data))

  # Random effects design matrices
  Z_list <- list()
  groups_list <- list()
  n_random_coefs <- integer(0)

  for (i in seq_along(parsed_formula$random_terms)) {
    rt <- parsed_formula$random_terms[[i]]

    # Get grouping variable(s)
    if (rt$nested) {
      # Create nested grouping variable
      group_vars <- rt$grouping
      group_factor <- interaction(data[, group_vars], drop = TRUE)
    } else {
      group_factor <- factor(data[[rt$grouping]])
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
    n_obs = length(y),
    n_fixed = ncol(X),
    n_random_terms = length(Z_list),
    n_random_coefs = n_random_coefs,
    n_groups = sapply(groups_list, function(g) length(unique(g)))
  )
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
