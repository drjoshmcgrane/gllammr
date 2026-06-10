#' Parse Random Effects Formula
#'
#' Internal function to parse lme4-style random effects formulas
#'
#' @param formula Random effects formula (e.g., ~ (1 | class) + (1 | school))
#' @param data Data frame containing grouping variables
#' @return List of parsed random effects terms
#' @keywords internal
parse_random_formula <- function(formula, data) {
  if (is.null(formula)) {
    return(NULL)
  }

  # Convert to formula if needed
  if (!inherits(formula, "formula")) {
    stop("random must be a formula (e.g., ~ (1 | group))")
  }

  # Extract random effects bars
  # We'll do this manually to avoid lme4 dependency
  terms_list <- extract_random_terms(formula)

  if (length(terms_list) == 0) {
    stop("No random effects found in formula. Use syntax like ~ (1 | group)")
  }

  # Parse each term
  result <- lapply(terms_list, function(term) {
    parse_single_random_term(term, data)
  })

  # Expand (x | g) into independent intercept + slope terms
  expanded <- list()
  for (t in result) {
    if (isTRUE(t$expand_both)) {
      t_int <- t; t_int$type <- "intercept"; t_int$re_covariate <- NULL
      t_int$expand_both <- FALSE
      t_slp <- t; t_slp$expand_both <- FALSE
      expanded <- c(expanded, list(t_int, t_slp))
    } else {
      expanded <- c(expanded, list(t))
    }
  }
  result <- expanded

  # Expand nested notation (e.g., school/class)
  result <- expand_nested_terms(result, data)

  # Check for duplicates (same grouping AND same design column)
  keys <- sapply(result, function(x) {
    paste(x$group_var, x$re_covariate %||% "(Intercept)")
  })
  if (any(duplicated(keys))) {
    stop("Duplicate random-effects terms in formula")
  }

  return(result)
}

#' Extract random effects terms from formula
#' @keywords internal
extract_random_terms <- function(formula) {
  # Get right-hand side
  rhs <- formula[[2]]

  # Find all bars (|)
  terms_list <- list()

  extract_bars <- function(expr) {
    if (length(expr) == 1) {
      return(NULL)
    }

    if (expr[[1]] == as.name("|")) {
      # Found a bar term
      return(list(expr))
    }

    if (expr[[1]] == as.name("+")) {
      # Addition, recurse on both sides
      left <- extract_bars(expr[[2]])
      right <- extract_bars(expr[[3]])
      return(c(left, right))
    }

    if (expr[[1]] == as.name("(")) {
      # Parentheses, recurse on content
      return(extract_bars(expr[[2]]))
    }

    return(NULL)
  }

  terms_list <- extract_bars(rhs)

  if (is.null(terms_list)) {
    return(list())
  }

  return(terms_list)
}

#' Parse a single random effects term
#' @keywords internal
parse_single_random_term <- function(term, data) {
  # term is like: (1 | class) or (1 | school/class)

  # Left side (before |)
  lhs <- term[[2]]

  # Right side (after |)
  rhs <- term[[3]]

  # Left side determines the random-effect design:
  #   1        -> random intercept
  #   0 + x    -> random slope only (person covariate x)
  #   x        -> intercept + independent slope (expanded by the caller)
  re_covariate <- NULL
  expand_both <- FALSE
  if (identical(lhs, 1)) {
    # intercept
  } else if (is.call(lhs) && identical(as.character(lhs[[1]]), "+") &&
             identical(lhs[[2]], 0)) {
    re_covariate <- as.character(lhs[[3]])
  } else if (is.name(lhs)) {
    re_covariate <- as.character(lhs)
    expand_both <- TRUE
  } else {
    stop("Person-level random effects support (1 | g), (0 + x | g), and ",
         "(x | g) [intercept + independent slope]; got: ", deparse(lhs))
  }
  if (!is.null(re_covariate) && !re_covariate %in% names(data)) {
    stop("Random-slope covariate '", re_covariate,
         "' not found in person_data")
  }

  # Parse grouping variable
  group_expr <- rhs

  # Check for nesting (/)
  if (length(group_expr) > 1 && group_expr[[1]] == as.name("/")) {
    # Nested notation: school/class
    # We'll handle this in expand_nested_terms
    outer_var <- as.character(group_expr[[2]])
    inner_var <- as.character(group_expr[[3]])

    return(list(
      type = if (is.null(re_covariate)) "intercept" else "slope",
      re_covariate = re_covariate,
      expand_both = expand_both,
      group_var = outer_var,
      nested_in = inner_var,
      is_nested = TRUE
    ))
  } else {
    # Simple grouping variable
    group_var <- as.character(group_expr)

    # Check if exists in data
    if (!group_var %in% names(data)) {
      stop("Grouping variable '", group_var, "' not found in person_data")
    }

    return(list(
      type = if (is.null(re_covariate)) "intercept" else "slope",
      re_covariate = re_covariate,
      expand_both = expand_both,
      group_var = group_var,
      nested_in = NULL,
      is_nested = FALSE
    ))
  }
}

#' Expand nested terms
#' @keywords internal
expand_nested_terms <- function(terms, data) {
  result <- list()

  for (term in terms) {
    if (term$is_nested) {
      # Expand school/class to school + school:class
      outer_var <- term$group_var
      inner_var <- term$nested_in

      # Check both exist
      if (!outer_var %in% names(data)) {
        stop("Outer grouping variable '", outer_var, "' not found in person_data")
      }
      if (!inner_var %in% names(data)) {
        stop("Inner grouping variable '", inner_var, "' not found in person_data")
      }

      # Add outer effect (carries the term's design)
      result[[length(result) + 1]] <- list(
        type = term$type,
        re_covariate = term$re_covariate,
        group_var = outer_var,
        nested_in = NULL,
        is_nested = FALSE
      )

      # Create interaction variable name
      interaction_var <- paste0(outer_var, ":", inner_var)

      # Add nested effect
      result[[length(result) + 1]] <- list(
        type = term$type,
        re_covariate = term$re_covariate,
        group_var = interaction_var,
        nested_in = outer_var,
        is_nested = FALSE,
        is_interaction = TRUE,
        interaction_components = c(outer_var, inner_var)
      )
    } else {
      result[[length(result) + 1]] <- term
    }
  }

  return(result)
}

#' Create grouping factor matrix
#'
#' @param terms Parsed random effects terms
#' @param data Person-level data
#' @return List with group_ids matrix and group information
#' @keywords internal
create_grouping_matrix <- function(terms, data) {
  n_persons <- nrow(data)
  n_re <- length(terms)

  # Initialize matrices (-1 indicates NA in group_ids)
  group_ids <- matrix(-1L, n_persons, n_re)
  re_design <- matrix(1, n_persons, n_re)
  n_groups <- integer(n_re)
  group_names <- character(n_re)

  for (i in seq_along(terms)) {
    term <- terms[[i]]
    group_var <- term$group_var

    # Get grouping variable
    if (term$is_interaction %||% FALSE) {
      # Create interaction on the fly
      components <- term$interaction_components
      group_factor <- interaction(data[[components[1]]],
                                   data[[components[2]]],
                                   drop = TRUE, sep = ":")
    } else {
      # Simple grouping variable
      if (!group_var %in% names(data)) {
        stop("Grouping variable '", group_var, "' not found in person_data")
      }
      group_factor <- data[[group_var]]
    }

    # Convert to factor
    group_factor <- as.factor(group_factor)

    # Store as 0-indexed integers, NA becomes -1
    group_int <- as.integer(group_factor) - 1L
    group_int[is.na(group_factor)] <- -1L

    group_ids[, i] <- group_int
    n_groups[i] <- nlevels(group_factor)
    cov_i <- term$re_covariate
    group_names[i] <- if (is.null(cov_i)) group_var
                      else paste0(group_var, ":", cov_i)

    # Design column: 1 for intercepts, the person covariate for slopes
    if (!is.null(cov_i)) {
      re_design[, i] <- as.numeric(data[[cov_i]])
    }
  }

  return(list(
    group_ids = group_ids,
    n_groups = n_groups,
    group_names = group_names,
    re_design = re_design,
    n_re = n_re
  ))
}

#' Null-coalescing operator
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
