#' Extract Variance Components from Multi-Level Models
#'
#' Generic function to extract variance components from fitted multi-level models
#'
#' @param x A fitted model object
#' @param ... Additional arguments passed to methods
#'
#' @return Variance components (method-specific format)
#'
#' @export
VarCorr <- function(x, ...) {
  UseMethod("VarCorr")
}


#' @export
VarCorr.default <- function(x, ...) {
  stop("VarCorr is only available for multi-level models. ",
       "Model does not contain random effects.")
}


#' Extract Variance Components from Multi-Level IRT Models
#'
#' Extract variance components from fitted multi-level IRT models
#'
#' @param x A fitted multi-level IRT model (class gllamm_irt_multilevel)
#' @param ... Additional arguments (not used)
#'
#' @return A data frame with variance components
#'
#' @examples
#' \dontrun{
#' # Fit multi-level model
#' fit <- fit_irt(responses, model = "2PL",
#'                person_data = data, random = ~ (1 | class))
#'
#' # Extract variance components
#' VarCorr(fit)
#' }
#'
#' @export
VarCorr.gllamm_irt_multilevel <- function(x, ...) {
  if (is.null(x$random_effects)) {
    stop("Model does not contain random effects")
  }

  # Build variance components table
  result <- data.frame(
    Groups = c(x$random_effects$group_names, "Person", "Residual"),
    Variance = c(x$random_effects$sigma_random^2,
                 x$ability_sd^2,
                 pi^2/3),  # Logistic variance for binary IRT
    Std.Dev = c(x$random_effects$sigma_random,
                x$ability_sd,
                sqrt(pi^2/3))
  )

  rownames(result) <- NULL

  class(result) <- c("VarCorr.gllamm", "data.frame")
  return(result)
}


#' @export
print.VarCorr.gllamm <- function(x, digits = 4, ...) {
  if (is.data.frame(x)) {
    # Data-frame form (multi-level IRT models)
    cat("Variance Components:\n")
    # Round only numeric columns
    x$Variance <- round(x$Variance, digits)
    x$Std.Dev <- round(x$Std.Dev, digits)
    # Use print.data.frame explicitly to avoid recursion
    print.data.frame(x, row.names = FALSE)
  } else {
    # List form (GLMM fits): one covariance matrix per random-effect term
    cat("Random effects variance components:\n")
    for (i in seq_along(x)) {
      cat("\n Group:", names(x)[i], "\n")
      if (is.matrix(x[[i]])) {
        print(round(x[[i]], digits))
      } else {
        cat("  Variance:", round(x[[i]], digits), "\n")
        cat("  Std.Dev.:", round(sqrt(x[[i]]), digits), "\n")
      }
    }
  }
  invisible(x)
}


#' Compute Intraclass Correlation Coefficients
#'
#' Calculate ICCs for multi-level IRT models
#'
#' @param x A fitted multi-level IRT model
#' @param level Optional: specific level to compute ICC for. If NULL, returns all levels.
#' @param ... Additional arguments (not used)
#'
#' @return A named vector of ICCs, or a single ICC if level specified
#'
#' @examples
#' \dontrun{
#' # All ICCs
#' icc(fit)
#'
#' # Specific level
#' icc(fit, level = "class")
#' }
#'
#' @export
icc <- function(x, ...) {
  UseMethod("icc")
}


#' @export
icc.default <- function(x, ...) {
  stop("icc is only available for multi-level models. ",
       "Model does not contain random effects.")
}


#' @export
icc.gllamm_irt_multilevel <- function(x, level = NULL, ...) {
  if (is.null(x$random_effects)) {
    stop("Model does not contain random effects")
  }

  if (is.null(level)) {
    # Return all ICCs
    return(x$random_effects$icc)
  } else {
    # Return specific level
    if (!level %in% names(x$random_effects$icc)) {
      stop("Level '", level, "' not found. Available levels: ",
           paste(names(x$random_effects$icc), collapse = ", "))
    }
    return(x$random_effects$icc[level])
  }
}


#' Extract Random Effects
#'
#' Generic function to extract random effects from fitted models
#'
#' @param object A fitted model object
#' @param ... Additional arguments passed to methods
#'
#' @return Random effects (method-specific format)
#'
#' @export
ranef <- function(object, ...) {
  UseMethod("ranef")
}


#' @export
ranef.default <- function(object, ...) {
  stop("ranef is only available for multi-level models. ",
       "Model does not contain random effects.")
}


#' Extract Random Effects from Multi-Level IRT Models
#'
#' Extract estimated random effects (group-level deviations) from fitted models
#'
#' @param object A fitted multi-level IRT model
#' @param level Which random effect level to extract. If NULL, returns all levels.
#' @param ... Additional arguments (not used)
#'
#' @return A named vector or matrix of random effects
#'
#' @examples
#' \dontrun{
#' # Extract class effects
#' ranef(fit, level = "class")
#'
#' # Extract all random effects
#' ranef(fit)
#' }
#'
#' @export
ranef.gllamm_irt_multilevel <- function(object, level = NULL, ...) {
  if (is.null(object$random_effects)) {
    stop("Model does not contain random effects")
  }

  n_re <- length(object$random_effects$group_names)

  if (is.null(level)) {
    # Return all random effects as a list
    result <- list()
    for (i in 1:n_re) {
      group_name <- object$random_effects$group_names[i]
      n_groups <- object$random_effects$n_groups[i]
      re_values <- object$random_effects$u_random[1:n_groups, i]
      names(re_values) <- paste0(group_name, "_", 1:n_groups)
      result[[group_name]] <- re_values
    }
    return(result)
  } else {
    # Return specific level
    level_idx <- which(object$random_effects$group_names == level)
    if (length(level_idx) == 0) {
      stop("Level '", level, "' not found. Available levels: ",
           paste(object$random_effects$group_names, collapse = ", "))
    }

    n_groups <- object$random_effects$n_groups[level_idx]
    re_values <- object$random_effects$u_random[1:n_groups, level_idx]
    names(re_values) <- paste0(level, "_", 1:n_groups)
    return(re_values)
  }
}


#' Extract Person Abilities from IRT Models
#'
#' Extract estimated person abilities (theta) from IRT models.
#' For multi-level models, this returns the person-level deviations (theta_0).
#' Use \code{composite_theta} to get total abilities including random effects.
#'
#' @param object A fitted IRT model
#' @param composite Logical. For multi-level models, return composite abilities
#'   (theta_0 + random effects) instead of just theta_0. Default FALSE.
#' @param ... Additional arguments (not used)
#'
#' @return A named vector of person abilities
#'
#' @examples
#' \dontrun{
#' # Person-level deviations
#' abilities(fit)
#'
#' # Total abilities (including class effects)
#' abilities(fit, composite = TRUE)
#' }
#'
#' @export
abilities <- function(object, ...) {
  UseMethod("abilities")
}


#' @export
abilities.gllamm_irt <- function(object, ...) {
  return(object$person_abilities)
}


#' @export
abilities.gllamm_irt_multilevel <- function(object, composite = FALSE, ...) {
  if (composite) {
    return(object$random_effects$composite_theta)
  } else {
    return(object$person_abilities)
  }
}


#' Extract Coefficients from IRT Models
#'
#' Extract item parameters from fitted IRT models
#'
#' @param object A fitted IRT model
#' @param type Type of coefficients: "item" for item parameters,
#'   "person" for person abilities, "random" for random effects (multi-level only)
#' @param ... Additional arguments passed to specific methods
#'
#' @return Item parameters (data frame), person abilities (vector),
#'   or random effects (list)
#'
#' @examples
#' \dontrun{
#' # Item parameters
#' coef(fit, type = "item")
#'
#' # Person abilities
#' coef(fit, type = "person")
#'
#' # Random effects (multi-level only)
#' coef(fit, type = "random")
#' }
#'
#' @export
coef.gllamm_irt <- function(object, type = c("item", "person"), ...) {
  type <- match.arg(type)

  if (type == "item") {
    return(object$item_parameters)
  } else if (type == "person") {
    return(object$person_abilities)
  }
}


#' @export
coef.gllamm_irt_multilevel <- function(object, type = c("item", "person", "random"), ...) {
  type <- match.arg(type)

  if (type == "item") {
    return(object$item_parameters)
  } else if (type == "person") {
    return(object$person_abilities)
  } else if (type == "random") {
    if (is.null(object$random_effects)) {
      stop("Model does not contain random effects")
    }
    return(ranef.gllamm_irt_multilevel(object))
  }
}
