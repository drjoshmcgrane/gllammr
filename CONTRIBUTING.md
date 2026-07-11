# Contributing to GLLAMMR

Thank you for your interest in contributing to GLLAMMR! This document
provides guidelines and instructions for contributing.

## Types of Contributions

We welcome several types of contributions:

1.  **Bug reports** - Report bugs via GitHub Issues
2.  **Bug fixes** - Submit pull requests fixing known bugs
3.  **New features** - Propose and implement new functionality
4.  **Documentation** - Improve documentation, vignettes, examples
5.  **Tests** - Add tests to improve coverage
6.  **Performance** - Optimize code for speed or memory

## Getting Started

### Development Setup

1.  **Fork and clone the repository**

``` bash
git clone https://github.com/drjoshmcgrane/gllammr.git
cd GLLAMMR
```

2.  **Install dependencies**

``` r

# In R
install.packages(c("TMB", "Matrix", "testthat", "roxygen2", "devtools"))
```

3.  **Compile TMB templates**

``` r

# In R
TMB::compile("src/gllamm_gaussian.cpp")
```

4.  **Load the package**

``` r

devtools::load_all()
```

### Development Workflow

1.  Create a new branch for your changes:

``` bash
git checkout -b feature/my-new-feature
```

2.  Make your changes

3.  Document your code with roxygen2 comments

4.  Add tests in `tests/testthat/`

5.  Run checks:

``` r

devtools::test()        # Run tests
devtools::check()       # Full package check
devtools::document()    # Update documentation
```

6.  Commit your changes:

``` bash
git add .
git commit -m "Add feature: description of changes"
```

7.  Push and create a pull request

## Code Style

### R Code

- Follow the [tidyverse style guide](https://style.tidyverse.org/)
- Use meaningful variable names
- Comment complex logic
- Keep functions focused and modular
- Maximum line length: 100 characters

Example:

``` r

#' Fit GLLAMM model
#'
#' @param formula Model formula
#' @param data Data frame
#' @return gllamm object
my_function <- function(formula, data) {
  # Validate inputs
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }

  # Main logic here
  result <- do_something(formula, data)

  return(result)
}
```

### C++ Code (TMB Templates)

- Follow TMB conventions
- Comment each major section
- Use clear variable names
- Indent with 2 spaces

Example:

``` cpp
// Compute log-likelihood for Gaussian response
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(y);
  DATA_MATRIX(X);

  // Parameters
  PARAMETER_VECTOR(beta);

  // Initialize negative log-likelihood
  Type nll = 0.0;

  // Likelihood contribution
  for (int i = 0; i < y.size(); i++) {
    Type mu = (X.row(i) * beta).sum();
    nll -= dnorm(y[i], mu, sigma, true);
  }

  return nll;
}
```

## Testing

All new features and bug fixes should include tests.

### Writing Tests

Tests use the `testthat` package:

``` r

test_that("function does what it should", {
  # Setup
  data <- data.frame(y = 1:10, x = 1:10, g = rep(1:2, 5))

  # Test
  result <- my_function(data)

  # Assertions
  expect_equal(result$n, 10)
  expect_true(result$converged)
  expect_s3_class(result, "gllamm")
})
```

### Test Guidelines

- Each test should test one thing
- Use descriptive test names
- Include tests for edge cases
- Test error handling
- Aim for \>90% code coverage

## Documentation

### Function Documentation

Use roxygen2 format:

``` r

#' One-line description
#'
#' Detailed description with multiple paragraphs if needed.
#'
#' @param param1 Description of parameter 1
#' @param param2 Description of parameter 2
#'
#' @return Description of return value
#'
#' @examples
#' \dontrun{
#' result <- my_function(x, y)
#' }
#'
#' @export
my_function <- function(param1, param2) {
  # Implementation
}
```

### Vignettes

Vignettes should: \* Start with a clear motivation/use case \* Include
runnable examples \* Explain both the “how” and “why” \* Use real or
realistic data \* Be concise but comprehensive

## Validation Requirements

New features should be validated against existing implementations:

### For GLMM Features

- Compare to lme4 where applicable
- Match within 1% for coefficients
- Match within 2% for standard errors

### For IRT Features

- Compare to mirt package
- Match within 2% for item parameters
- Use standard datasets (e.g., LSAT, Science)

### For Latent Class Features

- Compare to poLCA
- Match class probabilities within 2%

### For All Features

- Include simulation-recovery test
- Verify parameter recovery with known true values
- Test with various sample sizes

## Pull Request Process

1.  **Update documentation**: Run `devtools::document()`

2.  **Run all checks**: Ensure `devtools::check()` passes with no
    errors, warnings, or notes

3.  **Update NEWS.md**: Add entry describing your changes

4.  **Update tests**: Ensure all tests pass with `devtools::test()`

5.  **Write clear PR description**:

    - What problem does this solve?
    - What changes were made?
    - How was it tested?
    - Are there breaking changes?

6.  **Be responsive**: Address review comments promptly

## Reporting Bugs

Use GitHub Issues to report bugs. Good bug reports include:

1.  **Clear title**: Briefly describe the problem
2.  **Expected vs actual behavior**: What should happen vs what does
    happen
3.  **Minimal reproducible example**: Code that reproduces the issue
4.  **Session info**: Output of
    [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)
5.  **Error messages**: Full error messages and stack traces

Example:

``` markdown
## Bug: gllamm fails with crossed random effects

### Expected behavior
`gllamm(y ~ x + (1|g1) + (1|g2), data)` should fit a model with crossed random effects.

### Actual behavior
Error: "Currently only single random effects term supported"

### Reproducible example
\```r
library(GLLAMMR)
data <- data.frame(
  y = rnorm(100),
  x = rnorm(100),
  g1 = rep(1:10, each=10),
  g2 = rep(1:10, times=10)
)
fit <- gllamm(y ~ x + (1|g1) + (1|g2), data)
\```

### Session info
\```
R version 4.3.0 (2023-04-21)
GLLAMMR version 0.1.0
...
\```
```

## Feature Requests

Feature requests should include:

1.  **Use case**: Why is this feature needed?
2.  **Proposed API**: How should users interact with it?
3.  **Examples**: Show example usage
4.  **Alternatives**: What alternatives exist?

## Code of Conduct

Please be respectful and constructive in all interactions. We are
committed to providing a welcoming and inclusive environment for all
contributors.

## Questions?

If you have questions about contributing, please open a GitHub Issue
with the “question” label or contact the maintainers.

## Recognition

Contributors will be acknowledged in: \* Package DESCRIPTION file \*
Release notes \* Package documentation

Thank you for contributing to GLLAMMR!
