#' Package initialization
#'
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Check if TMB template is compiled
  dll_path <- file.path(libname, pkgname, "libs", paste0(pkgname, .Platform$dynlib.ext))

  if (!file.exists(dll_path)) {
    packageStartupMessage(
      "GLLAMMR: TMB templates not yet compiled.\n",
      "Run: TMB::compile(system.file('src/gllamm_gaussian.cpp', package = 'GLLAMMR'))\n",
      "Or use: GLLAMMR:::compile_gllamm_tmb()"
    )
  }
}


.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "GLLAMMR version ", utils::packageVersion("GLLAMMR"), "\n",
    "Generalized Linear Latent and Mixed Models\n",
    "Use citation('GLLAMMR') for citing this package in publications."
  )
}
