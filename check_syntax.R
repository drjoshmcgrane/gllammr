#!/usr/bin/env Rscript

# Check for syntax errors in all R files
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
errors <- c()

cat("Checking", length(r_files), "R files for syntax errors...\n\n")

for (f in r_files) {
  result <- tryCatch({
    parse(f)
    cat("✓", basename(f), "\n")
    NULL
  }, error = function(e) {
    paste(f, ":", conditionMessage(e))
  })
  if (!is.null(result)) {
    errors <- c(errors, result)
  }
}

cat("\n")
if (length(errors) > 0) {
  cat("❌ Syntax errors found:\n")
  cat(paste(errors, collapse = "\n"))
  cat("\n")
  quit(status = 1)
} else {
  cat("✅ All", length(r_files), "R files have valid syntax!\n")
  quit(status = 0)
}
