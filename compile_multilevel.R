#!/usr/bin/env Rscript

library(TMB)

cat("Compiling gllamm_irt_multilevel.hpp...\n")

# Remove any existing object files
if (file.exists("src/gllamm_irt_multilevel.o")) {
  file.remove("src/gllamm_irt_multilevel.o")
}

# Compile the template
result <- TMB::compile("src/gllamm_irt_multilevel.hpp")

if (result == 0) {
  cat("Compilation successful!\n")
  cat("Checking for .so file...\n")
  if (file.exists("src/gllamm_irt_multilevel.so")) {
    cat("✓ gllamm_irt_multilevel.so created successfully\n")
  } else {
    cat("✗ .so file not found\n")
  }
} else {
  cat("Compilation failed with code:", result, "\n")
}
