#!/usr/bin/env Rscript

# Installs required R packages for the R-backed GLMM analyzer.

required <- c(
  "jsonlite", "readxl", "dplyr", "glmmTMB", "DHARMa",
  "emmeans", "performance", "ggplot2", "broom.mixed",
  "sandwich", "lmtest", "MASS", "rmarkdown"
)

repos <- c(CRAN = "https://cloud.r-project.org")

# Prefer a writable user library when the default system library is not writable.
default_lib <- .libPaths()[1]
if (file.access(default_lib, 2) != 0) {
  user_lib <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(user_lib)) {
    user_lib <- file.path(Sys.getenv("USERPROFILE"), "Documents", "R", "win-library", as.character(getRversion()))
  }

  if (!dir.exists(user_lib)) {
    dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  }

  if (!dir.exists(user_lib) || file.access(user_lib, 2) != 0) {
    cat("ERROR: No writable R library path found.\n")
    cat("Tried user library:\n")
    cat(paste0("- ", user_lib, "\n"))
    cat("Please run R/RStudio once as your user, then retry.\n")
    quit(save = "no", status = 1)
  }

  .libPaths(c(user_lib, .libPaths()))
  cat("Using user R library for package installation:\n")
  cat(paste0("- ", user_lib, "\n\n"))
}

installed <- rownames(installed.packages())
missing <- setdiff(required, installed)

if (length(missing) == 0) {
  cat("All required packages are already installed.\n")
  quit(save = "no", status = 0)
}

cat("Installing missing packages:\n")
cat(paste0("- ", missing, collapse = "\n"), "\n")

install.packages(missing, repos = repos, dependencies = TRUE, lib = .libPaths()[1])

installed_after <- rownames(installed.packages())
still_missing <- setdiff(required, installed_after)

if (length(still_missing) > 0) {
  cat("\nSome packages are still missing:\n")
  cat(paste0("- ", still_missing, collapse = "\n"), "\n")
  quit(save = "no", status = 1)
}

cat("\nAll required packages installed successfully.\n")
