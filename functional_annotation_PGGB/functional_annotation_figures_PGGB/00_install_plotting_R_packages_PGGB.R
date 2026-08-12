#!/usr/bin/env Rscript

packages <- c("data.table", "ggplot2", "scales", "patchwork", "svglite")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

if (!length(missing)) {
  message("All plotting packages are already installed.")
  quit(status = 0)
}

message("Packages to install: ", paste(missing, collapse = ", "))
message("User library: ", Sys.getenv("R_LIBS_USER", unset = "R default user library"))
install.packages(missing, repos = "https://cloud.r-project.org", dependencies = TRUE)

still_missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) stop("Installation incomplete: ", paste(still_missing, collapse = ", "))
message("Plotting package installation completed.")
