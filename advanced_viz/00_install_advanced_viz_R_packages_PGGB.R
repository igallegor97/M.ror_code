#!/usr/bin/env Rscript
packages <- c("data.table", "ggplot2", "scales", "patchwork", "svglite",
              "ComplexUpset", "ggalluvial", "ggrepel", "ggdist")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (!length(missing)) {
  message("All advanced-viz packages are already installed.")
  quit(status = 0)
}
message("Installing: ", paste(missing, collapse = ", "))
install.packages(missing, repos = "https://cloud.r-project.org", dependencies = TRUE)
still_missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) stop("Installation incomplete: ", paste(still_missing, collapse = ", "))

