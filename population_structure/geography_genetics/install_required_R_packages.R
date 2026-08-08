#!/usr/bin/env Rscript
# Install required R packages into the active R library.

cran_packages <- c(
  "geosphere",
  "vegan",
  "ecodist",
  "ggplot2",
  "ggrepel",
  "sf",
  "rnaturalearth",
  "rnaturalearthdata",
  "ggspatial"
)

missing <- cran_packages[
  !vapply(cran_packages,requireNamespace,logical(1),quietly=TRUE)
]

if (length(missing)==0) {
  cat("All required R packages are already installed.\n")
} else {
  cat("Installing:",paste(missing,collapse=", "),"\n")
  install.packages(
    missing,
    repos="https://cloud.r-project.org",
    dependencies=TRUE
  )
}
