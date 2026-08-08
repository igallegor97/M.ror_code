#!/bin/bash

#$ -N comm_geo_fig_v2
#$ -q all.q
#$ -cwd
#$ -pe smp 2
#$ -V
#$ -o /Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/community_geography/logs/comm_geo_fig_v2_$JOB_ID.out
#$ -e /Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/community_geography/logs/comm_geo_fig_v2_$JOB_ID.err

set -euo pipefail

CODE_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/community_geography"

LOG_DIR="/Storage/data1/isabella.gallego/MAESTRIA/code/population_structure/community_geography/logs"

mkdir -p "${LOG_DIR}"

module purge
module load R/4.5.1

Rscript -e '
required <- c(
  "ggplot2",
  "scales",
  "gridExtra",
  "ggrepel"
)

missing <- required[
  !vapply(
    required,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing) > 0) {
  stop(
    paste(
      "Missing R packages:",
      paste(missing, collapse = ", ")
    )
  )
}
'

cd "${CODE_DIR}"

bash "${CODE_DIR}/run_community_geography_figures_v2.sh"
