###############################################################################
# OMAWorkshop | Packages needed for the IntegratedLearner tutorial
###############################################################################

## --- Installers --------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

options(repos = BiocManager::repositories())


## --- Packages used directly in the tutorial ----------------------------------

packages <- c(
  "MultiAssayExperiment",      # MAE container: experiments(), colData(), etc.
  "SummarizedExperiment",      # assay(), rowData(), colData()
  "mia",                       # subsetByPrevalent(), transformAssay()
  "matrixStats",               # rowSds()
  "tidyverse",                 # data wrangling, column_to_rownames()
  "miaViz"                     # plotLoadings()
)


## --- Install missing packages ------------------------------------------------

missing <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing)) {
  BiocManager::install(
    missing,
    update = FALSE,
    ask = FALSE
  )
}


## --- Install IntegratedLearner ------------------------------------------------

if (!requireNamespace("IntegratedLearner", quietly = TRUE)) {
  remotes::install_github(
    "himelmallick/IntegratedLearner",
    upgrade = "never",
    dependencies = TRUE
  )
}


## --- Load packages ------------------------------------------------------------

to_load <- c(
  "MultiAssayExperiment",
  "SummarizedExperiment",
  "mia",
  "matrixStats",
  "tidyverse",
  "IntegratedLearner",
  "miaViz"
)

invisible(
  lapply(
    to_load,
    function(pkg) {
      suppressPackageStartupMessages(
        library(pkg, character.only = TRUE)
      )
    }
  )
)

message(
  "Workshop libraries loaded: ",
  paste(to_load, collapse = ", ")
)
