
############################################################
# Preparing multi-omics data for IntegratedLearner
#
# Goal:
# We have two different omics datasets:
#   1. Microbiome data
#   2. Metabolomics data
#
# Before we can perform an integrated analysis using
# IntegratedLearner, we need to organize these datasets into
# one coordinated object.
#
# We use a MultiAssayExperiment (MAE) because it is designed
# specifically for storing multiple assays measured on the
# same biological samples.
#
# Importantly, MAE does NOT simply combine all of our features
# into one huge matrix. Each assay remains separate, while the
# sample information is coordinated across assays.
############################################################


############################################################
# 1. Load required packages
############################################################

library(dplyr)
library(tidyverse)
library(MultiAssayExperiment)

# dplyr:
# Provides tools for manipulating data frames, including
# filtering, selecting, joining, and working with columns.
#
# tidyverse:
# Loads a collection of packages for data manipulation.
# Here, one useful function is column_to_rownames().
#
# MultiAssayExperiment:
# Provides the MultiAssayExperiment class.
#
# This is the key package for this step because our microbiome
# and metabolomics measurements are stored as DIFFERENT assays.
#
# A MultiAssayExperiment gives us one object that contains:
#
#               MultiAssayExperiment
#                       |
#          -----------------------------
#          |                           |
#     microbiome                  metabolome
#          |                           |
#      features                     features
#          \                           /
#           \                         /
#             shared sample metadata
#
# This coordinated structure is what we want before moving
# into an integrated multi-omics analysis.


############################################################
# 2. Import the microbiome dataset
############################################################

tse <- readRDS("hmpibd_tse.rds")

# readRDS() restores an R object that was previously saved
# using saveRDS().
#
# 'tse' contains the microbiome measurements.
#
# It is likely a TreeSummarizedExperiment or another
# SummarizedExperiment-like object.
#
# Conceptually, it contains:
#
#   rows    = microbial features / taxa
#   columns = samples
#
# It also contains sample metadata in colData(tse).


############################################################
# 3. Import the metabolomics dataset
############################################################

metabo.se <- readRDS("hmpibd_met.rds")

# 'metabo.se' contains the metabolomics measurements.
#
# Conceptually:
#
#   rows    = metabolites
#   columns = samples
#
# The microbiome and metabolomics datasets therefore contain
# completely different TYPES and NUMBERS of features.
#
# This is one reason we do not simply cbind() the two datasets.
# Instead, MAE allows each assay to remain in its natural form.


############################################################
# 4. Extract the sample metadata from each assay
############################################################

tse.coldata <- as.data.frame(colData(tse))
metabo.coldata <- as.data.frame(colData(metabo.se))

# colData() contains information ABOUT the samples rather than
# the measured features themselves.
#
# Examples include:
#
#   subject_id
#   age
#   gender
#   antibiotic use
#   disease / study condition
#
# We convert the metadata to ordinary data.frames because we
# are about to compare and merge the metadata from the two
# assays.
#
# Why do this?
#
# Our microbiome dataset may contain samples that do NOT have
# metabolomics measurements, and vice versa.
#
# For an integrated analysis we usually want samples that have
# measurements available in BOTH omics layers.


############################################################
# 5. Remove assay-specific colData
############################################################

colData(tse) <- NULL
colData(metabo.se) <- NULL

# We temporarily remove the metadata stored inside the
# individual assay objects.
#
# Why?
#
# We are about to construct ONE shared sample metadata table
# for the entire MultiAssayExperiment.
#
# Keeping separate copies of the metadata inside each assay
# can create duplication or inconsistencies.
#
# For example:
#
# microbiome colData:  age = 25
# metabolome colData:  age = 26
#
# Which one should be used?
#
# By explicitly creating the metadata ourselves, we make the
# relationship between assays clear.
#
# The shared metadata will eventually be stored here:
#
#               mae
#                |
#              colData
#                |
#        shared sample information


############################################################
# 6. Define the variables that identify matching samples
############################################################

join.vars <- c(
  "row.names",
  "subject_id",
  "age",
  "gender",
  "antibiotics_current_use",
  "study_condition"
)

# These are the variables that we require to match between
# the microbiome and metabolomics metadata.
#
# "row.names" represents the sample identifier.
#
# The other variables provide additional confirmation that the
# records correspond to the same biological sample / subject
# information.
#
# This is important because we do NOT want to accidentally
# integrate microbiome information from one sample with
# metabolomics information from another sample.


############################################################
# 7. Find samples present in BOTH datasets
############################################################

meta <- tse.coldata |>
  merge(
    metabo.coldata,
    by = join.vars,
    all = FALSE
  ) |>
  column_to_rownames("Row.names")

# This is one of the most important steps in the script.
#
# merge(..., all = FALSE) performs an INNER JOIN.
#
# In other words:
#
# microbiome samples:
#
#     A   B   C   D
#
# metabolome samples:
#
#         B   C   D   E
#
# after the merge:
#
#         B   C   D
#
# Samples A and E are excluded because they do not have data
# from both assays.
#
# Why is this important?
#
# IntegratedLearner needs coordinated observations across the
# different omics layers. We therefore create a common set of
# samples before constructing the MAE.
#
# column_to_rownames("row.names") then restores the sample IDs
# as the row names of our metadata table.
#
# After this step:
#
#     rownames(meta)
#
# gives us the exact sample IDs that are shared between the
# microbiome and metabolomics datasets.


############################################################
# Optional live-coding check
############################################################

dim(meta)
head(meta)

# During a live tutorial, these are useful sanity checks.
#
# dim(meta)
# tells us how many shared samples and metadata variables
# remain.
#
# head(meta)
# lets us quickly inspect the resulting sample metadata.


############################################################
# 8. Subset both assays to the common samples
############################################################

tse <- tse[, rownames(meta)]
metabo.se <- metabo.se[, rownames(meta)]

# Now we subset BOTH omics datasets using exactly the same
# sample IDs.
#
# Remember:
#
# rows    = features
# columns = samples
#
# Therefore:
#
#     tse[, rownames(meta)]
#
# means:
#
# "keep all microbiome features, but only keep the samples
# that appear in our shared metadata."
#
# We do the same thing for the metabolomics assay.
#
# After this step we want:
#
# colnames(tse)
# colnames(metabo.se)
# rownames(meta)
#
# to refer to the same samples.


############################################################
# Optional sanity checks
############################################################

identical(colnames(tse), rownames(meta))
identical(colnames(metabo.se), rownames(meta))

# Ideally both commands return TRUE.
#
# This gives us confidence that the data are correctly aligned.
#
# Sample alignment is critical in multi-omics analysis.
#
# If the order were wrong, we could accidentally associate:
#
#     microbiome from patient A
#
# with
#
#     metabolome from patient B
#
# which would invalidate the downstream analysis.


############################################################
# 9. Build the MultiAssayExperiment
############################################################

mae <- MultiAssayExperiment(
  experiments = ExperimentList(
    microbiome = tse,
    metabolome = metabo.se
  ),
  colData = meta
)

# Now we create the actual MultiAssayExperiment.
#
# ExperimentList() tells R that these are our separate assays:
#
#   microbiome = tse
#   metabolome = metabo.se
#
# The important idea is:
#
# WE ARE NOT MERGING THE FEATURES.
#
# We are coordinating the datasets.
#
# The resulting structure is approximately:
#
# mae
# |
# |-- experiments
# |      |
# |      |-- microbiome
# |      |      |-- microbial features x samples
# |      |
# |      |-- metabolome
# |             |-- metabolite features x samples
# |
# |-- colData
#        |
#        |-- age
#        |-- gender
#        |-- antibiotic use
#        |-- study condition
#        |-- etc.
#
# This is exactly what makes MAE useful for multi-omics.
#
# Different assays can:
#
#   - contain different types of features
#   - contain different numbers of features
#   - use different underlying experiment classes
#
# while still being connected through their biological samples.


############################################################
# 10. Inspect the MAE
############################################################

mae

# Printing the object gives us a summary of its contents.
#
# During the tutorial, this is a useful point to stop and
# explain what we have achieved.
#
# We started with:
#
#       microbiome.rds       metabolome.rds
#              |                   |
#              |                   |
#           separate datasets
#
# and have transformed them into:
#
#                  MAE
#                   |
#        ------------------------
#        |                      |
#   microbiome             metabolome
#        |                      |
#        -------- samples -------
#                   |
#                metadata
#
# The datasets are still separate assays, but their samples
# are now coordinated.


############################################################
# Optional: show the assays stored inside the MAE
############################################################

experiments(mae)

# This lets us see the experiments contained inside the
# MultiAssayExperiment.


############################################################
# Optional: show the shared metadata
############################################################

head(colData(mae))

# This shows the metadata associated with the coordinated
# samples.


############################################################
# 11. Save the prepared MAE
############################################################

saveRDS(mae, "hmpibd_mae.rds")

# Finally, we save the complete MultiAssayExperiment.
#
# This means the data preparation only has to be done once.
#
# In the next stage of the analysis we can simply use:
#
#     mae <- readRDS("hmpibd_mae.rds")
#
# rather than repeating all of the sample matching and
# preprocessing steps.


############################################################
# TAKE-HOME MESSAGE
############################################################

# Why did we create a MultiAssayExperiment?
#
# Because IntegratedLearner is performing an INTEGRATED
# analysis across multiple data modalities.
#
# We therefore need a reliable way to represent:
#
#    microbiome measurements
#             +
#    metabolomics measurements
#             +
#    information about the samples
#
# while making sure that samples are correctly matched.
#
# MAE gives us that coordinated container.
#
# So the workflow is:
#
# Raw microbiome data
#          \
#           \
#            ---> Match samples ---> MultiAssayExperiment
#           /                           |
#          /                            |
# Raw metabolomics data                 v
#                                IntegratedLearner
#
# The MAE is therefore not the final biological analysis.
#
# It is the DATA INFRASTRUCTURE that makes the downstream
# integrated analysis reliable and reproducible.
############################################################

