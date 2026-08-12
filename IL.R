
###############################################################################
#
# Goal:
# Use microbiome + metabolomics data together to predict whether
# a sample comes from an individual with IBD.
#
# Overall workflow:
#
#   MultiAssayExperiment
#          |
#          |-- Microbiome
#          |
#          |-- Metabolome
#          |
#          v
#      Preprocess
#          |
#          v
#   Train / validation split
#          |
#          v
#   IntegratedLearner
#       /        \
#      /          \
# Early fusion   Late fusion
# (concatenate)  (stack predictions)
#      \          /
#       \        /
#        Compare AUC
#
###############################################################################


###############################################################################
# 0. LOAD THE MULTI-OMICS DATA
###############################################################################

# We created this MultiAssayExperiment in the previous session.
#
# A MultiAssayExperiment (MAE) is useful because it allows us to keep
# different omics datasets as separate assays while linking them through
# the same biological samples and sample metadata.

mae <- readRDS("hmpibd_mae.rds")


# What omics layers are inside the MAE?
experiments(mae)

# We should see something like:
#
#   microbiome
#   metabolome
#
# These datasets contain completely different types of features,
# but they represent the same set of samples.


# Look at the outcome we ultimately want to predict.
table(colData(mae)$study_condition)

# For example, we may have:
#
#   nonIBD / healthy
#   IBD
#
# This is going to become our binary machine-learning outcome.


###############################################################################
# 1. PREVALENCE FILTERING
###############################################################################

# First we remove features that are rarely observed.
#
# Why?
#
# Multi-omics datasets can contain thousands of features.
# Some microbial taxa or metabolites may occur in only a tiny number
# of samples.
#
# Very rare features:
#
#   - provide little information for prediction
#   - increase dimensionality
#   - can make models noisier
#   - increase computation
#
# Here we retain features with prevalence >= 10%.


# ---- Microbiome --------------------------------------------------------------

mae[["microbiome"]] <- subsetByPrevalent(
  mae[["microbiome"]],
  assay.type = "relative_abundance",
  prevalence = 0.1,
  include.lowest = TRUE
)


# ---- Metabolome --------------------------------------------------------------

mae[["metabolome"]] <- subsetByPrevalent(
  mae[["metabolome"]],
  assay.type = "nmr",
  prevalence = 0.1,
  include.lowest = TRUE
)


# A useful live-coding check:
dim(mae[["microbiome"]])
dim(mae[["metabolome"]])

# Remember:
#
# rows    = features
# columns = samples
#
# So the number of rows tells us how many features remain.


###############################################################################
# 2. REMOVE FEATURES WITH ALMOST NO VARIATION
###############################################################################

# Prevalence tells us whether a feature is PRESENT often enough.
#
# But that does not necessarily mean it varies between samples.
#
# Imagine a microbial feature with approximately the same abundance in
# every single person:
#
#       sample 1   0.050
#       sample 2   0.050
#       sample 3   0.051
#       sample 4   0.050
#
# That feature is unlikely to help distinguish IBD from healthy samples.
#
# So we calculate the standard deviation of every feature.


# ---- Microbiome feature variation --------------------------------------------

rowData(mae[["microbiome"]])$sd <- matrixStats::rowSds(
  assay(
    mae[["microbiome"]],
    "relative_abundance"
  )
)


# ---- Metabolomics feature variation ------------------------------------------

rowData(mae[["metabolome"]])$sd <- matrixStats::rowSds(
  assay(
    mae[["metabolome"]],
    "nmr"
  )
)


# Remove features with very little or no variation.
#
# The cutoff differs because the assays are measured on different scales.

mae[["microbiome"]] <- mae[["microbiome"]][
  rowData(mae[["microbiome"]])$sd > 0.001,
]

mae[["metabolome"]] <- mae[["metabolome"]][
  rowData(mae[["metabolome"]])$sd > 0,
]


# Again, inspect how many features remain.

dim(mae[["microbiome"]])
dim(mae[["metabolome"]])


###############################################################################
# 3. TRANSFORM THE DATA
###############################################################################

# Now we transform each omics layer.
#
# This is important because microbiome and metabolomics data have
# very different distributions and measurement properties.
#
# We should NOT assume that the same transformation is appropriate
# for both datasets.


###############################################################################
# 3A. MICROBIOME: CLR TRANSFORMATION
###############################################################################

# Microbiome relative abundances are COMPOSITIONAL.
#
# The values describe proportions of a total community.
#
# If one organism increases in relative abundance, the relative
# abundance of other organisms can appear to decrease simply because
# everything must share the same total.
#
# CLR = centered log-ratio transformation.
#
# The CLR transformation converts relative information into
# log-ratios, which are generally more appropriate for modelling
# compositional microbiome data.


mae[["microbiome"]] <- transformAssay(
  mae[["microbiome"]],
  assay.type = "relative_abundance",
  method = "clr",
  pseudocount = TRUE
)


###############################################################################
# 3B. METABOLOME: LOG10 TRANSFORMATION
###############################################################################

# Metabolite concentrations can be highly right-skewed.
#
# For example:
#
#   1
#   2
#   5
#   20
#   400
#
# A log transformation compresses very large values and often produces
# a distribution that is easier for downstream models to work with.
#
# We use a pseudocount because:
#
#   log10(0)
#
# is undefined.


mae[["metabolome"]] <- transformAssay(
  mae[["metabolome"]],
  assay.type = "nmr",
  method = "log10",
  pseudocount = 1
)


# Our transformed assay names should now include:
#
#   microbiome  -> "clr"
#   metabolome  -> "log10"

assayNames(mae[["microbiome"]])
assayNames(mae[["metabolome"]])


###############################################################################
# 4. PREPARE THE OUTCOME FOR INTEGRATEDLEARNER
###############################################################################

# IntegratedLearner needs to know two especially important things:
#
#   1. What are we trying to predict?
#   2. Which subject does each observation belong to?
#
# We create convenient columns called:
#
#   Y
#   subjectID


cd <- colData(mae)


# Check the available metadata variables.
colnames(cd)


# Convert disease status into a binary outcome:
#
#   healthy / non-IBD = 0
#   IBD               = 1

cd$Y <- as.numeric(cd$study_condition == "IBD")


# Check that the coding worked.
table(cd$study_condition, cd$Y)


# IntegratedLearner also needs a subject identifier.
#
# This is particularly important when multiple samples could belong
# to the same biological subject.

cd$subjectID <- cd$subject_id


# Put this metadata into each experiment.

colData(mae[["microbiome"]]) <- cd
colData(mae[["metabolome"]]) <- cd


###############################################################################
# 5. CREATE TRAINING AND VALIDATION DATA
###############################################################################

# We do NOT want to evaluate a predictive model only on the same data
# that were used to fit it.
#
# If we did that, performance could look artificially good.
#
# Instead, we hold out approximately 20% of samples as a validation set.
#
# Training set:
#       used to develop the model
#
# Validation set:

# Get metadata from one experiment
cd <- colData(mae[["microbiome"]])

# Check that our variables exist
table(cd$Y)
head(cd$subjectID)

# Get unique subjects
subjects <- unique(cd$subjectID)

length(subjects)

# Randomly assign ~20% of SUBJECTS to validation
set.seed(377)

valid_subjects <- sample(
  subjects,
  size = round(0.2 * length(subjects))
)

# Identify samples belonging to those subjects
in_valid <- cd$subjectID %in% valid_subjects

table(in_valid)

# Get sample IDs
train_samples <- rownames(cd)[!in_valid]
valid_samples <- rownames(cd)[in_valid]

length(train_samples)
length(valid_samples)

# Subset the MAE
mae_train <- mae[, train_samples]
mae_valid <- mae[, valid_samples]


# Number of samples in each omics layer
dim(mae_train[["microbiome"]])
dim(mae_train[["metabolome"]])

dim(mae_valid[["microbiome"]])
dim(mae_valid[["metabolome"]])


table(colData(mae_train[["microbiome"]])$Y)
table(colData(mae_valid[["microbiome"]])$Y)
###############################################################################
# 6. RUN INTEGRATEDLEARNER
###############################################################################

# Now we reach the main part of the tutorial.
#
# IntegratedLearner will use BOTH omics layers to predict IBD.
#
# We tell it:
#
#   WHICH experiments to use
#
#       microbiome
#       metabolome
#
#   WHICH transformed assay from each experiment to use
#
#       clr
#       log10
#
#   WHICH variable is our outcome
#
#       Y
#
#   WHICH variable identifies subjects
#
#       subjectID


fit <- IntegratedLearner(
  
  MAE_train = mae_train,
  
  # Completely separate validation dataset
  MAE_valid = mae_valid,
  
  
  # Which omics layers should be used?
  experiment = c(
    "microbiome",
    "metabolome"
  ),
  
  
  # Which assay within each layer should be modelled?
  #
  # The order corresponds to 'experiment':
  #
  # microbiome -> clr
  # metabolome -> log10
  
  assay.type = c(
    "clr",
    "log10"
  ),
  
  
  # Binary outcome:
  #
  # 0 = healthy
  # 1 = IBD
  
  outcome_col = "Y",
  
  
  # Identifies biological subjects
  
  subject_id_col = "subjectID",
  
  
  # Five-fold cross-validation within the training data
  
  folds = 5,
  
  
  # Machine-learning algorithm used separately
  # on each omics layer.
  
  base_learner = "SL.randomForest",
  
  
  # Model used to combine predictions from the different layers.
  
  meta_learner = "sl_nnls_auc",
  
  
  # TRUE = perform LATE FUSION
  
  run_stacked = TRUE,
  
  
  # TRUE = perform EARLY FUSION
  
  run_concat = TRUE,
  
  
  # We have a binary classification problem.
  
  family = binomial()
)


###############################################################################
# WHAT IS INTEGRATEDLEARNER ACTUALLY DOING?
###############################################################################

# This is the key conceptual part of the session.
#
#
# We have:
#
#             MICROBIOME             METABOLOME
#                 |                       |
#                 |                       |
#          Random Forest            Random Forest
#                 |                       |
#                 |                       |
#           probability              probability
#              of IBD                   of IBD
#                 |                       |
#                 ----------- ------------
#                            |
#                            v
#                      META-LEARNER
#                            |
#                            v
#                    FINAL PREDICTION
#
#
# This is LATE FUSION.
#
# Each data layer first learns independently.
# We then combine the predictions from those models.


###############################################################################
# EARLY FUSION vs LATE FUSION
###############################################################################

# IntegratedLearner lets us compare two approaches.


# ---- EARLY FUSION ------------------------------------------------------------
#
# Combine the FEATURES first:
#
#
#   microbial features + metabolite features
#                   |
#                   v
#             one big matrix
#                   |
#                   v
#              one model
#
#
# This is what:
#
#   run_concat = TRUE
#
# requests.


# ---- LATE FUSION -------------------------------------------------------------
#
# Build one model per omics layer:
#
#
#   Microbiome ---> Model -----\
#                               ---> Meta-model ---> prediction
#   Metabolome ---> Model -----/
#
#
# This is what:
#
#   run_stacked = TRUE
#
# requests.
#
# Late fusion lets us ask an interesting biological question:
#
#       How much predictive information is coming
#       from each omics layer?


###############################################################################
# 7. LOOK AT THE RESULTS
###############################################################################

# First look at the learned late-fusion weights.

fit$weights


# These tell us how strongly each omics layer contributes to
# the stacked prediction.
#
# Conceptually, we can ask:
#
#   Is microbiome carrying more predictive signal?
#
#              OR
#
#   Is metabolomics carrying more predictive signal?
#
#              OR
#
#   Do both contribute?


###############################################################################
# 8. COMPARE MODEL PERFORMANCE
###############################################################################

fit$AUC.train
fit$AUC.test


# AUC = Area Under the Receiver Operating Characteristic Curve.
#
# It measures how well the model can rank:
#
#        IBD
#
# versus
#
#        healthy
#
# samples.
#
#
# Rough intuition:
#
#   AUC = 0.50  -> approximately random ranking
#
#   AUC > 0.50  -> predictive information is present
#
#   AUC = 1.00  -> perfect discrimination
#
#
# More important than memorising a particular threshold is comparing:
#
#   microbiome alone
#
#   metabolome alone
#
#   early fusion
#
#   late fusion
#
#
# The key scientific question is:
#
#   "Does combining omics layers improve prediction compared
#    with using either layer alone?"


###############################################################################
# 9. PLOT ROC CURVES
###############################################################################

IntegratedLearner:::plot.learner(fit)


# The ROC plot allows us to visually compare the predictive models.
#
# We can compare:
#
#   microbiome model
#   metabolome model
#   early-fusion model
#   late-fusion model
#
# If integration adds useful information, we hope to see the integrated
# model outperforming at least some of the individual layers.


###############################################################################
# 10. FEATURE IMPORTANCE
###############################################################################

# Prediction performance tells us:
#
#     "Can the model predict IBD?"
#
# But we often also want to ask:
#
#     "Which biological features are helping the model?"
#
# Here we extract the fitted model from each omics layer.


models <- fit$model_fits$model_layers


# Each layer has feature-importance values.
#
# But remember: the late-fusion model also learned how much importance
# to assign to EACH OMICS LAYER.
#
# Therefore, we scale each layer's feature importance by its
# corresponding IntegratedLearner weight.


imp <- do.call(
  rbind,
  lapply(
    seq_along(models),
    function(i) {
      models[[i]]$importance * fit$weights[[i]]
    }
  )
)


###############################################################################
# 11. VISUALISE THE MOST IMPORTANT FEATURES
###############################################################################

plotLoadings(
  imp,
  ncomponents = 1,
  n = 20,
  show.color = FALSE
)


# This gives us a ranked view of features contributing most strongly
# to the integrated prediction.
#
# Now our interpretation moves from:
#
#       "Can we predict IBD?"
#
# to:
#
#       "Which microbial and metabolic features may be contributing
#        to that prediction?"


###############################################################################
#
# TAKE-HOME MESSAGE
#
###############################################################################

# The complete workflow was:
#
#
#         MultiAssayExperiment
#                 |
#                 v
#          Filter rare features
#                 |
#                 v
#       Remove invariant features
#                 |
#                 v
#            Transform data
#           /              \
#         CLR              log10
#     microbiome         metabolome
#           \              /
#            \            /
#             v          v
#          IntegratedLearner
#             /       \
#            /         \
#       Early fusion   Late fusion
#            \         /
#             \       /
#             Compare AUC
#                 |
#                 v
#           Layer weights
#                 |
#                 v
#         Feature importance
#
#
# Main idea:
#
# IntegratedLearner does not simply ask:
#
#     "Can microbiome predict IBD?"
#
# or
#
#     "Can metabolomics predict IBD?"
#
# It allows us to ask:
#
#     "Can these different biological layers work together
#      to improve prediction, and what does each layer contribute?"
#
###############################################################################

