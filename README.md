# Bioconductor Africa Multi-Omics Integration Workshop

This repository contains the practical material for a **Bioconductor Africa Seminar Series** session on multi-omics data integration using the Bioconductor ecosystem.

## Background


Researchers are increasingly interested in understanding how different data types — including transcriptomics, genomics, metabolomics, microbiome, and single-cell data — can be analysed together to answer a biological question, rather than being treated as completely separate analyses.

This workshop provides a practical demonstration of one such workflow using **microbiome and metabolomics data**. The session focuses on how multiple omics layers can be organised using `MultiAssayExperiment`, preprocessed appropriately, and integrated for prediction using `IntegratedLearner`.

## Learning objectives

By the end of the session, participants should be able to:

- understand why `MultiAssayExperiment` is useful for multi-omics studies;
- organise multiple omics assays and shared sample metadata;
- match and align samples across data types;
- apply assay-specific filtering and transformations;
- prepare repeated-measures data for predictive modelling;
- split data at the subject level to avoid information leakage;
- run `IntegratedLearner` on a `MultiAssayExperiment`;
- distinguish between **early fusion** and **late fusion** approaches;
- compare predictive performance using AUC and ROC curves; and
- inspect omics-layer weights and feature importance.

## Biological use case

The example uses paired:

- **microbiome data**, and
- **metabolomics data**

from an inflammatory bowel disease dataset (iHMP).

The prediction task is to distinguish **IBD** from non-IBD samples using information from both omics layers.


## Repository structure

```text
.
├── README.md
├── Install.R
├── build_mae.R
├── IL.R
├── data/hmpibd_tse.rds
├── data/hmpibd_met.rds
└── data/hmpibd_mae.rds
```

### `Install.R`

Installs and loads the R/Bioconductor packages required for the workshop.

Run this first when setting up the workshop environment.

### `build_mae.R`

Builds a `MultiAssayExperiment` from the microbiome and metabolomics datasets.

This script demonstrates:

- extracting sample metadata;
- identifying samples shared between assays;
- aligning the two omics datasets;
- storing them as separate experiments within one coordinated object; and
- saving the resulting `MultiAssayExperiment`.

### `IL.R`

Contains the main `IntegratedLearner` tutorial.

The workflow includes:

- prevalence filtering;
- near-zero-variance filtering;
- CLR transformation of microbiome data;
- log transformation of metabolomics data;
- preparation of the binary outcome;
- subject-level train/validation splitting;
- fitting `IntegratedLearner`;
- comparison of early and late fusion;
- AUC and ROC evaluation; and
- feature-importance visualisation.

## Data files

The repository includes the data required for the tutorial:

### `hmpibd_tse.rds`

Microbiome data stored as a Bioconductor experiment object.

### `hmpibd_met.rds`

Metabolomics data stored as a Bioconductor experiment object.

### `hmpibd_mae.rds`

A prepared `MultiAssayExperiment` containing the aligned microbiome and metabolomics datasets.

If you want to focus only on the `IntegratedLearner` analysis, you can start directly from `hmpibd_mae.rds`.

If you want to understand how the multi-omics object is constructed, run `build_mae.R` first.

## Setup

Clone or download this repository and open it as your working directory in R or RStudio.

Run:

```r
source("Install.R")
```

This will install missing dependencies and load the packages used in the workshop.

The main packages include:

- `MultiAssayExperiment`
- `SummarizedExperiment`
- `mia`
- `matrixStats`
- `tidyverse`
- `IntegratedLearner`
- `miaViz`

## Suggested workflow

For the full workshop, run the files in this order:

```text
1. Install.R
      |
      v
2. build_mae.R
      |
      v
3. IL.R
```

If the prepared MAE is already available and you want to go directly to the machine-learning section:

```text
1. Install.R
      |
      v
2. IL.R
```

