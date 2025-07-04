# Multi-Omics Clustering for Cancer Subtyping (MCLS)

This repository contains the complete implementation of the project:  
**"Integrating Multi-Omics Data for Cancer Disease Subtyping Using Machine Learning Approaches"**  
as part of the Project II course (BCSE498J) at VIT Vellore.

## 🧠 Project Overview

This project explores a method called **MCLS** (Multi-omics Clustering using Latent Subspace), which integrates multi-omics data (e.g., gene expression, methylation, miRNA) to cluster cancer patient samples into biologically meaningful subtypes. The key features of this approach include:

- Latent subspace learning using PCA and SVD
- Robust to **missing/incomplete omics views**
- Spectral clustering on shared latent space
- Evaluation using **survival analysis** and **clinical enrichment**

## 📁 Folder Structure

```plaintext
├── data/                 # Contains preprocessed multi-omics datasets (e.g., BIC, BREAST)
├── Complete_data.R       # Script for clustering on complete multi-omics data
├── Incomplete_data.R     # Script for clustering with missing modalities
├── run_data_col.R        # Subspace learning functions
├── specClust.R           # Spectral clustering implementation
├── enrichment.R          # Clinical enrichment analysis (age, gender, stage, TNM)
├── README.md             # Project documentation (this file)
```

## 🧪 Features

- Handles both **complete** and **incomplete** multi-omics datasets
- Avoids imputation by learning a projection from complete samples
- Uses **spectral clustering** with automatic selection of optimal cluster count
- Evaluates clusters using:
  - Kaplan-Meier survival curves
  - Clinical parameter enrichment (age, gender, stage)

## 🛠 Requirements

- R (version ≥ 4.0)
- Packages:
  - `survival`, `MASS`, `igraph`, `RSpectra`, `kknn`, `ggplot2`

## 🚀 How to Run

For complete data:
```r
source("Complete_data.R")
```

For incomplete data:
```r
source("Incomplete_data.R")
```
Example Dataset
The default scripts use the BIC-BREAST dataset from TCGA.
All realworld datasets are saved in "data/..." including the following seven multi-omics datasets:
BIC
COAD
GBM
KRCCC
LSCC
Bladder
Brain

## 📈 Results
Clusters are assigned to each sample

Survival p-values indicate if clusters have significant prognostic value

Clinical enrichment tells how well clusters correlate with metadata (age, stage, etc.)
