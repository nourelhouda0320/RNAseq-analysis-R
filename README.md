# RNA-seq Analysis in R

A hands-on RNA-seq data analysis project using R and Bioconductor, based on a mouse mammary gland RNA-seq dataset.

## Overview

This project implements an RNA-seq analysis workflow starting from a raw gene-level count matrix and sample metadata.

The analysis focuses on data preprocessing, quality control, expression filtering, exploratory analysis, identification of highly variable genes, and data visualization.

The project is based on publicly available RNA-seq training materials and has been reorganized into a reproducible R workflow.

## Objectives

The main objectives of this project are to:

- Import and preprocess RNA-seq count data
- Perform quality control of sequencing libraries
- Calculate counts per million (CPM)
- Filter lowly expressed genes
- Explore sample relationships using multidimensional scaling (MDS)
- Identify highly variable genes
- Visualize gene expression patterns using heatmaps
- Prepare the data for downstream differential expression analysis

## Dataset

The analysis uses the **GSE60450 mouse mammary gland RNA-seq dataset**.

The dataset contains gene-level RNA-seq count data together with sample metadata describing biological characteristics such as:

- Cell type
- Biological status
- Sample identity

The count matrix and corrected sample metadata are provided in the `data/` directory.

## Analysis Workflow

```text
RNA-seq Count Matrix
        │
        ▼
Data Import
        │
        ▼
Data Preprocessing
        │
        ▼
CPM Calculation
        │
        ▼
Low-expression Filtering
        │
        ▼
Library Size QC
        │
        ▼
Log2-CPM Transformation
        │
        ▼
MDS Analysis
        │
        ▼
Highly Variable Genes
        │
        ▼
Heatmap Visualization
```

## Methods

### 1. Data Preprocessing

The RNA-seq count matrix was imported into R and prepared for downstream analysis.

Gene identifiers were used as row names, and sample names were standardized to ensure consistency between the count matrix and the sample metadata.

### 2. Expression Filtering

Counts per million (CPM) were calculated using the `edgeR` package.

Genes with CPM greater than 0.5 in at least two samples were retained for downstream analysis. This step removes genes with very low expression across the dataset.

### 3. Quality Control

Quality control and exploratory analysis were performed using:

- Library size visualization
- Log2-CPM distribution
- Multidimensional scaling (MDS)

These analyses were used to assess sample quality and identify potential differences between biological conditions.

### 4. MDS Analysis

Multidimensional scaling was used to visualize relationships between RNA-seq samples.

Samples were colored according to:

- Cell type
- Biological status

This allows the main sources of variation between samples to be explored.

### 5. Highly Variable Genes

Gene-wise variance was calculated using the log2-CPM expression matrix.

The 500 genes with the highest variance across samples were selected for further visualization.

### 6. Heatmap Visualization

A heatmap was generated to visualize the expression patterns of the 500 most variable genes across the samples.

The heatmap provides an overview of similarities and differences in gene expression profiles between samples.

## Tools and Packages

- R
- Bioconductor
- edgeR
- limma
- Glimma
- gplots
- RColorBrewer
- org.Mm.eg.db

## Repository Structure

```text
RNAseq-analysis-R/
│
├── README.md
├── analysis.R
│
├── data/
│   ├── GSE60450_Lactation-GenewiseCounts.txt
│   └── SampleInfo_Corrected.txt
│
└── results/
    └── figures/


```
## Acknowledgements

This project is based on publicly available RNA-seq training materials from the Bioinformatics Core Shared Training.

Original training materials:

https://bioinformatics-core-shared-training.github.io/RNAseq-R/

The original training materials were developed for RNA-seq analysis using R and Bioconductor.
