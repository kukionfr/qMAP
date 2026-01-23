# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

qMAP (Quantitative Microanatomical Mapping) is a research codebase for analyzing human skin aging through morphometric features. The project performs statistical analysis correlating tissue features with patient age, including univariate/multivariate regression models and gender classification.

## Development Environment

- **MATLAB 2024a** (primary language)
- **Python 3.8.10** (supporting libraries)
- Tested on Windows 10

### Setup

```bash
# Python environment (conda)
conda env create -f environment.yml
conda activate qmap

# Or pip
pip install -r requirements.txt
```

## Running the Analysis

Open MATLAB, navigate to `source_code/`, and execute:

```matlab
hopkins_cohort_analysis.m
```

This generates the manuscript figures and outputs correlation metrics to Excel files.

## Code Architecture

### Main Script: `source_code/hopkins_cohort_analysis.m`

The analysis pipeline has sequential parts:

1. **Part 1 (Lines 1-47)**: Loads morphometric data from `input_data/data.xlsx`, aggregates section-level data to patient-level metrics (mean, median, std)
2. **Part 2 (Lines 48-72)**: Cohort selection - filters by body part (back), race (white), extracts demographic variables (gender, age)
3. **Figure 2C (Lines 79-182)**: Computes Spearman correlations and Cohen's d between features and age, identifies age-associated features (|ρ|>0.3, d>1)
4. **Figure 2F (Lines 183-270)**: Hierarchical clustering of age-associated features, Sankey diagram visualization
5. **Figure 3 (Lines 271-465)**: Univariate and bivariate GLM age prediction with leave-one-out cross-validation
6. **Figure 3H-I (Lines 516-685)**: PCA dimensionality reduction, multivariate regression (GLM, Lasso, SVM, Kernel) comparison
7. **Figure 4 (Lines 686-1000+)**: Gender analysis - sex differences in aging features, gender classification using discriminant analysis

### Helper Functions

- `cell2num.m`: Converts cell arrays with numeric content to double arrays
- `cohend.m`: Computes Cohen's d effect size between two groups
- `bjff3.p`: Figure formatting utility (compiled)

### Data Files

- `input_data/data.xlsx`: Patient morphometric measurements (sheet: "by section all")
- `input_data/feature_list.xlsx`: Feature metadata including "nonsense" flag for filtering
- `expected_output/plmed.mat`: Patient-level median feature values
- `expected_output/*.mat`: Cached demographic filters (ismale, iswhite, isback, age)

## Key Variables

- `plmed`: Patient-level median of morphometric features (patients × features matrix)
- `ccsample`: Logical index for selected cohort (white + back body part)
- `ccf`/`age_associated_features`: Logical index for age-correlated features
- `ccm`/`ismale`: Logical index for male patients
