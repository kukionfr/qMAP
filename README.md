# qMAP: Quantitative Microanatomical Mapping of Human Skin Aging

Code repository for *"Quantitative microanatomical analysis of human skin aging"*.

## Overview

qMAP is a computational pipeline for quantifying age-related morphological changes in human skin from H&E-stained histological sections. The pipeline combines deep learning-based tissue and nuclei segmentation with statistical analysis to identify and characterize tissue-level aging features.

**Key capabilities:**
- Automated tissue compartment labeling (12 classes) using DeepLabV3+
- Nuclei instance segmentation using a semi-supervised HoVerNet retraining pipeline
- Extraction of 109 morphometric features across tissue compartments
- Univariate and multivariate age prediction (MAE = 8.7 years with SVM)
- Sex-specific aging analysis and gender classification

## Repository Structure

```
qMAP/
├── source_code/
│   ├── hopkins_cohort_analysis.m    # Main analysis script (MATLAB)
│   ├── cell2num.m                   # Cell-to-numeric array conversion
│   ├── cohend.m                     # Cohen's d effect size computation
│   ├── bjff3.p                      # Figure formatting utility
│   ├── benchmark_hovernet.py        # HoVerNet benchmarking on standard datasets
│   ├── benchmark_stardist.py        # StarDist benchmarking on standard datasets
│   ├── benchmark_cellvit.py         # CellViT benchmarking on standard datasets
│   ├── benchmark_unified.py         # Unified benchmark runner and report generator
│   ├── evaluate_metrics.py          # Dice, AJI, PQ metric computation
│   ├── download_datasets.py         # Download NuInsSeg, MoNuSeg, CryoNuSeg
│   ├── prepare_nuinsseg.py          # NuInsSeg dataset preparation
│   ├── prepare_monuseg.py           # MoNuSeg dataset preparation
│   ├── prepare_cryonuseg.py         # CryoNuSeg dataset preparation
│   └── requirements_benchmark.txt   # Python dependencies for benchmarking
├── input_data/
│   ├── data.xlsx                    # Patient morphometric measurements
│   └── feature_list.xlsx            # Feature metadata
├── expected_output/                 # Cached intermediate results
├── environment.yml                  # Conda environment specification
└── requirements.txt                 # Python dependencies
```

## Requirements

- **MATLAB 2024a** (for main morphometric analysis)
- **Python 3.8+** (for benchmarking scripts)

### Setup

```bash
# Option 1: Conda (recommended)
conda env create -f environment.yml
conda activate qmap

# Option 2: pip
pip install -r requirements.txt

# For nuclei segmentation benchmarking
pip install -r source_code/requirements_benchmark.txt
```

## Usage

### Morphometric Analysis (Main Pipeline)

Open MATLAB, navigate to `source_code/`, and run:

```matlab
hopkins_cohort_analysis.m
```

This generates all manuscript figures and outputs correlation metrics. Expected run time: ~10 minutes.

### Nuclei Segmentation Benchmarking

To reproduce the cross-dataset benchmarking of pretrained nuclei segmentation models:

```bash
# Download benchmark datasets
python source_code/download_datasets.py

# Prepare datasets
python source_code/prepare_nuinsseg.py
python source_code/prepare_monuseg.py
python source_code/prepare_cryonuseg.py

# Run all benchmarks (HoVerNet, StarDist, CellViT on NuInsSeg, MoNuSeg, CryoNuSeg)
python source_code/benchmark_unified.py --all

# Or run individual models
python source_code/benchmark_hovernet.py
python source_code/benchmark_stardist.py
python source_code/benchmark_cellvit.py
```

Results are saved to `output/benchmark/` with a unified report at `output/benchmark/benchmark_report.md`.

## Cross-Dataset Benchmark Results

All models evaluated with pretrained weights only (no fine-tuning on evaluation data):

| Model | Training Data | NuInsSeg (n=665) | MoNuSeg (n=32) | CryoNuSeg (n=30) |
|-------|---------------|------------------|----------------|-------------------|
| HoVerNet | PanNuke | 0.497 / 0.313 / 0.275 | 0.790 / 0.442 / 0.410 | 0.777 / 0.524 / 0.424 |
| StarDist | MoNuSeg+TNBC | 0.452 / 0.284 / 0.278 | 0.753 / 0.425 / 0.411 | 0.733 / 0.501 / 0.416 |
| CellViT-SAM-H | PanNuke | 0.659 / 0.458 / 0.403 | 0.803 / 0.448 / 0.426 | 0.792 / 0.546 / 0.453 |

*Values: Dice / AJI / PQ. MoNuSeg restricted to 32 TCGA images (1000x1000).*

## Tested Environment

- Windows 10, MATLAB 2024a, Python 3.8.10
- Ubuntu 22.04, Python 3.8+ (benchmarking)

## License

See [LICENSE](LICENSE) file.
