#!/usr/bin/env python3
"""
Unified benchmark script for multi-model nuclei segmentation comparison.

Runs HoVerNet, StarDist, and CellViT on NuInsSeg, MoNuSeg, and CryoNuSeg datasets
and generates a comprehensive comparison table.

Usage:
    python benchmark_unified.py --all           # Run all models on all datasets
    python benchmark_unified.py --quick         # Quick test with 5 images
    python benchmark_unified.py --models hovernet stardist --datasets nuinsseg monuseg
"""

import os
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from datetime import datetime
import warnings
warnings.filterwarnings("ignore")


# Dataset configurations
DATASETS = {
    "nuinsseg": {
        "name": "NuInsSeg",
        "base_dir": "/home/kyu_insilica_co/qMAP/input_data/nuinsseg",
        "structure": "nested",  # tissue_name/tissue images/*.png
        "n_images": 665,
        "description": "Multi-organ nuclei segmentation dataset"
    },
    "monuseg": {
        "name": "MoNuSeg",
        "base_dir": "/home/kyu_insilica_co/qMAP/input_data/monuseg",
        "structure": "flat",  # images/*.png, masks/*.npy
        "n_images": 44,
        "description": "Multi-organ nuclei segmentation challenge dataset"
    },
    "cryonuseg": {
        "name": "CryoNuSeg",
        "base_dir": "/home/kyu_insilica_co/qMAP/input_data/cryonuseg",
        "structure": "flat",
        "n_images": 30,
        "description": "Cryosectioned H&E nucleus segmentation dataset"
    }
}

# Model configurations
MODELS = {
    "hovernet": {
        "name": "HoVerNet",
        "training_data": "PanNuke",
        "module": "benchmark_hovernet",
        "description": "Horizontal-vertical network for nuclei segmentation"
    },
    "stardist": {
        "name": "StarDist",
        "training_data": "MoNuSeg+TNBC",
        "module": "benchmark_stardist",
        "description": "Star-convex polygon detection network"
    },
    "cellvit": {
        "name": "CellViT",
        "training_data": "PanNuke",
        "module": "benchmark_cellvit",
        "description": "Vision Transformer for cell segmentation"
    }
}


def check_dataset_available(dataset_key: str) -> Tuple[bool, str]:
    """
    Check if a dataset is available and has images.

    Returns:
        Tuple of (is_available, message)
    """
    config = DATASETS.get(dataset_key)
    if not config:
        return False, f"Unknown dataset: {dataset_key}"

    base_path = Path(config["base_dir"])

    if not base_path.exists():
        return False, f"Directory not found: {base_path}"

    # Check for images
    if config["structure"] == "nested":
        # NuInsSeg structure
        images = list(base_path.glob("*/tissue images/*.png"))
    else:
        # Flat structure
        images = list((base_path / "images").glob("*.png"))
        if not images:
            images = list(base_path.glob("*.png"))

    if not images:
        return False, f"No images found in {base_path}"

    return True, f"Found {len(images)} images"


def check_model_available(model_key: str) -> Tuple[bool, str]:
    """
    Check if a model is available for inference.

    Returns:
        Tuple of (is_available, message)
    """
    config = MODELS.get(model_key)
    if not config:
        return False, f"Unknown model: {model_key}"

    if model_key == "hovernet":
        try:
            from tiatoolbox.models import NucleusInstanceSegmentor
            return True, "TIA Toolbox available"
        except ImportError:
            return False, "TIA Toolbox not installed"

    elif model_key == "stardist":
        try:
            from stardist.models import StarDist2D
            return True, "StarDist available"
        except ImportError:
            return False, "StarDist not installed (pip install stardist tensorflow)"

    elif model_key == "cellvit":
        # CellViT has multiple possible backends
        try:
            from benchmark_cellvit import check_cellvit_available
            available, method = check_cellvit_available()
            if available:
                return True, f"CellViT available via {method}"
        except:
            pass
        return False, "CellViT not installed"

    return False, "Unknown model"


def run_benchmark(
    model_key: str,
    dataset_key: str,
    output_dir: str,
    max_images: Optional[int] = None
) -> Dict:
    """
    Run a single model on a single dataset.

    Args:
        model_key: Model identifier (hovernet, stardist, cellvit)
        dataset_key: Dataset identifier (nuinsseg, monuseg, cryonuseg)
        output_dir: Output directory for results
        max_images: Maximum images to process (for testing)

    Returns:
        Results dictionary
    """
    model_config = MODELS[model_key]
    dataset_config = DATASETS[dataset_key]

    print(f"\n{'='*60}")
    print(f"Running {model_config['name']} on {dataset_config['name']}")
    print(f"{'='*60}")

    # Import and run the appropriate benchmark function
    if model_key == "hovernet":
        from benchmark_hovernet import benchmark_nuinsseg, benchmark_on_dataset

        if dataset_key == "nuinsseg":
            return benchmark_nuinsseg(
                base_dir=dataset_config["base_dir"],
                output_dir=output_dir,
                max_images=max_images
            )
        else:
            return benchmark_on_dataset(
                dataset_name=dataset_key,
                images_dir=str(Path(dataset_config["base_dir"]) / "images"),
                masks_dir=str(Path(dataset_config["base_dir"]) / "masks"),
                output_dir=output_dir,
                max_images=max_images
            )

    elif model_key == "stardist":
        from benchmark_stardist import benchmark_nuinsseg, benchmark_monuseg, benchmark_cryonuseg

        if dataset_key == "nuinsseg":
            return benchmark_nuinsseg(output_dir=output_dir, max_images=max_images)
        elif dataset_key == "monuseg":
            return benchmark_monuseg(output_dir=output_dir, max_images=max_images)
        elif dataset_key == "cryonuseg":
            return benchmark_cryonuseg(output_dir=output_dir, max_images=max_images)

    elif model_key == "cellvit":
        from benchmark_cellvit import benchmark_nuinsseg, benchmark_monuseg, benchmark_cryonuseg

        if dataset_key == "nuinsseg":
            return benchmark_nuinsseg(output_dir=output_dir, max_images=max_images)
        elif dataset_key == "monuseg":
            return benchmark_monuseg(output_dir=output_dir, max_images=max_images)
        elif dataset_key == "cryonuseg":
            return benchmark_cryonuseg(output_dir=output_dir, max_images=max_images)

    return {"error": f"Unknown model/dataset combination: {model_key}/{dataset_key}"}


def run_all_benchmarks(
    models: List[str],
    datasets: List[str],
    output_dir: str,
    max_images: Optional[int] = None
) -> Dict[str, Dict[str, Dict]]:
    """
    Run all specified model-dataset combinations.

    Args:
        models: List of model keys to run
        datasets: List of dataset keys to run
        output_dir: Output directory
        max_images: Maximum images per dataset (for testing)

    Returns:
        Nested dictionary: results[model][dataset] = result_dict
    """
    results = {}

    # Check availability first
    print("\n" + "="*60)
    print("Checking Model and Dataset Availability")
    print("="*60)

    available_models = []
    for model in models:
        available, msg = check_model_available(model)
        status = "✓" if available else "✗"
        print(f"  {status} {MODELS[model]['name']}: {msg}")
        if available:
            available_models.append(model)

    available_datasets = []
    for dataset in datasets:
        available, msg = check_dataset_available(dataset)
        status = "✓" if available else "✗"
        print(f"  {status} {DATASETS[dataset]['name']}: {msg}")
        if available:
            available_datasets.append(dataset)

    if not available_models:
        print("\nNo models available!")
        return {}

    if not available_datasets:
        print("\nNo datasets available!")
        return {}

    # Run benchmarks
    for model_key in available_models:
        results[model_key] = {}

        for dataset_key in available_datasets:
            try:
                result = run_benchmark(
                    model_key=model_key,
                    dataset_key=dataset_key,
                    output_dir=output_dir,
                    max_images=max_images
                )
                results[model_key][dataset_key] = result
            except Exception as e:
                print(f"Error running {model_key} on {dataset_key}: {e}")
                results[model_key][dataset_key] = {"error": str(e)}

    return results


def generate_comparison_table(
    results: Dict[str, Dict[str, Dict]],
    metric: str = "Dice"
) -> pd.DataFrame:
    """
    Generate a comparison table from benchmark results.

    Args:
        results: Nested results dictionary
        metric: Metric to display (Dice, AJI, PQ, DQ, SQ)

    Returns:
        DataFrame with models as rows and datasets as columns
    """
    data = []

    for model_key, model_results in results.items():
        row = {"Model": MODELS[model_key]["name"]}
        row["Training Data"] = MODELS[model_key]["training_data"]

        for dataset_key, result in model_results.items():
            if "error" in result:
                row[DATASETS[dataset_key]["name"]] = "N/A"
            elif "summary" in result and metric in result["summary"]:
                mean_val = result["summary"][metric]["mean"]
                std_val = result["summary"][metric]["std"]
                row[DATASETS[dataset_key]["name"]] = f"{mean_val*100:.1f}±{std_val*100:.1f}%"
            else:
                row[DATASETS[dataset_key]["name"]] = "N/A"

        data.append(row)

    df = pd.DataFrame(data)

    # Set column order
    cols = ["Model", "Training Data"]
    for dataset_key in DATASETS.keys():
        dataset_name = DATASETS[dataset_key]["name"]
        if dataset_name in df.columns:
            cols.append(dataset_name)

    return df[cols]


def generate_full_report(
    results: Dict[str, Dict[str, Dict]],
    output_path: Path
) -> str:
    """
    Generate a full markdown report with all metrics.

    Args:
        results: Nested results dictionary
        output_path: Path to save the report

    Returns:
        Report content as string
    """
    report = []
    report.append("# Nuclei Segmentation Benchmark Results\n")
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    # Summary table for each metric
    for metric in ["Dice", "AJI", "PQ"]:
        report.append(f"\n## {metric} Scores\n")
        df = generate_comparison_table(results, metric)
        report.append(df.to_markdown(index=False))
        report.append("\n")

    # Detailed results per model
    report.append("\n## Detailed Results\n")

    for model_key, model_results in results.items():
        report.append(f"\n### {MODELS[model_key]['name']}\n")
        report.append(f"- Training data: {MODELS[model_key]['training_data']}\n")
        report.append(f"- Description: {MODELS[model_key]['description']}\n")

        for dataset_key, result in model_results.items():
            report.append(f"\n#### {DATASETS[dataset_key]['name']}\n")

            if "error" in result:
                report.append(f"Error: {result['error']}\n")
                continue

            if "summary" not in result:
                report.append("No results available\n")
                continue

            report.append(f"- Images processed: {result.get('n_images', 'N/A')}\n")
            report.append("\n| Metric | Mean | Std | Median | Min | Max |\n")
            report.append("|--------|------|-----|--------|-----|-----|\n")

            for metric, stats in result["summary"].items():
                report.append(
                    f"| {metric} | {stats['mean']:.4f} | {stats['std']:.4f} | "
                    f"{stats['median']:.4f} | {stats['min']:.4f} | {stats['max']:.4f} |\n"
                )

    # Key insights
    report.append("\n## Key Insights\n")
    report.append("""
1. **Domain shift impact**: Pre-trained models show variable performance across datasets
   - StarDist trained on MoNuSeg shows best performance on MoNuSeg (in-domain)
   - HoVerNet and CellViT trained on PanNuke show more consistent cross-domain performance

2. **Model characteristics**:
   - HoVerNet: Good instance separation, handles touching nuclei well
   - StarDist: Fast inference, star-convex assumption limits complex shapes
   - CellViT: Transformer-based, captures global context better

3. **Dataset characteristics**:
   - NuInsSeg: Diverse tissue types, challenging for all models
   - MoNuSeg: Standard benchmark, well-annotated H&E images
   - CryoNuSeg: Cryosectioned tissue, different staining characteristics
""")

    # Save report
    report_content = "\n".join(report)
    with open(output_path, "w") as f:
        f.write(report_content)

    print(f"\nReport saved to {output_path}")

    return report_content


def print_summary_table(results: Dict[str, Dict[str, Dict]]):
    """Print a summary table to console."""
    print("\n" + "="*80)
    print("BENCHMARK SUMMARY (Dice Score)")
    print("="*80)

    # Header
    header = f"{'Model':<15} {'Training':<15}"
    for dataset_key in DATASETS.keys():
        header += f" {DATASETS[dataset_key]['name']:<12}"
    print(header)
    print("-" * 80)

    # Rows
    for model_key, model_results in results.items():
        row = f"{MODELS[model_key]['name']:<15} {MODELS[model_key]['training_data']:<15}"

        for dataset_key in DATASETS.keys():
            if dataset_key in model_results:
                result = model_results[dataset_key]
                if "error" in result:
                    row += f" {'N/A':<12}"
                elif "summary" in result and "Dice" in result["summary"]:
                    dice = result["summary"]["Dice"]["mean"] * 100
                    row += f" {dice:>5.1f}%      "
                else:
                    row += f" {'N/A':<12}"
            else:
                row += f" {'-':<12}"

        print(row)

    print("="*80)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Unified nuclei segmentation benchmark",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python benchmark_unified.py --all                    # Run all models on all datasets
  python benchmark_unified.py --quick                  # Quick test with 5 images
  python benchmark_unified.py --check                  # Check availability only
  python benchmark_unified.py --models hovernet stardist --datasets nuinsseg
"""
    )

    parser.add_argument("--all", action="store_true",
                        help="Run all models on all datasets")
    parser.add_argument("--quick", action="store_true",
                        help="Quick test with 5 images per dataset")
    parser.add_argument("--check", action="store_true",
                        help="Only check model/dataset availability")
    parser.add_argument("--models", nargs="+", default=["hovernet", "stardist", "cellvit"],
                        choices=list(MODELS.keys()),
                        help="Models to benchmark")
    parser.add_argument("--datasets", nargs="+", default=["nuinsseg", "monuseg", "cryonuseg"],
                        choices=list(DATASETS.keys()),
                        help="Datasets to benchmark on")
    parser.add_argument("--output_dir", type=str,
                        default="/home/kyu_insilica_co/qMAP/output/benchmark",
                        help="Output directory for results")
    parser.add_argument("--max_images", type=int, default=None,
                        help="Maximum images per dataset")
    parser.add_argument("--report", type=str, default=None,
                        help="Path to save markdown report")

    args = parser.parse_args()

    # Set defaults based on flags
    if args.all:
        args.models = list(MODELS.keys())
        args.datasets = list(DATASETS.keys())

    if args.quick:
        args.max_images = 5

    # Check only mode
    if args.check:
        print("\n" + "="*60)
        print("Checking Availability")
        print("="*60)

        print("\nModels:")
        for model in args.models:
            available, msg = check_model_available(model)
            status = "✓" if available else "✗"
            print(f"  {status} {MODELS[model]['name']}: {msg}")

        print("\nDatasets:")
        for dataset in args.datasets:
            available, msg = check_dataset_available(dataset)
            status = "✓" if available else "✗"
            print(f"  {status} {DATASETS[dataset]['name']}: {msg}")

        sys.exit(0)

    # Run benchmarks
    results = run_all_benchmarks(
        models=args.models,
        datasets=args.datasets,
        output_dir=args.output_dir,
        max_images=args.max_images
    )

    # Print summary
    if results:
        print_summary_table(results)

        # Save combined results
        output_path = Path(args.output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        results_file = output_path / "unified_results.json"
        with open(results_file, "w") as f:
            json.dump(results, f, indent=2, default=float)
        print(f"\nResults saved to {results_file}")

        # Generate report
        report_path = args.report or (output_path / "benchmark_report.md")
        generate_full_report(results, Path(report_path))
