#!/usr/bin/env python3
"""
Benchmark HoVerNet on nuclei segmentation datasets.

Uses TIA Toolbox for HoVerNet inference and computes standard metrics.
"""

import os
import sys
import json
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

# Import evaluation metrics
from evaluate_metrics import compute_all_metrics, compute_dataset_metrics


def run_hovernet_inference(
    image_paths: List[str],
    output_dir: str,
    model_name: str = "hovernet_fast-pannuke",
    batch_size: int = 8,
    num_workers: int = 4
) -> Dict[str, np.ndarray]:
    """
    Run HoVerNet inference using TIA Toolbox.

    Args:
        image_paths: List of paths to input images
        output_dir: Directory to save outputs
        model_name: Pretrained model name
        batch_size: Batch size for inference
        num_workers: Number of data loader workers

    Returns:
        Dictionary mapping image names to instance masks
    """
    from tiatoolbox.models import NucleusInstanceSegmentor

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    print(f"Running HoVerNet ({model_name}) on {len(image_paths)} images...")

    # Initialize segmentor
    segmentor = NucleusInstanceSegmentor(
        pretrained_model=model_name,
        batch_size=batch_size,
        num_postproc_workers=num_workers,
        verbose=False
    )

    # Run prediction (use CUDA if available)
    import torch
    device = "cuda" if torch.cuda.is_available() else "cpu"
    output = segmentor.predict(
        imgs=image_paths,
        save_dir=str(output_path),
        mode="tile",
        device=device,
        crash_on_exception=False
    )

    # Load predictions
    predictions = {}
    for img_path in image_paths:
        img_name = Path(img_path).stem

        # TIA Toolbox saves results as .dat files
        result_path = output_path / f"{img_name}.dat"
        if result_path.exists():
            import pickle
            with open(result_path, "rb") as f:
                result = pickle.load(f)

            # Extract instance map
            if "inst_map" in result:
                predictions[img_name] = result["inst_map"]
            else:
                # Reconstruct from nuclei dictionary
                inst_map = np.zeros_like(result.get("type_map", np.zeros((512, 512))), dtype=np.int32)
                for idx, (nuc_id, nuc_info) in enumerate(result.get("nuc", {}).items(), 1):
                    contour = nuc_info.get("contour", [])
                    if len(contour) > 0:
                        from skimage.draw import polygon
                        rr, cc = polygon(
                            [p[1] for p in contour],
                            [p[0] for p in contour],
                            inst_map.shape
                        )
                        inst_map[rr, cc] = idx
                predictions[img_name] = inst_map

    return predictions


def reconstruct_instance_map(nuclei_dict: Dict, image_shape: Tuple[int, int]) -> np.ndarray:
    """
    Reconstruct instance segmentation map from nuclei contours.

    Args:
        nuclei_dict: Dictionary with nucleus IDs as keys and info dicts as values
        image_shape: (height, width) of the output mask

    Returns:
        Instance segmentation mask where each nucleus has a unique integer ID
    """
    from skimage.draw import polygon

    inst_map = np.zeros(image_shape, dtype=np.int32)

    for idx, (nuc_id, nuc_info) in enumerate(nuclei_dict.items(), 1):
        contour = nuc_info.get("contour", [])
        if len(contour) > 0:
            # Contour format is [[x, y], ...], need to extract rows and cols
            contour = np.array(contour)
            rr, cc = polygon(contour[:, 1], contour[:, 0], image_shape)
            inst_map[rr, cc] = idx

    return inst_map


def run_hovernet_inference_simple(
    image_paths: List[str],
    output_dir: str,
    model_name: str = "hovernet_fast-pannuke",
    batch_size: int = 8
) -> Dict[str, np.ndarray]:
    """
    Run HoVerNet inference using TIA Toolbox with batch processing.

    Args:
        image_paths: List of paths to input images
        output_dir: Directory to save outputs
        model_name: Pretrained model name
        batch_size: Batch size for GPU inference

    Returns:
        Dictionary mapping image names to instance masks
    """
    from tiatoolbox.models.engine.nucleus_instance_segmentor import NucleusInstanceSegmentor
    from PIL import Image
    import torch
    import joblib

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    print(f"Running HoVerNet ({model_name}) on {len(image_paths)} images...")
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")

    # Initialize segmentor with larger batch size for A100
    segmentor = NucleusInstanceSegmentor(
        pretrained_model=model_name,
        batch_size=batch_size,
        verbose=False
    )

    predictions = {}
    failed = []

    for img_path in tqdm(image_paths, desc="HoVerNet inference"):
        img_name = Path(img_path).stem

        try:
            # Get image dimensions
            img = Image.open(img_path)
            img_shape = (img.height, img.width)

            # Run prediction on single image
            output = segmentor.predict(
                imgs=[img_path],
                save_dir=str(output_path / img_name),
                mode="tile",
                device=device,
                crash_on_exception=True
            )

            # Load result using joblib (TIA Toolbox uses joblib, not pickle)
            result_file = output_path / img_name / "0.dat"
            if result_file.exists():
                result = joblib.load(result_file)

                # TIA Toolbox returns dict of nuclei with contours
                if isinstance(result, dict):
                    # Check if it's nuclei dictionary format
                    first_key = list(result.keys())[0] if result else None
                    if first_key and isinstance(result.get(first_key), dict) and "contour" in result.get(first_key, {}):
                        # Reconstruct instance map from contours
                        inst_map = reconstruct_instance_map(result, img_shape)
                        predictions[img_name] = inst_map
                    elif "inst_map" in result:
                        predictions[img_name] = result["inst_map"]
                    else:
                        # Empty prediction
                        predictions[img_name] = np.zeros(img_shape, dtype=np.int32)
                elif isinstance(result, np.ndarray):
                    predictions[img_name] = result
                else:
                    predictions[img_name] = np.zeros(img_shape, dtype=np.int32)
            else:
                # Try alternative file pattern
                dat_files = [f for f in (output_path / img_name).glob("*.dat") if f.name != "file_map.dat"]
                if dat_files:
                    result = joblib.load(dat_files[0])
                    if isinstance(result, dict):
                        first_key = list(result.keys())[0] if result else None
                        if first_key and isinstance(result.get(first_key), dict) and "contour" in result.get(first_key, {}):
                            inst_map = reconstruct_instance_map(result, img_shape)
                            predictions[img_name] = inst_map

        except Exception as e:
            failed.append((img_name, str(e)))
            if len(failed) <= 3:
                print(f"Failed: {img_name} - {e}")

    if failed:
        print(f"\n{len(failed)} images failed (showing first 5):")
        for name, error in failed[:5]:
            print(f"  {name}: {error}")

    print(f"\nSuccessfully processed {len(predictions)}/{len(image_paths)} images")
    return predictions


def load_ground_truth_masks(manifest_path: str) -> Tuple[Dict[str, np.ndarray], List[str]]:
    """
    Load ground truth masks from manifest.

    Args:
        manifest_path: Path to manifest.json

    Returns:
        Dictionary mapping image IDs to masks, list of image paths
    """
    with open(manifest_path, "r") as f:
        manifest = json.load(f)

    gt_masks = {}
    image_paths = []

    for entry in manifest:
        image_id = entry["image_id"]
        mask_path = entry["mask_path"]
        image_paths.append(entry["image_path"])

        gt_masks[image_id] = np.load(mask_path)

    return gt_masks, image_paths


def benchmark_on_dataset(
    dataset_name: str,
    images_dir: str,
    masks_dir: str,
    output_dir: str,
    model_name: str = "hovernet_fast-pannuke",
    max_images: Optional[int] = None
) -> Dict:
    """
    Run benchmark on a dataset.

    Args:
        dataset_name: Name of the dataset
        images_dir: Directory containing images
        masks_dir: Directory containing ground truth masks
        output_dir: Directory for outputs
        model_name: HoVerNet model name
        max_images: Maximum number of images to process (for testing)

    Returns:
        Benchmark results dictionary
    """
    print(f"\n{'='*60}")
    print(f"Benchmarking on {dataset_name}")
    print(f"{'='*60}\n")

    images_path = Path(images_dir)
    masks_path = Path(masks_dir)
    output_path = Path(output_dir) / dataset_name

    # Find all images
    image_files = sorted(list(images_path.glob("*.png")) + list(images_path.glob("*.jpg")))
    if max_images:
        image_files = image_files[:max_images]

    print(f"Found {len(image_files)} images")

    # Run inference
    predictions = run_hovernet_inference_simple(
        image_paths=[str(f) for f in image_files],
        output_dir=str(output_path / "predictions"),
        model_name=model_name
    )

    # Load ground truth and compute metrics
    pred_masks = []
    gt_masks = []
    image_ids = []

    for img_file in image_files:
        img_name = img_file.stem

        # Find corresponding mask
        mask_file = masks_path / f"{img_name}.npy"
        if not mask_file.exists():
            mask_file = masks_path / f"{img_name}.tif"

        if not mask_file.exists():
            print(f"Warning: No mask for {img_name}")
            continue

        if img_name not in predictions:
            print(f"Warning: No prediction for {img_name}")
            continue

        # Load masks
        if mask_file.suffix == ".npy":
            gt_mask = np.load(mask_file)
        else:
            import tifffile
            gt_mask = tifffile.imread(mask_file)

        pred_masks.append(predictions[img_name])
        gt_masks.append(gt_mask)
        image_ids.append(img_name)

    print(f"\nComputing metrics for {len(pred_masks)} images...")

    # Compute metrics
    summary, per_image = compute_dataset_metrics(pred_masks, gt_masks, image_ids)

    # Print results
    print(f"\n{dataset_name} Results:")
    print("-" * 50)
    print(f"{'Metric':<10} {'Mean':>10} {'Std':>10} {'Median':>10}")
    print("-" * 50)
    for metric, stats in summary.items():
        print(f"{metric:<10} {stats['mean']:>10.4f} {stats['std']:>10.4f} {stats['median']:>10.4f}")

    # Save results
    results = {
        "dataset": dataset_name,
        "model": model_name,
        "n_images": len(pred_masks),
        "summary": summary,
        "per_image": per_image
    }

    results_file = output_path / "results.json"
    with open(results_file, "w") as f:
        json.dump(results, f, indent=2, default=float)

    print(f"\nResults saved to {results_file}")

    return results


def benchmark_nuinsseg(
    base_dir: str = "/home/kyu_insilica_co/qMAP/input_data/nuinsseg",
    output_dir: str = "/home/kyu_insilica_co/qMAP/output/benchmark",
    model_name: str = "hovernet_fast-pannuke",
    max_images: Optional[int] = None
):
    """
    Benchmark HoVerNet on NuInsSeg dataset.
    """
    print("Benchmarking on NuInsSeg dataset...")

    base_path = Path(base_dir)
    output_path = Path(output_dir) / "nuinsseg"
    output_path.mkdir(parents=True, exist_ok=True)

    # Collect all images and masks
    all_images = []
    all_masks = []
    all_ids = []

    for tissue_dir in sorted(base_path.iterdir()):
        if not tissue_dir.is_dir():
            continue

        image_dir = tissue_dir / "tissue images"
        mask_dir = tissue_dir / "label masks"

        if not image_dir.exists():
            continue

        for img_file in sorted(image_dir.glob("*.png")):
            mask_file = mask_dir / f"{img_file.stem}.tif"
            if mask_file.exists():
                all_images.append(str(img_file))
                all_masks.append(str(mask_file))
                all_ids.append(img_file.stem)

    print(f"Found {len(all_images)} image-mask pairs")

    if max_images:
        all_images = all_images[:max_images]
        all_masks = all_masks[:max_images]
        all_ids = all_ids[:max_images]

    # Run inference
    predictions = run_hovernet_inference_simple(
        image_paths=all_images,
        output_dir=str(output_path / "predictions"),
        model_name=model_name
    )

    # Compute metrics
    import tifffile

    pred_masks = []
    gt_masks = []
    valid_ids = []

    for img_path, mask_path, img_id in zip(all_images, all_masks, all_ids):
        if img_id not in predictions:
            continue

        gt_mask = tifffile.imread(mask_path)
        pred_mask = predictions[img_id]

        # Ensure same shape
        if pred_mask.shape != gt_mask.shape:
            print(f"Shape mismatch for {img_id}: pred {pred_mask.shape} vs gt {gt_mask.shape}")
            continue

        pred_masks.append(pred_mask)
        gt_masks.append(gt_mask)
        valid_ids.append(img_id)

    print(f"\nComputing metrics for {len(pred_masks)} images...")
    summary, per_image = compute_dataset_metrics(pred_masks, gt_masks, valid_ids)

    # Print and save results
    print(f"\nNuInsSeg Results:")
    print("-" * 50)
    print(f"{'Metric':<10} {'Mean':>10} {'Std':>10}")
    print("-" * 50)
    for metric, stats in summary.items():
        print(f"{metric:<10} {stats['mean']:>10.4f} {stats['std']:>10.4f}")

    results = {
        "dataset": "NuInsSeg",
        "model": model_name,
        "n_images": len(pred_masks),
        "summary": summary,
        "per_image": per_image
    }

    with open(output_path / "results.json", "w") as f:
        json.dump(results, f, indent=2, default=float)

    return results


def benchmark_monuseg(
    base_dir: str = "/home/kyu_insilica_co/qMAP/input_data/monuseg",
    output_dir: str = "/home/kyu_insilica_co/qMAP/output/benchmark",
    model_name: str = "hovernet_fast-pannuke",
    max_images: Optional[int] = None
):
    """
    Benchmark HoVerNet on MoNuSeg dataset.
    """
    return benchmark_on_dataset(
        dataset_name="monuseg",
        images_dir=str(Path(base_dir) / "images"),
        masks_dir=str(Path(base_dir) / "masks"),
        output_dir=output_dir,
        model_name=model_name,
        max_images=max_images
    )


def benchmark_cryonuseg(
    base_dir: str = "/home/kyu_insilica_co/qMAP/input_data/cryonuseg",
    output_dir: str = "/home/kyu_insilica_co/qMAP/output/benchmark",
    model_name: str = "hovernet_fast-pannuke",
    max_images: Optional[int] = None
):
    """
    Benchmark HoVerNet on CryoNuSeg dataset.
    """
    return benchmark_on_dataset(
        dataset_name="cryonuseg",
        images_dir=str(Path(base_dir) / "images"),
        masks_dir=str(Path(base_dir) / "masks"),
        output_dir=output_dir,
        model_name=model_name,
        max_images=max_images
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Benchmark HoVerNet on nuclei datasets")
    parser.add_argument("--dataset", type=str, default="nuinsseg",
                        choices=["nuinsseg", "monuseg", "cryonuseg", "all"],
                        help="Dataset to benchmark")
    parser.add_argument("--model", type=str, default="hovernet_fast-pannuke",
                        help="Pretrained model name")
    parser.add_argument("--output_dir", type=str,
                        default="/home/kyu_insilica_co/qMAP/output/benchmark",
                        help="Output directory")
    parser.add_argument("--max_images", type=int, default=None,
                        help="Max images to process (for testing)")

    args = parser.parse_args()

    if args.dataset == "nuinsseg" or args.dataset == "all":
        benchmark_nuinsseg(
            output_dir=args.output_dir,
            model_name=args.model,
            max_images=args.max_images
        )

    if args.dataset == "monuseg" or args.dataset == "all":
        benchmark_monuseg(
            output_dir=args.output_dir,
            model_name=args.model,
            max_images=args.max_images
        )

    if args.dataset == "cryonuseg" or args.dataset == "all":
        benchmark_cryonuseg(
            output_dir=args.output_dir,
            model_name=args.model,
            max_images=args.max_images
        )
