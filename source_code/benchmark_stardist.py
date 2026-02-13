#!/usr/bin/env python3
"""
Benchmark StarDist on nuclei segmentation datasets.

StarDist: Star-convex Polygons for Object Detection
- Uses 2D_versatile_he pretrained model (trained on MoNuSeg + TNBC)
- Expects H&E stained images
- Returns instance labels directly

Reference:
Schmidt et al., Cell Detection with Star-convex Polygons, MICCAI 2018
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


def normalize_image(img: np.ndarray) -> np.ndarray:
    """
    Normalize image for StarDist.

    StarDist expects images normalized per-channel to [0, 1] range.

    Args:
        img: Input image (H, W, C) or (H, W)

    Returns:
        Normalized image
    """
    from csbdeep.utils import normalize

    # StarDist's normalize function does percentile-based normalization
    # This is more robust to outliers than simple min-max
    return normalize(img, 1, 99.8)


def run_stardist_inference(
    image_paths: List[str],
    output_dir: str,
    model_name: str = "2D_versatile_he",
    prob_thresh: float = None,
    nms_thresh: float = None
) -> Dict[str, np.ndarray]:
    """
    Run StarDist inference on a list of images.

    Args:
        image_paths: List of paths to input images
        output_dir: Directory to save outputs
        model_name: Pretrained model name
            - "2D_versatile_he": Trained on MoNuSeg + TNBC (H&E images)
            - "2D_versatile_fluo": Trained on fluorescence images
            - "2D_paper_dsb2018": Trained on Data Science Bowl 2018
        prob_thresh: Probability threshold for detection (default: model default)
        nms_thresh: NMS threshold (default: model default)

    Returns:
        Dictionary mapping image names to instance masks
    """
    from stardist.models import StarDist2D
    from PIL import Image

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    print(f"Running StarDist ({model_name}) on {len(image_paths)} images...")

    # Load pretrained model
    model = StarDist2D.from_pretrained(model_name)

    # Get default thresholds if not specified
    if prob_thresh is None:
        prob_thresh = model.thresholds[0] if model.thresholds else 0.5
    if nms_thresh is None:
        nms_thresh = model.thresholds[1] if len(model.thresholds) > 1 else 0.4

    print(f"Using thresholds: prob={prob_thresh:.3f}, nms={nms_thresh:.3f}")

    predictions = {}
    failed = []

    for img_path in tqdm(image_paths, desc="StarDist inference"):
        img_name = Path(img_path).stem

        try:
            # Load image
            img = Image.open(img_path)
            if img.mode != 'RGB':
                img = img.convert('RGB')
            img_array = np.array(img)

            # Normalize
            img_normalized = normalize_image(img_array)

            # Run prediction
            labels, details = model.predict_instances(
                img_normalized,
                prob_thresh=prob_thresh,
                nms_thresh=nms_thresh
            )

            predictions[img_name] = labels.astype(np.int32)

            # Save prediction
            np.save(output_path / f"{img_name}.npy", labels)

        except Exception as e:
            failed.append((img_name, str(e)))
            if len(failed) <= 3:
                print(f"Failed: {img_name} - {e}")

    if failed:
        print(f"\n{len(failed)} images failed:")
        for name, error in failed[:5]:
            print(f"  {name}: {error}")

    print(f"\nSuccessfully processed {len(predictions)}/{len(image_paths)} images")
    return predictions


def run_stardist_inference_batch(
    image_paths: List[str],
    output_dir: str,
    model_name: str = "2D_versatile_he",
    batch_size: int = 8
) -> Dict[str, np.ndarray]:
    """
    Run StarDist inference with batching for efficiency.

    Note: StarDist doesn't natively support batched inference,
    but we can use GPU more efficiently by processing multiple
    images in sequence without reloading the model.

    Args:
        image_paths: List of paths to input images
        output_dir: Directory to save outputs
        model_name: Pretrained model name
        batch_size: Not used directly (StarDist processes one at a time)

    Returns:
        Dictionary mapping image names to instance masks
    """
    # For StarDist, batch processing is same as sequential
    # The model is loaded once and reused
    return run_stardist_inference(image_paths, output_dir, model_name)


def benchmark_stardist_on_dataset(
    dataset_name: str,
    images_dir: str,
    masks_dir: str,
    output_dir: str,
    model_name: str = "2D_versatile_he",
    max_images: Optional[int] = None
) -> Dict:
    """
    Run StarDist benchmark on a dataset.

    Args:
        dataset_name: Name of the dataset
        images_dir: Directory containing images
        masks_dir: Directory containing ground truth masks
        output_dir: Directory for outputs
        model_name: StarDist model name
        max_images: Maximum number of images to process (for testing)

    Returns:
        Benchmark results dictionary
    """
    print(f"\n{'='*60}")
    print(f"Benchmarking StarDist on {dataset_name}")
    print(f"{'='*60}\n")

    images_path = Path(images_dir)
    masks_path = Path(masks_dir)
    output_path = Path(output_dir) / dataset_name / "stardist"

    # Find all images
    image_files = sorted(
        list(images_path.glob("*.png")) +
        list(images_path.glob("*.jpg")) +
        list(images_path.glob("*.tif"))
    )

    if max_images:
        image_files = image_files[:max_images]

    print(f"Found {len(image_files)} images")

    # Run inference
    predictions = run_stardist_inference(
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
        mask_file = None
        for ext in ['.npy', '.tif', '.tiff', '.png']:
            candidate = masks_path / f"{img_name}{ext}"
            if candidate.exists():
                mask_file = candidate
                break

        if mask_file is None:
            print(f"Warning: No mask for {img_name}")
            continue

        if img_name not in predictions:
            print(f"Warning: No prediction for {img_name}")
            continue

        # Load masks
        if mask_file.suffix == '.npy':
            gt_mask = np.load(mask_file)
        elif mask_file.suffix in ['.tif', '.tiff']:
            import tifffile
            gt_mask = tifffile.imread(mask_file)
        else:
            from PIL import Image
            gt_mask = np.array(Image.open(mask_file))

        pred_mask = predictions[img_name]

        # Ensure same shape
        if pred_mask.shape != gt_mask.shape[:2]:
            print(f"Shape mismatch for {img_name}: pred {pred_mask.shape} vs gt {gt_mask.shape}")
            continue

        pred_masks.append(pred_mask)
        gt_masks.append(gt_mask)
        image_ids.append(img_name)

    print(f"\nComputing metrics for {len(pred_masks)} images...")

    # Compute metrics
    summary, per_image = compute_dataset_metrics(pred_masks, gt_masks, image_ids)

    # Print results
    print(f"\n{dataset_name} Results (StarDist):")
    print("-" * 50)
    print(f"{'Metric':<10} {'Mean':>10} {'Std':>10} {'Median':>10}")
    print("-" * 50)
    for metric, stats in summary.items():
        print(f"{metric:<10} {stats['mean']:>10.4f} {stats['std']:>10.4f} {stats['median']:>10.4f}")

    # Save results
    results = {
        "dataset": dataset_name,
        "model": f"stardist_{model_name}",
        "n_images": len(pred_masks),
        "summary": summary,
        "per_image": per_image
    }

    results_file = output_path / "results.json"
    output_path.mkdir(parents=True, exist_ok=True)
    with open(results_file, "w") as f:
        json.dump(results, f, indent=2, default=float)

    print(f"\nResults saved to {results_file}")

    return results


def benchmark_nuinsseg(
    base_dir: str = "/home/kyu_insilica_co/qMAP/input_data/nuinsseg",
    output_dir: str = "/home/kyu_insilica_co/qMAP/output/benchmark",
    model_name: str = "2D_versatile_he",
    max_images: Optional[int] = None
):
    """
    Benchmark StarDist on NuInsSeg dataset.
    """
    print("Benchmarking StarDist on NuInsSeg dataset...")

    base_path = Path(base_dir)
    output_path = Path(output_dir) / "nuinsseg" / "stardist"
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
    predictions = run_stardist_inference(
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
    print(f"\nNuInsSeg Results (StarDist):")
    print("-" * 50)
    print(f"{'Metric':<10} {'Mean':>10} {'Std':>10}")
    print("-" * 50)
    for metric, stats in summary.items():
        print(f"{metric:<10} {stats['mean']:>10.4f} {stats['std']:>10.4f}")

    results = {
        "dataset": "NuInsSeg",
        "model": f"stardist_{model_name}",
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
    model_name: str = "2D_versatile_he",
    max_images: Optional[int] = None
):
    """
    Benchmark StarDist on MoNuSeg dataset.

    Note: StarDist was trained on MoNuSeg, so expect good performance here.
    """
    return benchmark_stardist_on_dataset(
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
    model_name: str = "2D_versatile_he",
    max_images: Optional[int] = None
):
    """
    Benchmark StarDist on CryoNuSeg dataset.
    """
    return benchmark_stardist_on_dataset(
        dataset_name="cryonuseg",
        images_dir=str(Path(base_dir) / "images"),
        masks_dir=str(Path(base_dir) / "masks"),
        output_dir=output_dir,
        model_name=model_name,
        max_images=max_images
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Benchmark StarDist on nuclei datasets")
    parser.add_argument("--dataset", type=str, default="nuinsseg",
                        choices=["nuinsseg", "monuseg", "cryonuseg", "all"],
                        help="Dataset to benchmark")
    parser.add_argument("--model", type=str, default="2D_versatile_he",
                        choices=["2D_versatile_he", "2D_versatile_fluo", "2D_paper_dsb2018"],
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
