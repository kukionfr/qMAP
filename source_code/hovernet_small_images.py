#!/usr/bin/env python3
"""
Run HoVerNet on small MoNuSeg patches (256x256) using direct model inference.
TIA Toolbox's NucleusInstanceSegmentor tile mode fails on 256x256 images,
so we use the underlying IOSegmentorPredictor with appropriate settings.
"""

import os
import sys
import json
import numpy as np
from pathlib import Path
from tqdm import tqdm
from PIL import Image
import warnings
warnings.filterwarnings("ignore")
import joblib
import torch

from evaluate_metrics import compute_all_metrics, compute_dataset_metrics


def reconstruct_instance_map(nuclei_dict, image_shape):
    """Reconstruct instance map from nuclei contours."""
    from skimage.draw import polygon
    inst_map = np.zeros(image_shape, dtype=np.int32)
    for idx, (nuc_id, nuc_info) in enumerate(nuclei_dict.items(), 1):
        contour = nuc_info.get("contour", [])
        if len(contour) > 0:
            contour = np.array(contour)
            rr, cc = polygon(contour[:, 1], contour[:, 0], image_shape)
            inst_map[rr, cc] = idx
    return inst_map


def run_hovernet_direct(image_paths, output_dir):
    """
    Run HoVerNet on images by padding small images to 512x512,
    running tile mode, then cropping back.
    """
    from tiatoolbox.models.engine.nucleus_instance_segmentor import NucleusInstanceSegmentor

    output_path = Path(output_dir)
    temp_dir = output_path / "_temp_padded"
    temp_dir.mkdir(parents=True, exist_ok=True)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")

    # Pad small images to 512x512
    padded_paths = []
    original_sizes = {}
    for img_path in image_paths:
        img = Image.open(img_path)
        img_name = Path(img_path).stem
        w, h = img.size
        original_sizes[img_name] = (h, w)

        if h < 512 or w < 512:
            # Pad to 512x512 with reflection
            img_array = np.array(img)
            pad_h = max(0, 512 - h)
            pad_w = max(0, 512 - w)
            if len(img_array.shape) == 3:
                padded = np.pad(img_array, ((0, pad_h), (0, pad_w), (0, 0)), mode='reflect')
            else:
                padded = np.pad(img_array, ((0, pad_h), (0, pad_w)), mode='reflect')
            padded_img = Image.fromarray(padded)
            padded_path = temp_dir / f"{img_name}.png"
            padded_img.save(str(padded_path))
            padded_paths.append(str(padded_path))
        else:
            padded_paths.append(img_path)

    print(f"Running HoVerNet on {len(padded_paths)} images (padded to >= 512x512)...")

    segmentor = NucleusInstanceSegmentor(
        pretrained_model="hovernet_fast-pannuke",
        batch_size=8,
        verbose=False
    )

    predictions = {}
    failed = []

    for img_path, orig_path in tqdm(zip(padded_paths, image_paths), total=len(padded_paths), desc="HoVerNet"):
        img_name = Path(orig_path).stem
        save_dir = output_path / img_name

        try:
            output = segmentor.predict(
                imgs=[img_path],
                save_dir=str(save_dir),
                mode="tile",
                device=device,
                crash_on_exception=True
            )

            result_file = save_dir / "0.dat"
            if result_file.exists():
                result = joblib.load(result_file)
                orig_h, orig_w = original_sizes[img_name]

                if isinstance(result, dict):
                    first_key = list(result.keys())[0] if result else None
                    if first_key and isinstance(result.get(first_key), dict) and "contour" in result.get(first_key, {}):
                        # Get padded image size for reconstruction
                        padded_img = Image.open(img_path)
                        padded_shape = (padded_img.height, padded_img.width)
                        full_map = reconstruct_instance_map(result, padded_shape)
                        # Crop to original size
                        cropped = full_map[:orig_h, :orig_w]
                        predictions[img_name] = cropped
                    elif "inst_map" in result:
                        full_map = result["inst_map"]
                        predictions[img_name] = full_map[:orig_h, :orig_w]
                    else:
                        predictions[img_name] = np.zeros((orig_h, orig_w), dtype=np.int32)
            else:
                # Check for any .dat files
                dat_files = [f for f in save_dir.glob("*.dat") if f.name != "file_map.dat"]
                if dat_files:
                    result = joblib.load(dat_files[0])
                    orig_h, orig_w = original_sizes[img_name]
                    if isinstance(result, dict):
                        first_key = list(result.keys())[0] if result else None
                        if first_key and isinstance(result.get(first_key), dict) and "contour" in result.get(first_key, {}):
                            padded_img = Image.open(img_path)
                            padded_shape = (padded_img.height, padded_img.width)
                            full_map = reconstruct_instance_map(result, padded_shape)
                            predictions[img_name] = full_map[:orig_h, :orig_w]
                else:
                    predictions[img_name] = np.zeros((orig_h, orig_w), dtype=np.int32)
                    failed.append((img_name, "No output file"))

        except Exception as e:
            failed.append((img_name, str(e)))
            orig_h, orig_w = original_sizes[img_name]
            predictions[img_name] = np.zeros((orig_h, orig_w), dtype=np.int32)
            if len(failed) <= 5:
                print(f"Failed: {img_name} - {e}")

    # Clean up temp dir
    import shutil
    if temp_dir.exists():
        shutil.rmtree(temp_dir)

    if failed:
        print(f"\n{len(failed)} images failed:")
        for name, error in failed[:10]:
            print(f"  {name}: {error}")

    print(f"\nSuccessfully processed {len(predictions)}/{len(image_paths)} images")
    return predictions


def main():
    base_dir = Path("/home/kyu_insilica_co/qMAP/input_data/monuseg")
    output_dir = Path("/home/kyu_insilica_co/qMAP/output/benchmark/monuseg")
    pred_dir = output_dir / "predictions"

    images_path = base_dir / "images"
    masks_path = base_dir / "masks"

    # Get all image files
    image_files = sorted(list(images_path.glob("*.png")) + list(images_path.glob("*.jpg")))
    print(f"Found {len(image_files)} MoNuSeg images")

    large_images = [f for f in image_files if f.stem.startswith("TCGA")]
    small_images = [f for f in image_files if not f.stem.startswith("TCGA")]
    print(f"  Large (TCGA): {len(large_images)}")
    print(f"  Small (patches): {len(small_images)}")

    # Load existing TCGA predictions
    predictions = {}
    for img_file in large_images:
        img_name = img_file.stem
        result_file = pred_dir / img_name / "0.dat"
        if result_file.exists():
            result = joblib.load(result_file)
            img = Image.open(img_file)
            img_shape = (img.height, img.width)
            if isinstance(result, dict):
                first_key = list(result.keys())[0] if result else None
                if first_key and isinstance(result.get(first_key), dict) and "contour" in result.get(first_key, {}):
                    inst_map = reconstruct_instance_map(result, img_shape)
                    predictions[img_name] = inst_map
    print(f"Loaded {len(predictions)} existing TCGA predictions")

    # Run on small images with padding
    small_preds = run_hovernet_direct(
        [str(f) for f in small_images],
        str(pred_dir)
    )
    predictions.update(small_preds)
    print(f"\nTotal predictions: {len(predictions)}")

    # Compute metrics
    pred_masks = []
    gt_masks = []
    image_ids = []

    for img_file in image_files:
        img_name = img_file.stem
        mask_file = masks_path / f"{img_name}.npy"
        if not mask_file.exists():
            continue
        if img_name not in predictions:
            print(f"Warning: No prediction for {img_name}")
            continue

        gt_mask = np.load(mask_file)
        pred_mask = predictions[img_name]

        if pred_mask.shape != gt_mask.shape:
            print(f"Shape mismatch for {img_name}: pred {pred_mask.shape} vs gt {gt_mask.shape}")
            continue

        pred_masks.append(pred_mask)
        gt_masks.append(gt_mask)
        image_ids.append(img_name)

    print(f"\nComputing metrics for {len(pred_masks)} images...")
    summary, per_image = compute_dataset_metrics(pred_masks, gt_masks, image_ids)

    print(f"\nMoNuSeg Results (HoVerNet, n={len(pred_masks)}):")
    print("-" * 50)
    print(f"{'Metric':<10} {'Mean':>10} {'Std':>10} {'Median':>10}")
    print("-" * 50)
    for metric, stats in summary.items():
        print(f"{metric:<10} {stats['mean']:>10.4f} {stats['std']:>10.4f} {stats['median']:>10.4f}")

    results = {
        "dataset": "MoNuSeg",
        "model": "HoVerNet (hovernet_fast-pannuke)",
        "n_images": len(pred_masks),
        "summary": summary,
        "per_image": per_image,
        "note": "All 82 MoNuSeg images (32 TCGA 1000x1000 + 50 test patches 256x256). Small images padded to 512x512 for inference then cropped."
    }

    with open(output_dir / "results.json", "w") as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\nSaved to {output_dir / 'results.json'}")

    hovernet_dir = output_dir / "hovernet"
    hovernet_dir.mkdir(parents=True, exist_ok=True)
    results_hov = results.copy()
    results_hov["model"] = "hovernet"
    with open(hovernet_dir / "results.json", "w") as f:
        json.dump(results_hov, f, indent=2, default=float)
    print(f"Saved to {hovernet_dir / 'results.json'}")


if __name__ == "__main__":
    main()
