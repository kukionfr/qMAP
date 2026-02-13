#!/usr/bin/env python3
"""
Re-run HoVerNet on MoNuSeg 256x256 patches that failed in tile mode.
Then recompute metrics on all 82 images.
"""

import os
import sys
import json
import numpy as np
from pathlib import Path
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")
import joblib

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


def run_hovernet_on_small_images(image_paths, output_dir):
    """Run HoVerNet on small (256x256) images using wsi mode instead of tile mode."""
    from tiatoolbox.models.engine.nucleus_instance_segmentor import NucleusInstanceSegmentor
    from PIL import Image
    import torch

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")

    segmentor = NucleusInstanceSegmentor(
        pretrained_model="hovernet_fast-pannuke",
        batch_size=8,
        verbose=False
    )

    predictions = {}
    failed = []

    for img_path in tqdm(image_paths, desc="HoVerNet inference (small imgs)"):
        img_name = Path(img_path).stem
        result_dir = output_path / img_name

        # Check if prediction already exists
        result_file = result_dir / "0.dat"
        if result_file.exists():
            try:
                result = joblib.load(result_file)
                img = Image.open(img_path)
                img_shape = (img.height, img.width)
                if isinstance(result, dict):
                    first_key = list(result.keys())[0] if result else None
                    if first_key and isinstance(result.get(first_key), dict) and "contour" in result.get(first_key, {}):
                        inst_map = reconstruct_instance_map(result, img_shape)
                        predictions[img_name] = inst_map
                        continue
                    elif "inst_map" in result:
                        predictions[img_name] = result["inst_map"]
                        continue
            except Exception:
                pass

        try:
            img = Image.open(img_path)
            img_array = np.array(img)
            img_shape = (img.height, img.width)

            # For small images, pad to at least 270x270 to ensure tile mode works
            # HoVerNet uses 256x256 tiles, so a 256x256 image may fail in tile mode
            pad_h = max(0, 270 - img_shape[0])
            pad_w = max(0, 270 - img_shape[1])

            if pad_h > 0 or pad_w > 0:
                # Pad with reflected content
                padded = np.pad(img_array,
                    ((0, pad_h), (0, pad_w), (0, 0)),
                    mode='reflect')
                # Save temp padded image
                from PIL import Image as PILImage
                temp_path = output_path / f"_temp_{img_name}.png"
                PILImage.fromarray(padded).save(str(temp_path))
                input_path = str(temp_path)
            else:
                input_path = img_path

            # Run prediction
            output = segmentor.predict(
                imgs=[input_path],
                save_dir=str(result_dir),
                mode="tile",
                device=device,
                crash_on_exception=True
            )

            # Clean up temp file
            temp_path = output_path / f"_temp_{img_name}.png"
            if temp_path.exists():
                temp_path.unlink()

            # Load result
            if result_file.exists():
                result = joblib.load(result_file)
                if isinstance(result, dict):
                    first_key = list(result.keys())[0] if result else None
                    if first_key and isinstance(result.get(first_key), dict) and "contour" in result.get(first_key, {}):
                        # Reconstruct and crop to original size
                        full_inst_map = reconstruct_instance_map(result, (padded.shape[0] if pad_h > 0 else img_shape[0], padded.shape[1] if pad_w > 0 else img_shape[1]))
                        inst_map = full_inst_map[:img_shape[0], :img_shape[1]]
                        # Re-label to ensure contiguous IDs
                        from scipy.ndimage import label as scipy_label
                        inst_map_relabeled, _ = scipy_label(inst_map > 0)
                        # Actually preserve distinct nuclei
                        unique_ids = np.unique(inst_map)
                        new_map = np.zeros_like(inst_map)
                        new_id = 1
                        for uid in unique_ids:
                            if uid == 0:
                                continue
                            new_map[inst_map == uid] = new_id
                            new_id += 1
                        predictions[img_name] = new_map
                    elif "inst_map" in result:
                        full_map = result["inst_map"]
                        predictions[img_name] = full_map[:img_shape[0], :img_shape[1]]
                    else:
                        predictions[img_name] = np.zeros(img_shape, dtype=np.int32)
                else:
                    predictions[img_name] = np.zeros(img_shape, dtype=np.int32)
            else:
                # Try alternative file pattern
                dat_files = [f for f in result_dir.glob("*.dat") if f.name != "file_map.dat"]
                if dat_files:
                    result = joblib.load(dat_files[0])
                    if isinstance(result, dict):
                        first_key = list(result.keys())[0] if result else None
                        if first_key and isinstance(result.get(first_key), dict) and "contour" in result.get(first_key, {}):
                            full_inst_map = reconstruct_instance_map(result, (padded.shape[0] if pad_h > 0 else img_shape[0], padded.shape[1] if pad_w > 0 else img_shape[1]))
                            inst_map = full_inst_map[:img_shape[0], :img_shape[1]]
                            predictions[img_name] = inst_map
                else:
                    predictions[img_name] = np.zeros(img_shape, dtype=np.int32)

        except Exception as e:
            failed.append((img_name, str(e)))
            if len(failed) <= 5:
                print(f"Failed: {img_name} - {e}")

    if failed:
        print(f"\n{len(failed)} images failed (showing first 5):")
        for name, error in failed[:5]:
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

    # Separate large (TCGA) and small (numbered) images
    large_images = [f for f in image_files if f.stem.startswith("TCGA")]
    small_images = [f for f in image_files if not f.stem.startswith("TCGA")]
    print(f"  Large (TCGA): {len(large_images)}")
    print(f"  Small (patches): {len(small_images)}")

    # Load existing TCGA predictions
    from PIL import Image
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

    # Run inference on small images
    small_preds = run_hovernet_on_small_images(
        [str(f) for f in small_images],
        str(pred_dir)
    )
    predictions.update(small_preds)
    print(f"Total predictions: {len(predictions)}")

    # Load ground truth and compute metrics
    pred_masks = []
    gt_masks = []
    image_ids = []

    for img_file in image_files:
        img_name = img_file.stem
        mask_file = masks_path / f"{img_name}.npy"
        if not mask_file.exists():
            print(f"Warning: No mask for {img_name}")
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

    # Print results
    print(f"\nMoNuSeg Results (HoVerNet, n={len(pred_masks)}):")
    print("-" * 50)
    print(f"{'Metric':<10} {'Mean':>10} {'Std':>10} {'Median':>10}")
    print("-" * 50)
    for metric, stats in summary.items():
        print(f"{metric:<10} {stats['mean']:>10.4f} {stats['std']:>10.4f} {stats['median']:>10.4f}")

    # Save results to multiple locations
    results = {
        "dataset": "MoNuSeg",
        "model": "HoVerNet (hovernet_fast-pannuke)",
        "n_images": len(pred_masks),
        "summary": summary,
        "per_image": per_image,
        "note": "Re-evaluated on all 82 MoNuSeg images (32 TCGA + 50 test patches)"
    }

    # Save to monuseg/results.json (main)
    with open(output_dir / "results.json", "w") as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\nSaved to {output_dir / 'results.json'}")

    # Save to monuseg/hovernet/results.json
    hovernet_dir = output_dir / "hovernet"
    hovernet_dir.mkdir(parents=True, exist_ok=True)

    results_hovernet = results.copy()
    results_hovernet["model"] = "hovernet"
    with open(hovernet_dir / "results.json", "w") as f:
        json.dump(results_hovernet, f, indent=2, default=float)
    print(f"Saved to {hovernet_dir / 'results.json'}")


if __name__ == "__main__":
    main()
