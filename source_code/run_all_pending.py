#!/usr/bin/env python3
"""
Run all pending benchmark evaluations:
1. Re-evaluate CryoNuSeg CellViT from saved predictions (inference done, metrics needed)
2. Re-evaluate MoNuSeg HoVerNet with fixed masks + optimized metrics
3. Re-evaluate MoNuSeg StarDist with fixed masks + optimized metrics
4. Run CellViT on MoNuSeg (inference + metrics)
5. Run CellViT on NuInsSeg (inference + metrics)
"""

import os
import sys
import json
import numpy as np
from pathlib import Path
from tqdm import tqdm
import joblib
import time

sys.path.insert(0, '/home/kyu_insilica_co/qMAP/source_code')
os.chdir('/home/kyu_insilica_co/qMAP/source_code')

from evaluate_metrics import compute_all_metrics, compute_dataset_metrics

BASE = Path('/home/kyu_insilica_co/qMAP')


def save_results(results, output_path):
    """Save results JSON."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"Results saved to {output_path}")


def evaluate_from_predictions(pred_masks, gt_masks, image_ids, dataset_name, model_name, output_path):
    """Evaluate and save results."""
    print(f"\nComputing metrics for {len(pred_masks)} images ({model_name} on {dataset_name})...")
    start = time.time()
    summary, per_image = compute_dataset_metrics(pred_masks, gt_masks, image_ids)
    elapsed = time.time() - start
    print(f"Metrics computed in {elapsed:.1f}s")

    # Print
    print(f"\n{dataset_name} Results ({model_name}):")
    print("-" * 50)
    print(f"{'Metric':<10} {'Mean':>10} {'Std':>10}")
    print("-" * 50)
    for metric, stats in summary.items():
        print(f"{metric:<10} {stats['mean']:>10.4f} {stats['std']:>10.4f}")

    results = {
        "dataset": dataset_name,
        "model": model_name,
        "n_images": len(pred_masks),
        "summary": summary,
        "per_image": per_image
    }
    save_results(results, output_path)
    return results


def reconstruct_instance_map(nuclei_dict, image_shape):
    """Reconstruct instance map from HoVerNet contour dictionary."""
    from skimage.draw import polygon
    inst_map = np.zeros(image_shape, dtype=np.int32)
    for idx, (nuc_id, nuc_info) in enumerate(nuclei_dict.items(), 1):
        contour = nuc_info.get("contour", [])
        if len(contour) > 0:
            contour = np.array(contour)
            rr, cc = polygon(contour[:, 1], contour[:, 0], image_shape)
            inst_map[rr, cc] = idx
    return inst_map


# ============================================================
# 1. Re-evaluate CryoNuSeg CellViT
# ============================================================
print("=" * 60)
print("1. Re-evaluating CryoNuSeg CellViT (from saved predictions)")
print("=" * 60)

cryo_pred_dir = BASE / 'output/benchmark/cryonuseg/cellvit/predictions'
cryo_masks_dir = BASE / 'input_data/cryonuseg/masks'
cryo_images_dir = BASE / 'input_data/cryonuseg/images'

pred_masks, gt_masks, image_ids = [], [], []
for pred_file in sorted(cryo_pred_dir.glob('*.npy')):
    img_name = pred_file.stem
    mask_file = cryo_masks_dir / f"{img_name}.npy"
    if mask_file.exists():
        pred_masks.append(np.load(pred_file))
        gt_masks.append(np.load(mask_file))
        image_ids.append(img_name)

evaluate_from_predictions(
    pred_masks, gt_masks, image_ids,
    "cryonuseg", "cellvit",
    BASE / 'output/benchmark/cryonuseg/cellvit/results.json'
)

# ============================================================
# 2. Re-evaluate MoNuSeg HoVerNet
# ============================================================
print("\n" + "=" * 60)
print("2. Re-evaluating MoNuSeg HoVerNet (with fixed masks)")
print("=" * 60)

monuseg_masks_dir = BASE / 'input_data/monuseg/masks'
monuseg_hovernet_dir = BASE / 'output/benchmark/monuseg/predictions'

pred_masks, gt_masks, image_ids = [], [], []

for dat_dir in sorted(monuseg_hovernet_dir.iterdir()):
    if not dat_dir.is_dir():
        continue
    dat_file = dat_dir / '0.dat'
    if not dat_file.exists():
        continue

    img_name = dat_dir.name
    mask_file = monuseg_masks_dir / f"{img_name}.npy"
    if not mask_file.exists():
        continue

    # Load HoVerNet prediction (contour dictionary)
    nuclei_dict = joblib.load(dat_file)
    gt_mask = np.load(mask_file)

    # Reconstruct instance map
    pred_map = reconstruct_instance_map(nuclei_dict, gt_mask.shape[:2])

    pred_masks.append(pred_map)
    gt_masks.append(gt_mask)
    image_ids.append(img_name)

evaluate_from_predictions(
    pred_masks, gt_masks, image_ids,
    "monuseg", "hovernet",
    BASE / 'output/benchmark/monuseg/hovernet/results.json'
)

# ============================================================
# 3. Re-evaluate MoNuSeg StarDist
# ============================================================
print("\n" + "=" * 60)
print("3. Re-evaluating MoNuSeg StarDist (with fixed masks)")
print("=" * 60)

stardist_pred_dir = BASE / 'output/benchmark/monuseg/stardist/predictions'

pred_masks, gt_masks, image_ids = [], [], []
for pred_file in sorted(stardist_pred_dir.glob('*.npy')):
    img_name = pred_file.stem
    mask_file = monuseg_masks_dir / f"{img_name}.npy"
    if mask_file.exists():
        pred_map = np.load(pred_file)
        gt_mask = np.load(mask_file)
        if pred_map.shape == gt_mask.shape[:2]:
            pred_masks.append(pred_map)
            gt_masks.append(gt_mask)
            image_ids.append(img_name)

evaluate_from_predictions(
    pred_masks, gt_masks, image_ids,
    "monuseg", "stardist_2D_versatile_he",
    BASE / 'output/benchmark/monuseg/stardist/results.json'
)

# ============================================================
# 4. Run CellViT on MoNuSeg
# ============================================================
print("\n" + "=" * 60)
print("4. Running CellViT on MoNuSeg")
print("=" * 60)

from benchmark_cellvit import benchmark_monuseg
benchmark_monuseg(
    model_path=str(BASE / 'CellViT/models/pretrained/cellvit_sam_h.pth')
)

# ============================================================
# 5. Run CellViT on NuInsSeg
# ============================================================
print("\n" + "=" * 60)
print("5. Running CellViT on NuInsSeg")
print("=" * 60)

from benchmark_cellvit import benchmark_nuinsseg
benchmark_nuinsseg(
    model_path=str(BASE / 'CellViT/models/pretrained/cellvit_sam_h.pth')
)

print("\n" + "=" * 60)
print("ALL BENCHMARKS COMPLETE")
print("=" * 60)
