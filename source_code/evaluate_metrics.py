#!/usr/bin/env python3
"""
Compute segmentation evaluation metrics: Dice, AJI, PQ, DQ, SQ.

Implements standard nuclei segmentation metrics following:
- Kumar et al. (2017) for AJI
- Kirillov et al. (2019) for PQ

Optimized with bounding box early rejection and shared IoU matrix computation.
"""

import numpy as np
from scipy.ndimage import label as scipy_label
from scipy.optimize import linear_sum_assignment
from typing import Tuple, Dict, List
from collections import OrderedDict


def dice_coefficient(pred: np.ndarray, gt: np.ndarray) -> float:
    """
    Compute Dice coefficient (F1 score) for binary masks.

    Args:
        pred: Predicted binary mask
        gt: Ground truth binary mask

    Returns:
        Dice coefficient [0, 1]
    """
    pred_binary = pred > 0
    gt_binary = gt > 0

    intersection = np.logical_and(pred_binary, gt_binary).sum()
    total = pred_binary.sum() + gt_binary.sum()

    if total == 0:
        return 1.0  # Both empty

    return 2.0 * intersection / total


def get_instance_info_fast(instance_mask: np.ndarray) -> List[dict]:
    """
    Extract instance info with bounding boxes for fast IoU computation.

    Returns list of dicts with keys: id, bbox (rmin, rmax, cmin, cmax), area
    """
    instances = []
    unique_ids = np.unique(instance_mask)
    unique_ids = unique_ids[unique_ids > 0]

    for inst_id in unique_ids:
        coords = np.where(instance_mask == inst_id)
        if len(coords[0]) == 0:
            continue
        rmin, rmax = coords[0].min(), coords[0].max() + 1
        cmin, cmax = coords[1].min(), coords[1].max() + 1
        instances.append({
            'id': inst_id,
            'bbox': (rmin, rmax, cmin, cmax),
            'area': len(coords[0])
        })

    return instances


def compute_iou_fast(mask: np.ndarray, id1: int, bbox1: tuple,
                     id2: int, bbox2: tuple) -> float:
    """
    Compute IoU between two instances using bounding box intersection.
    Only examines the overlapping region of bounding boxes.
    """
    r1min, r1max, c1min, c1max = bbox1
    r2min, r2max, c2min, c2max = bbox2

    # Check bounding box overlap
    rmin = max(r1min, r2min)
    rmax = min(r1max, r2max)
    cmin = max(c1min, c2min)
    cmax = min(c1max, c2max)

    if rmin >= rmax or cmin >= cmax:
        return 0.0

    # Only examine overlapping region
    region = mask[rmin:rmax, cmin:cmax]
    m1 = region == id1
    m2 = region == id2

    intersection = np.logical_and(m1, m2).sum()
    if intersection == 0:
        return 0.0

    # For union, we need areas outside the overlap too
    # area1 in bbox1 that's not in overlap + area2 in bbox2 not in overlap + overlap region union
    # Simpler: use total areas
    # We pass areas separately to avoid recomputing
    return intersection  # Return raw intersection, compute IoU outside


def compute_iou_matrix_fast(pred: np.ndarray, gt: np.ndarray,
                            gt_info: List[dict], pred_info: List[dict]) -> np.ndarray:
    """
    Compute IoU matrix between GT and pred instances using bounding box optimization.
    """
    n_gt = len(gt_info)
    n_pred = len(pred_info)
    iou_matrix = np.zeros((n_gt, n_pred))

    for i, gi in enumerate(gt_info):
        gb = gi['bbox']
        gid = gi['id']

        for j, pi in enumerate(pred_info):
            pb = pi['bbox']
            pid = pi['id']

            # Bounding box overlap check
            rmin = max(gb[0], pb[0])
            rmax = min(gb[1], pb[1])
            cmin = max(gb[2], pb[2])
            cmax = min(gb[3], pb[3])

            if rmin >= rmax or cmin >= cmax:
                continue

            # Compute in overlapping region only
            region_gt = gt[rmin:rmax, cmin:cmax]
            region_pred = pred[rmin:rmax, cmin:cmax]
            m1 = region_gt == gid
            m2 = region_pred == pid

            intersection = np.logical_and(m1, m2).sum()
            if intersection == 0:
                continue

            union = gi['area'] + pi['area'] - intersection
            iou_matrix[i, j] = intersection / union if union > 0 else 0.0

    return iou_matrix


def aggregated_jaccard_index_from_iou(
    pred: np.ndarray, gt: np.ndarray,
    gt_info: List[dict], pred_info: List[dict],
    iou_matrix: np.ndarray
) -> float:
    """
    Compute AJI using pre-computed IoU matrix.
    """
    if len(gt_info) == 0 and len(pred_info) == 0:
        return 1.0
    if len(gt_info) == 0 or len(pred_info) == 0:
        return 0.0

    # Hungarian matching to maximize IoU
    row_ind, col_ind = linear_sum_assignment(-iou_matrix)

    # Compute AJI components
    intersection_total = 0
    union_total = 0
    matched_pred_idx = set()
    matched_row_set = set()

    for i, j in zip(row_ind, col_ind):
        if iou_matrix[i, j] > 0:
            gi = gt_info[i]
            pi = pred_info[j]

            # Recompute intersection from IoU: IoU = I / (A1 + A2 - I)
            # I = IoU * (A1 + A2) / (1 + IoU)
            iou_val = iou_matrix[i, j]
            inter = iou_val * (gi['area'] + pi['area']) / (1 + iou_val)
            union = gi['area'] + pi['area'] - inter

            intersection_total += inter
            union_total += union
            matched_pred_idx.add(j)
            matched_row_set.add(i)

    # Add unmatched GT areas to union
    for i, gi in enumerate(gt_info):
        if i not in matched_row_set:
            union_total += gi['area']

    # Add unmatched pred areas (false positives)
    for j, pi in enumerate(pred_info):
        if j not in matched_pred_idx:
            union_total += pi['area']

    if union_total == 0:
        return 0.0

    return intersection_total / union_total


def panoptic_quality_from_iou(
    gt_info: List[dict], pred_info: List[dict],
    iou_matrix: np.ndarray,
    iou_threshold: float = 0.5
) -> Tuple[float, float, float]:
    """
    Compute PQ, DQ, SQ using pre-computed IoU matrix.
    """
    if len(gt_info) == 0 and len(pred_info) == 0:
        return 1.0, 1.0, 1.0
    if len(gt_info) == 0:
        return 0.0, 0.0, 0.0
    if len(pred_info) == 0:
        return 0.0, 0.0, 0.0

    # Hungarian matching
    row_ind, col_ind = linear_sum_assignment(-iou_matrix)

    tp = 0
    matched_ious = []

    for i, j in zip(row_ind, col_ind):
        if iou_matrix[i, j] >= iou_threshold:
            tp += 1
            matched_ious.append(iou_matrix[i, j])

    fn = len(gt_info) - tp
    fp = len(pred_info) - tp

    if tp == 0:
        return 0.0, 0.0, 0.0

    dq = tp / (tp + 0.5 * fp + 0.5 * fn)
    sq = np.mean(matched_ious)
    pq = dq * sq

    return pq, dq, sq


def compute_all_metrics(pred: np.ndarray, gt: np.ndarray) -> Dict[str, float]:
    """
    Compute all segmentation metrics with optimized shared IoU matrix.

    Args:
        pred: Predicted instance mask
        gt: Ground truth instance mask

    Returns:
        Dictionary with Dice, AJI, PQ, DQ, SQ
    """
    dice = dice_coefficient(pred, gt)

    gt_info = get_instance_info_fast(gt)
    pred_info = get_instance_info_fast(pred)

    if len(gt_info) == 0 and len(pred_info) == 0:
        return OrderedDict([
            ("Dice", dice), ("AJI", 1.0),
            ("PQ", 1.0), ("DQ", 1.0), ("SQ", 1.0)
        ])
    if len(gt_info) == 0 or len(pred_info) == 0:
        return OrderedDict([
            ("Dice", dice), ("AJI", 0.0),
            ("PQ", 0.0), ("DQ", 0.0), ("SQ", 0.0)
        ])

    # Compute IoU matrix ONCE
    iou_matrix = compute_iou_matrix_fast(pred, gt, gt_info, pred_info)

    # Use shared IoU matrix for both metrics
    aji = aggregated_jaccard_index_from_iou(pred, gt, gt_info, pred_info, iou_matrix)
    pq, dq, sq = panoptic_quality_from_iou(gt_info, pred_info, iou_matrix)

    return OrderedDict([
        ("Dice", dice),
        ("AJI", aji),
        ("PQ", pq),
        ("DQ", dq),
        ("SQ", sq)
    ])


# Keep old API for backwards compatibility
def get_instance_info(instance_mask: np.ndarray) -> Dict[int, np.ndarray]:
    """Extract instance information from instance mask."""
    instances = {}
    unique_ids = np.unique(instance_mask)
    unique_ids = unique_ids[unique_ids > 0]
    for inst_id in unique_ids:
        instances[inst_id] = (instance_mask == inst_id)
    return instances


def compute_iou(mask1: np.ndarray, mask2: np.ndarray) -> float:
    """Compute IoU between two binary masks."""
    intersection = np.logical_and(mask1, mask2).sum()
    union = np.logical_or(mask1, mask2).sum()
    if union == 0:
        return 0.0
    return intersection / union


def aggregated_jaccard_index(pred: np.ndarray, gt: np.ndarray) -> float:
    """Compute AJI (legacy API, now uses optimized path internally)."""
    gt_info = get_instance_info_fast(gt)
    pred_info = get_instance_info_fast(pred)
    if len(gt_info) == 0 and len(pred_info) == 0:
        return 1.0
    if len(gt_info) == 0 or len(pred_info) == 0:
        return 0.0
    iou_matrix = compute_iou_matrix_fast(pred, gt, gt_info, pred_info)
    return aggregated_jaccard_index_from_iou(pred, gt, gt_info, pred_info, iou_matrix)


def panoptic_quality(pred: np.ndarray, gt: np.ndarray,
                     iou_threshold: float = 0.5) -> Tuple[float, float, float]:
    """Compute PQ (legacy API, now uses optimized path internally)."""
    gt_info = get_instance_info_fast(gt)
    pred_info = get_instance_info_fast(pred)
    if len(gt_info) == 0 and len(pred_info) == 0:
        return 1.0, 1.0, 1.0
    if len(gt_info) == 0 or len(pred_info) == 0:
        return 0.0, 0.0, 0.0
    iou_matrix = compute_iou_matrix_fast(pred, gt, gt_info, pred_info)
    return panoptic_quality_from_iou(gt_info, pred_info, iou_matrix, iou_threshold)


def compute_dataset_metrics(pred_masks: list, gt_masks: list,
                            image_ids: list = None) -> Tuple[Dict, list]:
    """
    Compute metrics across a dataset.

    Args:
        pred_masks: List of predicted instance masks
        gt_masks: List of ground truth instance masks
        image_ids: Optional list of image identifiers

    Returns:
        Tuple of (summary_stats, per_image_results)
    """
    from tqdm import tqdm

    n_samples = len(pred_masks)
    if image_ids is None:
        image_ids = [f"image_{i}" for i in range(n_samples)]

    all_results = []

    for idx in tqdm(range(n_samples), desc="Computing metrics"):
        metrics = compute_all_metrics(pred_masks[idx], gt_masks[idx])
        metrics["image_id"] = image_ids[idx]
        all_results.append(metrics)

    # Compute summary statistics
    metric_names = ["Dice", "AJI", "PQ", "DQ", "SQ"]
    summary = {}

    if len(all_results) == 0:
        for metric in metric_names:
            summary[metric] = {
                "mean": 0.0,
                "std": 0.0,
                "median": 0.0,
                "min": 0.0,
                "max": 0.0
            }
        return summary, all_results

    for metric in metric_names:
        values = [r[metric] for r in all_results]
        summary[metric] = {
            "mean": np.mean(values),
            "std": np.std(values),
            "median": np.median(values),
            "min": np.min(values),
            "max": np.max(values)
        }

    return summary, all_results


if __name__ == "__main__":
    # Test with synthetic data
    print("Testing evaluation metrics with synthetic data...")

    # Create synthetic masks
    np.random.seed(42)

    # Ground truth: 3 nuclei
    gt = np.zeros((256, 256), dtype=np.int32)
    gt[50:80, 50:80] = 1
    gt[100:150, 100:150] = 2
    gt[180:220, 180:220] = 3

    # Prediction: 3 nuclei with some offset
    pred = np.zeros((256, 256), dtype=np.int32)
    pred[55:85, 55:85] = 1  # Shifted
    pred[100:145, 105:145] = 2  # Slightly different
    pred[180:220, 180:220] = 3  # Perfect match

    metrics = compute_all_metrics(pred, gt)

    print("\nMetrics:")
    for name, value in metrics.items():
        print(f"  {name}: {value:.4f}")
