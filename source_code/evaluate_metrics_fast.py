#!/usr/bin/env python3
"""
Optimized evaluation metrics with bounding box pre-filtering.
"""

import numpy as np
from scipy.ndimage import label as scipy_label
from scipy.optimize import linear_sum_assignment
from typing import Tuple, Dict
from collections import OrderedDict
from tqdm import tqdm


def get_bounding_box(mask: np.ndarray) -> Tuple[int, int, int, int]:
    """Get bounding box (min_row, max_row, min_col, max_col) of a binary mask."""
    rows, cols = np.where(mask)
    if len(rows) == 0:
        return (0, 0, 0, 0)
    return (rows.min(), rows.max(), cols.min(), cols.max())


def boxes_overlap(box1: Tuple[int, int, int, int], box2: Tuple[int, int, int, int]) -> bool:
    """Check if two bounding boxes overlap."""
    min_r1, max_r1, min_c1, max_c1 = box1
    min_r2, max_r2, min_c2, max_c2 = box2
    
    # No overlap if one box is completely to the left/right/above/below the other
    if max_r1 < min_r2 or max_r2 < min_r1:
        return False
    if max_c1 < min_c2 or max_c2 < min_c1:
        return False
    
    return True


def dice_coefficient(pred: np.ndarray, gt: np.ndarray) -> float:
    """Compute Dice coefficient (F1 score) for binary masks."""
    pred_binary = pred > 0
    gt_binary = gt > 0
    
    intersection = np.logical_and(pred_binary, gt_binary).sum()
    total = pred_binary.sum() + gt_binary.sum()
    
    if total == 0:
        return 1.0
    
    return 2.0 * intersection / total


def get_instance_info(instance_mask: np.ndarray) -> Dict:
    """Extract instance information with bounding boxes."""
    instances = {}
    unique_ids = np.unique(instance_mask)
    unique_ids = unique_ids[unique_ids > 0]
    
    for inst_id in unique_ids:
        mask = (instance_mask == inst_id)
        bbox = get_bounding_box(mask)
        instances[inst_id] = {
            "mask": mask,
            "bbox": bbox
        }
    
    return instances


def compute_iou_fast(mask1: np.ndarray, mask2: np.ndarray) -> float:
    """Compute IoU between two binary masks (optimized)."""
    intersection = np.logical_and(mask1, mask2).sum()
    if intersection == 0:
        return 0.0
    
    union = np.logical_or(mask1, mask2).sum()
    if union == 0:
        return 0.0
    
    return intersection / union


def aggregated_jaccard_index(pred: np.ndarray, gt: np.ndarray) -> float:
    """Compute AJI with bounding box filtering."""
    gt_instances = get_instance_info(gt)
    pred_instances = get_instance_info(pred)
    
    if len(gt_instances) == 0 and len(pred_instances) == 0:
        return 1.0
    if len(gt_instances) == 0 or len(pred_instances) == 0:
        return 0.0
    
    gt_ids = list(gt_instances.keys())
    pred_ids = list(pred_instances.keys())
    
    # Compute IoU matrix with bounding box filtering
    iou_matrix = np.zeros((len(gt_ids), len(pred_ids)))
    for i, gt_id in enumerate(gt_ids):
        gt_bbox = gt_instances[gt_id]["bbox"]
        for j, pred_id in enumerate(pred_ids):
            pred_bbox = pred_instances[pred_id]["bbox"]
            
            # Skip if bounding boxes don't overlap
            if not boxes_overlap(gt_bbox, pred_bbox):
                continue
            
            # Compute IoU only for potentially overlapping instances
            iou_matrix[i, j] = compute_iou_fast(
                gt_instances[gt_id]["mask"],
                pred_instances[pred_id]["mask"]
            )
    
    # Hungarian matching
    row_ind, col_ind = linear_sum_assignment(-iou_matrix)
    
    # Compute AJI
    intersection_total = 0
    union_total = 0
    matched_pred = set()
    
    for i, j in zip(row_ind, col_ind):
        if iou_matrix[i, j] > 0:
            gt_mask = gt_instances[gt_ids[i]]["mask"]
            pred_mask = pred_instances[pred_ids[j]]["mask"]
            intersection_total += np.logical_and(gt_mask, pred_mask).sum()
            union_total += np.logical_or(gt_mask, pred_mask).sum()
            matched_pred.add(pred_ids[j])
    
    # Add unmatched GT
    for i, gt_id in enumerate(gt_ids):
        if i not in row_ind or iou_matrix[i, col_ind[list(row_ind).index(i)]] == 0:
            union_total += gt_instances[gt_id]["mask"].sum()
    
    # Add unmatched predictions
    for pred_id in pred_ids:
        if pred_id not in matched_pred:
            union_total += pred_instances[pred_id]["mask"].sum()
    
    if union_total == 0:
        return 0.0
    
    return intersection_total / union_total


def panoptic_quality(pred: np.ndarray, gt: np.ndarray,
                     iou_threshold: float = 0.5) -> Tuple[float, float, float]:
    """Compute PQ, DQ, SQ with bounding box filtering."""
    gt_instances = get_instance_info(gt)
    pred_instances = get_instance_info(pred)
    
    if len(gt_instances) == 0 and len(pred_instances) == 0:
        return 1.0, 1.0, 1.0
    if len(gt_instances) == 0:
        return 0.0, 0.0, 0.0
    if len(pred_instances) == 0:
        return 0.0, 0.0, 0.0
    
    gt_ids = list(gt_instances.keys())
    pred_ids = list(pred_instances.keys())
    
    # Compute IoU matrix with bounding box filtering
    iou_matrix = np.zeros((len(gt_ids), len(pred_ids)))
    for i, gt_id in enumerate(gt_ids):
        gt_bbox = gt_instances[gt_id]["bbox"]
        for j, pred_id in enumerate(pred_ids):
            pred_bbox = pred_instances[pred_id]["bbox"]
            
            if not boxes_overlap(gt_bbox, pred_bbox):
                continue
            
            iou_matrix[i, j] = compute_iou_fast(
                gt_instances[gt_id]["mask"],
                pred_instances[pred_id]["mask"]
            )
    
    # Hungarian matching
    row_ind, col_ind = linear_sum_assignment(-iou_matrix)
    
    # Count TP, FP, FN
    tp = 0
    matched_ious = []
    matched_gt = set()
    matched_pred = set()
    
    for i, j in zip(row_ind, col_ind):
        if iou_matrix[i, j] >= iou_threshold:
            tp += 1
            matched_ious.append(iou_matrix[i, j])
            matched_gt.add(gt_ids[i])
            matched_pred.add(pred_ids[j])
    
    fn = len(gt_ids) - len(matched_gt)
    fp = len(pred_ids) - len(matched_pred)
    
    if tp == 0:
        return 0.0, 0.0, 0.0
    
    dq = tp / (tp + 0.5 * fp + 0.5 * fn)
    sq = np.mean(matched_ious)
    pq = dq * sq
    
    return pq, dq, sq


def compute_all_metrics(pred: np.ndarray, gt: np.ndarray) -> Dict[str, float]:
    """Compute all segmentation metrics."""
    dice = dice_coefficient(pred, gt)
    aji = aggregated_jaccard_index(pred, gt)
    pq, dq, sq = panoptic_quality(pred, gt)
    
    return OrderedDict([
        ("Dice", dice),
        ("AJI", aji),
        ("PQ", pq),
        ("DQ", dq),
        ("SQ", sq)
    ])


def compute_dataset_metrics(pred_masks: list, gt_masks: list,
                            image_ids: list = None) -> Tuple[Dict, list]:
    """Compute metrics across a dataset."""
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
