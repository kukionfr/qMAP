#!/usr/bin/env python3
"""
Prepare CryoNuSeg dataset for nuclei segmentation benchmarking.

CryoNuSeg: Cryosectioned H&E Stained Nucleus Segmentation Dataset
- 30 image patches (512x512)
- Instance segmentation masks

This script standardizes the dataset format to match our benchmark pipeline:
- PNG images in images/
- NPY instance masks in masks/
"""

import os
import glob
import numpy as np
from pathlib import Path
from PIL import Image
import json
from tqdm import tqdm
from typing import Dict, Tuple, Optional
import re


def load_mask_from_various_formats(mask_path: Path) -> Optional[np.ndarray]:
    """
    Load instance mask from various file formats.

    Supports:
    - .npy: NumPy array
    - .tif/.tiff: TIFF image
    - .png: PNG image (may need relabeling)
    - .mat: MATLAB file

    Args:
        mask_path: Path to mask file

    Returns:
        Instance mask as numpy array, or None if failed
    """
    try:
        suffix = mask_path.suffix.lower()

        if suffix == '.npy':
            return np.load(mask_path).astype(np.int32)

        elif suffix in ['.tif', '.tiff']:
            import tifffile
            mask = tifffile.imread(mask_path)
            return mask.astype(np.int32)

        elif suffix == '.png':
            # PNG masks may be RGB or grayscale
            img = Image.open(mask_path)
            mask = np.array(img)

            if len(mask.shape) == 3:
                # RGB mask - need to convert to instance IDs
                # Often each unique RGB color is a different instance
                mask = rgb_to_instance_mask(mask)
            return mask.astype(np.int32)

        elif suffix == '.mat':
            from scipy.io import loadmat
            data = loadmat(mask_path)
            # Find the mask array (usually named 'inst_map', 'mask', or similar)
            for key in ['inst_map', 'mask', 'instance_map', 'label']:
                if key in data:
                    return data[key].astype(np.int32)
            # Try first non-metadata key
            for key, value in data.items():
                if not key.startswith('__') and isinstance(value, np.ndarray):
                    return value.astype(np.int32)

        print(f"Unknown mask format: {suffix}")
        return None

    except Exception as e:
        print(f"Error loading mask {mask_path}: {e}")
        return None


def rgb_to_instance_mask(rgb_mask: np.ndarray) -> np.ndarray:
    """
    Convert RGB mask to instance mask.

    Each unique RGB color becomes a unique instance ID.

    Args:
        rgb_mask: RGB image (H, W, 3)

    Returns:
        Instance mask (H, W) with unique integer IDs
    """
    h, w = rgb_mask.shape[:2]
    instance_mask = np.zeros((h, w), dtype=np.int32)

    # Get unique colors
    pixels = rgb_mask.reshape(-1, rgb_mask.shape[-1])
    unique_colors = np.unique(pixels, axis=0)

    instance_id = 0
    for color in unique_colors:
        # Skip background (typically black or white)
        if np.all(color == 0) or np.all(color == 255):
            continue

        mask = np.all(rgb_mask == color, axis=-1)
        if mask.sum() > 0:
            instance_id += 1
            instance_mask[mask] = instance_id

    return instance_mask


def relabel_mask(mask: np.ndarray) -> np.ndarray:
    """
    Relabel mask to have contiguous instance IDs starting from 1.

    Args:
        mask: Instance mask with possibly non-contiguous IDs

    Returns:
        Relabeled mask with IDs 0, 1, 2, ..., N
    """
    from scipy.ndimage import label

    binary_mask = mask > 0
    labeled, n_objects = label(binary_mask)

    return labeled.astype(np.int32)


def process_cryonuseg(
    input_dir: str,
    output_dir: str
) -> Dict:
    """
    Process CryoNuSeg dataset.

    Expected input structures:

    Structure 1 (Kaggle format):
    input_dir/
      tissue_images/
        *.png
      label_masks/
        *.tif

    Structure 2:
    input_dir/
      images/
        *.png
      masks/
        *.npy or *.tif or *.png

    Structure 3:
    input_dir/
      *.png (images)
      *_mask.png or *_label.tif (masks)

    Args:
        input_dir: Path to raw CryoNuSeg data
        output_dir: Path to output directory

    Returns:
        Manifest dictionary
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)

    # Create output directories
    images_dir = output_path / "images"
    masks_dir = output_path / "masks"
    images_dir.mkdir(parents=True, exist_ok=True)
    masks_dir.mkdir(parents=True, exist_ok=True)

    # Try to find the data structure
    image_files = []
    mask_map = {}  # Maps image stem to mask path

    # Structure 1: tissue_images/ and label_masks/
    if (input_path / "tissue_images").exists():
        image_dir = input_path / "tissue_images"
        mask_dir = input_path / "label_masks"

        for img_file in sorted(image_dir.glob("*.png")):
            image_files.append(img_file)
            # Try various mask naming conventions
            for ext in ['.tif', '.tiff', '.npy', '.png']:
                mask_file = mask_dir / f"{img_file.stem}{ext}"
                if mask_file.exists():
                    mask_map[img_file.stem] = mask_file
                    break

    # Structure 2: images/ and masks/
    elif (input_path / "images").exists():
        image_dir = input_path / "images"
        mask_dir = input_path / "masks" if (input_path / "masks").exists() else input_path / "labels"

        for img_file in sorted(list(image_dir.glob("*.png")) + list(image_dir.glob("*.jpg"))):
            image_files.append(img_file)
            for ext in ['.npy', '.tif', '.tiff', '.png', '.mat']:
                mask_file = mask_dir / f"{img_file.stem}{ext}"
                if mask_file.exists():
                    mask_map[img_file.stem] = mask_file
                    break

    # Structure 3: Mixed in same directory
    else:
        # Find all potential images (not masks)
        all_files = list(input_path.glob("*.png")) + list(input_path.glob("*.jpg"))

        for f in sorted(all_files):
            name = f.stem.lower()
            # Skip mask files
            if any(x in name for x in ['mask', 'label', 'seg', 'gt']):
                continue
            image_files.append(f)

            # Find corresponding mask
            for pattern in [f"{f.stem}_mask", f"{f.stem}_label", f"{f.stem}_seg", f"{f.stem}mask"]:
                for ext in ['.png', '.tif', '.tiff', '.npy']:
                    mask_file = input_path / f"{pattern}{ext}"
                    if mask_file.exists():
                        mask_map[f.stem] = mask_file
                        break
                if f.stem in mask_map:
                    break

    print(f"Found {len(image_files)} images")
    print(f"Found {len(mask_map)} masks")

    if not image_files:
        print(f"No images found in {input_path}")
        print("Expected directory structure:")
        print("  tissue_images/*.png + label_masks/*.tif")
        print("  OR images/*.png + masks/*.npy")
        return {"images": []}

    manifest = {"images": []}

    for img_file in tqdm(image_files, desc="Processing CryoNuSeg"):
        img_name = img_file.stem

        if img_name not in mask_map:
            print(f"Warning: No mask for {img_name}")
            continue

        # Load image
        img = Image.open(img_file)
        if img.mode != 'RGB':
            img = img.convert('RGB')
        img_array = np.array(img)
        img_shape = (img_array.shape[0], img_array.shape[1])

        # Load mask
        mask = load_mask_from_various_formats(mask_map[img_name])

        if mask is None:
            print(f"Warning: Could not load mask for {img_name}")
            continue

        # Ensure mask matches image size
        if mask.shape[:2] != img_shape:
            print(f"Warning: Size mismatch for {img_name}: image {img_shape} vs mask {mask.shape[:2]}")
            # Try to resize mask
            mask_img = Image.fromarray(mask.astype(np.uint16))
            mask_img = mask_img.resize((img_shape[1], img_shape[0]), Image.NEAREST)
            mask = np.array(mask_img).astype(np.int32)

        # Relabel to ensure contiguous IDs
        mask = relabel_mask(mask)

        # Count nuclei
        n_nuclei = mask.max()

        # Save outputs
        output_image_path = images_dir / f"{img_name}.png"
        output_mask_path = masks_dir / f"{img_name}.npy"

        img.save(output_image_path)
        np.save(output_mask_path, mask)

        manifest["images"].append({
            "image_id": img_name,
            "image_path": str(output_image_path),
            "mask_path": str(output_mask_path),
            "width": img_shape[1],
            "height": img_shape[0],
            "n_nuclei": int(n_nuclei)
        })

    # Save manifest
    manifest_path = output_path / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"\nProcessed {len(manifest['images'])} images")
    if manifest['images']:
        print(f"Total nuclei: {sum(m['n_nuclei'] for m in manifest['images'])}")
    print(f"Manifest saved to {manifest_path}")

    return manifest


def get_cryonuseg_stats(data_dir: str):
    """
    Print statistics about processed CryoNuSeg dataset.

    Args:
        data_dir: Path to processed CryoNuSeg directory
    """
    data_path = Path(data_dir)
    manifest_path = data_path / "manifest.json"

    if not manifest_path.exists():
        print(f"Manifest not found at {manifest_path}")
        return

    with open(manifest_path, "r") as f:
        manifest = json.load(f)

    images = manifest.get("images", [])

    if not images:
        print("No images in manifest")
        return

    print("\n" + "="*60)
    print("CryoNuSeg Dataset Statistics")
    print("="*60)

    print(f"\nTotal images: {len(images)}")

    # Image size stats
    widths = [m['width'] for m in images]
    heights = [m['height'] for m in images]
    nuclei_counts = [m['n_nuclei'] for m in images]

    print(f"\nImage dimensions:")
    print(f"  Width:  min={min(widths)}, max={max(widths)}, mean={np.mean(widths):.0f}")
    print(f"  Height: min={min(heights)}, max={max(heights)}, mean={np.mean(heights):.0f}")

    print(f"\nNuclei per image:")
    print(f"  Min: {min(nuclei_counts)}")
    print(f"  Max: {max(nuclei_counts)}")
    print(f"  Mean: {np.mean(nuclei_counts):.1f}")
    print(f"  Total: {sum(nuclei_counts)}")

    print("\nPer-image details:")
    print(f"{'Image':<40} {'Size':>15} {'Nuclei':>8}")
    print("-" * 65)
    for m in images[:10]:  # Show first 10
        size_str = f"{m['width']}x{m['height']}"
        print(f"{m['image_id']:<40} {size_str:>15} {m['n_nuclei']:>8}")

    if len(images) > 10:
        print(f"... and {len(images) - 10} more images")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Prepare CryoNuSeg dataset")
    parser.add_argument("--input_dir", type=str,
                        default="/home/kyu_insilica_co/qMAP/input_data/cryonuseg/raw",
                        help="Path to raw CryoNuSeg data")
    parser.add_argument("--output_dir", type=str,
                        default="/home/kyu_insilica_co/qMAP/input_data/cryonuseg",
                        help="Path to output directory")
    parser.add_argument("--stats_only", action="store_true",
                        help="Only print statistics of processed data")

    args = parser.parse_args()

    if args.stats_only:
        get_cryonuseg_stats(args.output_dir)
    else:
        process_cryonuseg(args.input_dir, args.output_dir)
        get_cryonuseg_stats(args.output_dir)
