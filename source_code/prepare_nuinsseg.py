#!/usr/bin/env python3
"""
Prepare NuInsSeg dataset for HoVerNet benchmarking.

Collects tissue images and label masks, converts 16-bit TIF instance masks
to format suitable for evaluation.
"""

import os
import glob
import numpy as np
from pathlib import Path
from PIL import Image
import tifffile
import json
from tqdm import tqdm


def collect_nuinsseg_data(base_dir: str, output_dir: str):
    """
    Collect and organize NuInsSeg dataset.

    Args:
        base_dir: Path to NuInsSeg dataset root
        output_dir: Path to output directory
    """
    base_path = Path(base_dir)
    output_path = Path(output_dir)

    # Create output directories
    images_dir = output_path / "images"
    masks_dir = output_path / "masks"
    images_dir.mkdir(parents=True, exist_ok=True)
    masks_dir.mkdir(parents=True, exist_ok=True)

    # Find all tissue subdirectories
    tissue_dirs = [d for d in base_path.iterdir() if d.is_dir()]

    data_manifest = []

    for tissue_dir in tqdm(tissue_dirs, desc="Processing tissues"):
        tissue_name = tissue_dir.name

        # Get tissue images
        image_dir = tissue_dir / "tissue images"
        mask_dir = tissue_dir / "label masks"

        if not image_dir.exists() or not mask_dir.exists():
            print(f"Skipping {tissue_name}: missing images or masks")
            continue

        # Process each image
        for image_path in sorted(image_dir.glob("*.png")):
            image_name = image_path.stem
            mask_path = mask_dir / f"{image_name}.tif"

            if not mask_path.exists():
                print(f"Warning: No mask for {image_name}")
                continue

            # Copy image
            output_image_path = images_dir / f"{image_name}.png"
            if not output_image_path.exists():
                img = Image.open(image_path)
                img.save(output_image_path)

            # Load and save mask as numpy array
            output_mask_path = masks_dir / f"{image_name}.npy"
            if not output_mask_path.exists():
                mask = tifffile.imread(mask_path)
                # Instance mask: each unique value is a different nucleus
                np.save(output_mask_path, mask.astype(np.int32))

            data_manifest.append({
                "image_id": image_name,
                "tissue": tissue_name,
                "image_path": str(output_image_path),
                "mask_path": str(output_mask_path)
            })

    # Save manifest
    manifest_path = output_path / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(data_manifest, f, indent=2)

    print(f"\nProcessed {len(data_manifest)} images")
    print(f"Manifest saved to {manifest_path}")

    return data_manifest


def get_tissue_stats(base_dir: str):
    """
    Get statistics about the NuInsSeg dataset.
    """
    base_path = Path(base_dir)
    tissue_dirs = [d for d in base_path.iterdir() if d.is_dir()]

    stats = {}
    total_images = 0

    for tissue_dir in tissue_dirs:
        tissue_name = tissue_dir.name
        image_dir = tissue_dir / "tissue images"

        if image_dir.exists():
            n_images = len(list(image_dir.glob("*.png")))
            stats[tissue_name] = n_images
            total_images += n_images

    print(f"\nNuInsSeg Dataset Statistics:")
    print(f"{'Tissue':<40} {'Images':>8}")
    print("-" * 50)
    for tissue, count in sorted(stats.items()):
        print(f"{tissue:<40} {count:>8}")
    print("-" * 50)
    print(f"{'Total':<40} {total_images:>8}")

    return stats


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Prepare NuInsSeg dataset")
    parser.add_argument("--input_dir", type=str,
                        default="/home/kyu_insilica_co/qMAP/input_data/nuinsseg",
                        help="Path to NuInsSeg dataset")
    parser.add_argument("--output_dir", type=str,
                        default="/home/kyu_insilica_co/qMAP/input_data/nuinsseg_processed",
                        help="Path to output directory")
    parser.add_argument("--stats_only", action="store_true",
                        help="Only print statistics")

    args = parser.parse_args()

    if args.stats_only:
        get_tissue_stats(args.input_dir)
    else:
        get_tissue_stats(args.input_dir)
        collect_nuinsseg_data(args.input_dir, args.output_dir)
