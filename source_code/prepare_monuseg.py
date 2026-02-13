#!/usr/bin/env python3
"""
Prepare MoNuSeg dataset for nuclei segmentation benchmarking.

MoNuSeg uses XML annotations with polygon vertices. This script:
1. Parses XML annotation files to extract nuclei polygons
2. Converts polygon vertices to instance segmentation masks
3. Outputs standardized PNG images and NPY masks
"""

import os
import glob
import numpy as np
from pathlib import Path
from PIL import Image
import xml.etree.ElementTree as ET
from skimage.draw import polygon
import json
from tqdm import tqdm
from typing import List, Tuple, Dict, Optional


def parse_xml_annotation(xml_path: str) -> List[np.ndarray]:
    """
    Parse MoNuSeg XML annotation file to extract nuclei polygons.

    MoNuSeg XML format:
    <Annotations>
      <Annotation>
        <Regions>
          <Region>
            <Vertices>
              <Vertex X="..." Y="..."/>
              ...
            </Vertices>
          </Region>
          ...
        </Regions>
      </Annotation>
    </Annotations>

    Args:
        xml_path: Path to XML annotation file

    Returns:
        List of polygon coordinates (each as Nx2 numpy array)
    """
    polygons = []

    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()

        # Handle different XML structures
        # Structure 1: Annotations/Annotation/Regions/Region/Vertices/Vertex
        for annotation in root.findall('.//Annotation'):
            for region in annotation.findall('.//Region'):
                vertices = region.find('Vertices')
                if vertices is not None:
                    coords = []
                    for vertex in vertices.findall('Vertex'):
                        x = float(vertex.get('X', 0))
                        y = float(vertex.get('Y', 0))
                        coords.append([x, y])

                    if len(coords) >= 3:  # Valid polygon
                        polygons.append(np.array(coords))

        # Structure 2: Direct Region elements
        if not polygons:
            for region in root.findall('.//Region'):
                vertices = region.find('Vertices')
                if vertices is not None:
                    coords = []
                    for vertex in vertices.findall('Vertex'):
                        x = float(vertex.get('X', 0))
                        y = float(vertex.get('Y', 0))
                        coords.append([x, y])

                    if len(coords) >= 3:
                        polygons.append(np.array(coords))

        # Structure 3: Contours with Points
        if not polygons:
            for contour in root.findall('.//Contour'):
                points = contour.findall('Point')
                if points:
                    coords = []
                    for point in points:
                        x = float(point.get('X', point.get('x', 0)))
                        y = float(point.get('Y', point.get('y', 0)))
                        coords.append([x, y])

                    if len(coords) >= 3:
                        polygons.append(np.array(coords))

    except ET.ParseError as e:
        print(f"XML parse error in {xml_path}: {e}")
    except Exception as e:
        print(f"Error processing {xml_path}: {e}")

    return polygons


def polygons_to_instance_mask(
    polygons: List[np.ndarray],
    image_shape: Tuple[int, int]
) -> np.ndarray:
    """
    Convert list of polygons to instance segmentation mask.

    Args:
        polygons: List of polygon coordinates (each Nx2 array, [x, y] format)
        image_shape: (height, width) of output mask

    Returns:
        Instance segmentation mask (H, W) where each nucleus has unique ID
    """
    instance_mask = np.zeros(image_shape, dtype=np.int32)

    for idx, poly in enumerate(polygons, start=1):
        if len(poly) < 3:
            continue

        # Clip coordinates to image bounds
        poly_clipped = poly.copy()
        poly_clipped[:, 0] = np.clip(poly_clipped[:, 0], 0, image_shape[1] - 1)  # x -> width
        poly_clipped[:, 1] = np.clip(poly_clipped[:, 1], 0, image_shape[0] - 1)  # y -> height

        # polygon() expects (row, col) = (y, x)
        rr, cc = polygon(poly_clipped[:, 1], poly_clipped[:, 0], image_shape)

        if len(rr) > 0:
            instance_mask[rr, cc] = idx

    return instance_mask


def process_monuseg_training(
    input_dir: str,
    output_dir: str
) -> Dict:
    """
    Process MoNuSeg training data.

    Expected input structure (Kaggle format):
    input_dir/
      Tissue Images/
        TCGA-*.png
      Annotations/
        TCGA-*.xml

    Or alternative structure:
    input_dir/
      MoNuSegTrainingData/
        Tissue Images/
        Annotations/

    Args:
        input_dir: Path to raw MoNuSeg data
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

    # Find data directories (handle various structures)
    possible_image_dirs = [
        input_path / "Tissue Images",
        input_path / "Tissue images",
        input_path / "tissue_images",
        input_path / "images",
        input_path / "MoNuSegTrainingData" / "Tissue Images",
        input_path / "MoNuSegTrainingData" / "Tissue images",
        input_path / "Train" / "Images",
        input_path,
    ]

    possible_annotation_dirs = [
        input_path / "Annotations",
        input_path / "annotations",
        input_path / "MoNuSegTrainingData" / "Annotations",
        input_path / "Train" / "Annotations",
        input_path,
    ]

    image_dir = None
    annotation_dir = None

    for d in possible_image_dirs:
        if d.exists() and list(d.glob("*.png")) + list(d.glob("*.tif")):
            image_dir = d
            break

    for d in possible_annotation_dirs:
        if d.exists() and list(d.glob("*.xml")):
            annotation_dir = d
            break

    if image_dir is None:
        print(f"Could not find images directory in {input_path}")
        print("Expected: 'Tissue Images' or 'images' subdirectory")
        return {"images": []}

    if annotation_dir is None:
        print(f"Could not find annotations directory in {input_path}")
        print("Expected: 'Annotations' subdirectory")
        return {"images": []}

    print(f"Image directory: {image_dir}")
    print(f"Annotation directory: {annotation_dir}")

    # Find all images
    image_files = sorted(
        list(image_dir.glob("*.png")) +
        list(image_dir.glob("*.tif")) +
        list(image_dir.glob("*.tiff"))
    )

    print(f"Found {len(image_files)} images")

    manifest = {"images": []}

    for img_file in tqdm(image_files, desc="Processing MoNuSeg"):
        img_name = img_file.stem

        # Find corresponding XML annotation
        xml_file = annotation_dir / f"{img_name}.xml"
        if not xml_file.exists():
            # Try without extension variations
            xml_candidates = list(annotation_dir.glob(f"{img_name}*.xml"))
            if xml_candidates:
                xml_file = xml_candidates[0]
            else:
                print(f"Warning: No annotation for {img_name}")
                continue

        # Load image
        img = Image.open(img_file)
        if img.mode != 'RGB':
            img = img.convert('RGB')
        img_array = np.array(img)
        img_shape = (img_array.shape[0], img_array.shape[1])

        # Parse annotations
        polygons = parse_xml_annotation(str(xml_file))

        if len(polygons) == 0:
            print(f"Warning: No valid polygons in {xml_file.name}")
            continue

        # Convert to instance mask
        instance_mask = polygons_to_instance_mask(polygons, img_shape)

        # Count nuclei
        n_nuclei = len(np.unique(instance_mask)) - 1  # Exclude background

        # Save outputs
        output_image_path = images_dir / f"{img_name}.png"
        output_mask_path = masks_dir / f"{img_name}.npy"

        img.save(output_image_path)
        np.save(output_mask_path, instance_mask)

        manifest["images"].append({
            "image_id": img_name,
            "image_path": str(output_image_path),
            "mask_path": str(output_mask_path),
            "width": img_shape[1],
            "height": img_shape[0],
            "n_nuclei": n_nuclei
        })

    # Save manifest
    manifest_path = output_path / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"\nProcessed {len(manifest['images'])} images")
    print(f"Total nuclei: {sum(m['n_nuclei'] for m in manifest['images'])}")
    print(f"Manifest saved to {manifest_path}")

    return manifest


def process_monuseg_test(
    input_dir: str,
    output_dir: str
) -> Dict:
    """
    Process MoNuSeg test data (if available).

    Args:
        input_dir: Path to raw MoNuSeg test data
        output_dir: Path to output directory

    Returns:
        Manifest dictionary
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)

    # Test data structure may differ
    # Often: MoNuSegTestData/Tissue Images + Annotations

    return process_monuseg_training(input_dir, output_dir)


def get_monuseg_stats(data_dir: str):
    """
    Print statistics about processed MoNuSeg dataset.

    Args:
        data_dir: Path to processed MoNuSeg directory
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
    print("MoNuSeg Dataset Statistics")
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

    parser = argparse.ArgumentParser(description="Prepare MoNuSeg dataset")
    parser.add_argument("--input_dir", type=str,
                        default="/home/kyu_insilica_co/qMAP/input_data/monuseg/raw",
                        help="Path to raw MoNuSeg data")
    parser.add_argument("--output_dir", type=str,
                        default="/home/kyu_insilica_co/qMAP/input_data/monuseg",
                        help="Path to output directory")
    parser.add_argument("--stats_only", action="store_true",
                        help="Only print statistics of processed data")

    args = parser.parse_args()

    if args.stats_only:
        get_monuseg_stats(args.output_dir)
    else:
        process_monuseg_training(args.input_dir, args.output_dir)
        get_monuseg_stats(args.output_dir)
