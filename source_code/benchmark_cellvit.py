#!/usr/bin/env python3
"""
Benchmark CellViT on nuclei segmentation datasets.

CellViT: Vision Transformer for Precise Cell Segmentation and Classification
- Uses Vision Transformer architecture
- Trained on PanNuke dataset
- Expects 1024x1024 patches (will pad smaller images)

Reference:
Hörst et al., CellViT: Vision Transformers for Precise Cell Segmentation and Classification, 2023

Note: CellViT has more complex setup requirements. This script provides
two implementation paths:
1. Using CellViT repository directly (recommended)
2. Using TIA Toolbox if CellViT is integrated

Installation:
    git clone https://github.com/TIO-IKIM/CellViT.git
    cd CellViT
    pip install -e .
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

# Add CellViT to Python path
CELLVIT_PATH = '/home/kyu_insilica_co/qMAP/CellViT'
if CELLVIT_PATH not in sys.path:
    sys.path.insert(0, CELLVIT_PATH)


def check_cellvit_available() -> Tuple[bool, str]:
    """
    Check if CellViT is available for import.

    Returns:
        Tuple of (is_available, installation_method)
    """
    # Try direct CellViT import (our working CPU implementation)
    try:
        from models.segmentation.cell_segmentation.cellvit import CellViTSAM
        from models.segmentation.cell_segmentation.cellvit_shared import CellViTSAMShared
        from utils.tools import unflatten_dict
        import torch
        import albumentations as A
        return True, "cellvit_direct"
    except ImportError as e:
        print(f"CellViT import failed: {e}")
        pass

    return False, "not_found"


def pad_to_size(img: np.ndarray, target_size: int = 1024) -> Tuple[np.ndarray, Tuple[int, int]]:
    """
    Pad image to target size for CellViT.

    CellViT expects 1024x1024 patches. This function pads smaller images
    with reflection padding to minimize edge artifacts.

    Args:
        img: Input image (H, W, C)
        target_size: Target size (default 1024)

    Returns:
        Tuple of (padded_image, original_size)
    """
    h, w = img.shape[:2]
    original_size = (h, w)

    if h >= target_size and w >= target_size:
        return img, original_size

    # Calculate padding
    pad_h = max(0, target_size - h)
    pad_w = max(0, target_size - w)

    # Symmetric padding
    pad_top = pad_h // 2
    pad_bottom = pad_h - pad_top
    pad_left = pad_w // 2
    pad_right = pad_w - pad_left

    # Reflection padding
    if len(img.shape) == 3:
        padded = np.pad(
            img,
            ((pad_top, pad_bottom), (pad_left, pad_right), (0, 0)),
            mode='reflect'
        )
    else:
        padded = np.pad(
            img,
            ((pad_top, pad_bottom), (pad_left, pad_right)),
            mode='reflect'
        )

    return padded, (pad_top, pad_left, h, w)


def unpad_mask(mask: np.ndarray, padding_info: Tuple[int, int, int, int]) -> np.ndarray:
    """
    Remove padding from mask to match original image size.

    Args:
        mask: Padded mask
        padding_info: (pad_top, pad_left, orig_h, orig_w)

    Returns:
        Unpadded mask
    """
    pad_top, pad_left, orig_h, orig_w = padding_info
    return mask[pad_top:pad_top + orig_h, pad_left:pad_left + orig_w]


def run_cellvit_inference_direct(
    image_paths: List[str],
    output_dir: str,
    model_path: str,
    batch_size: int = 1
) -> Dict[str, np.ndarray]:
    """
    Run CellViT inference using direct model loading (CPU-compatible).

    This function uses the CellViT repository structure and works on CPU.

    Args:
        image_paths: List of paths to input images
        output_dir: Directory to save outputs
        model_path: Path to CellViT model checkpoint
        batch_size: Batch size for inference (ignored, always 1 for simplicity)

    Returns:
        Dictionary mapping image names to instance masks
    """
    import torch
    import torch.nn.functional as F
    from PIL import Image
    import albumentations as A
    from models.segmentation.cell_segmentation.cellvit import CellViTSAM
    from models.segmentation.cell_segmentation.cellvit_shared import CellViTSAMShared
    from utils.tools import unflatten_dict

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    print(f"Running CellViT on {len(image_paths)} images...")
    device = "cpu"  # Force CPU
    print(f"Using device: {device}")

    # Load model
    print(f"Loading CellViT model from {model_path}")
    checkpoint = torch.load(model_path, map_location=device)
    
    # Extract config
    config = unflatten_dict(checkpoint['config'], '.')
    arch = checkpoint['arch']
    
    print(f"Model architecture: {arch}")
    print(f"Backbone: {config['model']['backbone']}")
    
    # Create model
    if arch == "CellViTSAM":
        model_class = CellViTSAM
    elif arch == "CellViTSAMShared":
        model_class = CellViTSAMShared
    else:
        raise NotImplementedError(f"Architecture {arch} not supported")
    
    model = model_class(
        model_path=None,
        num_nuclei_classes=config['data']['num_nuclei_classes'],
        num_tissue_classes=config['data']['num_tissue_classes'],
        vit_structure=config['model']['backbone'],
        regression_loss=False
    )
    
    # Load weights
    model.load_state_dict(checkpoint['model_state_dict'])
    model.to(device)
    model.eval()
    print("Model loaded successfully")
    
    # Setup transforms
    transform_settings = config.get('transformations', {})
    if 'normalize' in transform_settings:
        mean = transform_settings['normalize'].get('mean', [0.5, 0.5, 0.5])
        std = transform_settings['normalize'].get('std', [0.5, 0.5, 0.5])
    else:
        mean = [0.5, 0.5, 0.5]
        std = [0.5, 0.5, 0.5]
    
    transforms = A.Compose([A.Normalize(mean=mean, std=std)])

    predictions = {}
    failed = []

    for img_path in tqdm(image_paths, desc="CellViT inference"):
        img_name = Path(img_path).stem

        try:
            # Load image
            img = Image.open(img_path)
            if img.mode != 'RGB':
                img = img.convert('RGB')
            img_array = np.array(img).astype(np.float32)
            original_shape = img_array.shape[:2]

            # Pad to be divisible by 16 (patch size)
            h, w = img_array.shape[:2]
            patch_size = 16
            pad_h = (patch_size - h % patch_size) % patch_size
            pad_w = (patch_size - w % patch_size) % patch_size
            
            if pad_h > 0 or pad_w > 0:
                img_array = np.pad(
                    img_array,
                    ((0, pad_h), (0, pad_w), (0, 0)),
                    mode='reflect'
                )
                padded = True
            else:
                padded = False
            
            # Normalize
            transformed = transforms(image=img_array)
            img_normalized = transformed['image']
            
            # Convert to tensor
            img_tensor = torch.from_numpy(img_normalized).permute(2, 0, 1).unsqueeze(0)
            img_tensor = img_tensor.float().to(device)

            # Run inference
            with torch.no_grad():
                result = model(img_tensor)
            
            # Postprocess
            result['nuclei_binary_map'] = F.softmax(result['nuclei_binary_map'], dim=1)
            result['nuclei_type_map'] = F.softmax(result['nuclei_type_map'], dim=1)
            
            instance_maps, instance_types = model.calculate_instance_map(result, magnification=40)
            inst_map = instance_maps[0].cpu().numpy().astype(np.int32)
            
            # Remove padding
            if padded:
                inst_map = inst_map[:original_shape[0], :original_shape[1]]

            predictions[img_name] = inst_map

            # Save prediction
            np.save(output_path / f"{img_name}.npy", inst_map)

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


def run_cellvit_inference_tiatoolbox(
    image_paths: List[str],
    output_dir: str,
    model_name: str = "cellvit-pannuke",
    batch_size: int = 8
) -> Dict[str, np.ndarray]:
    """
    Run CellViT inference using TIA Toolbox (if available).

    Args:
        image_paths: List of paths to input images
        output_dir: Directory to save outputs
        model_name: Model name in TIA Toolbox
        batch_size: Batch size for inference

    Returns:
        Dictionary mapping image names to instance masks
    """
    try:
        from tiatoolbox.models.engine.nucleus_instance_segmentor import NucleusInstanceSegmentor
        from PIL import Image
        import torch
        import joblib

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        print(f"Running CellViT via TIA Toolbox on {len(image_paths)} images...")
        device = "cuda" if torch.cuda.is_available() else "cpu"
        print(f"Using device: {device}")

        # Initialize segmentor
        segmentor = NucleusInstanceSegmentor(
            pretrained_model=model_name,
            batch_size=batch_size,
            verbose=False
        )

        predictions = {}
        failed = []

        for img_path in tqdm(image_paths, desc="CellViT inference"):
            img_name = Path(img_path).stem

            try:
                # Get image dimensions
                img = Image.open(img_path)
                img_shape = (img.height, img.width)

                # Run prediction
                output = segmentor.predict(
                    imgs=[img_path],
                    save_dir=str(output_path / img_name),
                    mode="tile",
                    device=device,
                    crash_on_exception=True
                )

                # Load result
                result_file = output_path / img_name / "0.dat"
                if result_file.exists():
                    result = joblib.load(result_file)

                    if isinstance(result, dict):
                        # Reconstruct from contours if needed
                        from benchmark_hovernet import reconstruct_instance_map
                        first_key = list(result.keys())[0] if result else None
                        if first_key and isinstance(result.get(first_key), dict):
                            inst_map = reconstruct_instance_map(result, img_shape)
                            predictions[img_name] = inst_map
                        elif "inst_map" in result:
                            predictions[img_name] = result["inst_map"]
                    elif isinstance(result, np.ndarray):
                        predictions[img_name] = result

            except Exception as e:
                failed.append((img_name, str(e)))
                if len(failed) <= 3:
                    print(f"Failed: {img_name} - {e}")

        if failed:
            print(f"\n{len(failed)} images failed")

        print(f"\nSuccessfully processed {len(predictions)}/{len(image_paths)} images")
        return predictions

    except ImportError:
        print("TIA Toolbox CellViT not available")
        return {}


def run_cellvit_inference_huggingface(
    image_paths: List[str],
    output_dir: str,
    model_name: str = "TIO-IKIM/CellViT-SAM-H-x40",
    batch_size: int = 1
) -> Dict[str, np.ndarray]:
    """
    Run CellViT inference using HuggingFace model.

    Args:
        image_paths: List of paths to input images
        output_dir: Directory to save outputs
        model_name: HuggingFace model identifier
        batch_size: Batch size for inference

    Returns:
        Dictionary mapping image names to instance masks
    """
    try:
        import torch
        from PIL import Image
        from transformers import AutoModel, AutoProcessor

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        print(f"Running CellViT via HuggingFace on {len(image_paths)} images...")
        device = "cuda" if torch.cuda.is_available() else "cpu"
        print(f"Using device: {device}")

        # Load model and processor
        # Note: This is a placeholder - actual HF integration may vary
        processor = AutoProcessor.from_pretrained(model_name)
        model = AutoModel.from_pretrained(model_name).to(device)
        model.eval()

        predictions = {}
        failed = []

        for img_path in tqdm(image_paths, desc="CellViT inference"):
            img_name = Path(img_path).stem

            try:
                # Load and preprocess image
                img = Image.open(img_path).convert('RGB')
                img_array = np.array(img)
                original_size = img_array.shape[:2]

                # Pad if needed
                img_padded, padding_info = pad_to_size(img_array, target_size=1024)

                # Process
                inputs = processor(images=Image.fromarray(img_padded), return_tensors="pt")
                inputs = {k: v.to(device) for k, v in inputs.items()}

                with torch.no_grad():
                    outputs = model(**inputs)

                # Extract instance mask
                if hasattr(outputs, 'instance_map'):
                    inst_map = outputs.instance_map.cpu().numpy()
                else:
                    # Fallback: use semantic segmentation as binary mask
                    inst_map = outputs.logits.argmax(dim=1).cpu().numpy()[0]

                # Remove padding
                if padding_info != original_size:
                    inst_map = unpad_mask(inst_map, padding_info)

                predictions[img_name] = inst_map.astype(np.int32)

                # Save
                np.save(output_path / f"{img_name}.npy", inst_map)

            except Exception as e:
                failed.append((img_name, str(e)))
                if len(failed) <= 3:
                    print(f"Failed: {img_name} - {e}")

        if failed:
            print(f"\n{len(failed)} images failed")

        print(f"\nSuccessfully processed {len(predictions)}/{len(image_paths)} images")
        return predictions

    except ImportError as e:
        print(f"HuggingFace transformers not available: {e}")
        return {}


def run_cellvit_inference(
    image_paths: List[str],
    output_dir: str,
    model_path: Optional[str] = None,
    batch_size: int = 1
) -> Dict[str, np.ndarray]:
    """
    Run CellViT inference with automatic backend selection.

    Tries multiple backends in order:
    1. Direct CellViT
    2. TIA Toolbox
    3. HuggingFace

    Args:
        image_paths: List of paths to input images
        output_dir: Directory to save outputs
        model_path: Path to model checkpoint (optional)
        batch_size: Batch size for inference

    Returns:
        Dictionary mapping image names to instance masks
    """
    available, method = check_cellvit_available()

    if method == "cellvit_direct":
        return run_cellvit_inference_direct(image_paths, output_dir, model_path, batch_size)
    elif method == "tiatoolbox":
        return run_cellvit_inference_tiatoolbox(image_paths, output_dir, batch_size=batch_size)
    else:
        print("\nCellViT is not available. Installation options:")
        print("\n1. Install from GitHub:")
        print("   git clone https://github.com/TIO-IKIM/CellViT.git")
        print("   cd CellViT && pip install -e .")
        print("\n2. Download pretrained model from:")
        print("   https://github.com/TIO-IKIM/CellViT/releases")
        print("\n3. Or try HuggingFace (if available):")
        print("   pip install transformers")
        return {}


def benchmark_cellvit_on_dataset(
    dataset_name: str,
    images_dir: str,
    masks_dir: str,
    output_dir: str,
    model_path: Optional[str] = None,
    max_images: Optional[int] = None
) -> Dict:
    """
    Run CellViT benchmark on a dataset.

    Args:
        dataset_name: Name of the dataset
        images_dir: Directory containing images
        masks_dir: Directory containing ground truth masks
        output_dir: Directory for outputs
        model_path: Path to model checkpoint
        max_images: Maximum number of images to process

    Returns:
        Benchmark results dictionary
    """
    print(f"\n{'='*60}")
    print(f"Benchmarking CellViT on {dataset_name}")
    print(f"{'='*60}\n")

    images_path = Path(images_dir)
    masks_path = Path(masks_dir)
    output_path = Path(output_dir) / dataset_name / "cellvit"

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
    predictions = run_cellvit_inference(
        image_paths=[str(f) for f in image_files],
        output_dir=str(output_path / "predictions"),
        model_path=model_path
    )

    if not predictions:
        print("No predictions generated. CellViT may not be installed.")
        return {
            "dataset": dataset_name,
            "model": "cellvit",
            "n_images": 0,
            "summary": {},
            "per_image": [],
            "error": "CellViT not available"
        }

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
            continue

        if img_name not in predictions:
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
            print(f"Shape mismatch for {img_name}")
            continue

        pred_masks.append(pred_mask)
        gt_masks.append(gt_mask)
        image_ids.append(img_name)

    print(f"\nComputing metrics for {len(pred_masks)} images...")

    # Compute metrics
    summary, per_image = compute_dataset_metrics(pred_masks, gt_masks, image_ids)

    # Print results
    print(f"\n{dataset_name} Results (CellViT):")
    print("-" * 50)
    print(f"{'Metric':<10} {'Mean':>10} {'Std':>10} {'Median':>10}")
    print("-" * 50)
    for metric, stats in summary.items():
        print(f"{metric:<10} {stats['mean']:>10.4f} {stats['std']:>10.4f} {stats['median']:>10.4f}")

    # Save results
    results = {
        "dataset": dataset_name,
        "model": "cellvit",
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
    model_path: Optional[str] = None,
    max_images: Optional[int] = None
):
    """Benchmark CellViT on NuInsSeg dataset."""
    print("Benchmarking CellViT on NuInsSeg dataset...")

    base_path = Path(base_dir)
    output_path = Path(output_dir) / "nuinsseg" / "cellvit"
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
    predictions = run_cellvit_inference(
        image_paths=all_images,
        output_dir=str(output_path / "predictions"),
        model_path=model_path
    )

    if not predictions:
        return {"dataset": "NuInsSeg", "model": "cellvit", "error": "Not available"}

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

        if pred_mask.shape != gt_mask.shape:
            continue

        pred_masks.append(pred_mask)
        gt_masks.append(gt_mask)
        valid_ids.append(img_id)

    print(f"\nComputing metrics for {len(pred_masks)} images...")
    summary, per_image = compute_dataset_metrics(pred_masks, gt_masks, valid_ids)

    # Print and save results
    print(f"\nNuInsSeg Results (CellViT):")
    print("-" * 50)
    for metric, stats in summary.items():
        print(f"{metric:<10} {stats['mean']:>10.4f} {stats['std']:>10.4f}")

    results = {
        "dataset": "NuInsSeg",
        "model": "cellvit",
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
    model_path: Optional[str] = None,
    max_images: Optional[int] = None
):
    """Benchmark CellViT on MoNuSeg dataset."""
    return benchmark_cellvit_on_dataset(
        dataset_name="monuseg",
        images_dir=str(Path(base_dir) / "images"),
        masks_dir=str(Path(base_dir) / "masks"),
        output_dir=output_dir,
        model_path=model_path,
        max_images=max_images
    )


def benchmark_cryonuseg(
    base_dir: str = "/home/kyu_insilica_co/qMAP/input_data/cryonuseg",
    output_dir: str = "/home/kyu_insilica_co/qMAP/output/benchmark",
    model_path: Optional[str] = None,
    max_images: Optional[int] = None
):
    """Benchmark CellViT on CryoNuSeg dataset."""
    return benchmark_cellvit_on_dataset(
        dataset_name="cryonuseg",
        images_dir=str(Path(base_dir) / "images"),
        masks_dir=str(Path(base_dir) / "masks"),
        output_dir=output_dir,
        model_path=model_path,
        max_images=max_images
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Benchmark CellViT on nuclei datasets")
    parser.add_argument("--dataset", type=str, default="nuinsseg",
                        choices=["nuinsseg", "monuseg", "cryonuseg", "all"],
                        help="Dataset to benchmark")
    parser.add_argument("--model_path", type=str, default=None,
                        help="Path to CellViT model checkpoint")
    parser.add_argument("--output_dir", type=str,
                        default="/home/kyu_insilica_co/qMAP/output/benchmark",
                        help="Output directory")
    parser.add_argument("--max_images", type=int, default=None,
                        help="Max images to process (for testing)")
    parser.add_argument("--check", action="store_true",
                        help="Check if CellViT is available")

    args = parser.parse_args()

    if args.check:
        available, method = check_cellvit_available()
        if available:
            print(f"CellViT is available via: {method}")
        else:
            print("CellViT is not available")
        sys.exit(0)

    if args.dataset == "nuinsseg" or args.dataset == "all":
        benchmark_nuinsseg(
            output_dir=args.output_dir,
            model_path=args.model_path,
            max_images=args.max_images
        )

    if args.dataset == "monuseg" or args.dataset == "all":
        benchmark_monuseg(
            output_dir=args.output_dir,
            model_path=args.model_path,
            max_images=args.max_images
        )

    if args.dataset == "cryonuseg" or args.dataset == "all":
        benchmark_cryonuseg(
            output_dir=args.output_dir,
            model_path=args.model_path,
            max_images=args.max_images
        )
