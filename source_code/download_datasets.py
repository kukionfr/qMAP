#!/usr/bin/env python3
"""
Download MoNuSeg and CryoNuSeg datasets for nuclei segmentation benchmarking.

Supports multiple download sources:
- Kaggle API (requires ~/.kaggle/kaggle.json)
- HuggingFace datasets
- Direct URLs
"""

import os
import sys
import zipfile
import tarfile
import shutil
from pathlib import Path
from typing import Optional
import urllib.request
import subprocess


def check_kaggle_auth() -> bool:
    """Check if Kaggle API credentials are configured."""
    kaggle_json = Path.home() / ".kaggle" / "kaggle.json"
    if kaggle_json.exists():
        # Ensure proper permissions
        os.chmod(kaggle_json, 0o600)
        return True
    return False


def download_with_kaggle(dataset: str, output_dir: Path) -> bool:
    """
    Download dataset using Kaggle API.

    Args:
        dataset: Kaggle dataset identifier (e.g., 'ipateam/nuinsseg')
        output_dir: Directory to save downloaded files

    Returns:
        True if successful, False otherwise
    """
    try:
        import kaggle
        output_dir.mkdir(parents=True, exist_ok=True)
        kaggle.api.dataset_download_files(dataset, path=str(output_dir), unzip=True)
        print(f"Successfully downloaded {dataset} via Kaggle API")
        return True
    except Exception as e:
        print(f"Kaggle download failed: {e}")
        return False


def download_with_huggingface(dataset: str, output_dir: Path) -> bool:
    """
    Download dataset using HuggingFace datasets library.

    Args:
        dataset: HuggingFace dataset identifier
        output_dir: Directory to save downloaded files

    Returns:
        True if successful, False otherwise
    """
    try:
        from datasets import load_dataset
        output_dir.mkdir(parents=True, exist_ok=True)
        ds = load_dataset(dataset, cache_dir=str(output_dir))
        print(f"Successfully downloaded {dataset} via HuggingFace")
        return True
    except Exception as e:
        print(f"HuggingFace download failed: {e}")
        return False


def download_file(url: str, output_path: Path, description: str = "Downloading") -> bool:
    """
    Download a file from URL with progress bar.

    Args:
        url: URL to download from
        output_path: Path to save the file
        description: Description for progress bar

    Returns:
        True if successful, False otherwise
    """
    try:
        from tqdm import tqdm

        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Get file size if available
        response = urllib.request.urlopen(url)
        total_size = int(response.headers.get('content-length', 0))

        with tqdm(total=total_size, unit='B', unit_scale=True, desc=description) as pbar:
            def reporthook(count, block_size, total_size):
                pbar.update(block_size)

            urllib.request.urlretrieve(url, str(output_path), reporthook)

        print(f"Downloaded to {output_path}")
        return True
    except Exception as e:
        print(f"Download failed: {e}")
        return False


def extract_archive(archive_path: Path, output_dir: Path) -> bool:
    """Extract zip or tar archive."""
    try:
        output_dir.mkdir(parents=True, exist_ok=True)

        if archive_path.suffix == '.zip':
            with zipfile.ZipFile(archive_path, 'r') as zf:
                zf.extractall(output_dir)
        elif archive_path.suffix in ['.tar', '.gz', '.tgz']:
            with tarfile.open(archive_path, 'r:*') as tf:
                tf.extractall(output_dir)
        else:
            print(f"Unknown archive format: {archive_path.suffix}")
            return False

        print(f"Extracted to {output_dir}")
        return True
    except Exception as e:
        print(f"Extraction failed: {e}")
        return False


def download_monuseg(output_dir: Path) -> bool:
    """
    Download MoNuSeg dataset.

    MoNuSeg: Multi-organ Nuclei Segmentation
    - 30 training images + 14 test images (1000x1000)
    - XML annotations with polygon vertices

    Sources:
    1. Kaggle: andrewmvd/monuseg-2018-training-data
    2. Grand Challenge: https://monuseg.grand-challenge.org/
    3. HuggingFace: RationAI/MoNuSeg
    """
    print("\n" + "="*60)
    print("Downloading MoNuSeg Dataset")
    print("="*60)

    monuseg_dir = output_dir / "monuseg"

    # Check if already downloaded
    if (monuseg_dir / "images").exists() and len(list((monuseg_dir / "images").glob("*"))) > 0:
        print(f"MoNuSeg already exists at {monuseg_dir}")
        return True

    # Try Kaggle first (training data)
    if check_kaggle_auth():
        print("Trying Kaggle API...")
        kaggle_datasets = [
            "andrewmvd/monuseg-2018-training-data",
            "ipateam/monuseg"
        ]
        for dataset in kaggle_datasets:
            try:
                if download_with_kaggle(dataset, monuseg_dir / "raw"):
                    return True
            except:
                continue

    # Try direct download from Grand Challenge mirrors
    print("Trying direct download...")

    # Alternative: Download from data hosting
    urls = [
        # Training data
        ("https://drive.google.com/uc?export=download&id=1NKkSQ5T0ZNQ8aUhh0a8Dt2YKYCQXIViw",
         "MoNuSegTrainingData.zip"),
        # Test data
        ("https://drive.google.com/uc?export=download&id=1G54vsOdxWY1hG7dzmkeK3r0xz9s-heyQ",
         "MoNuSegTestData.zip"),
    ]

    raw_dir = monuseg_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    print("\nNote: MoNuSeg may need manual download from Grand Challenge:")
    print("  https://monuseg.grand-challenge.org/Data/")
    print("\nAlternatively, use Kaggle:")
    print("  kaggle datasets download -d andrewmvd/monuseg-2018-training-data")

    # Create placeholder directory structure
    (monuseg_dir / "images").mkdir(parents=True, exist_ok=True)
    (monuseg_dir / "masks").mkdir(parents=True, exist_ok=True)
    (monuseg_dir / "annotations").mkdir(parents=True, exist_ok=True)

    print(f"\nCreated directory structure at {monuseg_dir}")
    print("Please download data manually and place in appropriate directories.")

    return False


def download_cryonuseg(output_dir: Path) -> bool:
    """
    Download CryoNuSeg dataset.

    CryoNuSeg: Cryosectioned Nucleus Segmentation
    - 30 image patches (512x512)
    - Instance segmentation masks

    Sources:
    1. Kaggle: ipateam/cryonuseg
    2. Original paper supplementary materials
    """
    print("\n" + "="*60)
    print("Downloading CryoNuSeg Dataset")
    print("="*60)

    cryonuseg_dir = output_dir / "cryonuseg"

    # Check if already downloaded
    if (cryonuseg_dir / "images").exists() and len(list((cryonuseg_dir / "images").glob("*"))) > 0:
        print(f"CryoNuSeg already exists at {cryonuseg_dir}")
        return True

    # Try Kaggle first
    if check_kaggle_auth():
        print("Trying Kaggle API...")
        try:
            if download_with_kaggle("ipateam/cryonuseg", cryonuseg_dir / "raw"):
                return True
        except:
            pass

    # Try alternative sources
    print("Trying direct download...")

    # CryoNuSeg is available on Zenodo and through the authors
    urls = [
        # Zenodo DOI link (if available)
        # GitHub release (if available)
    ]

    raw_dir = cryonuseg_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    print("\nNote: CryoNuSeg may need manual download:")
    print("  Kaggle: kaggle datasets download -d ipateam/cryonuseg")
    print("  Or from the original paper supplementary materials")

    # Create placeholder directory structure
    (cryonuseg_dir / "images").mkdir(parents=True, exist_ok=True)
    (cryonuseg_dir / "masks").mkdir(parents=True, exist_ok=True)

    print(f"\nCreated directory structure at {cryonuseg_dir}")
    print("Please download data manually and place in appropriate directories.")

    return False


def setup_kaggle_auth():
    """Interactive Kaggle authentication setup."""
    kaggle_dir = Path.home() / ".kaggle"
    kaggle_json = kaggle_dir / "kaggle.json"

    if kaggle_json.exists():
        print(f"Kaggle credentials found at {kaggle_json}")
        return True

    print("\n" + "="*60)
    print("Kaggle API Setup")
    print("="*60)
    print("\nTo use Kaggle API, you need to:")
    print("1. Create a Kaggle account at https://www.kaggle.com")
    print("2. Go to Account Settings -> API -> Create New API Token")
    print("3. Download kaggle.json and place it at ~/.kaggle/kaggle.json")
    print("\nOr run:")
    print("  mkdir -p ~/.kaggle")
    print("  echo '{\"username\":\"YOUR_USERNAME\",\"key\":\"YOUR_KEY\"}' > ~/.kaggle/kaggle.json")
    print("  chmod 600 ~/.kaggle/kaggle.json")

    kaggle_dir.mkdir(parents=True, exist_ok=True)

    return False


def download_all(output_dir: str):
    """Download all datasets."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    print("\n" + "="*60)
    print("Nuclei Segmentation Dataset Downloader")
    print("="*60)

    # Check Kaggle auth
    has_kaggle = check_kaggle_auth()
    if not has_kaggle:
        print("\nKaggle API not configured. Attempting alternative download methods.")
        setup_kaggle_auth()
    else:
        print("\nKaggle API credentials found.")

    # Download datasets
    results = {}

    results['monuseg'] = download_monuseg(output_path)
    results['cryonuseg'] = download_cryonuseg(output_path)

    # Summary
    print("\n" + "="*60)
    print("Download Summary")
    print("="*60)
    for dataset, success in results.items():
        status = "SUCCESS" if success else "NEEDS MANUAL DOWNLOAD"
        print(f"  {dataset}: {status}")

    return results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Download nuclei segmentation datasets")
    parser.add_argument("--output_dir", type=str,
                        default="/home/kyu_insilica_co/qMAP/input_data",
                        help="Output directory for datasets")
    parser.add_argument("--dataset", type=str, default="all",
                        choices=["monuseg", "cryonuseg", "all"],
                        help="Dataset to download")
    parser.add_argument("--setup_kaggle", action="store_true",
                        help="Show Kaggle setup instructions")

    args = parser.parse_args()

    if args.setup_kaggle:
        setup_kaggle_auth()
    elif args.dataset == "all":
        download_all(args.output_dir)
    elif args.dataset == "monuseg":
        download_monuseg(Path(args.output_dir))
    elif args.dataset == "cryonuseg":
        download_cryonuseg(Path(args.output_dir))
