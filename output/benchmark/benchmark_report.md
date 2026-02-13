# Nuclei Segmentation Benchmark Report

## Models Evaluated

| Model | Architecture | Training Data | Checkpoint |
|-------|-------------|---------------|------------|
| HoVerNet | hovernet_fast | PanNuke | pannuke |
| StarDist | 2D_versatile_he | MoNuSeg + TNBC | pretrained |
| CellViT | CellViT-SAM-H | PanNuke | SAM ViT-H backbone |

## Datasets

| Dataset | Images | Resolution | Tissue Types |
|---------|--------|------------|--------------|
| NuInsSeg | 665 | Variable | 31 human/mouse organs |
| MoNuSeg | 32 | 1000x1000 | Multiple organs (H&E) |
| CryoNuSeg | 30 | 512x512 | 10 human organs (frozen sections) |

---

## Results Summary

### Dice Coefficient (Mean +/- Std)

| Model | NuInsSeg (n=665) | MoNuSeg (n=32) | CryoNuSeg (n=30) |
|-------|-----------------|---------|-------------------|
| HoVerNet | 0.497 +/- 0.245 | 0.790 +/- 0.101 | 0.777 +/- 0.073 |
| StarDist | 0.452 +/- 0.236 | 0.753 +/- 0.102 | 0.733 +/- 0.065 |
| **CellViT** | **0.659 +/- 0.258** | **0.803 +/- 0.103** | **0.792 +/- 0.073** |

### Aggregated Jaccard Index (AJI) (Mean +/- Std)

| Model | NuInsSeg (n=665) | MoNuSeg (n=32) | CryoNuSeg (n=30) |
|-------|-----------------|---------|-------------------|
| HoVerNet | 0.313 +/- 0.179 | 0.442 +/- 0.131 | 0.524 +/- 0.097 |
| StarDist | 0.284 +/- 0.173 | 0.425 +/- 0.131 | 0.501 +/- 0.088 |
| **CellViT** | **0.458 +/- 0.207** | **0.448 +/- 0.137** | **0.546 +/- 0.103** |

### Panoptic Quality (PQ) (Mean +/- Std)

| Model | NuInsSeg (n=665) | MoNuSeg (n=32) | CryoNuSeg (n=30) |
|-------|-----------------|---------|-------------------|
| HoVerNet | 0.275 +/- 0.162 | 0.410 +/- 0.130 | 0.424 +/- 0.104 |
| StarDist | 0.278 +/- 0.159 | 0.411 +/- 0.132 | 0.416 +/- 0.101 |
| **CellViT** | **0.403 +/- 0.192** | **0.426 +/- 0.134** | **0.453 +/- 0.112** |

### Detection Quality (DQ) (Mean +/- Std)

| Model | NuInsSeg (n=665) | MoNuSeg (n=32) | CryoNuSeg (n=30) |
|-------|-----------------|---------|-------------------|
| HoVerNet | 0.369 +/- 0.205 | 0.563 +/- 0.171 | 0.551 +/- 0.122 |
| StarDist | 0.384 +/- 0.212 | 0.564 +/- 0.173 | 0.547 +/- 0.120 |
| **CellViT** | **0.534 +/- 0.242** | **0.581 +/- 0.176** | **0.590 +/- 0.131** |

### Segmentation Quality (SQ) (Mean +/- Std)

| Model | NuInsSeg (n=665) | MoNuSeg (n=32) | CryoNuSeg (n=30) |
|-------|-----------------|---------|-------------------|
| HoVerNet | 0.668 +/- 0.217 | 0.722 +/- 0.038 | 0.764 +/- 0.024 |
| StarDist | 0.670 +/- 0.185 | 0.722 +/- 0.037 | 0.757 +/- 0.021 |
| **CellViT** | **0.685 +/- 0.211** | **0.727 +/- 0.041** | **0.764 +/- 0.025** |

---

## Compact Comparison (Dice / AJI / PQ)

| Model | Training Data | NuInsSeg (n=665) | MoNuSeg (n=32) | CryoNuSeg (n=30) |
|-------|---------------|----------|---------|-----------|
| HoVerNet | PanNuke | 0.497 / 0.313 / 0.275 | 0.790 / 0.442 / 0.410 | 0.777 / 0.524 / 0.424 |
| StarDist | MoNuSeg+TNBC | 0.452 / 0.284 / 0.278 | 0.753 / 0.425 / 0.411 | 0.733 / 0.501 / 0.416 |
| CellViT-SAM-H | PanNuke | **0.659** / **0.458** / **0.403** | **0.803** / **0.448** / **0.426** | **0.792** / **0.546** / **0.453** |

---

## Key Findings

1. **CellViT-SAM-H achieves the best overall performance** across all three datasets on instance-level metrics (AJI, PQ), with the largest advantage on NuInsSeg (+33% Dice, +46% AJI over HoVerNet).

2. **All models show significant domain shift on NuInsSeg**: Dice scores drop to 0.45-0.66 compared to 0.73-0.79 on MoNuSeg/CryoNuSeg. NuInsSeg's 31 diverse tissue types challenge pretrained models.

3. **CellViT achieves the highest Dice on MoNuSeg** (Dice 0.803) on the 32 TCGA images (1000x1000). HoVerNet (Dice 0.790) and StarDist (Dice 0.753) are close behind. Instance-level metrics (AJI, PQ) are similar across all three models on this dataset.

4. **StarDist shows no in-domain advantage on MoNuSeg** despite being trained on MoNuSeg+TNBC data (Dice 0.753 vs. HoVerNet's 0.790). This may indicate the pretrained model's capacity limitations or differences in the evaluation protocol.

5. **CryoNuSeg (frozen sections)** shows similar trends with CellViT leading, despite frozen tissue presenting different morphology than standard FFPE sections.

6. **Instance-level metrics tell a different story than pixel-level**: While Dice scores are relatively close between models on MoNuSeg/CryoNuSeg, AJI and PQ on NuInsSeg and CryoNuSeg reveal CellViT's superior instance detection and segmentation.

## Notes

- MoNuSeg: All three models evaluated on 32 TCGA images (1000x1000) only. The 50 additional 256x256 test patches were excluded to ensure fair comparison, as the models require different padding strategies for small images (CellViT pads to 1024, HoVerNet to 512, StarDist no padding), which introduces asymmetric artifacts.
- All inference was performed on CPU (no GPU available).
- Metrics: Dice (pixel-level), AJI (instance-level, Kumar et al. 2017), PQ/DQ/SQ (panoptic, Kirillov et al. 2019).
- Updated 2026-02-08: MoNuSeg restricted to 32 TCGA images for fair cross-model comparison (padding fairness audit).
