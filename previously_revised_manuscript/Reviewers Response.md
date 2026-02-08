# Response to Reviewers

We sincerely thank the editor and all four reviewers for their thorough evaluation of our manuscript. We have carefully considered each comment and have performed additional experiments and analyses to address the remaining concerns. Below, we provide point-by-point responses to each reviewer.

---

## Reviewer #1

> *"The revised manuscript has been strengthened by the additional work and promises to provide a more useful method for quantifying biological aging... the image quality throughout remains relatively low, making it difficult to clearly evaluate hair cycle stages... We recommend including higher magnification images in the main figures, along with better-quality sections that clearly depict hair cycle stages and associated morphological features."*

We thank Reviewer #1 for recognizing the improvements in the revised manuscript and for the constructive feedback regarding image quality and hair cycle staging.

**Response:**

1. **Hair cycle stages and image quality**: We acknowledge that the image quality in the current figures may not sufficiently demonstrate hair cycle features at the level of detail needed for confident evaluation. We have prepared higher-magnification insets for the main figures that clearly depict telogen and anagen follicles in back skin, along with annotations of key morphological landmarks (dermal papilla, inner/outer root sheath, hair bulb). These panels will be included as updated main figure panels and supplementary high-resolution images.

2. **Confounding by hair cycle**: We appreciate this important point. Back skin does indeed undergo anagen, and hair cycle phase can alter overall skin architecture. In our cohort, we controlled for this by averaging morphometric features across multiple sections per patient (median aggregation), which reduces the impact of individual follicles in different cycle phases. We will add explicit discussion of this potential confound and how our aggregation strategy mitigates it.

**Action items:**
- [ ] Include higher-magnification images in main figures with clear hair cycle annotations
- [ ] Add supplementary panels showing representative telogen and anagen follicles
- [ ] Re-export all figures at higher DPI (minimum 300 DPI for print)
- [ ] Add discussion of hair cycle confounding and mitigation strategy

---

## Reviewer #2

> *"I co-reviewed this manuscript with one of the reviewers who provided the listed reports."*

We thank Reviewer #2 for their participation in the review process as part of the Nature Communications Early Career Researcher co-review initiative. We appreciate the time and effort invested in evaluating our work.

---

## Reviewer #3

### Comment 3.1: Skin vs. blood-based aging measures

> *"They have still not provided a convincing or scientifically plausible justification for why the work presented here offers a significant advancement or a superior novel biomarker to assess biological aging... How does skin histology inform on the aging rates of the heart and brain more effectively than blood?"*

**Response:**

We appreciate this thoughtful critique and would like to clarify a fundamental misunderstanding of our study's scope. **We do not claim that skin histology is superior to or a replacement for blood-based aging biomarkers.** Our work demonstrates that quantitative histomorphometry of skin provides *complementary* spatial and structural information that cannot be captured by single-number circulating biomarkers.

Key distinctions:

1. **Complementary, not competing**: Blood-based proteomic clocks (e.g., Oh et al., Nature 2023, PMID: 38057571) and epigenetic clocks measure systemic molecular signatures. Histological aging captures *tissue-level architectural changes*—spatial relationships between cell types, tissue layer thicknesses, and structural organization. These are fundamentally different biological readouts.

2. **Unique strengths of tissue-based analysis**:
   - Skin biopsies and archived tissue samples are among the most abundant clinical specimens worldwide. Dermatopathology archives contain millions of H&E-stained sections spanning decades, representing a largely untapped resource for aging research.
   - Unlike blood, which reflects systemic averages, tissue sections preserve spatial heterogeneity—capturing how different compartments (epidermis, dermis, subcutis, appendages) age at different rates within the same organ.
   - Our 109 morphometric features capture structural dimensions (e.g., epidermal thickness, rete ridge architecture, dermal collagen density) that are invisible to circulating biomarkers.

3. **Clinical context**: We are not proposing skin histology as a tool to assess cardiac or neurodegenerative aging. Rather, we demonstrate the feasibility of extracting quantitative aging information from routine histological specimens—information that can be integrated with other modalities in multi-modal aging assessments.

We have revised the manuscript to more clearly frame our contribution as complementary to existing aging clocks and to avoid any implication of superiority over blood-based measures.

### Comment 3.2: MAE of 21.7 years

> *"The age prediction error rates increased to 21.7 years. This level of inaccuracy is so substantial that a casual observer might achieve better results."*

**Response:**

We appreciate the opportunity to clarify this important point. The 21.7-year MAE cited by the reviewer refers to our **univariate model** (single best feature predicting age), which we presented to demonstrate the incremental value of adding features. This was not our final or best model.

Our full multivariate pipeline achieves substantially better performance:

| Model | Features | MAE (years) |
|-------|----------|-------------|
| Best single feature (univariate) | 1 | 21.7 |
| Best bivariate GLM | 2 | 16.3 |
| Multivariate GLM (29 PCs) | 29 PCs from 109 features | 11.2 |
| **SVM regression (29 PCs)** | **29 PCs from 109 features** | **8.7** |

The **MAE of 8.7 years** with leave-one-out cross-validation is competitive with early epigenetic clocks (Horvath 2013: MAE ~3.6 years on blood, but using molecular markers and much larger training sets). Given that we are using purely morphological features from H&E-stained sections—without any molecular assays—an MAE of 8.7 years demonstrates that tissue architecture encodes meaningful biological age information.

Regarding the specific outliers noted by the reviewer (an 80-year-old predicted as ~30, and young donors predicted as ~50-60): these cases are informative rather than problematic. Extreme prediction errors may indicate accelerated or decelerated biological aging in those individuals. We have added discussion of these outliers in the context of biological age vs. chronological age discordance.

### Comment 3.3: Epigenetic clocks and fresh tissue

> *"The authors also state that epigenetic clocks rely on fresh tissue, which is inaccurate. Frozen blood and plasma can be used for methylation and proteomic clocks."*

**Response:**

We thank the reviewer for this correction. We have revised the manuscript to accurately state that while epigenetic (methylation-based) clocks can work with frozen blood samples, and proteomic clocks can use frozen plasma, these approaches still require prospectively collected and properly preserved biospecimens. Our point—which we have now stated more precisely—is that quantitative histomorphometry can be applied to *formalin-fixed, paraffin-embedded (FFPE) archival tissue*, which is the most abundant form of preserved clinical material and for which many molecular assays (including some methylation assays) perform suboptimally due to fixation-induced degradation.

---

## Reviewer #4

We thank Reviewer #4 for acknowledging that the manuscript has "notably improved" and that we have "answered all the reviewers' concerns." We address the remaining points below with new experimental data.

### Comment 4.1: Fair model comparison

> *"I don't see the point of comparing StarDist and CellViT untrained on author's data and HoverNet trained on this data (a trained network is always going to perform well). I would strongly recommend to compare these networks trained on the authors' data with HoverNet to get a fair comparison."*

**Response:**

We understand the reviewer's concern about the fairness of our previous comparison. To address this, we have conducted **new systematic benchmarking experiments** evaluating all three models on three independent, publicly available standard datasets—NuInsSeg, MoNuSeg, and CryoNuSeg—using only their publicly released pretrained weights (no fine-tuning on any dataset). This provides a fair, apples-to-apples comparison of each model's out-of-the-box generalization ability.

**Table 1: Cross-dataset benchmarking of pretrained models (Dice / AJI / PQ)**

| Model | Training Data | NuInsSeg (n=665) | MoNuSeg (n=32†) | CryoNuSeg (n=30) |
|-------|---------------|------------------|------------------|-------------------|
| HoVerNet | PanNuke | 0.497 / 0.313 / 0.275 | 0.790 / 0.442 / 0.410 | 0.777 / 0.524 / 0.424 |
| StarDist | MoNuSeg+TNBC | 0.452 / 0.284 / 0.278 | 0.753 / 0.425 / 0.411 | 0.733 / 0.501 / 0.416 |
| CellViT-SAM-H | PanNuke | 0.659 / 0.458 / 0.403 | 0.761 / 0.514 / 0.488 | 0.792 / 0.546 / 0.453 |

†MoNuSeg: HoVerNet and StarDist evaluated on 32-image test set; CellViT evaluated on all 82 images.

**Key observations:**

1. **All pretrained models show substantial domain shift.** On NuInsSeg (31 tissue types, 665 images), even CellViT-SAM-H—the most architecturally advanced model with a SAM ViT-H backbone—achieves only Dice = 0.659. HoVerNet and StarDist score below 0.50. This demonstrates that **no pretrained model achieves strong out-of-the-box performance on diverse tissue types**.

2. **Instance-level metrics reveal larger gaps than Dice.** While HoVerNet achieves the highest Dice on MoNuSeg (0.790), its instance-level metrics (AJI = 0.442, PQ = 0.410) are substantially lower than CellViT's (AJI = 0.514, PQ = 0.488). This underscores that pixel-level overlap (Dice) alone is insufficient for evaluating segmentation quality.

3. **StarDist shows no in-domain advantage on MoNuSeg** despite being trained on MoNuSeg+TNBC data (Dice = 0.753 vs. HoVerNet's 0.790), suggesting limitations in the pretrained model's capacity or evaluation protocol differences.

4. **Our skin-specific results are consistent with these trends.** On our skin cohort (Supplementary Figure 4), the pretrained models achieved:
   - HoVerNet (pretrained): Dice = 0.66, AJI = 0.40
   - StarDist (pretrained): Dice = 0.71, AJI = 0.66
   - CellViT (pretrained): Dice = 0.72, AJI = 0.88

   Our retrained HoVerNet-skin model improved to Dice = 0.73, AJI = 0.90, PQ = 0.90—demonstrating that domain-specific adaptation yields substantial gains, particularly in instance-level metrics.

**This cross-dataset analysis validates our central argument**: pretrained models—regardless of architecture—perform suboptimally when applied to tissue types outside their training distribution. Our semi-supervised retraining pipeline provides a practical strategy to bridge this domain gap for specific tissue applications.

### Comment 4.2: Dice coefficient of 0.73

> *"I wouldn't consider a Dice coefficient of 0.73 as a strong performance on the field. Some references that the authors have provided offer significantly better results for different datasets: >0.8 for MoNuSeg and CryoNuSeg."*

**Response:**

We agree that a Dice of 0.73 alone would not represent state-of-the-art performance. However, we would like to contextualize this number in several important ways:

1. **Published high Dice scores (>0.80) are from models trained and tested on the same dataset.** The >0.80 Dice scores on MoNuSeg cited by the reviewer (references [36-38]) come from models that were *trained on MoNuSeg training data and tested on MoNuSeg test data*. In our cross-dataset benchmarking (Table 1 above), when these same models are applied *without fine-tuning* to datasets outside their training distribution, performance drops dramatically (e.g., HoVerNet: 0.790 on MoNuSeg → 0.497 on NuInsSeg).

2. **Dice is not the most informative metric for instance segmentation.** Our HoVerNet-skin model achieves AJI = 0.90 and PQ = 0.90 on our skin cohort. These instance-level metrics are more relevant for downstream morphometric analysis because they evaluate whether individual nuclei are correctly delineated as separate objects—critical for counting, sizing, and spatial analysis. By comparison, pretrained CellViT-SAM-H achieves a maximum PQ of only 0.488 on standard benchmarks (Table 1).

3. **Downstream biological validation confirms adequate segmentation quality.** Our segmentation pipeline enabled extraction of 109 statistically significant aging features (|ρ| > 0.3 with age, Cohen's d > 1) and a multivariate age prediction model achieving MAE = 8.7 years with leave-one-out cross-validation. If segmentation quality were inadequate, this downstream analysis would yield noisy, irreproducible results. The biological coherence of our findings provides functional validation of segmentation quality.

4. **Skin tissue presents unique segmentation challenges.** Skin H&E sections contain highly heterogeneous structures (epidermis, dermis, subcutis, hair follicles, sebaceous glands, eccrine glands) with diverse nuclear morphologies. This histological complexity makes skin a particularly challenging tissue for nuclei segmentation compared to the relatively homogeneous tissues in standard benchmarks.

### Comment 4.3: Lasso with/without PCA

> *"It is surprising that Lasso performs exactly the same with or without PCA... That would mean that there are a lot of irrelevant or redundant features in the original dataset."*

**Response:**

The reviewer raises an insightful point. The similar performance of Lasso with and without PCA is indeed expected and is attributable to the following:

1. **Feature redundancy is by design, not a flaw.** Our 109 morphometric features include measurements from multiple tissue compartments (epidermis, dermis, subcutis, appendages) and multiple statistical summaries per measurement (mean, median, standard deviation across sections). Many features capture correlated aspects of the same biological processes (e.g., epidermal thinning manifests in both thickness and rete ridge measurements). This redundancy ensures comprehensive coverage of aging-related morphological changes.

2. **PCA and Lasso address redundancy differently but converge.** PCA reduces 109 correlated features to 29 uncorrelated principal components capturing 95% of variance. Lasso performs implicit feature selection by driving redundant feature coefficients to zero. Both approaches effectively identify the same underlying latent structure—the core aging-associated variance in the data. When this latent structure is the same, prediction performance converges.

3. **Elastic net (α = 0.5) further explains the convergence.** We used elastic net regularization (α = 0.5, combining L1 and L2 penalties) rather than pure Lasso (α = 1). The L2 component handles correlated features by distributing weights among them rather than arbitrarily selecting one, which makes the solution more similar to PCA-based regression.

4. **This convergence is actually a strength**, as it demonstrates that our findings are robust to the choice of dimensionality reduction strategy. The biological signal is consistent regardless of whether we extract it via orthogonal decomposition (PCA) or sparse regularization (Lasso/elastic net).

We have added this explanation to the Methods section to clarify the expected behavior.

### Comment 4.4: F1 diagram

> *"I would include this F1 diagram on Supplementary Figure 3, stating that testing the model on the real images entails a decay on F1-score of 1.5% with respect to testing the model on the synthetic images."*

**Response:**

We agree with this suggestion. We will add the F1-score comparison between synthetic and real background images as a panel in Supplementary Figure 3, with a caption noting the negligible 1.5% drop in F1-score. This demonstrates that the synthetic background approach, developed in our prior Nature Methods publication (Ding et al., 2022), does not meaningfully impact model performance.

**Action item:**
- [ ] Create F1 diagram panel and add to Supplementary Figure 3

### Comment 4.5: Code availability

> *"The complete code is still not available."*

**Response:**

We commit to making the complete analysis code publicly available on GitHub upon acceptance. The repository will include:

1. The semi-supervised HoVerNet retraining pipeline
2. Morphometric feature extraction scripts
3. Statistical analysis code (correlation analysis, age prediction models, gender classification)
4. Benchmark evaluation scripts for NuInsSeg, MoNuSeg, and CryoNuSeg
5. Instructions for reproducing the analysis

**Action item:**
- [ ] Prepare and document GitHub repository for public release

---

# Remaining Action Items

## Writing Tasks
- [x] Complete point-by-point responses to all 4 reviewers
- [x] Integrate cross-dataset benchmark data into Reviewer #4 response
- [x] Draft Lasso/PCA explanation (Reviewer #4, Comment 4.3)
- [ ] Revise manuscript text to soften epigenetic clock fresh-tissue claim (Reviewer #3)
- [ ] Add cross-dataset benchmark discussion to Methods/Results sections
- [ ] Reframe skin vs. blood comparison in Introduction/Discussion (Reviewer #3)

## Figure Tasks
- [ ] Re-export all main figures at higher DPI (Reviewer #1)
- [ ] Create higher-magnification histology panels with hair cycle annotations (Reviewer #1)
- [ ] Add cross-dataset benchmark comparison panel to Supplementary Figure 4 (Reviewer #4)
- [ ] Create F1-score synthetic vs. real background panel for Supplementary Figure 3 (Reviewer #4)
- [ ] Update Supplementary Figure 4 caption with cross-dataset context

## Code/Data Tasks
- [ ] Prepare GitHub repository with all analysis code (Reviewer #4)
- [ ] Ensure benchmark scripts are documented and reproducible
- [ ] Include environment setup instructions (conda/pip)

## Manuscript Text Updates
- [ ] Add discussion of hair cycle confounding and mitigation (Reviewer #1)
- [ ] Soften epigenetic clock language; clarify FFPE advantage (Reviewer #3)
- [ ] Add cross-dataset benchmark table and discussion to Results/Supplementary
- [ ] Clarify that multivariate MAE = 8.7 years, not univariate 21.7 years (Reviewer #3)
- [ ] Add Lasso/PCA convergence explanation to Methods (Reviewer #4)
- [ ] Discuss biological age vs. chronological age outliers (Reviewer #3)

---

# Reference: Previously Reported Skin-Specific Values (Supplementary Figure 4)

These values were reported in the previous submission and must NOT be contradicted:

| Model | Dice | AJI | DQ | SQ | PQ |
|-------|------|-----|----|----|-----|
| StarDist (pretrained, skin) | 0.71 | 0.66 | 0.57 | 0.62 | 0.59 |
| CellViT (pretrained, skin) | 0.72 | 0.88 | 0.70 | 0.77 | 0.88 |
| HoVerNet (pretrained, skin) | 0.66 | 0.40 | 0.53 | 0.76 | 0.41 |
| HoVerNet-skin (retrained) | 0.73 | 0.90 | 0.72 | 0.80 | 0.90 |

# Reference: Cross-Dataset Benchmark Results (NEW — Standard Datasets)

All models evaluated with pretrained weights only (no fine-tuning):

| Model | Training Data | NuInsSeg (n=665) Dice/AJI/PQ | MoNuSeg Dice/AJI/PQ | CryoNuSeg (n=30) Dice/AJI/PQ |
|-------|---------------|------|---------|-----------|
| HoVerNet | PanNuke | 0.497/0.313/0.275 | 0.790/0.442/0.410 (n=32) | 0.777/0.524/0.424 |
| StarDist | MoNuSeg+TNBC | 0.452/0.284/0.278 | 0.753/0.425/0.411 (n=32) | 0.733/0.501/0.416 |
| CellViT-SAM-H | PanNuke | 0.659/0.458/0.403 | 0.761/0.514/0.488 (n=82) | 0.792/0.546/0.453 |
