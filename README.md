# Breast Cancer Classification

**PCA, logistic regression, and threshold analysis** on the Wisconsin Breast Cancer dataset.

![R](https://img.shields.io/badge/R-analysis-276DC3?logo=r&logoColor=white) ![PCA](https://img.shields.io/badge/PCA-dimension%20reduction-167D8D) ![Logistic regression](https://img.shields.io/badge/Model-logistic%20regression-BC5A36)

> Compare an interpretable nine-feature logistic model with a six-component PCA model, then examine how classification thresholds change malignant-case sensitivity.

## Findings

| Analysis | Reported result |
| --- | --- |
| Full-feature model, threshold 0.3 | 95.10% accuracy; 95.35% sensitivity; 4 false negatives |
| Full-feature model, threshold 0.5 | 93.63% accuracy; 91.86% sensitivity; 7 false negatives |
| Full-feature model, threshold 0.7 | 91.67% accuracy; 86.05% sensitivity; 12 false negatives |
| Five-fold logistic cross-validation | 96.78% accuracy; 95.40% sensitivity; 97.52% specificity |

The report identifies false negatives as the most concerning error. It reports PC1 at 65.53% of variance and the first six components at 93.29% cumulatively. For the PCA model, the narrative says 94.12% test accuracy; its displayed output says 93.63%. Both values are preserved as they appear in the original report.

## Method

The R analysis cleans the `mlbench::BreastCancer` data, removes incomplete rows, converts nine ordinal features to numeric, and makes a seeded 70/30 stratified split. It explores scaled training-set PCA, fits logistic regression on the original features and on six principal components, compares thresholds, and includes a five-fold cross-validation analysis.

The report's stated rationale is to standardize before PCA so features contribute on a common variance scale and to stratify the split to preserve class proportions.

## Files

- [`analysis.R`](analysis.R): sequential R script containing the report's code.
- [`Breast-Cancer-Analysis.R.html`](Breast-Cancer-Analysis.R.html): original project report.
- [`report.pdf`](report.pdf): concise methods, rationale, and findings report.
- [`report.html`](report.html): HTML source for the PDF report.

## Run

Requires R and the `mlbench`, `caret`, and `ggplot2` packages.

```r
install.packages(c("mlbench", "caret", "ggplot2"))
```

```sh
Rscript analysis.R
```

*Methods, rationale, interpretations, and metrics are summarized from the original HTML report.*
