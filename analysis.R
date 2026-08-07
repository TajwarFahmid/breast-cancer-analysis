# Sequential R script containing the analysis code chunks from
# Breast-Cancer-Analysis.R.html. Section order follows the original report.

## Step 1: Load and prepare the Wisconsin Breast Cancer data

#########################
##### Mini Project 4 ####
#### by Tajwar Fahmid ###
#########################

### Part A: Preprocessing the dataset

## Data Preparation

# Load data from mlbench
library(mlbench)

data(BreastCancer, package = "mlbench")
mydata <- BreastCancer

# Rename mlbench columns to match your script's variable names
names(mydata)[names(mydata) == "Id"]              <- "Case"
names(mydata)[names(mydata) == "Cl.thickness"]    <- "Thick"
names(mydata)[names(mydata) == "Cell.size"]       <- "Size"
names(mydata)[names(mydata) == "Cell.shape"]      <- "Shape"
names(mydata)[names(mydata) == "Marg.adhesion"]   <- "Adhesion"
names(mydata)[names(mydata) == "Bare.nuclei"]     <- "Nuclei"
names(mydata)[names(mydata) == "Class"]           <- "Tumor"
names(mydata)[names(mydata) == "Turmor"]          <- "Tumor" # In case already misspelled

# Clean missing values and convert types
mydata[mydata == "?"] <- NA
mydata$Nuclei <- as.numeric(as.character(mydata$Nuclei))
mydata <- na.omit(mydata)

# Convert all ordinal factor features to numeric
cols_to_scale <- c("Thick", "Size", "Nuclei", "Shape", "Adhesion", 
                   "Epith.c.size", "Bl.cromatin", "Normal.nucleoli", "Mitoses")
mydata[cols_to_scale] <- lapply(mydata[cols_to_scale], function(x) as.numeric(as.character(x)))

# Convert target to binary (0/1)
mydata$Tumor <- ifelse(mydata$Tumor == "malignant", 1, 0)

# Why standardization is important:
# Standardization is important before performing PCA because it ensures each variable contributes with unit variance and so no single variable can dominate the principal components.

## Train-Test Split

## Step 2: Create the seeded stratified training and test sets

# 70-30 split
library(caret)

set.seed(7435)
train_indices <- createDataPartition(mydata$Tumor, p = 0.70, list = FALSE)

train_data <- mydata[train_indices, ]
test_data  <- mydata[-train_indices, ]

# Why stratified sampling is important: 
# This method preserves the original target class distribution and is helpful because without it, our model can have extreme sampling bias where an important minority class may not be properly represented in the sample.

## PCA for Training

## Step 3: Explore training-feature principal components
# Correlation matrix
train_features <- train_data[, cols_to_scale]
R_train <- cor(train_features, use = "complete.obs")
round(R_train, 3)

# Eigenvalues and eigenvectors for R_train
decomp <- eigen(R_train)
eigenvalues <- decomp$values
eigenvectors <- decomp$vectors

eigenvalues

eigenvectors

# Verification with prcomp()
pca_out <- prcomp(train_features, center = TRUE, scale. = TRUE)
pca_out

# Table of eigenvalues, proportion of variance explained, and cumulative variance
eigenvalues <- pca_out$sdev^2
prop_var <- eigenvalues / sum(eigenvalues)
cum_var  <- cumsum(prop_var)

pca_variance_table <- data.frame(
  PC = paste0("PC", seq_along(eigenvalues)),
  Eigenvalue = round(eigenvalues, 4),
  Proportion_Variance = round(prop_var, 4),
  Cumulative_Variance = round(cum_var, 4)
)

print(pca_variance_table, row.names = FALSE)

# How many components for at least 90% explainability of variance:
# We need to select at least 6 components to explain at least 90% of the variance in the data.

## Scree Plot and Component Selection
library(ggplot2)

# Build scree dataframe directly from the computed table
scree_df <- data.frame(
  PC = pca_variance_table$PC,
  Eigenvalue = pca_variance_table$Eigenvalue,
  PropVar = pca_variance_table$Proportion_Variance * 100,
  CumVar  = pca_variance_table$Cumulative_Variance * 100
)

# Plot individual and cumulative variance explained
ggplot(scree_df, aes(x = PC)) +
  geom_col(aes(y = PropVar), fill = "#2b5c8f", width = 0.6, alpha = 0.85) +
  geom_line(aes(y = CumVar, group = 1), color = "#d95f02", linewidth = 1) +
  geom_point(aes(y = CumVar), color = "#d95f02", size = 2.5) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "darkred", linewidth = 0.8) +
  annotate("text", x = 3.5, y = 93, label = "90% Variance Threshold (PC6)", color = "darkred", size = 3.5) +
  scale_y_continuous(
    name = "Individual Variance Explained (%)",
    sec.axis = sec_axis(~., name = "Cumulative Variance Explained (%)")
  ) +
  labs(
    title = "Scree Plot: Wisconsin Breast Cancer Features",
    subtitle = "Bars: Individual Variance | Line: Cumulative Variance",
    x = "Principal Component"
  ) +
  theme_minimal(base_size = 12)

# Components based on elbow:
# The scree elbow selects 1 or 2 components, where the variance curve sharply levels off. This contradicts the 90% variance rule, which requires 6 components to retain subtle secondary signals.

## PC Loadings and Interpretations

# Loading vectors for 4 PCs
loadings_4pc <- pca_out$rotation[, 1:4]
round(loadings_4pc, 4)

# Part (b): Plain Language Interpretation
# PC1 (Global Malignancy Index): Captures overall structural deformity and severity across all 9 features; large negative values indicate aggressive, highly abnormal cells.
# PC2 (Proliferation Rate): Dominated almost entirely by Mitoses; measures active cell division independent of static cell shape.
# PC3 (Clump Thickness): Driven primarily by Thick opposing Adhesion; contrasts dense, multilayered cell clusters against flatter tissue sheets.
# PC4 (Cellular Integrity): Contrasts Nuclei against Epith.c.size; separates stripped, bare nuclei from intact, cytoplasm-rich epithelial cells.

# Part (c): Dominant Feature Contributions
# PC1: Size (14.7%) and Shape (14.4%) contribute most, though variance is shared evenly across all features except Mitoses.
# PC2: Mitoses (86.4%) overwhelmingly drives this component.
# PC3: Thick (68.9%) dominates, followed by Adhesion (17.9%).
# PC4: Nuclei (30.1%) and Epith.c.size (26.6%) form the primary opposing weights.

## Visualization
# 1. Combine PCA scores with class labels
pca_plot_df <- data.frame(
  PC1 = pca_out$x[, 1],
  PC2 = pca_out$x[, 2],
  Tumor = factor(train_data$Tumor, levels = c(0, 1), labels = c("Benign", "Malignant"))
)

# 2. Generate the 2D projection scatterplot
ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = Tumor)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("Benign" = "#2b5c8f", "Malignant" = "#d95f02")) +
  labs(
    title = "PCA Projection: PC1 vs. PC2",
    subtitle = "Wisconsin Breast Cancer Training Set",
    x = "PC1 (65.53% Variance Explained)",
    y = "PC2 (8.74% Variance Explained)",
    color = "Tumor Class"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

# Interpretation - The plot shows distinct separation along PC1, where benign cases cluster tightly at positive values while malignant cases disperse widely into negative values. PC2 primarily captures variance within the malignant cohort, where aggressive tumors with high mitotic rates spread vertically.

### Part B: Logistic Regression Without PCA

## Step 4: Fit and evaluate the full-feature logistic model
## Fitting log reg model
# Fit logistic regression on the unscaled training split
glm_full <- glm(
  Tumor ~ Thick + Size + Shape + Adhesion + Epith.c.size + Nuclei + Bl.cromatin + Normal.nucleoli + Mitoses,
  data = train_data,
  family = binomial(link = "logit")
)

summary(glm_full)

# Model Summary, coefficient, p-values
# Extract coefficients, standard errors, z-values, and p-values
coef_table <- summary(glm_full)$coefficients
round(coef_table, 4)

# Compute Odds Ratios (OR) and 95% Confidence Intervals
odds_ratios <- exp(cbind(OR = coef(glm_full), confint.default(glm_full)))
round(odds_ratios, 4)

## Evaluating log reg model
# Confusion matrix
# Predict probabilities on the unseen test set
test_probs <- predict(glm_full, newdata = test_data, type = "response")

# Classify using standard threshold 0.5
test_preds_50 <- ifelse(test_probs >= 0.5, 1, 0)

# Build confusion matrix (Positive Class = Malignant / 1)
conf_50 <- confusionMatrix(
  data = factor(test_preds_50, levels = c(0, 1), labels = c("Benign", "Malignant")),
  reference = factor(test_data$Tumor, levels = c(0, 1), labels = c("Benign", "Malignant")),
  positive = "Malignant"
)

conf_50$table

# Most concerning misclassification for this problem:
# A false negative is the most concerning misclassification for this case because predicting a malignant tumor to be benign will cause treatment to be delayed and witheld, this will allow cancer cells to multiply and metastasize, which directly reduces patient survival.

# Recalculating tables at different thresholds:
# Function to evaluate classification metrics across custom cutoffs
eval_cutoff <- function(probs, actual, cutoff) {
  preds <- ifelse(probs >= cutoff, 1, 0)
  cm <- confusionMatrix(
    data = factor(preds, levels = c(0, 1), labels = c("Benign", "Malignant")),
    reference = factor(actual, levels = c(0, 1), labels = c("Benign", "Malignant")),
    positive = "Malignant"
  )
  data.frame(
    Cutoff      = cutoff,
    TP          = cm$table["Malignant", "Malignant"],
    FP          = cm$table["Malignant", "Benign"],
    TN          = cm$table["Benign", "Benign"],
    FN          = cm$table["Benign", "Malignant"],
    Sensitivity = round(unname(cm$byClass["Sensitivity"]), 4),
    Specificity = round(unname(cm$byClass["Specificity"]), 4),
    Accuracy    = round(unname(cm$overall["Accuracy"]), 4)
  )
}

# Tabulate across 0.3, 0.5, and 0.7 cutoffs
threshold_table <- do.call(rbind, lapply(c(0.3, 0.5, 0.7), function(c) {
  eval_cutoff(test_probs, test_data$Tumor, c)
}))

print(threshold_table, row.names = FALSE)

# The logistic regression model achieves solid diagnostic performance across all cutoffs, peaking at 95.10% accuracy and 95.35% sensitivity at the 0.3 threshold. 
# As the cutoff rises toward 0.7, false negatives jump from 4 to 12, trading away vital cancer detection capability for only a 0.84% gain in specificity. 
# Consequently, the 0.3 cutoff is optimal for clinical screening, providing the lowest missed-cancer rate alongside the highest overall test accuracy.

### Part C: Logistic Regression With PCA

## Step 5: Fit and evaluate the six-component logistic model
## Model
# 1. Extract the 6 selected PC training scores and pair with Tumor target
train_pca_data <- data.frame(
  Tumor = train_data$Tumor,
  pca_out$x[, 1:6]
)

# 2. Fit Logistic Regression using the first 6 Principal Components
glm_pca <- glm(
  Tumor ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6,
  data = train_pca_data,
  family = binomial(link = "logit")
)

summary(glm_pca)

## Comparing models and confusion matrix
# Extract coefficient table rounded to 4 decimal places
pca_coef_table <- summary(glm_pca)$coefficients
round(pca_coef_table, 4)

# Odds Ratios and 95% Confidence Intervals
pca_or <- exp(cbind(OR = coef(glm_pca), confint.default(glm_pca)))
round(pca_or, 4)

# Predict on test set using the first 6 components
test_features <- test_data[, cols_to_scale]
test_pca_scores <- predict(pca_out, newdata = test_features)[, 1:6]
pca_test_probs <- predict(glm_pca, newdata = as.data.frame(test_pca_scores), type = "response")
pca_test_preds <- ifelse(pca_test_probs >= 0.5, 1, 0)

conf_pca <- confusionMatrix(
  data = factor(pca_test_preds, levels = c(0, 1), labels = c("Benign", "Malignant")),
  reference = factor(test_data$Tumor, levels = c(0, 1), labels = c("Benign", "Malignant")),
  positive = "Malignant"
)

conf_pca$table

round(conf_pca$overall["Accuracy"], 4)

round(conf_pca$byClass[c("Sensitivity", "Specificity", "Precision")], 4)

## Better model?
# While the 6-component PCA logistic model achieves slightly higher test accuracy (94.12% vs. 93.63%) by resolving severe feature collinearity, the original logistic model should be preferred in a medical setting because its coefficients map directly to tangible cytological markersâ€”such as bare nuclei and clump thicknessâ€”which pathologists and tumor boards can transparently interpret, explain to patients, and align with established clinical grading protocols without sacrificing diagnostic power.

### Bonus

## Step 6: Run the report's five-fold cross-validation analysis
# Ensure target is a properly formatted factor for caret classification
cv_data <- mydata[, c("Tumor", cols_to_scale)]
cv_data$Tumor <- factor(cv_data$Tumor, levels = c(0, 1), labels = c("Benign", "Malignant"))

# Set up stratified 5-fold cross-validation
set.seed(7435)
train_control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary, # Computes ROC, Sens, Spec
  savePredictions = "final"
)

# Fit Logistic Regression via 5-Fold CV
cv_model <- train(
  Tumor ~ .,
  data = cv_data,
  method = "glm",
  family = binomial(link = "logit"),
  trControl = train_control,
  metric = "ROC"
)

# Extract pooled out-of-fold confusion matrix
cv_conf_matrix <- confusionMatrix(
  data = cv_model$pred$pred,
  reference = cv_model$pred$obs,
  positive = "Malignant"
)

print(cv_conf_matrix)

# The 5-fold cross-validated logistic model demonstrates strictly superior and more reliable performance than the single 70/30 test split across every evaluation metric, raising accuracy from 93.63% to 96.78%, sensitivity from 91.86% to 95.40%, and specificity from 94.92% to 97.52%.
# By evaluating out-of-fold predictions across all 683 observations rather than relying on a single 204-sample partition, cross-validation eliminates split-dependent variance, recovers borderline cases that were disproportionately concentrated in the test holdout, and provides a robust, statistically sound estimate of the model's true clinical generalization ability.

