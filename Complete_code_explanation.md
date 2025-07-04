## 📌 Complete_data.R Explanation

### 📦 1. **Load Required Libraries and Source Scripts**

```r
library(survival)
library("MASS")
source("run_data_col.R")
source("specClust.R")
source("enrichment.R")
```

* `survival`: For survival analysis (Kaplan-Meier, log-rank test).
* `MASS`: Contains multivariate functions (used internally, e.g., `ginv`).
* `run_data_col.R`: Loads the function `data_col` for latent space extraction.
* `specClust.R`: Loads the custom spectral clustering function.
* `enrichment.R`: Performs clinical enrichment analysis on clusters.

---

### 🧩 2. **Namespace Reassignment Hack**

```r
assignInNamespace('specClust', specClust, ns='kknn')
environment(specClust) <- asNamespace('kknn')
```

* Injects the local `specClust` function into the `kknn` package namespace.
* This is likely a workaround to avoid namespace conflicts or reuse internal helpers from `kknn`.

---

### 📁 3. **Set File Paths**

```r
file1 <- "BIC"    # The folder name (e.g., Breast Invasive Carcinoma)
file2 <- "BREAST" # The actual cancer type prefix used in file names
```

---

### 📂 4. **Read Input Files (Gene, Methylation, miRNA)**

```r
file1_name <- paste0("data/", file1, "/", file2, "_Gene_Expression.txt")
file2_name <- paste0("data/", file1, "/", file2, "_Methy_Expression.txt")
file3_name <- paste0("data/", file1, "/", file2, "_Mirna_Expression.txt")
```

Creates paths like:
`data/BIC/BREAST_Gene_Expression.txt`
`data/BIC/BREAST_Methy_Expression.txt`
`data/BIC/BREAST_Mirna_Expression.txt`

```r
v1 <- read.table(file1_name, header = T, sep = "\t", row.names = 1)
v2 <- read.table(file2_name, header = T, sep = "\t", row.names = 1)
v3 <- read.table(file3_name, header = T, sep = "\t", row.names = 1)
```

Reads the gene expression, methylation, and miRNA data as matrices.
Each row = gene/feature, each column = sample (with matching names).

---

### 🧬 5. **Combine Data into List**

```r
data <- list(mRNA = v1, Methy = v2, miRNA = v3, clinical = NULL)
```

Combines the three omics views into a single list structure.

---

### 🔍 6. **Latent Subspace Learning**

```r
simul_result <- data_col(
    data,
    incomplete_data = F,
    incomplete_sample_name,
    remain_view = 1,
    dim1 = 13,
    dim2 = 4,
    pca_scale = F
)
emb <- simul_result$l_space
```

This is the **core of MCLS**:

* `data_col()`:

  * Applies PCA and SVD on the full sample set.
  * Learns a **latent subspace** of dimensionality `dim1 = 13`.
  * `remain_view = 1` means the first view (mRNA) is retained fully.
  * Returns a common low-dimensional embedding for all samples.
  * `emb`: The latent space matrix (samples × latent features).

---

### 🧭 7. **Spectral Clustering**

```r
cl_s <- specClust(emb, 5, nn = 20)
```

* Clusters the latent representations using **spectral clustering**.
* `5`: Number of clusters (subtypes).
* `nn = 20`: Builds similarity graph using 20 nearest neighbors.

---

### 🧪 8. **Survival Analysis**

```r
labels <- cl_s$cluster
names(labels) <- names(v1)

survfile <- paste0("data/", file1, "/", file2, "_Survival.txt")
surv <- read.table(survfile, header = T, sep = "\t")
survresult <- survdiff(Surv(Survival, Death) ~ labels, data = surv)
p.val <- 1 - pchisq(survresult$chisq, length(survresult$n) - 1)
```

* Retrieves **cluster labels** and names them according to sample IDs.
* Loads survival data: columns like `Survival` (time), `Death` (0/1).
* Performs **log-rank test** (via `survdiff`) to see if **survival curves differ across clusters**.
* Computes the **p-value** from the Chi-square distribution.

---

### 🩺 9. **Clinical Enrichment Analysis**

```r
para_num <- enrich(folder = file1, sign = T, label = labels)
```

* Runs clinical association tests using `enrich.R`.
* Looks for significant relationships between **cluster labels** and **clinical features** (e.g., stage, subtype).
* `para_num`: Number of clinical parameters significantly enriched across clusters.

---

### 📤 10. **Print Results**

```r
print(p.val)
print(para_num)
```

* Prints the survival analysis **p-value** (significance of subtype survival separation).
* Prints number of clinically **enriched parameters**.

---

## ✅ Summary of What the Code Does

| Step | Description                                                          |
| ---- | -------------------------------------------------------------------- |
| 1    | Loads gene, methylation, and miRNA data for one cancer type (BREAST) |
| 2    | Learns latent subspace from complete multi-omics data                |
| 3    | Applies spectral clustering to group patients into subtypes          |
| 4    | Tests if subtypes differ significantly in survival (Kaplan-Meier)    |
| 5    | Tests if clinical features are significantly enriched in any subtype |

---


## 📌 InComplete_data.R Explanation

### 📦 1. **Load Libraries and Source Scripts**

```r
library(survival)
library("MASS")
source("run_data_col.R")
source("specClust.R")
source("enrichment.R")
```

* **`survival`**: For survival curve analysis (e.g., Kaplan-Meier, log-rank test).
* **`MASS`**: For matrix operations (like `ginv()` – generalized inverse).
* The `source(...)` commands load custom MCLS functions:

  * `run_data_col.R` (not used directly here, but contains helper functions)
  * `specClust.R`: Spectral clustering function.
  * `enrichment.R`: Clinical enrichment analysis.

---

### ⚠️ 2. **Namespace Hack**

```r
assignInNamespace('specClust', specClust, ns='kknn')
environment(specClust) <- asNamespace('kknn')
```

* This manually injects your `specClust` function into the **kknn** namespace, likely to fix scoping issues when `specClust` calls other package functions.

---

### 🔧 3. **Define `get_emb()` for Latent Subspace Projection**

```r
get_emb <- function(data, complete_names, incomplete_names, all_names, dim=6, dim_pca=3, pca_scale=F)
```

#### 🔬 Purpose:

Projects incomplete samples into the **latent subspace** learned from complete samples using **PCA and SVD**.

#### ✅ Main Steps Inside:

1. **PCA on Each Omics**:

   ```r
   mRNA_pca <- prcomp(t(data[["mRNA"]]), ...)
   ```

   * Transpose to samples × features.
   * Extract `dim_pca` components per omics (15 default).

2. **Separate Complete vs. Incomplete Samples**:

   ```r
   complete_mRNA_pca <- mRNA_pca$x[complete_names, 1:dim_pca]
   incomplete_Methy_pca <- Methy_pca$x[incomplete_names, 1:dim_pca]
   ```

3. **Stack Complete Omics (Horizontal Join)**:

   ```r
   completed_data_pca <- cbind(complete_mRNA_pca, complete_Methy_pca, complete_miRNA_pca)
   ```

4. **SVD on Stacked Complete Data**:

   ```r
   complete_data_svd <- svd(completed_data_pca)
   yc <- complete_data_svd$u[, 1:dim]
   ```

   * `yc`: Latent representation of complete samples (samples × `dim`)

5. **Project Incomplete Samples**:

   ```r
   yw_data1 <- t(yc) %*% t(ginv(complete_Methy_pca)) %*% t(incomplete_Methy_pca)
   yw_data2 <- t(yc) %*% t(ginv(complete_miRNA_pca)) %*% t(incomplete_miRNA_pca)
   yw <- t((yw_data1 + yw_data2)/2)
   ```

6. **Return Combined Latent Space**:

   ```r
   y <- rbind(yc, yw)
   return(y)
   ```

---

### 🧬 4. **Load Incomplete Omics Data**

```r
file1 <- "BIC"
file2 <- "BREAST"
```

These are used to build file paths:

```r
file1_name <- paste0("data/", file1, "/partial_data/", file2, "_Gene_0.1.txt")
file2_name <- paste0("data/", file1, "/", file2, "_Methy_Expression.txt")
file3_name <- paste0("data/", file1, "/", file2, "_Mirna_Expression.txt")
```

* Gene expression: Only **90% genes are missing** for some samples (0.1 means 10% present).
* Methylation and miRNA: assumed complete.

```r
v1 <- read.table(file1_name, ...)
v2 <- read.table(file2_name, ...)
v3 <- read.table(file3_name, ...)
```

---

### 📋 5. **Split Samples into Complete vs Incomplete**

```r
v1_colsum <- colSums(v1)
...
if ((v1_colsum[i] != 0) & (v2_colsum[i] != 0) & (v3_colsum[i] != 0)) {
  complete_names[...] <- ...
} else {
  incomplete_names[...] <- ...
}
```

* **Complete samples**: Have all omics data non-zero.
* **Incomplete samples**: One or more modalities are missing (i.e., sum = 0).

---

### 📐 6. **Build Latent Subspace Representation**

```r
data <- list(mRNA=v1, Methy=v2, miRNA=v3)
emb <- get_emb(data, complete_names, incomplete_names, all_names, dim=4, dim_pca=15, pca_scale=F)
```

* `dim = 4`: Subspace dimension (low-rank representation).
* `dim_pca = 15`: Top 15 PCA components per omics.

---

### 📊 7. **Cluster the Embedded Samples**

```r
cl_s <- specClust(emb, 5, nn = 20)
labels <- cl_s$cluster
names(labels) <- names(v1)
```

* Clustering is done on the latent subspace.
* `5`: Number of clusters (cancer subtypes).
* `nn = 20`: Used to construct similarity graph (k-nearest neighbors).

---

### ⚰️ 8. **Survival Analysis**

```r
survfile <- paste("data/", file1, "/", file2, "_Survival.txt", sep="")
surv <- read.table(survfile, header = T)
survresult <- survdiff(Surv(Survival, Death) ~ labels, data = surv)
p.val <- 1 - pchisq(survresult$chisq, length(survresult$n) - 1)
```

* `surv`: Contains `Survival` time and `Death` status for each patient.
* `survdiff()`: Performs **log-rank test** to test survival difference between clusters.
* `p.val`: Chi-square p-value (smaller = stronger cluster-survival association).

---

### 🧪 9. **Clinical Enrichment Analysis**

```r
para_num <- enrich(folder=file1, sign=T, label=labels)
```

* Tests for **clinical features** (e.g., stage, subtype) being enriched in any cluster.
* Returns number of significant clinical features (`para_num`).

---

### 📤 10. **Output Results**

```r
print(p.val)
print(para_num)
```

* Outputs:

  * `p.val`: How significant are the survival differences between subtypes.
  * `para_num`: How many clinical parameters are significantly enriched.

---

## ✅ Summary Table

| Step | Description                                             |
| ---- | ------------------------------------------------------- |
| 1–2  | Load libraries and helper functions                     |
| 3    | Define `get_emb()` to handle incomplete data projection |
| 4    | Load omics datasets with missing modalities             |
| 5    | Identify complete vs. incomplete samples                |
| 6    | Apply PCA and SVD to learn shared latent space          |
| 7    | Apply spectral clustering to latent space               |
| 8    | Run survival analysis to test cluster relevance         |
| 9    | Perform enrichment analysis on clinical parameters      |
| 10   | Output p-values and enrichment result                   |

---

## 📌 enrichment.R Explanation
---

## 🧠 Goal of These Functions

After clustering patients into subtypes, we want to know:

* Do **clinical features** like age, gender, or tumor stage **vary significantly** between subtypes?
* If yes, this implies that the subtypes might have **biological or clinical relevance**.

---

## 📄 Function 1: `cal.age.enrichment(age, cluster)`

```r
cal.age.enrichment <- function(age, cluster){
  age.test.result = kruskal.test(as.numeric(unlist(age)), as.numeric(unlist(cluster)))
  return(age.test.result)
}
```

### 🔍 Step-by-Step:

| Line                | What it does                                   |
| ------------------- | ---------------------------------------------- |
| `unlist(age)`       | Converts the age table to a numeric vector.    |
| `unlist(cluster)`   | Flattens cluster labels to match sample count. |
| `kruskal.test(...)` | Runs a **non-parametric Kruskal-Wallis test**. |
| `return(...)`       | Returns the test result (includes p-value).    |

### ✅ Purpose:

Tests if **age distributions differ significantly** across clusters.

---

## 📄 Function 2: `cal.discrete.enrichment(file_dir)`

```r
cal.discrete.enrichment <- function(file_dir){
  para_tb = table(as.data.frame(file_dir[,2:3]))
  enrichment.result <- chisq.test(para_tb)
  return(enrichment.result)
}
```

### 🔍 Step-by-Step:

| Line              | What it does                                      |
| ----------------- | ------------------------------------------------- |
| `file_dir[,2:3]`  | Takes clinical feature and cluster label columns. |
| `table(...)`      | Creates a contingency table.                      |
| `chisq.test(...)` | Performs a **Chi-square test of independence**.   |
| `return(...)`     | Returns the test result (includes p-value).       |

### ✅ Purpose:

Used for categorical (discrete) variables like **gender**, **stage**, **pathologic M/N/T**, etc.

---

## 📄 Function 3: `enrich(folder, sign=F, label)`

This is the **main function** that calls the above two and returns the number of **clinically enriched features** (i.e., features whose distribution varies significantly across clusters).

---

### 🔍 Full Step-by-Step Breakdown:

```r
enrich <- function(folder, sign=F, label) {
```

* `folder`: Name of cancer dataset folder (e.g., "BIC")
* `sign`: Whether to also test tumor-related features (T/N/M/stage)
* `label`: Vector of cluster assignments for each sample

---

### 🧪 Step 1: AGE (continuous feature)

```r
age_na <- paste0("data/", folder, "/age.txt")
age <- read.table(age_na, row.names = 1)
age_test_res = cal.age.enrichment(age, labels)
p_value = as.data.frame(as.matrix(age_test_res)[3,])
if (p_value < 0.05) para_num <- para_num + 1
```

* Loads patient age
* Tests if age varies by cluster
* If significant (`p < 0.05`), increments `para_num`

---

### 🚻 Step 2: GENDER (categorical)

```r
para_gen_na <- paste0("data/", folder, "/gender.txt")
para_gen <- read.table(para_gen_na, sep = "\t", header = T)
para_gen.pval = cal.discrete.enrichment(cbind(para_gen, labels))
...
if (p_value_gen < 0.05) para_num <- para_num + 1
```

* Reads gender info
* Tests gender-cluster independence using Chi-square test

---

### ✅ Step 3: (Conditional on `sign = TRUE`)

If you set `sign = TRUE`, the function **also checks tumor characteristics**:

| Feature             | File Name              |
| ------------------- | ---------------------- |
| **M** (metastasis)  | `pathologic_M.txt`     |
| **N** (lymph nodes) | `pathologic_N.txt`     |
| **T** (tumor size)  | `pathologic_T.txt`     |
| **Stage**           | `pathologic_stage.txt` |

Each is tested using `cal.discrete.enrichment()`:

```r
para_M_na <- paste0("data/", folder, "/pathologic_M.txt")
para_M <- read.table(para_M_na, sep = "\t", header = T)
...
if(p_value_M < 0.05) para_num <- para_num +1
```

This block is repeated for **N**, **T**, and **Stage**.

---

### 🔚 Final Step: Return the number of significant features

```r
return(para_num)
```

This value is used in the main clustering script to measure how **clinically meaningful** the clusters are.

---

## ✅ Summary Table

| Function                    | Description                                                              | Test Used      |
| --------------------------- | ------------------------------------------------------------------------ | -------------- |
| `cal.age.enrichment()`      | Tests if age differs across clusters                                     | Kruskal-Wallis |
| `cal.discrete.enrichment()` | Tests if a categorical feature (e.g. gender) is associated with clusters | Chi-square     |
| `enrich()`                  | Counts how many features are significantly enriched                      | Uses above two |

---

## 📦 Input File Expectations

Each clinical parameter is expected to be stored as a separate `.txt` file inside:
`data/<folder>/`

| File Name              | Example Values        |
| ---------------------- | --------------------- |
| `age.txt`              | Continuous age values |
| `gender.txt`           | "Male", "Female"      |
| `pathologic_M.txt`     | M0, M1                |
| `pathologic_N.txt`     | N0, N1, N2            |
| `pathologic_T.txt`     | T1, T2, T3            |
| `pathologic_stage.txt` | Stage I, II, III, IV  |

---

## 📌 run_data_col.R Explanation
The `data_col` function is a **core part** of the MCLS pipeline (Multi-omics Clustering based on Latent Subspace). Its purpose is to integrate **multi-omics data** (e.g., gene expression, methylation, miRNA) into a shared **low-dimensional latent space**, handling both **complete and incomplete data**.

---

## ✅ High-Level Purpose

This function:

* Accepts multi-omics data (e.g., mRNA, methylation, miRNA).
* Applies **PCA and SVD** to reduce dimensionality.
* Returns a **shared latent representation (l\_space)** of all samples.
* Supports:

  * ✅ **Complete data** (all views available for all samples)
  * ✅ **Incomplete data** (some views missing for some samples)

---

## 📄 Function Signature

```r
data_col <- function(data, incomplete_data = FALSE, incomplete_sample_name = FALSE, 
                     remain_view = FALSE, dim1 = 2, dim2 = 2, pca_scale = FALSE, seed = 0)
```

| Parameter                | Description                                        |
| ------------------------ | -------------------------------------------------- |
| `data`                   | A list of data matrices (`mRNA`, `Methy`, `miRNA`) |
| `incomplete_data`        | Set to `TRUE` if some samples have missing views   |
| `incomplete_sample_name` | Vector of sample names that are incomplete         |
| `remain_view`            | Number of available views per incomplete sample    |
| `dim1`                   | PCA components to keep per omics                   |
| `dim2`                   | Final dimension of latent space                    |
| `pca_scale`              | Whether to standardize data for PCA                |
| `seed`                   | For reproducibility of random selection            |

---

## 🧪 Case 1: **Complete Data**

```r
if (incomplete_data == F)
```

### 🔍 What happens:

1. **PCA** on each omics:

   ```r
   mRNA_pca <- prcomp(t(data[["mRNA"]]), ...)
   ```

2. **Concatenate PCA-reduced features**:

   ```r
   completed_data <- cbind(mRNA_pca$x[,1:dim1], Methy_pca$x[,1:dim1], miRNA_pca$x[,1:dim1])
   ```

3. **SVD on the concatenated data**:

   ```r
   data_svd <- svd(completed_data)
   U1 <- data_svd$u[, 1:dim2]
   ```

   → `U1`: Latent representation of all samples (rows = samples, cols = `dim2`)

4. **Projection matrices (g)**:
   For each omics:

   ```r
   g <- list(ginv(mRNA_pca$x[,1:dim1]) %*% U1, ...)
   ```

   These project omics features → latent space.

5. **Return**:

   ```r
   return(list(l_space = U1, f = f, g = g))
   ```

---

## 🧪 Case 2: **Incomplete Data**

```r
else {  # incomplete_data == TRUE
```

### 🔍 What happens:

#### 1. Split complete and incomplete samples

```r
sample_names <- colnames(data[[1]])
complete_sample_name <- setdiff(sample_names, incomplete_sample_name)
```

#### 2. Subset complete samples for PCA:

```r
completed_data[[i]] <- subset(data[[i]], select = complete_sample_name)
```

#### 3. PCA and SVD on complete samples (same as before)

```r
mRNA_pca, Methy_pca, miRNA_pca → combined → svd → U1
```

#### 4. Randomly choose available views for incomplete samples

```r
for (i in incomplete_sample_name) {
  random_selected_view[[i]] <- sample(c(1,2,3), size = remain_view)
}
```

* This randomly assigns 1 or 2 of the 3 omics for each incomplete sample.

#### 5. Project incomplete samples to latent space

##### If `remain_view == 2`:

* For each sample:

  ```r
  a = (data_from_view1 × PCA rotation) × projection g
  b = (data_from_view2 × PCA rotation) × projection g
  converted_incomplete_data[[i]] = (a + b) / 2
  ```

##### If `remain_view == 1`:

* Same as above, but with only one view (`a` only)

#### 6. Append incomplete samples to `U1`

```r
U <- rbind(U1, converted_incomplete_data[[i]])
```

#### 7. Return

```r
return(list(l_space = U, f = f, g = g, views = random_selected_view))
```

---

## 🧠 Summary Table

| Stage                  | Complete Data          | Incomplete Data                    |
| ---------------------- | ---------------------- | ---------------------------------- |
| PCA                    | ✔️ (on all samples)    | ✔️ (on complete samples only)      |
| SVD                    | ✔️ on combined PCA     | ✔️ on combined PCA (complete only) |
| Projection `g`         | ✔️ from omics → latent | ✔️ from complete data              |
| Latent space `l_space` | All samples            | Complete + Projected Incomplete    |
| Views used             | All 3                  | Random 1 or 2 views                |

---

## 📌 What Are `f` and `g`?

* `f`: List of PCA rotation matrices (for each omics)

  * Projects **original omics → PCA**
* `g`: Projection matrices

  * Projects **PCA → latent space**

Together, they help reconstruct how an individual sample is represented in the final low-dimensional shared space.

---

## 📤 Output of `data_col`

```r
list(
  l_space = U,      # Final latent representation
  f = f,            # PCA rotations (one for each omics)
  g = g,            # Projection matrices to latent space
  views = ...       # (only for incomplete) random views per sample
)
```

---

## ✅ Example Usage

```r
result <- data_col(data, incomplete_data = TRUE, incomplete_sample_name = c("S1", "S5"), remain_view = 2, dim1 = 13, dim2 = 4)
embedding <- result$l_space
```
## 📌 specClust.R Explanation
The function `specClust()` in this code implements a **spectral clustering algorithm** with **automatic neighbor selection**, **graph connectivity checking**, and **eigenvector analysis**. It's an essential part of the MCLS project, used to cluster samples based on their **low-dimensional latent representation** (`l_space`) produced by previous steps like PCA and SVD.

Let’s break it down **step-by-step** and explain each block.

---

## 🔁 Function Signature

```r
specClust <- function (data, centers = NULL, nn = 7, method = "symmetric", 
                       gmax = NULL, max.iter = 10000, ...) 
```

| Parameter  | Meaning                                                                           |
| ---------- | --------------------------------------------------------------------------------- |
| `data`     | Matrix of samples (rows) and features (columns), often the latent space `l_space` |
| `centers`  | Number of clusters `k` (if `NULL`, it will estimate `k`)                          |
| `nn`       | Number of nearest neighbors used to build graph                                   |
| `method`   | Type of Laplacian ("symmetric" or "random walk")                                  |
| `gmax`     | Max number of connected components allowed                                        |
| `max.iter` | Max iterations for ARPACK eigen decomposition                                     |
| `...`      | Additional arguments for `kmeans()`                                               |

---

## 🔍 Step-by-Step Explanation

### 1. 📦 Preprocessing Data

```r
if (is.data.frame(data)) data = as.matrix(data)
da = apply(data, 1, paste, collapse = "#")
indUnique = which(!duplicated(da))
indAll = match(da, da[indUnique])
data2 = data
data = data[indUnique, ]
```

* Ensures `data` is a matrix.
* Removes duplicate rows.
* Stores mappings to recover original sample positions later.

---

### 2. ⚖️ Normalize the data

```r
data = scale(data, FALSE, TRUE)
```

* Standardizes data across **columns (features)** to have unit variance (but not mean-centered).

---

### 3. 🔁 Adjust `gmax` based on `centers`

```r
if (is.null(gmax)) {
  if (!is.null(centers)) gmax = centers - 1L
  else gmax = 1L
}
```

* Sets `gmax` to limit how fragmented the graph is allowed to be.

---

### 4. 🔗 Construct Graph and Check Connectivity

```r
while (test) {
  DC = mydist(data, nn)
  sif <- rbind(1:n, as.vector(DC[[2]]))
  g <- graph(sif, directed = FALSE)
  g <- decompose(g, min.vertices = 4)
  ...
}
```

* Uses `mydist()` to find nearest neighbors (`nn`).
* Constructs an undirected graph where nodes = samples, edges = similarity.
* If the graph breaks into too many components (`> gmax`), increases `nn` and rebuilds until it's connected enough.

---

### 5. 🧠 Compute Affinity/Weight Matrix `W`

```r
W <- DC[[1]]  # distances
wi <- W[, nn] # distance to farthest neighbor
SC <- matrix(1, nrow(W), nn)
SC[] <- wi[DC[[2]]] * wi
W = W^2/SC
```

* Rescales distances using a similarity metric.
* Applies a Gaussian kernel-like transformation:

```r
alpha = 1/(2 * (nn + 1))
qua = abs(qnorm(alpha))
W = W * qua
W = dnorm(W, sd = 1)
```

* Applies a normal distribution to convert distances into affinities.
* This is a type of **adaptive kernel similarity**.

---

### 6. 🧮 Construct Laplacian Matrix

```r
L = Laplacian(DC, nn, method)
```

* Creates the **graph Laplacian** matrix used for spectral clustering.

  * `method="symmetric"` uses normalized symmetric Laplacian.
  * `method="random"` would use a random walk Laplacian.

---

### 7. 🔍 Eigen Decomposition via ARPACK

```r
U <- arpack(f, extra = L, options = list(...), sym = TRUE)
```

* Computes the **eigenvectors of the Laplacian** using `arpack()`, a fast method for large matrices.
* `U[[1]]`: eigenvalues
  `U[[2]]`: eigenvectors

---

### 8. 🔎 Determine `k` Automatically (if `centers` is `NULL`)

```r
tmp = which.max(diff(U[[1]])) + 1
centers = which.min(AUC(U[[1]][1:tmp]))
```

* Looks for the largest **gap** between eigenvalues (scree plot elbow).
* Then chooses the best `k` using **area under curve (AUC)** minimization.

---

### 9. 🧽 Normalize Eigenvectors (if `method == symmetric`)

```r
rs = sqrt(rowSums(U[[2]]^2))
U[[2]] = U[[2]] / rs
```

* Required for symmetric Laplacian to ensure unit-length rows.

---

### 10. 📊 Apply K-Means Clustering

```r
result = kmeans(U[[2]], centers = centers, nstart = 20, ...)
```

* Clusters the samples in the eigenvector space.

---

### 11. 🔁 Assign Archetypes & Restore Mapping

```r
archeType = getClosest(U[[2]][indAll, ], result$centers)
```

* `getClosest()` finds the closest point to each cluster center.
* Remaps the data back to include any duplicated rows that were removed earlier.

---

### 12. 🧾 Output Object

```r
result$eigenvalue = U[[1]]
result$eigenvector = U[[2]]
result$data = data2
result$indAll = indAll
result$indUnique = indUnique
result$L = L
result$archetype = archeType
result$call = call
class(result) = c("specClust", "kmeans")
```

* Adds extra useful info to the result, like eigenvalues, Laplacian, etc.

---

## ✅ Final Output: `specClust` Object

Returns an object similar to `kmeans()` output with extra fields:

| Field          | Description                          |
| -------------- | ------------------------------------ |
| `$cluster`     | Cluster assignments                  |
| `$centers`     | Cluster centroids                    |
| `$eigenvector` | Eigenvectors of Laplacian            |
| `$eigenvalue`  | Eigenvalues                          |
| `$L`           | Laplacian matrix                     |
| `$archetype`   | Closest data points to each centroid |
| `$call`        | Original function call               |

---

## 🔁 Dependencies

You must have:

* `mydist()`: a function that computes nearest neighbors
* `Laplacian()`: creates Laplacian matrix from neighbor data
* `arpack()`: computes eigenvectors (from `RSpectra` or similar)
* `getClosest()`: matches cluster centers to data points
* `graph()` and `decompose()`: from the **`igraph`** package

---

## 💡 Summary

| Step | Description                                  |
| ---- | -------------------------------------------- |
| 1    | Remove duplicates and scale data             |
| 2    | Build nearest-neighbor graph                 |
| 3    | Compute similarity matrix                    |
| 4    | Create Laplacian matrix                      |
| 5    | Extract eigenvectors                         |
| 6    | Cluster samples in spectral space            |
| 7    | Return result with clustering and graph info |

---
