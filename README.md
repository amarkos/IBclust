# IBclust Package

**IBclust** is an R package for clustering datasets using the Information Bottleneck method and its variants. This package supports datasets with mixed-type variables (nominal, ordinal, and continuous), as well as datasets that are purely continuous or categorical. The IB approach preserves the most relevant information while forming concise and interpretable clusters, guided by principles from information theory.
<!-- , as introduced in [Costa, Papatsouma, and Markos (2024)](https://arxiv.org/abs/2407.03389) -->
## Installation

You can install the latest version of the package directly from GitHub using `devtools`:

```r
install.packages("devtools")  # Install devtools if not already installed
devtools::install_github("amarkos/IBclust")  # Install IBclust from GitHub
```

## Getting Started
Below is a comprehensive example demonstrating how to use the package for clustering mixed-type, continuous, and categorical datasets, and displaying the results. The examples make use of the Deterministic Information Bottleneck (DIB) method for clustering; other options include the Agglomerative IB for hierarchical clustering, the Generalised IB and the standard IB for fuzzy clustering.

```r
library(IBclust)

# Example Mixed-Type Data
data <- data.frame(
  cat_var = factor(sample(letters[1:3], 100, replace = TRUE)),      # Nominal categorical variable
  ord_var = factor(sample(c("low", "medium", "high"), 100, replace = TRUE),
                   levels = c("low", "medium", "high"),
                   ordered = TRUE),                                # Ordinal variable
  cont_var1 = rnorm(100),                                          # Continuous variable 1
  cont_var2 = runif(100)                                           # Continuous variable 2
)

# Perform Mixed-Type Clustering using the Deterministic variant and automatic bandwidth selection
result_mix <- DIBmix(X = data, ncl = 3)
cat("Mixed-Type Clustering Results:\n")
print(result_mix$Cluster)
print(result_mix$Entropy)
print(result_mix$MutualInfo)

# Example Continuous Data
X_cont <- as.data.frame(matrix(rnorm(1000), ncol = 5))  # 200 observations, 5 features

# Perform Continuous Data Clustering 
result_cont <- DIBmix(X = X_cont, ncl = 3, s = -1, nstart = 50)
cat("Continuous Clustering Results:\n")
print(result_cont$Cluster)
print(result_cont$Entropy)
print(result_cont$MutualInfo)

# Example Categorical Data
X_cat <- data.frame(
  Var1 = factor(sample(letters[1:3], 200, replace = TRUE)),  # Nominal variable
  Var2 = factor(sample(letters[4:6], 200, replace = TRUE)),  # Nominal variable
  Var3 = factor(sample(c("low", "medium", "high"), 200, replace = TRUE),
                levels = c("low", "medium", "high"), ordered = TRUE)  # Ordinal variable
)

# Perform Categorical Data Clustering
result_cat <- DIBmix(X = X_cat, ncl = 3, lambda = -1, nstart = 50)
cat("Categorical Clustering Results:\n")
print(result_cat$Cluster)
print(result_cat$Entropy)
print(result_cat$MutualInfo)
```

## Contributing
Contributions are welcome! If you encounter issues, have suggestions, or would like to enhance the package, please feel free to submit an issue or a pull request on the GitHub repository.

## License
This package is distributed under the GPL-3 License. See the [GNU General Public License version 3](https://www.gnu.org/licenses/gpl-3.0.html) for details.

