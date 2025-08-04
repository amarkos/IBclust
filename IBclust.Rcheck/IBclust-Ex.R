pkgname <- "IBclust"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('IBclust')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("AIBcat")
### * AIBcat

flush(stderr()); flush(stdout())

### Name: AIBcat
### Title: Cluster Categorical Data Using the Agglomerative Information
###   Bottleneck Algorithm
### Aliases: AIBcat
### Keywords: clustering

### ** Examples

# Simulated categorical data
set.seed(123)
X <- data.frame(
  Var1 = as.factor(sample(letters[1:3], 200, replace = TRUE)),  # Nominal variable
  Var2 = as.factor(sample(letters[4:6], 200, replace = TRUE)),  # Nominal variable
  Var3 = factor(sample(c("low", "medium", "high"), 200, replace = TRUE),
                levels = c("low", "medium", "high"), ordered = TRUE)  # Ordinal variable
)

# Run AIBcat with automatic lambda selection 
result <- AIBcat(X = X, lambda = -1)

# Print clustering results
plot(result$dendrogram, xlab = "", sub = "")  # Plot dendrogram



cleanEx()
nameEx("AIBcont")
### * AIBcont

flush(stderr()); flush(stdout())

### Name: AIBcont
### Title: Cluster Continuous Data Using the Agglomerative Information
###   Bottleneck Algorithm
### Aliases: AIBcont
### Keywords: clustering

### ** Examples

# Generate simulated continuous data
set.seed(123)
X <- matrix(rnorm(1000), ncol = 5)  # 200 observations, 5 features

# Run AIBcont with automatic bandwidth selection 
result <- AIBcont(X = X, s = -1, scale = TRUE)

# Print clustering results
plot(result$dendrogram, xlab = "", sub = "")  # Plot dendrogram



cleanEx()
nameEx("AIBmix")
### * AIBmix

flush(stderr()); flush(stdout())

### Name: AIBmix
### Title: Agglomerative Information Bottleneck Clustering for Mixed-Type
###   Data
### Aliases: AIBmix
### Keywords: clustering

### ** Examples

# Example dataset with categorical, ordinal, and continuous variables
set.seed(123)
data <- data.frame(
  cat_var = factor(sample(letters[1:3], 100, replace = TRUE)),      # Nominal categorical variable
  ord_var = factor(sample(c("low", "medium", "high"), 100, replace = TRUE),
                   levels = c("low", "medium", "high"),
                   ordered = TRUE),                                # Ordinal variable
  cont_var1 = rnorm(100),                                          # Continuous variable 1
  cont_var2 = runif(100)                                           # Continuous variable 2
)

# Perform Mixed-Type Hierarchical Clustering with Agglomerative IB
result <- AIBmix(X = data, catcols = 1:2, contcols = 3:4, lambda = -1, s = -1, scale = TRUE)

# Print clustering results
plot(result$dendrogram, xlab = "", sub = "")  # Plot dendrogram



cleanEx()
nameEx("DIBcat")
### * DIBcat

flush(stderr()); flush(stdout())

### Name: DIBcat
### Title: Cluster Categorical Data Using the Deterministic Information
###   Bottleneck Algorithm
### Aliases: DIBcat
### Keywords: clustering

### ** Examples

# Simulated categorical data
set.seed(123)
X <- data.frame(
  Var1 = as.factor(sample(letters[1:3], 200, replace = TRUE)),  # Nominal variable
  Var2 = as.factor(sample(letters[4:6], 200, replace = TRUE)),  # Nominal variable
  Var3 = factor(sample(c("low", "medium", "high"), 200, replace = TRUE),
                levels = c("low", "medium", "high"), ordered = TRUE)  # Ordinal variable
)

# Run DIBcat with automatic lambda selection and multiple initializations
result <- DIBcat(X = X, ncl = 3, lambda = -1, nstart = 10)

# Print clustering results
print(result$Cluster)       # Cluster assignments
print(result$Entropy)       # Final entropy
print(result$MutualInfo)    # Mutual information



cleanEx()
nameEx("DIBcont")
### * DIBcont

flush(stderr()); flush(stdout())

### Name: DIBcont
### Title: Cluster Continuous Data Using the Deterministic Information
###   Bottleneck Algorithm
### Aliases: DIBcont
### Keywords: clustering

### ** Examples

# Generate simulated continuous data
set.seed(123)
X <- matrix(rnorm(200), ncol = 5)  # 200 observations, 5 features

# Run DIBcont with automatic bandwidth selection and multiple initializations
result <- DIBcont(X = X, ncl = 3, s = -1, nstart = 10)

# Print clustering results
print(result$Cluster)       # Cluster assignments
print(result$Entropy)       # Final entropy
print(result$MutualInfo)    # Mutual information



cleanEx()
nameEx("DIBmix")
### * DIBmix

flush(stderr()); flush(stdout())

### Name: DIBmix
### Title: Deterministic Information Bottleneck Clustering for Mixed-Type
###   Data
### Aliases: DIBmix
### Keywords: clustering

### ** Examples

# Example 1: Basic Mixed-Type Clustering
set.seed(123)

# Create a more realistic dataset with mixed variable types
data <- data.frame(
  # Categorical variables
  education = factor(sample(c("High School", "Bachelor", "Master", "PhD"), 150, 
                           replace = TRUE, prob = c(0.4, 0.3, 0.2, 0.1))),
  employment = factor(sample(c("Full-time", "Part-time", "Unemployed"), 150, 
                            replace = TRUE, prob = c(0.6, 0.25, 0.15))),
  
  # Ordinal variable
  satisfaction = factor(sample(c("Low", "Medium", "High"), 150, replace = TRUE),
                       levels = c("Low", "Medium", "High"), ordered = TRUE),
  
  # Continuous variables  
  income = rlnorm(150, meanlog = 10, sdlog = 0.5),  # Log-normal income
  age = rnorm(150, mean = 35, sd = 10),             # Normally distributed age
  experience = rpois(150, lambda = 8)               # Years of experience
)

# Perform Mixed-Type Clustering
result <- DIBmix(X = data, ncl = 3, 
                catcols = c(1, 2, 3),    # education, employment, satisfaction
                contcols = c(4, 5, 6),   # income, age, experience
                nstart = 50)

# View results
print(paste("Number of clusters found:", length(unique(result$Cluster))))
print(paste("Mutual Information:", round(result$MutualInfo, 3)))
table(result$Cluster)

# Example 2: Comparing cat_first parameter
# When categorical variables are more informative
result_cat_first <- DIBmix(X = data, ncl = 3,
                          catcols = c(1, 2, 3), 
                          contcols = c(4, 5, 6),
                          cat_first = TRUE,  # Prioritize categorical variables
                          nstart = 50)

# When continuous variables are more informative (default)
result_cont_first <- DIBmix(X = data, ncl = 3,
                           catcols = c(1, 2, 3), 
                           contcols = c(4, 5, 6),
                           cat_first = FALSE,
                           nstart = 50)

# Compare clustering performance
if (requireNamespace("mclust", quietly = TRUE)){  # For adjustedRandIndex
  print(paste("Agreement between approaches:", 
              round(mclust::adjustedRandIndex(result_cat_first$Cluster, 
                    result_cont_first$Cluster), 3)))
  }



cleanEx()
nameEx("GIBcat")
### * GIBcat

flush(stderr()); flush(stdout())

### Name: GIBcat
### Title: Cluster Categorical Data Using the Generalised Information
###   Bottleneck Algorithm
### Aliases: GIBcat
### Keywords: clustering

### ** Examples

# Simulated categorical data
set.seed(123)
X <- data.frame(
  Var1 = as.factor(sample(letters[1:3], 200, replace = TRUE)),  # Nominal variable
  Var2 = as.factor(sample(letters[4:6], 200, replace = TRUE)),  # Nominal variable
  Var3 = factor(sample(c("low", "medium", "high"), 200, replace = TRUE),
                levels = c("low", "medium", "high"), ordered = TRUE)  # Ordinal variable
)

# Run GIBcat with automatic lambda selection and multiple initializations
result <- GIBcat(X = X, ncl = 2, beta = 25, alpha = 0.75, lambda = -1, nstart = 10)

# Print clustering results
print(result$Cluster)       # Cluster membership matrix
print(result$Entropy)       # Entropy of final clustering
print(result$CondEntropy)   # Conditional entropy of final clustering
print(result$MutualInfo)    # Mutual information between Y and T



cleanEx()
nameEx("GIBcont")
### * GIBcont

flush(stderr()); flush(stdout())

### Name: GIBcont
### Title: Cluster Continuous Data Using the Generalised Information
###   Bottleneck Algorithm
### Aliases: GIBcont
### Keywords: clustering

### ** Examples

# Generate simulated continuous data
set.seed(123)
X <- matrix(rnorm(200), ncol = 5)  # 200 observations, 5 features

# Run GIBcont with automatic bandwidth selection and multiple initializations
result <- GIBcont(X = X, ncl = 2, beta = 50, alpha = 0.75, s = -1, nstart = 10)

# Print clustering results
print(result$Cluster)       # Cluster membership matrix
print(result$Entropy)       # Entropy of final clustering
print(result$CondEntropy)   # Conditional entropy of final clustering
print(result$MutualInfo)    # Mutual information between Y and T



cleanEx()
nameEx("GIBmix")
### * GIBmix

flush(stderr()); flush(stdout())

### Name: GIBmix
### Title: Generalised Information Bottleneck Clustering for Mixed-Type
###   Data
### Aliases: GIBmix
### Keywords: clustering

### ** Examples

# Example dataset with categorical, ordinal, and continuous variables
set.seed(123)
data <- data.frame(
  cat_var = factor(sample(letters[1:3], 100, replace = TRUE)),      # Nominal categorical variable
  ord_var = factor(sample(c("low", "medium", "high"), 100, replace = TRUE),
                   levels = c("low", "medium", "high"),
                   ordered = TRUE),                                # Ordinal variable
  cont_var1 = rnorm(100),                                          # Continuous variable 1
  cont_var2 = runif(100)                                           # Continuous variable 2
)

# Perform Mixed-Type Fuzzy Clustering with Generalised IB
result <- GIBmix(X = data, ncl = 3, beta = 2, alpha = 0.5, 
catcols = 1:2, contcols = 3:4, nstart = 20)

# Print clustering results
print(result$Cluster)       # Cluster membership matrix
print(result$Entropy)       # Entropy of final clustering
print(result$CondEntropy)   # Conditional entropy of final clustering
print(result$MutualInfo)    # Mutual information between Y and T



cleanEx()
nameEx("IBcat")
### * IBcat

flush(stderr()); flush(stdout())

### Name: IBcat
### Title: Cluster Categorical Data Using the Information Bottleneck
###   Algorithm
### Aliases: IBcat
### Keywords: clustering

### ** Examples

# Simulated categorical data
set.seed(123)
X <- data.frame(
  Var1 = as.factor(sample(letters[1:3], 200, replace = TRUE)),  # Nominal variable
  Var2 = as.factor(sample(letters[4:6], 200, replace = TRUE)),  # Nominal variable
  Var3 = factor(sample(c("low", "medium", "high"), 200, replace = TRUE),
                levels = c("low", "medium", "high"), ordered = TRUE)  # Ordinal variable
)

# Run IBcat with automatic lambda selection and multiple initializations
result <- IBcat(X = X, ncl = 3, beta = 15, lambda = -1, nstart = 10)

# Print clustering results
print(result$Cluster)       # Cluster membership matrix
print(result$InfoXT)       # Mutual information between X and T
print(result$MutualInfo)    # Mutual information between Y and T



cleanEx()
nameEx("IBcont")
### * IBcont

flush(stderr()); flush(stdout())

### Name: IBcont
### Title: Cluster Continuous Data Using the Information Bottleneck
###   Algorithm
### Aliases: IBcont
### Keywords: clustering

### ** Examples

# Generate simulated continuous data
set.seed(123)
X <- matrix(rnorm(200), ncol = 5)  # 200 observations, 5 features

# Run IBcont with automatic bandwidth selection and multiple initializations
result <- IBcont(X = X, ncl = 3, beta = 50, s = -1, nstart = 10)

# Print clustering results
print(result$Cluster)       # Cluster membership matrix
print(result$InfoXT)       # Mutual information between X and T
print(result$MutualInfo)    # Mutual information between Y and T



cleanEx()
nameEx("IBmix")
### * IBmix

flush(stderr()); flush(stdout())

### Name: IBmix
### Title: Information Bottleneck Clustering for Mixed-Type Data
### Aliases: IBmix
### Keywords: clustering

### ** Examples

# Example dataset with categorical, ordinal, and continuous variables
set.seed(123)
data <- data.frame(
  cat_var = factor(sample(letters[1:3], 100, replace = TRUE)),      # Nominal categorical variable
  ord_var = factor(sample(c("low", "medium", "high"), 100, replace = TRUE),
                   levels = c("low", "medium", "high"),
                   ordered = TRUE),                                # Ordinal variable
  cont_var1 = rnorm(100),                                          # Continuous variable 1
  cont_var2 = runif(100)                                           # Continuous variable 2
)

# Perform Mixed-Type Fuzzy Clustering
result <- IBmix(X = data, ncl = 3, beta = 2, catcols = 1:2, contcols = 3:4, nstart = 10)

# Print clustering results
print(result$Cluster)       # Cluster membership matrix
print(result$InfoXT)       # Mutual information between X and T
print(result$MutualInfo)    # Mutual information between Y and T



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
