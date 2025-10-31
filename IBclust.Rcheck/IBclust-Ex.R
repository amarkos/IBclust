pkgname <- "IBclust"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('IBclust')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
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
data_mix <- data.frame(
  cat_var = factor(sample(letters[1:3], 100, replace = TRUE)),      # Nominal categorical variable
  ord_var = factor(sample(c("low", "medium", "high"), 100, replace = TRUE),
                   levels = c("low", "medium", "high"),
                   ordered = TRUE),                                # Ordinal variable
  cont_var1 = rnorm(100),                                          # Continuous variable 1
  cont_var2 = runif(100)                                           # Continuous variable 2
)

# Perform Mixed-Type Hierarchical Clustering with Agglomerative IB
result_mix <- AIBmix(X = data_mix, lambda = -1, s = -1, scale = TRUE)

# Print clustering results
plot(result_mix, type = "dendrogram", xlab = "", sub = "", cex = 0.5)  # Plot dendrogram
plot(result_mix, type = "info", col = "black", pch = 16)  # Plot dendrogram

# Simulated categorical data example
set.seed(123)
data_cat <- data.frame(
  Var1 = as.factor(sample(letters[1:3], 200, replace = TRUE)),  # Nominal variable
  Var2 = as.factor(sample(letters[4:6], 200, replace = TRUE)),  # Nominal variable
  Var3 = factor(sample(c("low", "medium", "high"), 200, replace = TRUE),
                levels = c("low", "medium", "high"), ordered = TRUE)  # Ordinal variable
)

# Run AIBmix with automatic lambda selection 
result_cat <- AIBmix(X = data_cat, lambda = -1)

# Print clustering results
plot(result_cat, type = "dendrogram", xlab = "", sub = "", cex = 0.5)  # Plot dendrogram

# Results summary
summary(result_cat)

# Simulated continuous data example
set.seed(123)
# Continuous data with 200 observations, 5 features
data_cont <- as.data.frame(matrix(rnorm(1000), ncol = 5))

# Run AIBmix with automatic bandwidth selection 
result_cont <- AIBmix(X = data_cont, s = -1, scale = TRUE)

# Print concise summary ofoutput
print(result_cont)

# Print clustering results
plot(result_cont, type = "dendrogram", xlab = "", sub = "", cex = 0.5)  # Plot dendrogram



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
data_mix <- data.frame(
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
result_mix <- DIBmix(X = data_mix, ncl = 3, nstart = 5)

# View results
print(paste("Number of clusters found:", length(unique(result_mix$Cluster))))
print(paste("Mutual Information:", round(result_mix$MutualInfo, 3)))
table(result_mix$Cluster)

# Example 2: Comparing cat_first parameter
# When categorical variables are more informative
result_cat_first <- DIBmix(X = data_mix, ncl = 3,
                           cat_first = TRUE,  # Prioritize categorical variables
                           nstart = 5)

# When continuous variables are more informative (default)
result_cont_first <- DIBmix(X = data_mix, ncl = 3,
                            cat_first = FALSE,
                            nstart = 5)

# Compare clustering performance
if (requireNamespace("mclust", quietly = TRUE)){  # For adjustedRandIndex
  print(paste("Agreement between approaches:", 
              round(mclust::adjustedRandIndex(result_cat_first$Cluster, 
                    result_cont_first$Cluster), 3)))
  }

plot(result_cat_first, type = "sizes") # Bar plot of cluster sizes
plot(result_cat_first, type = "info")  # Information-theoretic quantities plot
plot(result_cat_first, type = "beta")  # Plot of evolution of beta

# Simulated categorical data example
data_cat <- data.frame(
  Var1 = as.factor(sample(letters[1:3], 200, replace = TRUE)),  # Nominal variable
  Var2 = as.factor(sample(letters[4:6], 200, replace = TRUE)),  # Nominal variable
  Var3 = factor(sample(c("low", "medium", "high"), 200, replace = TRUE),
                levels = c("low", "medium", "high"), ordered = TRUE)  # Ordinal variable
)

# Perform hard clustering on categorical data with Deterministic IB
result_cat <- DIBmix(X = data_cat, ncl = 3, lambda = -1, nstart = 5)

# Print clustering results
print(result_cat$Cluster)       # Cluster assignments
print(result_cat$Entropy)       # Final entropy
print(result_cat$MutualInfo)    # Mutual information

# Simulated continuous data example
set.seed(123)
# Continuous data with 200 observations, 5 features
data_cont <- as.data.frame(matrix(rnorm(1000), ncol = 5))

# Perform hard clustering on continuous data with Deterministic IB
result_cont <- DIBmix(X = data_cont, ncl = 3, s = -1, nstart = 5)

# Print clustering results
print(result_cont$Cluster)       # Cluster assignments
print(result_cont$Entropy)       # Final entropy
print(result_cont$MutualInfo)    # Mutual information

# Summary of output
print(result_cont)
summary(result_cont)



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
data_mix <- data.frame(
  cat_var = factor(sample(letters[1:3], 100, replace = TRUE)),      # Nominal categorical variable
  ord_var = factor(sample(c("low", "medium", "high"), 100, replace = TRUE),
                   levels = c("low", "medium", "high"),
                   ordered = TRUE),                                # Ordinal variable
  cont_var1 = rnorm(100),                                          # Continuous variable 1
  cont_var2 = runif(100)                                           # Continuous variable 2
)

# Perform Mixed-Type Fuzzy Clustering with Generalised IB
result_mix <- GIBmix(X = data_mix, ncl = 3, beta = 2, alpha = 0.5, nstart = 5)

# Print clustering results
print(result_mix$Cluster)       # Cluster membership matrix
print(result_mix$Entropy)       # Entropy of final clustering
print(result_mix$CondEntropy)   # Conditional entropy of final clustering
print(result_mix$MutualInfo)    # Mutual information between Y and T

# Summary of output
summary(result_mix)

# Simulated categorical data example
set.seed(123)
data_cat <- data.frame(
  Var1 = as.factor(sample(letters[1:3], 200, replace = TRUE)),  # Nominal variable
  Var2 = as.factor(sample(letters[4:6], 200, replace = TRUE)),  # Nominal variable
  Var3 = factor(sample(c("low", "medium", "high"), 200, replace = TRUE),
                levels = c("low", "medium", "high"), ordered = TRUE)  # Ordinal variable
)

# Perform Fuzzy Clustering on categorical data with Generalised IB
result_cat <- GIBmix(X = data_cat, ncl = 2, beta = 25, alpha = 0.75, lambda = -1, nstart = 5)

# Print clustering results
print(result_cat$Cluster)       # Cluster membership matrix
print(result_cat$Entropy)       # Entropy of final clustering
print(result_cat$CondEntropy)   # Conditional entropy of final clustering
print(result_cat$MutualInfo)    # Mutual information between Y and T

# Simulated continuous data example
set.seed(123)
# Continuous data with 200 observations, 5 features
data_cont <- as.data.frame(matrix(rnorm(1000), ncol = 5))

# Perform Fuzzy Clustering on continuous data with Generalised IB
result_cont <- GIBmix(X = data_cont, ncl = 2, beta = 50, alpha = 0.75, s = -1, nstart = 5)

# Print clustering results
print(result_cont$Cluster)       # Cluster membership matrix
print(result_cont$Entropy)       # Entropy of final clustering
print(result_cont$CondEntropy)   # Conditional entropy of final clustering
print(result_cont$MutualInfo)    # Mutual information between Y and T

plot(result_cont, type = "sizes") # Bar plot of cluster sizes (hardened assignments)
plot(result_cont, type = "info")  # Information-theoretic quantities plot



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
data_mix <- data.frame(
  cat_var = factor(sample(letters[1:3], 100, replace = TRUE)),      # Nominal categorical variable
  ord_var = factor(sample(c("low", "medium", "high"), 100, replace = TRUE),
                   levels = c("low", "medium", "high"),
                   ordered = TRUE),                                # Ordinal variable
  cont_var1 = rnorm(100),                                          # Continuous variable 1
  cont_var2 = runif(100)                                           # Continuous variable 2
)

# Perform Mixed-Type Fuzzy Clustering
result_mix <- IBmix(X = data_mix, ncl = 3, beta = 2, nstart = 1)

# Print clustering results
print(result_mix$Cluster)       # Cluster membership matrix
print(result_mix$InfoXT)        # Mutual information between X and T
print(result_mix$MutualInfo)    # Mutual information between Y and T

# Summary of output
summary(result_mix)

# Simulated categorical data example
set.seed(123)
data_cat <- data.frame(
  Var1 = as.factor(sample(letters[1:3], 100, replace = TRUE)),  # Nominal variable
  Var2 = as.factor(sample(letters[4:6], 100, replace = TRUE)),  # Nominal variable
  Var3 = factor(sample(c("low", "medium", "high"), 100, replace = TRUE),
                levels = c("low", "medium", "high"), ordered = TRUE)  # Ordinal variable
)

# Perform fuzzy clustering on categorical data with standard IB
result_cat <- IBmix(X = data_cat, ncl = 3, beta = 15, lambda = -1, nstart = 2, maxiter = 200)

# Print clustering results
print(result_cat$Cluster)       # Cluster membership matrix
print(result_cat$InfoXT)        # Mutual information between X and T
print(result_cat$MutualInfo)    # Mutual information between Y and T

plot(result_cat, type = "sizes") # Bar plot of cluster sizes (hardened assignments)
plot(result_cat, type = "info")  # Information-theoretic quantities plot

# Simulated continuous data example
set.seed(123)
# Continuous data with 100 observations, 5 features
data_cont <- as.data.frame(matrix(rnorm(500), ncol = 5))

# Perform fuzzy clustering on continuous data with standard IB
result_cont <- IBmix(X = data_cont, ncl = 3, beta = 50, s = -1, nstart = 2)

# Print clustering results
print(result_cont$Cluster)       # Cluster membership matrix
print(result_cont$InfoXT)        # Mutual information between X and T
print(result_cont$MutualInfo)    # Mutual information between Y and T



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
