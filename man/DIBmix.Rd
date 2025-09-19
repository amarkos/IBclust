\name{DIBmix}
\alias{DIBmix}
\title{Deterministic Information Bottleneck Clustering for Mixed-Type Data}
\description{
The \code{DIBmix} function implements the Deterministic Information Bottleneck (DIB) algorithm
for clustering datasets containing continuous, categorical (nominal and ordinal), and mixed-type variables.
This method optimizes an information-theoretic objective to preserve
relevant information in the cluster assignments while achieving effective data compression
\insertCite{costa_dib_2025}{IBclust}.
}
\usage{
DIBmix(X, ncl, randinit = NULL,
       s = -1, lambda = -1, scale = TRUE,
       maxiter = 100, nstart = 100,
       contkernel = "gaussian",
       nomkernel = "aitchisonaitken", ordkernel = "liracine",
       cat_first = FALSE, verbose = FALSE)
}
\arguments{
\item{X}{A data frame containing the input data to be clustered. It should include categorical variables
(\code{factor} for nominal and \code{ordered} for ordinal) and continuous variables (\code{numeric}).}

\item{ncl}{An integer specifying the number of clusters.}

\item{randinit}{An optional vector specifying the initial cluster assignments. If \code{NULL}, cluster assignments are initialized randomly.}

\item{s}{A numeric value or vector specifying the bandwidth parameter(s) for continuous variables. The values must be greater than \eqn{0}. The default value is \eqn{-1}, which enables the automatic selection of optimal bandwidth(s). Argument is ignored when no variables are continuous.}

\item{lambda}{A numeric value or vector specifying the bandwidth parameter for categorical variables. The default value is \eqn{-1}, which enables automatic determination of the optimal bandwidth. For nominal variables and \code{nomkernel = 'aitchisonaitken'}, the maximum allowable value of \code{lambda} is \eqn{(l - 1)/l}, where \eqn{l} represents the number of categories, whereas for \code{nomkernel = 'liracine'} the maximum allowable value is \eqn{1}. For ordinal variables, the maximum allowable value of \code{lambda} is \eqn{1}, regardless of what \code{ordkernel} is being used. Argument is ignored when all variables are continuous.}

\item{scale}{A logical value indicating whether the continuous variables should be scaled to have unit variance before clustering. Defaults to \code{TRUE}. Argument is ignored when all variables are categorical.}

\item{maxiter}{The maximum number of iterations allowed for the clustering algorithm. Defaults to \eqn{100}.}

\item{nstart}{The number of random initializations to run. The best clustering solution is returned. Defaults to \eqn{100}.}

\item{contkernel}{Kernel used for continuous variables. Can be one of \code{gaussian} (default) or \code{epanechnikov}. Argument is ignored when no variables are continuous.}

\item{nomkernel}{Kernel used for nominal (unordered categorical) variables. Can be one of \code{aitchisonaitken} (default) or \code{liracine}. Argument is ignored when no variables are nominal.}

\item{ordkernel}{Kernel used for ordinal (ordered categorical) variables. Can be one of \code{liracine} (default) or \code{wangvanryzin}. Argument is ignored when no variables are ordinal.}

\item{cat_first}{A logical value indicating whether bandwidth selection is prioritised for the categorical variables, instead of the continuous. Defaults to \code{FALSE}. Set to \code{TRUE} if you suspect that the continuous variables are not informative of the cluster structure. Can only be \code{TRUE} when data is of mixed-type and all bandwidths are selected automatically (i.e. \code{s = -1}, \code{lambda = -1}).}

\item{verbose}{Logical. Defaults to \code{FALSE} to suppress progress messages. Change to \code{TRUE} to print.}
}
\value{
An object of class \code{"gibclust"} representing the final clustering result. The returned object is a list with the following components:
  \item{Cluster}{An integer vector giving the cluster assignments for each data point.}
  \item{Entropy}{A numeric value representing the entropy of the cluster assignments at convergence.}
  \item{CondEntropy}{A numeric value representing the conditional entropy of cluster assignment, given the observation weights \eqn{H(T \mid X)}.}
  \item{MutualInfo}{A numeric value representing the mutual information, \eqn{I(Y;T)}, between the original labels (\eqn{Y}) and the cluster assignments (\eqn{T}).}
  \item{InfoXT}{A numeric value representing the mutual information, \eqn{I(X;T)}, between the original observations weights (\eqn{X}) and the cluster assignments (\eqn{T}).}
  \item{beta}{A numeric vector of the final beta values used in the iterative procedure.}
  \item{alpha}{A numeric value of the strength of conditional entropy used, controlling fuzziness of the solution. This is by default equal to \eqn{0} for \code{DIBmix}.}
  \item{s}{A numeric vector of bandwidth parameters used for the continuous variables. A value of \eqn{-1} is returned if all variables are categorical.}
  \item{lambda}{A numeric vector of bandwidth parameters used for the categorical variables. A value of \eqn{-1} is returned if all variables are continuous.}
  \item{call}{The matched call.}
  \item{ncl}{Number of clusters.}
  \item{n}{Number of observations.}
  \item{iters}{Number of iterations used to obtain the returned solution.}
  \item{converged}{Logical indicating whether convergence was reached before \code{maxiter}.}
  \item{conv_tol}{Numeric convergence tolerance; by default \eqn{0} for \code{DIBmix}.}
  \item{contcols}{Indices of continuous columns in \code{X}.}
  \item{catcols}{Indices of categorical columns in \code{X}.}
  \item{kernels}{List with names of kernels used for continuous, nominal, and ordinal features.}
  
Objects of class \code{"gibclust"} support the following methods:
  \itemize{
    \item \code{print.gibclust}: Display a concise description of the cluster assignment.
    \item \code{summary.gibclust}: Show detailed information including cluster sizes, 
          information-theoretic metrics, hyperparameters, and convergence details.
    \item \code{plot.gibclust}: Produce diagnostic plots:
      \itemize{
        \item \code{type = "sizes"}: barplot of cluster sizes or hardened sizes (IB/GIB).
        \item \code{type = "info"}: barplot of entropy, conditional entropy, and mutual information.
        \item \code{type = "beta"}: trajectory of \eqn{\log \beta} over iterations (only available for hard clustering outputs obtained using \code{DIBmix}).
    }
  }
}


\details{
The \code{DIBmix} function clusters data while retaining maximal information about the original variable
distributions. The Deterministic Information Bottleneck algorithm optimizes an information-theoretic
objective that balances information preservation and compression. Bandwidth parameters for categorical
(nominal, ordinal) and continuous variables are adaptively determined if not provided. This iterative
process identifies stable and interpretable cluster assignments by maximizing mutual information while
controlling complexity. The method is well-suited for datasets with mixed-type variables and integrates
information from all variable types effectively.

The following kernel functions can be used to estimate densities for the clustering procedure. For continuous variables:

\itemize{
  \item \emph{Gaussian (RBF) kernel \insertCite{silverman_density_1998}{IBclust}:}
  \deqn{K_c\left(\frac{x - x'}{s}\right) = \frac{1}{\sqrt{2\pi}} \exp\left\{-\frac{\left(x - x'\right)^2}{2s^2}\right\}, \quad s > 0.}
  \item \emph{Epanechnikov kernel \insertCite{epanechnikov1969non}{IBclust}:}
  \deqn{K_c(x - x'; s) = \begin{cases}
    \frac{3}{4\sqrt{5}}\left(1 - \frac{(x-x')^2}{5s^2} \right), & \text{if } \frac{(x - x')^2}{s^2} < 5 \\
    0, & \text{otherwise}
\end{cases}, \quad s > 0.}
}

For nominal (unordered categorical variables):

\itemize{
\item \emph{Aitchison & Aitken kernel \insertCite{aitchison_kernel_1976}{IBclust}:}
\deqn{K_u(x = x'; \lambda) = \begin{cases}
    1 - \lambda, & \text{if } x = x' \\
    \frac{\lambda}{\ell - 1}, & \text{otherwise}
\end{cases}, \quad 0 \leq \lambda \leq \frac{\ell - 1}{\ell}.}
\item \emph{Li & Racine kernel \insertCite{ouyang2006cross}{IBclust}:}
\deqn{K_u(x = x'; \lambda) = \begin{cases}
    1, & \text{if } x = x' \\
    \lambda, & \text{otherwise}
\end{cases}, \quad 0 \leq \lambda \leq 1.}
}

For ordinal (ordered categorical) variables:

\itemize{
\item \emph{Li & Racine kernel \insertCite{li_nonparametric_2003}{IBclust}:}
\deqn{K_o(x = x'; \nu) = \begin{cases}
    1, & \text{if } x = x' \\
    \nu^{|x - x'|}, & \text{otherwise}
\end{cases}, \quad 0 \leq \nu \leq 1.}
\item \emph{Wang & van Ryzin kernel \insertCite{wang1981class}{IBclust}:}
\deqn{K_o(x = x'; \nu) = \begin{cases}
    1 - \nu, & \text{if } x = x' \\
    \frac{1-\nu}{2}\nu^{|x - x'|}, & \text{otherwise}
\end{cases}, \quad 0 \leq \nu \leq 1.}
}

The bandwidth parameters \eqn{s}, \eqn{\lambda}, and \eqn{\nu} control the smoothness of the density estimate and are automatically determined by the algorithm if not provided by the user using the approach in \insertCite{costa_dib_2025;textual}{IBclust}. \eqn{\ell} is the number of levels of the categorical variable. For ordinal variables, the lambda parameter of the function is used to define \eqn{\nu}.
}


\examples{
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
}

\author{
Efthymios Costa, Ioanna Papatsouma, Angelos Markos
}
\references{
{
\insertRef{costa_dib_2025}{IBclust}

\insertRef{aitchison_kernel_1976}{IBclust}

\insertRef{li_nonparametric_2003}{IBclust}

\insertRef{silverman_density_1998}{IBclust}

\insertRef{ouyang2006cross}{IBclust}

\insertRef{wang1981class}{IBclust}

\insertRef{epanechnikov1969non}{IBclust}
}
}
\keyword{clustering}