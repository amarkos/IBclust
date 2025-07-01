\name{AIBmix}
\alias{AIBmix}
\title{Agglomerative Information Bottleneck Clustering for Mixed-Type Data}
\description{
The \code{AIBmix} function implements the Agglomerative Information Bottleneck (AIB) algorithm
for hierarchical clustering of datasets containing mixed-type variables, including categorical (nominal and ordinal)
and continuous variables. This method merges clusters so that information retention is maximised at each step to create meaningful clusters,
leveraging bandwidth parameters to handle different categorical data types (nominal and ordinal) effectively \insertCite{slonim_aib_1999}{IBclust}.
}
\usage{
AIBmix(X, catcols, contcols, lambda = -1, s = -1, scale = TRUE)
}
\arguments{
  \item{X}{
    A data frame containing the categorical data to be clustered. All variables should be categorical,
    either \code{factor} (for nominal variables) or \code{ordered} (for ordinal variables).
  }
  \item{catcols}{A vector indicating the indices of the categorical variables in \code{X}.}
  \item{contcols}{A vector indicating the indices of the continuous variables in \code{X}.}
  \item{lambda}{
    A numeric value or vector specifying the bandwidth parameter for categorical variables. The default value is \eqn{-1}, which enables automatic determination of the optimal bandwidth. For nominal variables, the maximum allowable value of \code{lambda} is \eqn{(l - 1)/l}, where \eqn{l} represents the number of categories. For ordinal variables, the maximum allowable value of \code{lambda} is \eqn{1}.
  }
  \item{s}{
  A numeric value or vector specifying the bandwidth parameter(s) for continuous variables. The values must be greater than \eqn{0}. The default value is \eqn{-1}, which enables the automatic selection of optimal bandwidth(s).
  }
  \item{scale}{
    A logical value indicating whether the continuous variables should be scaled to have unit variance before clustering. Defaults to \code{TRUE}.
  }
}

\value{
A list containing the following elements:
  \item{merges}{A data frame with 2 columns and \eqn{n} rows, showing which observations are merged at each step.}
  \item{merge_costs}{A numeric vector tracking the cost incurred by each merge \eqn{I(Z_{m} ; Y) - I(Z_{m-1} ; Y)}.}
  \item{partitions}{A list containing \eqn{n} sub-lists. Each sub-list includes the cluster partition at each step.}
  \item{I_Z_Y}{A numeric vector including the mutual information \eqn{I(Z_{m}; Y)} as the number of clusters \eqn{m} increases.}
  \item{I_X_Y}{A numeric value of the mutual information \eqn{I(X; Y)} between observation indices and location.}
  \item{info_ret}{A numeric vector of length \eqn{n} including the fraction of the original information retained after each merge.}
  \item{dendrogram}{A dendrogram visualising the cluster hierarchy. The height is determined by the cost of cluster merges.}
}

\details{
The \code{AIBmix} function produces a hierarchical agglomerative clustering of the data while retaining maximal information about the original variable
distributions. The Agglomerative Information Bottleneck algorithm uses an information-theoretic criterion to merge clusters so that information retention is maximised at each step,
hence creating meaningful clusters with maximal information about the original distribution. Bandwidth parameters for categorical
(nominal, ordinal) and continuous variables are adaptively determined if not provided. This process identifies stable and interpretable cluster assignments by maximizing mutual information while
controlling complexity. The method is well-suited for datasets with mixed-type variables and integrates
information from all variable types effectively. 

The following kernel functions are used to estimate densities for the clustering procedure:
\itemize{
    \item \emph{Continuous variables: Gaussian kernel}
    \deqn{K_c\left(\frac{x-x'}{s}\right) = \frac{1}{\sqrt{2\pi}} \exp\left\{ - \frac{\left(x-x'\right)^2}{2s^2} \right\}, \quad s > 0.}
    \item \emph{Nominal categorical variables: Aitchison & Aitken kernel}
    \deqn{K_u\left(x = x' ; \lambda\right) = \begin{cases}
          1-\lambda & \text{if } x = x' \\
          \frac{\lambda}{\ell-1} & \text{otherwise}
        \end{cases}, \quad 0 \leq \lambda \leq \frac{\ell-1}{\ell}.}
    \item \emph{Ordinal categorical variables: Li & Racine kernel}
    \deqn{K_o\left(x = x' ; \nu\right) = \begin{cases}
          1 & \text{if } x = x' \\
          \nu^{|x - x'|} & \text{otherwise}
        \end{cases}, \quad 0 \leq \nu \leq 1.}
}
}


\examples{
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
}
\seealso{
\code{\link{AIBcat}}, \code{\link{AIBcont}}
}
\author{
Efthymios Costa, Ioanna Papatsouma, Angelos Markos
}
\references{
{
\insertRef{slonim_aib_1999}{IBclust}

\insertRef{aitchison_kernel_1976}{IBclust}

\insertRef{li_nonparametric_2003}{IBclust}

\insertRef{silverman_density_1998}{IBclust}
}
}
\keyword{clustering}
