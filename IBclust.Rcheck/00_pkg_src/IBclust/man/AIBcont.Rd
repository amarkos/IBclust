\name{AIBcont}
\alias{AIBcont}
\title{Cluster Continuous Data Using the Agglomerative Information Bottleneck Algorithm}
\description{
The \code{AIBcont} function implements the Agglomerative Information Bottleneck (AIB) algorithm
for hierarchical clustering of datasets containing categorical variables. This method merges clusters
so that information retention is maximised at each step to create meaningful clusters,
leveraging bandwidth parameters to handle
different categorical data types (nominal and ordinal) effectively \insertCite{slonim_aib_1999}{IBclust}.
}
\usage{
AIBcont(X, s = -1, scale = TRUE, contkernel = "gaussian")
}
\arguments{
  \item{X}{
    A data frame containing the categorical data to be clustered. All variables should be categorical,
    either \code{factor} (for nominal variables) or \code{ordered} (for ordinal variables).
  }
  \item{s}{
  A numeric value or vector specifying the bandwidth parameter(s) for continuous variables. The values must be greater than \eqn{0}. The default value is \eqn{-1}, which enables the automatic selection of optimal bandwidth(s).
  }
  \item{scale}{
    A logical value indicating whether the continuous variables should be scaled to have unit variance before clustering. Defaults to \code{TRUE}.
  }
  \item{contkernel}{
    Kernel used for continuous variables. Can be one of \code{gaussian} (default) or \code{epanechnikov}.
  }
}

\value{
A list containing the following elements:
  \item{merges}{A data frame with 2 columns and \eqn{n} rows, showing which observations are merged at each step.}
  \item{merge_costs}{A numeric vector tracking the cost incurred by each merge \eqn{I(Z_m ; Y) - I(Z_{m-1} ; Y)}.}
  \item{partitions}{A list containing \eqn{n} sub-lists. Each sub-list includes the cluster partition at each step.}
  \item{I_Z_Y}{A numeric vector including the mutual information \eqn{I(Z_{m}; Y)} as the number of clusters \eqn{m} increases.}
  \item{I_X_Y}{A numeric value of the mutual information \eqn{I(X; Y)} between observation indices and location.}
  \item{info_ret}{A numeric vector of length \eqn{n} including the fraction of the original information retained after each merge.}
  \item{dendrogram}{A dendrogram visualising the cluster hierarchy. The height is determined by the cost of cluster merges.}
}

\details{
The \code{AIBcont} function applies the Agglomerative Information Bottleneck algorithm to do hierarchical clustering of datasets containing only continuous variables, both nominal and ordinal. The algorithm uses an information-theoretic criterion to merge clusters so that information retention is maximised at each step to create meaningful clusters with maximal information about the original distribution.

To estimate the distributions of continuous features, the function utilizes specialized kernel functions:

\itemize{
  \item \emph{Gaussian (RBF) kernel \insertCite{silverman_density_1998}{IBclust}:}
  \deqn{K_c\left(\frac{x - x'}{s}\right) = \frac{1}{\sqrt{2\pi}} \exp\left\{-\frac{\left(x - x'\right)^2}{2s^2}\right\}, \quad s > 0.}
  \item \emph{Epanechnikov kernel \insertCite{epanechnikov1969non}{IBclust}:}
  \deqn{K_c(x - x'; s) = \begin{cases}
    \frac{3}{4\sqrt{5}}\left(1 - \frac{(x-x')^2}{5s^2} \right), & \text{if } \frac{(x - x')^2}{s^2} < 5 \\
    0, & \text{otherwise}
\end{cases}, \quad s > 0.}
}

The bandwidth parameter \eqn{s}, which controls the smoothness of the density estimate, is automatically determined by the algorithm if not provided by the user.
}


\examples{
# Generate simulated continuous data
set.seed(123)
X <- matrix(rnorm(1000), ncol = 5)  # 200 observations, 5 features

# Run AIBcont with automatic bandwidth selection 
result <- AIBcont(X = X, s = -1, scale = TRUE)

# Print clustering results
plot(result$dendrogram, xlab = "", sub = "")  # Plot dendrogram
}
\seealso{
\code{\link{AIBmix}}, \code{\link{AIBcat}}
}
\author{
Efthymios Costa, Ioanna Papatsouma, Angelos Markos
}
\references{
{
\insertRef{slonim_aib_1999}{IBclust}

\insertRef{silverman_density_1998}{IBclust}

\insertRef{epanechnikov1969non}{IBclust}
}
}
\keyword{clustering}
